! Copyright (c) 2013-2014, 2019-2022 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2014, 2018 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016-2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2022 Jana Gericke <jana.gericke@dlr.de>
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! **************************************************************************** !

!> This module provides the facility to control the execution of the solver.
!!
!! It is configured by the `sim_control` table in the configuration.
!! From this table three settings will be read:
!!
!! * `time_control` defines the control of time stepping see
!!   [[tem_timeControl_module]].
!! * `abort_criteria` defines additional abort criteria to be heeded see
!!   [[tem_abortCriteria_module]].
!! * `delay_check` defines whether checks should be delayed (default is false).
!!
!! If there is no `sim_control` found in the configuration, an attempt will be
!! made to read the `time_control` table instead, while an empty
!! `abort_criteria` table is assumed.
!!
!! Thus, the general structure of the `sim_control` configuration has the
!! following form:
!!
!!```lua
!!  sim_control = {
!!    time_control = {},
!!    abort_criteria = {},
!!    delay_check = false
!!  }
!!```
!!
!! Alternatively, if there are no abort criteria to be specified, and
!! `delay_check` should not be activated in general, just the
!! `time_control` table may be specified:
!!
!!```lua
!!    time_control = {}
!!```
!!
!! See the [[tem_timeControl_module]] for details on the `time_control` table.
!! See the [[tem_abortCriteria_module]] for details on the `time_control` table.
!!
!! The `delay_check` flag indicates whether to use a nonblocking allreduce
!! for the synchronization of status flags, to decide on triggered events during
!! the computation. It relaxes the synchronization requirements but introduces
!! a delay by one check interval until all processes get set status bits.
!!
!! A simple complete example without checks for steady state would be:
!!
!!```lua
!!  sim_control = {
!!    time_control = {
!!      min = 0,
!!      max = 10.0,
!!      interval = {iter = 5}
!!    },
!!    abort_criteria = {
!!      stop_file = 'stop',
!!    },
!!    delay_check = false
!!  }
!!```
!!
module tem_simControl_module
  use env_module, only: rk, labelLen

  use tem_comm_env_module, only: tem_comm_env_type
  use tem_time_module, only: tem_time_type,     &
    &                        tem_time_reset,    &
    &                        tem_time_advance,  &
    &                        tem_time_dump,     &
    &                        tem_time_set_clock,&
    &                        tem_time_out,      &
    &                        tem_time_sim_id,   &
    &                        tem_time_iter_id,  &
    &                        tem_time_clock_id, &
    &                        tem_time_n_ids

  use tem_timeControl_module, only: tem_timeControl_type,         &
    &                               tem_timeControl_load,         &
    &                               tem_timeControl_out,          &
    &                               tem_timeControl_dump,         &
    &                               tem_timeControl_start_at_sim, &
    &                               tem_timeControl_reached_max,  &
    &                               tem_timeControl_triggered,    &
    &                               tem_timeControl_update,       &
    &                               tem_timeControl_reset_trigger

  use tem_status_module, only: tem_status_type,        &
    &                          tem_status_clear,       &
    &                          tem_status_dump,        &
    &                          tem_status_communicate, &
    &                          tem_status_communicate_delayed, &
    &                          tem_stat_nFlags,        &
    &                          tem_stat_max_sim,       &
    &                          tem_stat_max_iter,      &
    &                          tem_stat_max_clock,     &
    &                          tem_stat_interval,      &
    &                          tem_stat_stop_file

  use tem_abortCriteria_module, only: tem_abortCriteria_type, &
    &                                 tem_abortCriteria_new,  &
    &                                 tem_abortCriteria_load, &
    &                                 tem_abortCriteria_out,  &
    &                                 tem_abortCriteria_dump, &
    &                                 tem_stop_file_exists,   &
    &                                 tem_solverAborts_type

  use tem_convergence_module,       only: tem_convergence_reset

  use tem_timer_module,       only: tem_addTimer,   &
    &                               tem_startTimer, &
    &                               tem_stopTimer

  use tem_logging_module,    only: logUnit

  use aotus_module,     only: flu_State, aot_get_val
  use aot_table_module, only: aot_table_open, &
    &                         aot_table_close

  use aot_out_module,   only: aot_out_type,       &
    &                         aot_out_val,        &
    &                         aot_out_open_table, &
    &                         aot_out_close_table

  use mpi

  implicit none

  private

  public :: tem_simControl_type
  public :: tem_simControl_start
  public :: tem_simControl_load
  public :: tem_simControl_out
  public :: tem_simControl_dump
  public :: tem_simControl_dump_now
  public :: tem_simControl_syncUpdate
  public :: tem_simControl_clearStat
  public :: tem_simControl_steadyState_reset


  !> Data structure to describe the overall control of a simulation.
  !!
  !! This comprises the current time in all available definitions.
  type tem_simControl_type
    !> Representation of the current time.
    type(tem_time_type) :: now

    !> Time control, when the simulation should end, and definition of
    !! special interval, at which regular actions should take place.
    !!
    !! The minimum setting has no significance here and is always set to
    !! the time, provided when loading the sim control.
    type(tem_timeControl_type) :: timeControl

    !> Further abort criteria.
    type(tem_abortCriteria_type) :: abortCriteria

    !> Flag collection to describe the status of the simulation.
    type(tem_status_type) :: status

    !> Use nonblocking operations for gobal checks and delay evaluation by
    !! one check interval (see timeControl%check_iter)
    logical :: delay_check = .false.

    !> Handle for the syncUpdate timer to measure the time spent on syncUpdate
    !! calls.
    integer :: syncUpdate_timer
  end type tem_simControl_type


contains


  ! ************************************************************************ !
  !> Start a sim control by resetting its time object.
  !!
  !! Note, that the actual control needs to be filled afterwards with
  !! tem_simControl_load.
  subroutine tem_simControl_start(me)
    ! -------------------------------------------------------------------- !
    !> The simulation control structure to start.
    type(tem_simControl_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !

    call tem_time_reset(me%now)

    call tem_status_clear(me = me%status)

  end subroutine tem_simControl_start
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Load sim control settings from a configuration script.
  !!
  !! The sim control should be started beforhand to ensure a sane setting of
  !! the current time.
  !! The main setting here, is the time_control, which is also attempted to
  !! be read directly, if there is no sim_control table provided.
  !! Solvers may pass solverAborts to load additional abort criteria that
  !! are to be loaded from the configuration.
  subroutine tem_simControl_load(me, conf, parent, key, solverAborts)
    ! -------------------------------------------------------------------- !
    !> Simulation control parameters to set.
    type(tem_simControl_type), intent(inout) :: me

    !> Handle to the configuration script to load the settings from.
    type(flu_state) :: conf

    !> Potential parent table, in which the simulation control table is to be
    !! found.
    integer, intent(in), optional :: parent

    !> Name for the simulation control table. Default is 'sim_control'.
    character(len=*), optional :: key

    !> Solver specific abort criteria to load.
    class(tem_solverAborts_type), intent(inout), optional :: solverAborts
    ! -------------------------------------------------------------------- !
    character(len=labelLen) :: loc_key
    integer :: thandle
    integer :: iErr
    ! -------------------------------------------------------------------- !

    loc_key = 'sim_control'
    if (present(key)) loc_key = key

    call aot_table_open( L       = conf,    &
      &                  parent  = parent,  &
      &                  thandle = thandle, &
      &                  key     = loc_key  )

    call aot_get_val(L       = conf,           &
      &              thandle = thandle,        &
      &              val     = me%delay_check, &
      &              key     = 'delay_check',  &
      &              default = .false.,        &
      &              ErrCode = iErr            )

    if (thandle /= 0) then
      ! There is a sim control table in the configuration.
      call tem_timeControl_load(me          = me%timeControl, &
        &                       conf        = conf,           &
        &                       parent      = thandle,        &
        &                       delay_check = me%delay_check  )

      call tem_abortCriteria_load( me     = me%abortCriteria,  &
        &                          conf   = conf,              &
        &                          parent = thandle,           &
        &                          solverAborts = solverAborts )
    else
      ! No sim control table found, try to load the time control table itself.
      call tem_timeControl_load(me%timeControl, conf, parent)
      me%abortCriteria = tem_abortCriteria_new()
    end if

    if (me%delay_check) then
      write(logUnit(1),*) 'Delaying status checks by one interval' &
        // '(delay_check=True).'
    end if

    call aot_table_close(L = conf, thandle = thandle)

    ! Overwrite the min of the simulation time range to point to the current
    ! simulation time.
    call tem_timeControl_start_at_sim(me = me%timeControl, now = me%now)

    call tem_addTimer( timerHandle = me%syncUpdate_timer, &
      &                timerName = 'tem_syncUpdate'       )

  end subroutine tem_simControl_load
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Write sim control settings to a configuration script.
  subroutine tem_simControl_out(me, conf, key)
    ! -------------------------------------------------------------------- !
    !> The simulation control settings to write to a Lua table.
    type(tem_simControl_type), intent(inout) :: me

    !> Handle for the Lua script to write to.
    type(aot_out_type) :: conf

    !> Name for the simulation control table. Default is sim_control.
    character(len=*), optional :: key
    ! -------------------------------------------------------------------- !
    character(len=labelLen) :: loc_key
    ! -------------------------------------------------------------------- !

    loc_key = 'sim_control'
    if (present(key)) loc_key = key

    call aot_out_open_table( put_conf = conf,   &
      &                      tname    = loc_key )

    call tem_time_out( me   = me%now, &
      &                conf = conf    )

    call tem_timeControl_out( me   = me%timeControl, &
      &                       conf = conf            )

    call tem_abortCriteria_out( me   = me%abortCriteria, &
      &                         conf = conf              )

    call aot_out_val( put_conf = conf,           &
      &               val      = me%delay_check, &
      &               vname    = 'delay_check'   )

    call aot_out_close_table( put_conf = conf )

  end subroutine tem_simControl_out
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Dump the current time (now) of the simControl to outUnit.
  !!
  !! This also updates the clock in now to really show the current time.
  subroutine tem_simControl_dump_now(me, outUnit)
    ! -------------------------------------------------------------------- !
    !> Simulation control settings to write to outUnit.
    type(tem_simControl_type), intent(inout) :: me

    !> File unit to write to.
    integer, intent(in) :: outUnit
    ! -------------------------------------------------------------------- !

    call tem_time_set_clock(me%now)
    call tem_time_dump(me%now, outUnit)

  end subroutine tem_simControl_dump_now
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Dump simcontrol information to the specified outUnit.
  subroutine tem_simControl_dump(me, outUnit)
    ! -------------------------------------------------------------------- !
    !> Simulation control settings to write to outUnit.
    type(tem_simControl_type), intent(inout) :: me

    !> File unit to write to.
    integer, intent(in) :: outUnit
    ! -------------------------------------------------------------------- !

    write(outUnit,*) '+-------------------------------------------+'
    write(outUnit,*) ' - Now:'
    call tem_simControl_dump_now(me, outUnit)
    write(outUnit,*) ''
    write(outUnit,*) ' - Time Control:'
    call tem_timeControl_dump(me%timeControl, outUnit)
    write(outUnit,*) ''
    write(outUnit,*) ' - Abort Criteria:'
    call tem_abortCriteria_dump(me%abortCriteria, outUnit)
    write(outUnit,*) ''
    write(outUnit,*) ' - Status:'
    call tem_status_dump(me%status, outUnit)
    write(outUnit,*) '+-------------------------------------------+'

  end subroutine tem_simControl_dump
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Synchronize the status bits across all processes and update the time.
  subroutine tem_simControl_syncUpdate(me, proc, dt, d_iter, outUnit)
    ! -------------------------------------------------------------------- !
    !> Simulation control information.
    type(tem_simControl_type), intent(inout) :: me

    !> Unit to write messages to.
    !!
    !! If this argument is present, the current time will be printed whenever
    !! the interval of the simControl is triggered.
    integer, intent(in), optional :: outUnit

    !> Communicator to use for the communication of status flags.
    type(tem_comm_env_type), intent(in) :: proc

    !> Time step to use for updating the simulation time.
    !!
    !! If this is not given, no advance of the time will be done.
    real(kind=rk), intent(in), optional :: dt

    !> Number of iterations to add to the current number of iterations.
    !! (Default: 1)
    integer, intent(in), optional :: d_iter
    ! -------------------------------------------------------------------- !
    logical :: max_reached(tem_time_n_ids)
    logical :: stat_interval
    logical :: out_interval
    ! -------------------------------------------------------------------- !

    call tem_startTimer(timerHandle = me%syncUpdate_timer)

    if (present(dt)) then
      call tem_time_advance( me     = me%now, &
        &                    sim_dt = dt,     &
        &                    iter   = d_iter  )
    end if

    max_reached = tem_timeControl_reached_max(me%timeControl, me%now)

    me%status%bits(tem_stat_max_sim)   = max_reached(tem_time_sim_id)
    me%status%bits(tem_stat_max_iter)  = max_reached(tem_time_iter_id)
    if (mod(me%now%iter, me%timeControl%check_iter) == 0) then
      me%status%bits(tem_stat_max_clock) = max_reached(tem_time_clock_id)

      me%status%bits(tem_stat_interval) = tem_timeControl_triggered( &
        &                                   me  = me%timeControl,    &
        &                                   now = me%now             )

      me%status%bits(tem_stat_stop_file) &
        &  = tem_stop_file_exists( abortCriteria = me%abortCriteria, &
        &                          rank          = proc%rank         )

      stat_interval = me%status%bits(tem_stat_interval)

      if (me%delay_check) then
        call tem_status_communicate_delayed(me = me%status, comm = proc%comm)
        me%status%bits(tem_stat_max_sim)  = max_reached(tem_time_sim_id)
        me%status%bits(tem_stat_max_iter) = max_reached(tem_time_iter_id)
        out_interval = stat_interval
      else
        call tem_status_communicate(me = me%status, comm = proc%comm)
        out_interval = me%status%bits(tem_stat_interval)
      end if

      if (present(outUnit) .and. out_interval) then
        call tem_time_set_clock(me%now)
        call tem_time_dump(me%now, outUnit)
      end if

      call tem_timeControl_update( me             = me%timeControl, &
        &                          now            = me%now,         &
        &                          hasTriggered   = out_interval,   &
        &                          localTriggered = stat_interval   )

    end if

    call tem_stopTimer(timerHandle = me%syncUpdate_timer)

  end subroutine tem_simControl_syncUpdate
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Clear the status flags in the simcontrol.
  subroutine tem_simControl_clearStat(me)
    ! -------------------------------------------------------------------- !
    !> Simulation control information.
    type(tem_simControl_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !

    if (mod(me%now%iter, me%timeControl%check_iter) == 0) then
      call tem_status_clear(me%status)
    end if

  end subroutine tem_simControl_clearStat
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Reset trigger, status bit and current time for steady state
  subroutine tem_simControl_steadyState_reset(me)
    ! -------------------------------------------------------------------- !
    type(tem_simControl_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !
    ! clear status bit
    call tem_status_clear(me%status)

    ! reset current time
    call tem_time_reset(me%now)

    ! Run steady state solver untill the solution convergences
    me%timeControl%max%sim = huge(me%timeControl%max%sim)
    me%timeControl%max%iter = huge(me%timeControl%max%iter)

    ! reset simcontrol trigger
    call tem_timeControl_reset_trigger(me%timeControl)

    ! reset convergence to check for new steady state
    call tem_convergence_reset(me%abortCriteria%convergence )

  end subroutine tem_simControl_steadyState_reset
  ! ************************************************************************ !


end module tem_simControl_module
