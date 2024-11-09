! Copyright (c) 2012-2014,2017,2019-2020,2022 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2014, 2018 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012 Laura Didinger <l.didinger@grs-sim.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012, 2014 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2013-2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2013 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2014, 2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
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
!> This module provides functions to control explicit time stepping solvers.
!!
!! It makes use of the [[tem_time_module]] with its definition of time in terms
!! of various measurements, allowing the control of events based on these units.
!! A time control is always active within a certain range, given by its min and
!! max value, and emits actions to take based on the given interval.
!! If multiple measurements are given for any of the measures, this time is
!! always taken by that unit which occurs first.
!! All three settings (min, max and interval) might use different time
!! measures independently.
!!
!! The `time_control` table takes the following form:
!!
!!```lua
!!     time_control = {
!!       min = {},
!!       max = {},
!!       interval = {},
!!       check_iter = 123,
!!       delay_check = false
!!     }
!!```
!!
!! Here `min`, `max` and `interval` are each time definitions, see
!! [[tem_time_module]]. The time definition is providing, a measure of time
!! in terms of simulation time, number of iterations and passed running time
!! (clock): `{sim = 1.0, iter = 1, clock = 60.0}`.
!!
!! `max` describes the end of the time span. If this `time_control` is for the
!! overall simulation it defines the end of the run, after any of the time
!! measures is reached, the simulation will stop.
!!
!! `min` describes the beginning of the timespan for this control object.
!! This may for example be used to start writing tracking data only after
!! a certain point in time. Again, if multiple time definitions
!! (sim, iter, clock) are defined, whichever first will be encountered will
!! start the timespan of this control.
!!
!! `interval` describes intervals within in the `min` and `max` points in time
!! where something is to happen, for example that a restart file is to be
!! written.
!!
!! If a time definition is not provided as a table the given value is
!! interpreted as specifying a time in terms of simulated time (`sim`).
!! See [[tem_time_module]] for details on the time definition.
!!
!! `check_iter` allows you to control, how often the trigger status of the
!! time control is to be checked in terms of iterations.
!! These checks involve communication and may have a performance impact if
!! performed in every iteration and the iterations are very short.
!! By setting the `check_iter` to some larger value the communications can
!! be decreased as they will only be performed every `check_iter` iteration.
!! Be aware that increasing the `check_iter` setting also decreases the
!! accuracy for the time control. There might be performed multiple iterations
!! beyond the intended specified trigger in this case.
!! Typically it is not necessary to specify the `check_iter` setting. Without
!! providing it a default of one will be used.
!!
!! `delay_check` provides the option to delay the evaluation for the clock time
!! trigger by one check interval by using a nonblocking allreduce.
!! It relaxes the synchronization requirements, but introduces a delay in the
!! actual triggering of the clock event.
!! The default of this setting depends on the corresponding `delay_check`
!! setting in simControl. Thus, if it is set in the simControl table itself,
!! that setting will also be used here unless overwritten by explicitly
!! providing it here.
!!
module tem_timeControl_module

  use mpi

  use env_module, only: rk, labelLen
  use tem_tools_module, only: tem_horizontalSpacer
  use tem_time_module, only: tem_time_type,          &
    &                        tem_time_load,          &
    &                        tem_time_out,           &
    &                        tem_time_dump,          &
    &                        tem_time_never,         &
    &                        tem_time_last_interval, &
    &                        tem_time_default_zero,  &
    &                        operator(+),            &
    &                        operator(-),            &
    &                        max,                    &
    &                        tem_time_ge_trigger,    &
    &                        tem_time_gt_trigger,    &
    &                        tem_time_ge_each,       &
    &                        tem_time_needs_reduce,  &
    &                        tem_time_n_ids
  use tem_logging_module,    only: logUnit

  use aotus_module, only: flu_State
  use aot_table_module, only: aot_table_open,  &
    &                         aot_table_close, &
    &                         aot_get_val

  use aot_out_module, only: aot_out_type,       &
    &                       aot_out_val,        &
    &                       aot_out_open_table, &
    &                       aot_out_close_table

  implicit none

  private


  public :: tem_timeControl_type
  public :: tem_timeControl_load
  public :: tem_timeControl_dump
  public :: tem_timeControl_out
  public :: tem_timeControl_triggered
  public :: tem_timeControl_globalTriggered
  public :: tem_timeControl_update
  public :: tem_timeControl_check
  public :: tem_timeControl_start_at_sim
  public :: tem_timeControl_reached_max
  public :: tem_timeControl_align_trigger
  public :: tem_timeControl_reset_trigger


  !> Definition of a time control.
  !!
  !! The control is active in the range of time between min and max.
  !! It will trigger its action after a time interval specified in interval.
  !! For all time definitions always the one that occurs first is being used.
  type tem_timeControl_type
    !> Minimal point in time, from where on, this control should be active.
    !! Whichever time definition happens first will be used.
    type(tem_time_type) :: min

    !> Maximal point in time, after which the control should not be active
    !! anymore. Whichever time definition happens first will be used.
    type(tem_time_type) :: max

    !> A regular interval at which an action should be triggered between
    !! min and max.
    type(tem_time_type) :: interval

    !> Keep track of the next point in time, at which an action should be
    !! triggered by this control.
    type(tem_time_type) :: trigger

    !> Configuration flag, whether to use nonblocking allreduces to determine
    !! whether an event has triggered across all processes.
    !!
    !! Note that this delays the reaction on the trigger by one check interval
    !! (given by check_iter in number of iterations).
    logical :: delay_check = .false.

    !> Trigger checking can involve communication and is potentially hurting
    !! the performance.
    !!
    !! With this setting, the iteration interval at which these trigger
    !! updates should be done, can be controlled.
    !! Per default each iteration a check is done, but if this is too
    !! frequent, it can be increased here.
    !! However, it should be noted, that all trigger checks are only done
    !! every check_iter iteration.
    integer :: check_iter = 1

    !> IAllreduce request handle
    integer :: check_request = MPI_REQUEST_NULL

    !> Flag to indicate whether trigger was globally activated
    logical :: globally_triggered

    !> Flag to indicate if this control object needs a MPI_reduce to determine
    !! trigger status.
    logical :: needs_reduce

    !> Flag that indicates whether the minimal point in time specified in min
    !! has already been reached.
    logical :: min_reached = .false.
  end type tem_timeControl_type

  integer, parameter :: sim = 1
  integer, parameter :: iter = 2
  integer, parameter :: clock = 3


contains


  ! ************************************************************************ !
  !> Load a time control definition from a Lua script.
  !!
  !! The time control description me is loaded from conf within parent and
  !! under the name key.
  !! If no key is provided the name is assumed to be 'time_control'.
  !! If the table is not found at all, all components of the control are set
  !! to never.
  subroutine tem_timeControl_load(me, conf, parent, key, delay_check)
    ! -------------------------------------------------------------------- !
    !> Time control definition to load from a Lua config.
    type(tem_timeControl_type), intent(out) :: me

    !> Handle for the Lua script.
    type(flu_state) :: conf

    !> Parent table to read from.
    integer, intent(in), optional :: parent

    !> Name of the time control table. Default: 'time_control'
    character(len=*), intent(in), optional :: key

    !> Default setting for the delay_check.
    !!
    !! If set to true the check will use a nonblocking iAllreduce and delay
    !! the evaluation by one check_iter block.
    !! This setting may be overwritten by the user in the timecontrol block.
    logical, intent(in), optional :: delay_check
    ! -------------------------------------------------------------------- !
    type(tem_time_type) :: usermin
    integer :: thandle, iErr
    character(len=labelLen) :: localKey
    logical :: def_delay
    ! -------------------------------------------------------------------- !

    if (present(key)) then
      localKey = key
    else
      localKey = 'time_control'
    end if

    if (present(delay_check)) then
      def_delay = delay_check
    else
      def_delay = .false.
    end if

    call aot_table_open( L       = conf,    &
      &                  parent  = parent,  &
      &                  thandle = thandle, &
      &                  key     = localKey )

    if (thandle /= 0) then

      call tem_time_load(me     = usermin, &
        &                conf   = conf,    &
        &                key    = 'min',   &
        &                parent = thandle  )

      call tem_time_load(me     = me%max, &
        &                conf   = conf,   &
        &                key    = 'max',  &
        &                parent = thandle )

      call tem_time_load(me     = me%interval, &
        &                conf   = conf,        &
        &                key    = 'interval',  &
        &                parent = thandle      )

      call aot_get_val(L       = conf,          &
        &              thandle = thandle,       &
        &              val     = me%check_iter, &
        &              key     = 'check_iter',  &
        &              default = 1,             &
        &              ErrCode = iErr           )

      call aot_get_val(L       = conf,           &
        &              thandle = thandle,        &
        &              val     = me%delay_check, &
        &              key     = 'delay_check',  &
        &              default = def_delay,      &
        &              ErrCode = iErr            )

      if (me%delay_check) then
        write(logUnit(1),*) 'Delaying clock checks by one interval' &
          // '(delay_check=True).'
      end if

      ! check iter can not be smaller than 1
      if ( me%check_iter < 1 ) then
        me%check_iter = 1
      end if

    else

      ! No table defined, set all times to never.
      usermin     = tem_time_never()
      me%max      = tem_time_never()
      me%interval = tem_time_never()
      me%check_iter = 1

    end if

    call aot_table_close(L = conf, thandle = thandle)

    me%min = tem_time_default_zero( time       = usermin,    &
      &                             dependency = me%interval )

    me%trigger = me%min
    me%min_reached = .false.

    me%needs_reduce = (     tem_time_needs_reduce(me%min)      &
      &                .or. tem_time_needs_reduce(me%max)      &
      &                .or. tem_time_needs_reduce(me%interval) )

  end subroutine tem_timeControl_load
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Align the trigger to intervals since min.
  !!
  !! Only the time components given in the configuration will be considered
  !! for the alignment the other components remain untouched.
  subroutine tem_timeControl_align_trigger(me, conf, now, parent, key)
    ! -------------------------------------------------------------------- !
    !> Time control definition to load from a Lua config.
    type(tem_timeControl_type), intent(inout) :: me

    !> Handle for the Lua script.
    type(flu_state) :: conf

    !> Current point in time to find alignement of trigger.
    type(tem_time_type), intent(in) :: now

    !> Parent table to read from.
    integer, intent(in), optional :: parent

    !> Name of the time control table. Default: 'time_control'
    character(len=*), intent(in), optional :: key
    ! -------------------------------------------------------------------- !
    integer :: thandle
    character(len=labelLen) :: localKey
    logical :: alignmask(3)
    type(tem_time_type) :: align_interval
    ! -------------------------------------------------------------------- !

    if (present(key)) then
      localKey = key
    else
      localKey = 'time_control'
    endif

    call aot_table_open( L       = conf,    &
      &                  parent  = parent,  &
      &                  thandle = thandle, &
      &                  key     = localKey )

    if (thandle /= 0) then
      call load_alignmask( mask   = alignmask,       &
        &                  conf   = conf,            &
        &                  key    = 'align_trigger', &
        &                  parent = thandle          )
    else
      alignmask = .false.
    end if

    align_interval%sim = 0.0_rk
    align_interval%iter = 0
    align_interval%clock = 0.0_rk

    if (alignmask(sim)) align_interval%sim = me%interval%sim
    if (alignmask(iter)) align_interval%iter = me%interval%iter
    if (alignmask(clock)) align_interval%clock = me%interval%clock

    me%trigger = me%min + tem_time_last_interval( now      = now-me%min,    &
      &                                           interval = align_interval )

    me%min_reached = tem_time_ge_trigger(now, me%min)

  end subroutine tem_timeControl_align_trigger
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine load_alignmask(mask, conf, key, parent)
    !> Time to be read from the Lua script
    logical, intent(out) :: mask(3)

    !> Handle to the Lua script.
    type(flu_state), intent(inout) :: conf

    !> Name of the table containing the time definition. Default: 'time'.
    character(len=*), intent(in) :: key

    !> Handle to the parent table.
    integer, intent(in), optional :: parent
    ! -------------------------------------------------------------------- !
    integer :: iErr
    integer :: thandle
    ! -------------------------------------------------------------------- !

    call aot_table_open(L       = conf,     &
      &                 parent  = parent,   &
      &                 thandle = thandle,  &
      &                 key     = trim(key) )

    if (thandle /= 0) then
      ! The mask is given as a table load its components accordingly.
      call aot_get_val(L       = conf,      &
        &              thandle = thandle,   &
        &              val     = mask(sim), &
        &              key     = 'sim',     &
        &              default = .false.,   &
        &              ErrCode = iErr       )

      call aot_get_val(L       = conf,       &
        &              thandle = thandle,    &
        &              val     = mask(iter), &
        &              key     = 'iter',     &
        &              default = .false.,    &
        &              ErrCode = iErr        )

      call aot_get_val(L       = conf,        &
        &              thandle = thandle,     &
        &              val     = mask(clock), &
        &              key     = 'clock',     &
        &              default = .false.,     &
        &              ErrCode = iErr         )
    else
      ! The mask is not given as a table, try to interpret it as a setting for
      ! the simtime.
      mask = .false.
      call aot_get_val(L       = conf,      &
        &              thandle = parent,    &
        &              key     = trim(key), &
        &              val     = mask(sim), &
        &              default = .false.,   &
        &              ErrCode = iErr       )
    end if

    call aot_table_close(conf, thandle)

  end subroutine load_alignmask
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Write a time control definition to a Lua script.
  subroutine tem_timeControl_out(me, conf, key)
    ! -------------------------------------------------------------------- !
    !> Time control definition to write to a Lua config.
    type(tem_timeControl_type), intent(in) :: me

    !> Handle for the Lua script.
    type(aot_out_type) :: conf

    !> Name of the time control table. Default: 'time_control'
    character(len=*), intent(in), optional :: key
    ! -------------------------------------------------------------------- !
    character(len=labelLen) :: localKey
    ! -------------------------------------------------------------------- !

    if (present(key)) then
      localKey = key
    else
      localKey = 'time_control'
    endif

    call aot_out_open_table( put_conf = conf,    &
      &                      tname    = localkey )

    call tem_time_out(me     = me%min, &
      &               conf   = conf,   &
      &               key    = 'min'   )

    call tem_time_out(me     = me%max, &
      &               conf   = conf,   &
      &               key    = 'max'   )

    call tem_time_out(me     = me%interval, &
      &               conf   = conf,        &
      &               key    = 'interval'   )

    call aot_out_val( put_conf = conf,          &
      &               val      = me%check_iter, &
      &               vname    = 'check_iter'   )

    call aot_out_val( put_conf = conf,           &
      &               val      = me%delay_check, &
      &               vname    = 'delay_check'   )

    call aot_out_close_table(put_conf = conf)

  end subroutine tem_timeControl_out
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Dump timecontrol information to the specified outUnit.
  subroutine tem_timeControl_dump(me, outUnit)
    ! -------------------------------------------------------------------- !
    !> Time control to write on outUnit.
    type(tem_timeControl_type), intent(inout) :: me

    !> The file unit to write the time control to.
    integer, intent(in) :: outUnit
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    write(outUnit,*) '+-------------------------------------------+'
    write(outUnit,'(A,I0)') 'Iterations between trigger checks: ', me%check_iter
    write(outUnit,*) ' - Min:'
    call tem_time_dump(me%min, outUnit)
    write(outUnit,*) ''
    write(outUnit,*) ' - Max:'
    call tem_time_dump(me%max, outUnit)
    write(outUnit,*) ''
    write(outUnit,*) ' - Interval:'
    call tem_time_dump(me%interval, outUnit)
    write(outUnit,*) '+-------------------------------------------+'

  end subroutine tem_timeControl_dump
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Set the begin of the control interval in me to now.
  !!
  !! Setting only the simulation time, while putting all other counters to
  !! never avoids the introduction of new dependencies, that might result in
  !! the need for communication.
  subroutine tem_timeControl_start_at_sim(me, now)
    ! -------------------------------------------------------------------- !
    !> Time control that should be started at now.
    type(tem_timeControl_type), intent(inout) :: me

    !> Time that should be used as starting point for the time control.
    type(tem_time_type), intent(in) :: now
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    me%min = tem_time_never()
    if ( me%interval%clock < huge(me%interval%clock) ) then
      me%min%clock = now%clock
    end if
    me%min%iter = now%iter
    me%min%sim = now%sim
    me%trigger = me%min

    ! As me%min is set to now, it is also reached.
    me%min_reached = .true.

  end subroutine tem_timeControl_start_at_sim
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Returns if the timeControl has triggered since last update.
  !!
  !! This is true if now >= me%trigger and with in the bounds of min and max.
  !! Please note that, to allow arbitrary settings of min and interval, this
  !! routine might change the timeControl data given in me, by setting the
  !! trigger to now, when min is reached for the first time.
  !! This is required to allow independent time definitions for min, max and
  !! interval.
  function tem_timeControl_triggered(me, now) result(hasTriggered)
    ! -------------------------------------------------------------------- !
    !> Time control to check if it was triggered.
    type(tem_timeControl_type), intent(inout) :: me

    !> Current time that is to be used as comparison for the trigger check.
    type(tem_time_type), intent(in) :: now

    !> Result indicating if the time control has triggered.
    logical :: hasTriggered
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    hasTriggered = .false.

    ! Time control intervals *only* are checked every check_iter iteration,
    ! in between the intervals, the triggered status remains false.
    if (mod(now%iter, me%check_iter) == 0) then

      ! As long as the min was not reached yet, we need to do some extra checks,
      ! to ensure that the trigger can be set correctly when min is reached for
      ! first time.
      if (.not. me%min_reached) then

        me%min_reached = tem_time_ge_trigger(now, me%min)

        ! 'Now' has reached min, set the trigger accordingly for all those
        ! entries that are relevant in the interval configuration.
        if (me%min_reached) then
          me%trigger = tem_time_never()
          if (me%interval%sim   < huge(me%interval%sim)) then
            me%trigger%sim = now%sim
          end if
          if (me%interval%iter  < huge(me%interval%iter)) then
            me%trigger%iter = now%iter
          end if
          if (me%interval%clock < huge(me%interval%clock)) then
            me%trigger%clock = now%clock
          end if
        end if

      end if

      if (me%min_reached .and. (.not. tem_time_gt_trigger(now, me%max)) ) then
        hasTriggered = tem_time_ge_trigger(now, me%trigger)
      end if

    end if

  end function tem_timeControl_triggered
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine checks globally if the control has triggered.
  !!
  !! It takes care of communication as well.
  !! A reduction of the trigger status might be needed, depending on the time
  !! definitions in the trigger.
  !! This communication is done with the MPI communicator comm, and all
  !! processes calling this routine should be members of comm.
  !! The communication is only done, if necessary.
  !! If the trigger became active since the last check or update, the triggered
  !! argument will be set to true.
  !! If delay_check is true, the communication is done with a nonblocking
  !! allreduce, which performs the check essentially on the previous check
  !! interval, rather than the current one.
  !!
  !! If this should be done in combination with other status communications to
  !! avoid unnecessary synchronisation points, the separate routines
  !! tem_timeControl_triggered and tem_timeControl_update have to be used
  !! instead.
  function tem_timeControl_globalTriggered(me, now, comm) result (hasTriggered)
    ! -------------------------------------------------------------------- !
    !> Time control to check if it was triggered across all processes.
    type(tem_timeControl_type), intent(inout) :: me

    !> Current time to use for the comparison.
    type(tem_time_type), intent(in) :: now

    !> Communicator to use for the global reduction.
    integer, intent(in) :: comm

    !> Result indicating if the time control has triggered.
    logical :: hasTriggered
    ! -------------------------------------------------------------------- !
    logical :: local_triggered
    integer :: iError
    integer :: sync_status(MPI_STATUS_SIZE)
    ! -------------------------------------------------------------------- !

    local_triggered = tem_timeControl_triggered(me, now)
    hasTriggered = local_triggered

    if (me%needs_reduce .and. me%delay_check) then
      if (me%check_request /= MPI_REQUEST_NULL) then
        call MPI_WAIT(me%check_request, sync_status, iError)
        hasTriggered = me%globally_triggered
      end if
    end if

    if (me%needs_reduce .and. (mod(now%iter, me%check_iter) == 0)) then
      if (me%delay_check) then
        me%globally_triggered = local_triggered
        call MPI_IAllreduce(MPI_IN_PLACE, me%globally_triggered,             &
          &                 1, MPI_LOGICAL, MPI_LOR, comm, me%check_request, &
          &                 iError                                           )
        ! has_triggered is set by the check after waiting on the completion
        ! of the previous allreduce above.
      else
        call MPI_Allreduce(local_triggered, me%globally_triggered, 1, &
          &                MPI_LOGICAL, MPI_LOR, comm, iError         )
        hasTriggered = me%globally_triggered
      end if
    end if

  end function tem_timeControl_globalTriggered
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Returns for each of the time measurements, if the max specification was
  !! reached.
  !!
  !! An array of logicals of the length tem_time_n_ids is returned, indicating
  !! for each measurement, if the max time of the timeControl was reached.
  pure function tem_timeControl_reached_max(me, now) result(at_max)
    ! -------------------------------------------------------------------- !
    !> Time control to compare agains its max settings.
    type(tem_timeControl_type), intent(in) :: me

    !> Current time to compare the max settings to.
    type(tem_time_type), intent(in) :: now

    !> Resulting array indicating for each time definition, if its max setting
    !! was reached.
    logical :: at_max(tem_time_n_ids)
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    at_max = tem_time_ge_each(now, me%max)

  end function tem_timeControl_reached_max
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Update the given timeControl if it just triggered.
  !!
  !! The timeControl will be updated to trigger after the next complete
  !! interval, or at least after now.
  !! The update is only done if the timeControl actually triggered since the
  !! last update, or the optional argument hasTriggered is true.
  !!
  !! Usually, this routine should be called right after checking the
  !! status of the time control with tem_timeControl_triggered.
  !! However, due to the fact, that the hasTriggered might need to be reduced
  !! across all processes in between (if the clock measurement is involved
  !! in the trigger), this communication has to be done first, and the
  !! hasTriggered argument can be used to pass the result from the allreduce.
  !! To avoid unnecessary communication, this allreduce might be used for
  !! other flags as well, therefore it is not included in these routines.
  !!
  !! If such a separation is not desirable, use tem_timeControl_check, which
  !! probe the trigger and update it if needed, including communication if
  !! the clock time setting is used in the trigger definition.
  subroutine tem_timeControl_update(me, now, hasTriggered, localTriggered)
    ! -------------------------------------------------------------------- !
    !> Time control object to update.
    type(tem_timeControl_type), intent(inout) :: me

    !> Current time to use for the update.
    type(tem_time_type), intent(in) :: now

    !> Flag to indicate if the time control already has triggered.
    !!
    !! If this argument is not present, the check for the trigger status of
    !! the time control will be done internally.
    logical, intent(in), optional :: hasTriggered

    !> Flag to indicate if the local time control already has triggered.
    !!
    !! This will be used in place of the hasTriggered, if delayCheck is
    !! true to avoid subsequent multiple checks of the same interval.
    logical, intent(in), optional :: localTriggered
    ! -------------------------------------------------------------------- !
    logical :: triggered
    ! -------------------------------------------------------------------- !

    if (present(hasTriggered)) then
      triggered = hasTriggered
      if (present(localTriggered) .and. me%delay_check) then
        triggered = localTriggered
      end if
    else
      triggered = tem_timeControl_triggered(me, now)
    end if

    if (triggered) then
      me%trigger = max(me%trigger + me%interval, now)
    end if

  end subroutine tem_timeControl_update
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine checks if the control has triggered, and if so updates it.
  !!
  !! It takes care of communication as well.
  !! A reduction of the trigger status might be needed, depending on the time
  !! definitions in the trigger.
  !! This communication is done with the MPI communicator comm, and all
  !! processes calling this routine should be members of comm.
  !! The communication is only done, if necessary.
  !! If the trigger became active since the last check or update, the triggered
  !! argument will be set to true.
  !!
  !! If this should be done in combination with other status communications to
  !! avoid unnecessary synchronisation points, the separate routines
  !! tem_timeControl_triggered and tem_timeControl_update have to be used
  !! instead.
  subroutine tem_timeControl_check(me, now, comm, triggered)
    ! -------------------------------------------------------------------- !
    !> Time control settings to check.
    type(tem_timeControl_type), intent(inout) :: me

    !> Current time to check the control against.
    type(tem_time_type), intent(in) :: now

    !> Communicator to use in the reduction.
    integer, intent(in) :: comm

    !> Result of the check, indicating if the time control was triggered now.
    logical, intent(out) :: triggered
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    triggered = tem_timeControl_globalTriggered(me, now, comm)

    call tem_timeControl_update(me, now, triggered)

  end subroutine tem_timeControl_check
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine resets trigger to min and sets min_reached to false
  subroutine tem_timeControl_reset_trigger(me)
    ! -------------------------------------------------------------------- !
    !> Time control settings to check.
    type(tem_timeControl_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    me%trigger = me%min
    me%min_reached = .false.

  end subroutine tem_timeControl_reset_trigger
  ! ************************************************************************ !

end module tem_timeControl_module
