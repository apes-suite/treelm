! Copyright (c) 2013-2017, 2019, 2021 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2015, 2018 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
! ****************************************************************************** !
!> author: Harald Klimach
!! Initialize general infrastructure settings of treelm.
!!
!! This module provides a central point to initialize the
!! various general settings in treelm.
module tem_general_module

  ! include aotus modules
  use aotus_module, only: flu_State, aot_get_val

  ! include treelm modules
  use env_module,             only: rk, tem_load_env_params, io_buffer_size,   &
    &                               printRuntimeInfo, init_env, fin_env,       &
    &                               print_self_status, null_device, pathLen,   &
    &                               stdOutUnit, labelLen
  use tem_aux_module,         only: tem_abort
  use tem_abortCriteria_module, only: tem_solverAborts_type
  use tem_comm_module,        only: tem_commpattern_type, tem_load_commpattern
  use tem_logging_module,     only: logUnit
  use tem_solveHead_module,   only: tem_solveHead_type, tem_init_solveHead
  use tem_restart_module,     only: tem_restart_type
  use tem_timer_module,       only: tem_addTimer, tem_startTimer,              &
    &                               tem_stopTimer, tem_getTimerVal,            &
    &                               tem_getTimerName, tem_timer_dump_glob,     &
    &                               tem_timer_loadconfig_glob
  use tem_balance_module,     only: tem_balance_type, tem_balance_load
  use tem_status_module,      only: tem_status_type, tem_status_run_terminate
  use tem_comm_env_module,    only: tem_comm_env_init, tem_comm_env_fin,       &
    &                               tem_comm_env_type
  use tem_tools_module,       only: tem_horizontalSpacer
  use tem_simControl_module,  only: tem_simControl_type, tem_simControl_load,  &
    &                               tem_simControl_dump, tem_simControl_start


  use tem_precice_module, only: precice_available, tem_precice_load
  use tem_sparse_comm_module, only: tem_sparse_comm_load

  implicit none

  private

  public :: tem_start
  public :: tem_finalize
  public :: tem_load_general
  public :: tem_general_type

  !> Global parameter type contains all general information
  !! needed for all solvers
  type tem_general_type

    !> General description of the deployed solver.
    type(tem_solveHead_type) :: solver

    !> contains current simulation time, timeControl, abortCriteria and
    !! simulation status
    type(tem_simControl_type) :: simControl

    !> MPI communication enviroment including MPI communicator.
    !!
    !!@todo HK: not sure, if this should be here!
    type(tem_comm_env_type) :: proc

    !> MPI communication pattern type.
    type(tem_commPattern_type) :: commPattern

    !> Global restart type
    !!
    !!@todo HK: not sure, if this should be here!
    type(tem_restart_type) :: restart

    !> Load balancing information.
    !!
    !!@todo HK: not sure, if this should be here!
    type(tem_balance_type) :: balance

    !> Filename for solver timing output
    character(len=pathLen) :: timingFile

  end type tem_general_type


contains


  ! ************************************************************************** !
  !> Load general treelm settings from the Lua script in conf.
  !!
  subroutine tem_load_general( me, conf, timingFile, solverAborts )
    ! ----------------------------------------------------------------------
    !> global general parameter
    type( tem_general_type ), intent(inout) :: me
    !> Handle to the Lua script containing the configuration.
    type(flu_state) :: conf
    !> Default timing filename provided by the caller, overwritten by config
    !! file.
    character(len=*), optional, intent(in) :: timingFile

    !> Solver specific abort criteria to load.
    !!
    !! See [[tem_abortCriteria_module]] for details on this additional
    !! abortCriteria parameters, that the solver may define to be loaded
    !! from the configuration.
    class(tem_solverAborts_type), intent(inout), optional :: solverAborts
    ! ----------------------------------------------------------------------
    integer :: iError
    character(len=pathLen) :: def_timingFile
    ! ----------------------------------------------------------------------

    if ( me%proc%isRoot ) then
      call tem_horizontalSpacer(fUnit = logUnit(1))
      write(logUnit(1),"(A)") 'Loading general parameters:'
    end if

    ! load global enviromental parameters from config file
    call tem_load_env_params(conf)

    if ( me%proc%isRoot ) then
      write(logUnit(1),"(A)") 'Using '//trim(null_device)//' as null device.'

      if (printRuntimeInfo) then
        write(logUnit(1),"(A)") 'Will print run time info in the end'
        write(logUnit(1),"(A)") '(/proc/self/status).'
      else
        write(logUnit(1),"(A)") 'Will NOT print run time info.'
      end if

      write(logUnit(1),"(A,I0)") 'Size of the IO Buffer (MB): ', &
        &                          (io_buffer_size*8/1024/1024)
    end if

    ! load simulation name
    call aot_get_val( L       = conf,              &
      &               key     = 'simulation_name', &
      &               val     = me%solver%simName, &
      &               ErrCode = iError,            &
      &               default = 'simulation'       )

    if ( me%proc%isRoot ) then
      write(logUnit(1),"(A)") 'Simulation Name: ' // trim( me%solver%simName )
    end if

    ! load simulation time control
    call tem_simControl_load( me           = me%simControl, &
      &                       conf         = conf,          &
      &                       solverAborts = solverAborts   )

    if ( me%proc%isRoot ) then
      call tem_simControl_dump(me = me%simControl, outUnit = logUnit(1))
    end if

    ! Get the setting, whether to use sparse communication patterns or not.
    call tem_sparse_comm_load(conf)

    ! load communication pattern and initialize commuication
    ! infrastructure
    call tem_load_commpattern(conf = conf, me = me%commpattern)

    ! load timing filename
    if (present(timingFile)) then
      def_timingFile = timingFile
    else
      def_timingFile = 'timing.res'
    end if

    call aot_get_val( L       = conf,           &
      &               key     = 'timing_file',  &
      &               val     = me%timingFile,  &
      &               default = def_timingFile, &
      &               ErrCode = iError          )
    write(logUnit(1),"(A)") 'Write timings to: ' // trim( me%timingFile )

    ! Get load balancing config
    call tem_balance_load( conf = conf,      &
      &                    me   = me%balance )

    ! Load configuration of timer output from timer table.
    call tem_timer_loadconfig_glob( conf )

    ! load precice if available
    if (precice_available) then
      write(logUnit(1),"(A)") 'Loading precice data'
      call tem_precice_load(conf = conf)
    end if
    call tem_horizontalSpacer(fUnit = logUnit(1))

  end subroutine tem_load_general
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Initialize the environment. Should be the very first call in the program.
  !!
  subroutine tem_start(codeName, general, comm, simControl)
    ! ----------------------------------------------------------------------
    !> name of code
    character(len=*), intent(in) :: codeName
    !> encapsulates global parameters which are common for all solvers
    type(tem_general_type), intent(out) :: general
    !> mpi communicator if it is predefined as in apesmate
    integer, intent(in), optional :: comm
    !> simulation control to initialize
    type(tem_simControl_type), intent(out), optional :: simControl
    ! ----------------------------------------------------------------------
    integer :: nProcs, nThreads
    ! ----------------------------------------------------------------------

    ! Initialize all logunits to point to the stdout unit.
    logunit = stdoutunit

    ! if comm is present initialize environment already called
    ! so should not be called again
    if(.not. present(comm)) call init_env()

    ! initialize mpi environment
    call tem_comm_env_init(general%proc, comm)
    nProcs   = general%proc%comm_size
    nThreads = general%proc%nThreads

    ! initialize solverHead
    call tem_init_solveHead( me      = general%solver, &
      &                      solName = codeName        )

    if (present(simControl)) call tem_simControl_start(simControl)

    if ( general%proc%isRoot ) then
      write(logUnit(1),*) "Starting up " // trim(codeName) &
        &                 // " with nprocs: ", nProcs
      !$ write(logUnit(1),*)"               and nThreads pp: ", nThreads
    end if

    call tem_addTimer( timerHandle = general%solver%timerHandle, &
      &                timerName   = trim(codeName)              )
    call tem_startTimer(timerHandle = general%solver%timerHandle )

  end subroutine tem_start
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Finalize the environment, should be the very last call in the program.
  subroutine tem_finalize(general)
    ! ----------------------------------------------------------------------
    !> encapsulates global parameters which are common for all solvers
    type(tem_general_type), intent(in) :: general
    ! ----------------------------------------------------------------------
    character(len=labelLen) :: timerName
    real(kind=rk) :: timerValue
    ! ----------------------------------------------------------------------

    if ( tem_status_run_terminate(general%simControl%status) ) then
      call tem_abort()
    end if

    call tem_timer_dump_glob( comm   = general%proc%comm,     &
      &                       myrank = general%proc%rank,     &
      &                       nProcs = general%proc%comm_size )

    if ( general%proc%isRoot ) then
      if ( printRuntimeInfo ) call print_self_status()

      timerName  = tem_getTimerName(timerHandle = general%solver%timerHandle )
      timerValue = tem_getTimerVal( timerHandle = general%solver%timerHandle )

      write(logUnit(1),*)
      write(logUnit(1),"(A,F10.2,A)") " Done with " // trim(timerName) &
        &                             // " in ", timerValue, ' s'
    end if

    ! finialize mpi
    call fin_env()

  end subroutine tem_finalize
  ! ************************************************************************ !

end module tem_general_module
! **************************************************************************** !
