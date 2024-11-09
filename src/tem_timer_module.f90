! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2013,2015,2017,2021 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2014 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016-2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2018 Robin Weihe <robin.weihe@student.uni-siegen.de>
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
!> Timer infrastructure methods for the TreElM module
!!
!! This module provides a convenience timing functionality for coarse grain
!! measurements by encapsulating mpi_wtime and providing a set of timers to
!! track.
!!
module tem_timer_module

  ! include treelm modules
  use mpi

  use flu_binding, only: flu_next
  use aotus_module, only: aot_get_val, flu_state
  use aot_table_ops_module, only: aot_table_open, aot_table_close, &
    &                             aot_table_top, aot_table_first

  use env_module,            only: rk, rk_mpi, labelLen, pathLen
  use tem_aux_module,        only: tem_open, tem_abort
  use tem_grow_array_module, only: grw_realArray_type, &
    &                              grw_intArray_type, &
    &                              grw_logicalArray_type, &
    &                              append
  use tem_dyn_array_module, only: dyn_labelArray_type, &
    &                             append, positionofval
  use tem_logging_module, only: logUnit
  use tem_tools_module, only: upper_to_lower

  implicit none
  private

  public :: tem_timer_type
  public :: tem_labeledtimer_type
  public :: tem_timerconfig_type

  public :: tem_appendTimers
  public :: tem_addTimer
  public :: tem_startTimer
  public :: tem_stopTimer
  public :: tem_resetTimer
  public :: tem_getTimerName
  public :: tem_getTimerVal
  public :: tem_getMaxTimerVal
  public :: tem_getMinTimerVal
  public :: tem_getAveTimerVal
  public :: tem_getSumTimerVal
  public :: tem_getNTimers
  public :: tem_timer_dumplabeled
  public :: tem_timer_loadconfig
  public :: tem_timer_dump_glob
  public :: tem_timer_loadconfig_glob
  public :: tem_set_timerConfig
  public :: tem_get_timerConfig

  integer, parameter :: tem_timer_ignored = 0
  integer, parameter :: tem_timer_summary = 1
  integer, parameter :: tem_timer_details = 2

  ! Collection of (unnamed) timers
  type tem_timer_type
    !> Number of timers in this collection
    integer :: nTimers = 0

    !> Start timing values
    type(grw_realArray_type) :: tStart

    !> timer running?
    type(grw_logicalArray_type) :: running

    !> timing value
    type(grw_realArray_type) :: duration
  end type tem_timer_type


  !> Configuration of the output for the timers.
  type tem_timerconfig_type
    !> Name of the file to write the timings into.
    !!
    !! If this is an empty string (the default), no timings will be written.
    character(len=PathLen) :: filename = ''

    !> Label of the timer to apply the verbosity setting to.
    type(dyn_labelArray_type) :: label

    !> Defines to which detail the corresponding timer with this label
    !! should be printed:
    !!
    !! - time_timer_ignored: do not print the timer information at all
    !! - time_timer_summary: print min, max and sum of the timer over all
    !!                       processes
    !! - time_timer_details: gather the timings from all processes and print
    !!                       them all for this timer in a separate file
    !!
    !! By default a summary will be printed for each timer.
    type(grw_intArray_type) :: verbosity
  end type tem_timerconfig_type


  ! Timers with names, each label can only be present once in this datatype.
  type tem_labeledtimer_type
    !> Actual timer data
    type(tem_timer_type) :: timedat

    !> Output configuration
    type(tem_timerconfig_type) :: config

    !> Label to use for each timer (unique).
    type(dyn_labelArray_type) :: label
  end type tem_labeledtimer_type

  type(tem_labeledtimer_type), save :: timer


contains


  !> Append nVals new timers to the timer collection provided in 'me'.
  subroutine tem_appendTimers(me, nVals)
    ! -------------------------------------------------------------------- !
    !> Timer object to extend by nVals timers
    type(tem_timer_type), intent(inout) :: me

    !> Number of timers to append to me.
    integer, intent(in) :: nVals
    ! -------------------------------------------------------------------- !
    logical :: running(nVals)
    real(kind=rk) :: zeroes(nVals)
    ! -------------------------------------------------------------------- !

    running = .false.
    zeroes = 0.0_rk

    call append( me  = me%running, &
      &          val = running     )
    call append( me  = me%tStart, &
      &          val = zeroes     )
    call append( me  = me%duration, &
      &          val = zeroes       )

    me%nTimers = me%nTimers + nVals

  end subroutine tem_appendTimers
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Setup a new timer object, reset the values and give it a label
  !! for later identification
  subroutine tem_addTimer( me, timerName, timerHandle )
    ! -------------------------------------------------------------------- !
    !> timer object
    type(tem_labeledtimer_type), intent(inout), optional :: me
    !> Name for this timer
    character(len=*), intent(in) :: timerName
    !> A handle to reference this timer
    integer, intent(out) :: timerHandle
    ! -------------------------------------------------------------------- !
    logical :: is_newtimer
    ! -------------------------------------------------------------------- !

    if ( present(me) ) then
      call append( me       = me%label,        &
        &          val      = trim(timerName), &
        &          pos      = timerHandle,     &
        &          wasAdded = is_newtimer      )
      if (is_newtimer) then
        call tem_appendTimers( me    = me%timedat, &
          &                    nVals = 1           )
      end if
    else
      call append( me       = timer%label,     &
        &          val      = trim(timerName), &
        &          pos      = timerHandle,     &
        &          wasAdded = is_newtimer      )
      if (is_newtimer) then
        call tem_appendTimers( me    = timer%timedat, &
          &                    nVals = 1              )
      end if
    end if

  end subroutine tem_addTimer
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Start the timer for the given timer handle
  subroutine tem_startTimer( me, timerHandle )
    ! -------------------------------------------------------------------- !
    !> timer object
    type(tem_timer_type), intent(inout), optional :: me
    !> Handle of the timer to start
    integer, intent(in) :: timerHandle
    ! -------------------------------------------------------------------- !

    if ( present(me) ) then
      if ( .not. me%running%val(timerHandle) ) then
        me%running%val(timerHandle) = .true.
        me%tStart%val(timerHandle)  = mpi_wtime()
      end if
    else
      if ( .not.timer%timedat%running%val(timerHandle) ) then
        timer%timedat%running%val(timerHandle) = .true.
        timer%timedat%tStart%val(timerHandle)  = mpi_wtime()
      end if
    end if

  end subroutine tem_startTimer
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Reset the timer to zero for the given timer handle
  subroutine tem_resetTimer( me, timerHandle )
    ! -------------------------------------------------------------------- !
    !> timer object
    type(tem_timer_type), intent(inout), optional :: me
    !> Handle of the timer to start
    integer, intent(in) :: timerHandle
    ! -------------------------------------------------------------------- !

    if ( present(me) ) then
      me%tStart%val(timerHandle) = mpi_wtime()
      me%duration%val(timerHandle)  = 0._rk
    else
      timer%timedat%tStart%val(timerHandle) = mpi_wtime()
      timer%timedat%duration%val(timerHandle)  = 0._rk
    end if

  end subroutine tem_resetTimer
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Stop the timer for the given timer handle and update the timings
  subroutine tem_stopTimer( me, timerHandle )
    ! -------------------------------------------------------------------- !
    !> timer object
    type(tem_timer_type), intent(inout), optional :: me
    !> Handle of the timer to stop
    integer, intent(in) :: timerHandle
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: endTime
    ! -------------------------------------------------------------------- !

    endTime = mpi_wtime()

    if ( present(me) ) then
      if ( me%running%val(timerHandle) ) then
        me%running%val(timerHandle) = .false.
        me%duration%val(timerHandle) = me%duration%val(timerHandle) &
          &                            + endTime                    &
          &                            - me%tStart%val(timerHandle)
      end if
    else
      if ( timer%timedat%running%val(timerHandle) ) then
        timer%timedat%running%val(timerHandle) = .false.
        timer%timedat%duration%val(timerHandle) &
          &  = timer%timedat%duration%val(timerHandle) &
          &    + endTime                               &
          &    - timer%timedat%tStart%val(timerHandle)
      end if
    end if

  end subroutine tem_stopTimer
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Write out the timer name and its value
  pure function tem_getTimerName( me, timerHandle ) result( timerName )
    ! -------------------------------------------------------------------- !
    !> timer object
    type(tem_labeledtimer_type), intent(in), optional :: me
    !> timer handle
    integer, intent(in) :: timerHandle
    !> name of the timer
    character(len=labelLen) :: timerName
    ! -------------------------------------------------------------------- !

    if ( present(me) ) then
      timerName = trim( me%label%val(timerHandle) )
    else
      timerName = trim( timer%label%val(timerHandle) )
    end if

  end function tem_getTimerName
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Write out the timer name and its value
  function tem_getTimerVal( me, timerHandle ) result( retValue )
    ! -------------------------------------------------------------------- !
    !> timer object
    type(tem_timer_type), intent(inout), optional :: me
    !> timer handle
    integer, intent(in) :: timerHandle
    !> timer value
    real(kind=rk) :: retValue
    ! -------------------------------------------------------------------- !

    if ( present(me) ) then
      if ( me%running%val(timerHandle) ) then
        call tem_stoptimer(me, timerhandle)
      end if
      retValue = me%duration%val(timerHandle)
    else
      if ( timer%timedat%running%val(timerHandle) ) then
        call tem_stoptimer(timerhandle = timerhandle)
      end if
      retValue = timer%timedat%duration%val(timerHandle)
    end if

  end function tem_getTimerVal
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Get the maximum timer duration across all partitions.
  !!
  !!@Note This assumes the same timerhandle to be used across all partitions
  !!      for all timers.
  function tem_getMaxTimerVal( me, timerHandle, comm ) result( retValue )
    ! -------------------------------------------------------------------- !
    !> timer object
    type(tem_timer_type), intent(inout), optional :: me
    !> timer handle
    integer, intent(in) :: timerHandle
    !> communicator handle
    integer, intent(in) :: comm
    !> timer value
    real(kind=rk) :: retValue
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: send, recv
    integer :: iError
    ! -------------------------------------------------------------------- !

    send = tem_getTimerVal( me = me, timerHandle = timerHandle )

    call MPI_Allreduce( send, recv, 1, rk_mpi, MPI_MAX, comm, iError )

    retValue = recv

  end function tem_getMaxTimerVal
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Get the minimum timer duration across all partitions.
  !!
  !!@Note This assumes the same timerhandle to be used across all partitions
  !!      for all timers.
  function tem_getMinTimerVal( me, timerHandle, comm ) result( retValue )
    ! -------------------------------------------------------------------- !
    !> timer object
    type(tem_timer_type), intent(inout), optional :: me
    !> timer handle
    integer, intent(in) :: timerHandle
    !> communicator handle
    integer, intent(in) :: comm
    !> timer value
    real(kind=rk) :: retValue
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: send, recv
    integer :: iError
    ! -------------------------------------------------------------------- !

    send = tem_getTimerVal( me = me, timerHandle = timerHandle )

    call MPI_Allreduce( send, recv, 1, rk_mpi, MPI_MIN, comm, iError )

    retValue = recv

  end function tem_getMinTimerVal
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Get the sum of timer durations across all partitions.
  !!
  !!@Note This assumes the same timerhandle to be used across all partitions
  !!      for all timers.
  function tem_getSumTimerVal( me, timerHandle, comm ) result( retValue )
    ! -------------------------------------------------------------------- !
    !> timer object
    type(tem_timer_type), intent(inout), optional :: me
    !> timer handle
    integer, intent(in) :: timerHandle
    !> communicator handle
    integer, intent(in) :: comm
    !> timer value
    real(kind=rk) :: retValue
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: send, recv
    integer :: iError
    ! -------------------------------------------------------------------- !

    send = tem_getTimerVal( me = me, timerHandle = timerHandle )

    call MPI_Allreduce( send, recv, 1, rk_mpi, MPI_SUM, comm, iError )

    retValue = recv

  end function tem_getSumTimerVal
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Get the average of timer durations across all partitions.
  !!
  !!@Note This assumes the same timerhandle to be used across all partitions
  !!      for all timers.
  function tem_getAveTimerVal( me, timerHandle, comm, nProcs ) result(retValue)
    ! -------------------------------------------------------------------- !
    !> timer object
    type(tem_timer_type), intent(inout), optional :: me
    !> timer handle
    integer, intent(in) :: timerHandle
    !> communicator handle
    integer, intent(in) :: comm
    !> Number of processes in the communicator.
    integer, intent(in) :: nProcs
    !> timer value
    real(kind=rk) :: retValue
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: send, recv
    integer :: iError
    ! -------------------------------------------------------------------- !

    send = tem_getTimerVal( me = me, timerHandle = timerHandle )

    call MPI_Allreduce( send, recv, 1, rk_mpi, MPI_SUM, comm, iError )

    retValue = recv / nProcs

  end function tem_getAveTimerVal
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Get the number of timers.
  !!
  pure function tem_getNTimers( me ) result ( n )
    ! -------------------------------------------------------------------- !
    !> timer object
    type(tem_timer_type), intent(in), optional :: me
    !> number of timer handles
    integer :: n
    ! -------------------------------------------------------------------- !

    if ( present(me) ) then
      n = me%nTimers
    else
      n = timer%timedat%nTimers
    end if

  end function tem_getNTimers
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Load the configuration for the timer output.
  !!
  !! The user may specify how which timers are to be written.
  !! We use a timer table to describe the file to write to and the level of
  !! detail for each timer.
  !!
  !!    :::lua
  !!    timer = {
  !!      file = 'timeinfo', -- this is the default
  !!      details = {
  !!        {'overall', 'details'}, -- will be written to
  !!                                -- timeinfo_overall.details
  !!        {'init', 'summary'},    -- summary is the default
  !!        {'output', 'ignored'}   -- will not be printed
  !!      }
  !!    }
  !!
  !! You might list an arbitrary number of timers with their level of detail,
  !! those which are not stated will be printed in summary form.
  !! If no timer table is given, or the file name is an empty string, no
  !! timer information will be written.
  !!
  !!@note Names of the timers will only be checked during dumping, the names
  !!      and verbosity settings are case insensitive.
  !!      If the same name is provided multiple times, the first occurence will
  !!      take precedence.
  subroutine tem_timer_loadconfig(timer_config, conf, parent)
    ! -------------------------------------------------------------------- !
    !> Timer configuration to load from the Lua configuration script.
    type(tem_timerconfig_type), intent(out) :: timer_config
    !> Handle to the Lua configuration script.
    type(flu_state) :: conf
    !> Handle of the table containing the requested table.
    integer, intent(in), optional :: parent
    ! -------------------------------------------------------------------- !
    integer :: timertab
    integer :: settingtab
    integer :: detailtab
    character(len=labelLen) :: label
    character(len=labelLen) :: verbosity
    character(len=labelLen) :: tmp
    integer :: verbconf
    integer :: pos
    logical :: wasAdded
    integer :: iError
    ! -------------------------------------------------------------------- !

    call aot_table_open( L       = conf,     &
      &                  parent  = parent,   &
      &                  thandle = timertab, &
      &                  key     = 'timer'   )

    if (timertab /= 0) then
      call aot_get_val( L       = conf,                 &
        &               thandle = timertab,             &
        &               key     = 'file',               &
        &               default = 'timerinfo',          &
        &               errcode = iError,               &
        &               val     = timer_config%filename )
      write(logUnit(3),*) 'Timer info is written to file: ',&
      &  trim(timer_config%filename)

      call aot_table_open( L       = conf,      &
        &                  parent  = timertab,  &
        &                  thandle = detailtab, &
        &                  key     = 'details'  )

      if (aot_table_first(L = conf, thandle = detailtab)) then
        write(logunit(3),*) 'Reading timer configuration'
        do
          settingtab = aot_table_top(conf)
          call aot_get_val( L       = conf,       &
            &               thandle = settingtab, &
            &               pos     = 1,          &
            &               errcode = iError,     &
            &               val     = label       )
          call aot_get_val( L       = conf,       &
            &               thandle = settingtab, &
            &               pos     = 2,          &
            &               errcode = iError,     &
            &               val     = verbosity   )
          call aot_table_close(L = conf, thandle = settingtab)
          call append( me       = timer_config%label,    &
            &          val      = upper_to_lower(label), &
            &          pos      = pos,                   &
            &          wasAdded = wasAdded               )
          if (wasAdded) then
            tmp = upper_to_lower(verbosity)
            select case(trim(tmp))
            case ('ignore')
              verbconf = tem_timer_ignored
            case ('summary')
              verbconf = tem_timer_summary
            case ('details')
              verbconf = tem_timer_details
            case default
              write(logunit(1),*) 'WARNING: Unknown verbosity setting ' &
                &                 // trim(tmp) // ' for timer '         &
                &                 // trim(label) // ' !'
              write(logunit(1),*) 'Please set it to one of the following:'
              write(logunit(1),*) '  * ignore'
              write(logunit(1),*) '  * summary'
              write(logunit(1),*) '  * details'
              write(logunit(1),*)
              write(logunit(1),*) 'Using the default setting summary ' &
                &                 // 'for timer ' // trim(label)
              verbconf = tem_timer_summary
            end select
            call append( me       = timer_config%verbosity, &
              &          val      = verbconf                )
          end if
          ! Get the next entry from the timer table or leave the loop, if there
          ! is none anymore.
          if ( .not. flu_next(L = conf, index = detailtab) ) EXIT
        end do
      end if

    end if

    call aot_table_close(L = conf, thandle = timertab)

  end subroutine tem_timer_loadconfig
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Load the configuration for the global (module) timer.
  subroutine tem_timer_loadconfig_glob(conf, parent)
    ! -------------------------------------------------------------------- !
    !> Timer configuration to load from the Lua configuration script.
    !> Handle to the Lua configuration script.
    type(flu_state) :: conf
    !> Handle of the table containing the requested table.
    integer, intent(in), optional :: parent
    ! -------------------------------------------------------------------- !

    call tem_timer_loadconfig( timer_config = timer%config, &
      &                        conf         = conf,         &
      &                        parent       = parent        )

   end subroutine tem_timer_loadconfig_glob
  ! ************************************************************************ !



  ! ************************************************************************ !
  subroutine tem_timer_dumplabeled(me, comm, myrank, nProcs)
    ! -------------------------------------------------------------------- !
    !> timer object
    type(tem_labeledtimer_type), intent(inout) :: me
    !> communicator handle
    integer, intent(in) :: comm
    !> MPI rank of the calling process.
    integer, intent(in) :: myrank
    !> Number of processes in the communicator.
    integer, intent(in) :: nProcs
    ! -------------------------------------------------------------------- !
    character(len=labelLen) :: timerlabel
    character(len=pathLen) :: detail_file
    integer :: iProc
    integer :: funit
    integer :: dunit
    integer :: timerpos
    integer :: itimer
    integer :: verbosity
    integer :: iError
    logical :: unmatched(me%config%label%nVals)
    real(kind=rk) :: mintime, maxtime, sumtime
    real(kind=rk), allocatable :: proctimer(:)
    ! -------------------------------------------------------------------- !

    unmatched = .true.

    if (me%config%filename /= '') then

      if (myrank == 0) then
        call tem_open( file    = trim(me%config%filename), &
          &            newunit = funit,                    &
          &            action  = 'write',                  &
          &            status  = 'replace',                &
          &            form    = 'formatted'               )
        write(funit,*) 'Timings for a run on ', nProcs, 'processes'
        write(funit,'(a24,3(1x,a16))') 'timer', 'min', 'max', 'sum'
        allocate(proctimer(nProcs))
      else
        allocate(proctimer(0))
      end if
      do itimer=1,me%label%nVals
        verbosity = tem_timer_summary
        timerlabel = upper_to_lower(me%label%val(itimer))
        timerpos = PositionOfVal( me  = me%config%label, &
          &                       val = trim(timerlabel) )
        if (timerpos > 0) then
          verbosity = me%config%verbosity%val(timerpos)
          unmatched(timerpos) = .false.
        end if

        select case(verbosity)
        case (tem_timer_summary)
          mintime = tem_getMinTimerVal(me%timedat, itimer, comm)
          maxtime = tem_getMaxTimerVal(me%timedat, itimer, comm)
          sumtime = tem_getSumTimerVal(me%timedat, itimer, comm)
          if (myrank == 0) then
            write(funit,'(a24,3(1x,en16.6))') trim(me%label%val(itimer)), &
              &                               mintime, maxtime, sumtime
          end if

        case (tem_timer_details)
          call MPI_Gather( me%timedat%duration%val(iTimer), 1, rk_mpi, &
            &              proctimer, 1, rk_mpi, 0, comm, iError       )
          if (myrank == 0) then
            detail_file = trim(me%config%filename) // '_' &
              &           // trim(me%label%val(itimer)) // '.details'
            write(funit,'(a1,a22,a1,1x,a)') 'D', trim(me%label%val(itimer)), &
              &                             ':', trim(detail_file)
            call tem_open( file    = trim(detail_file), &
              &            newunit = dunit,             &
              &            action  = 'write',           &
              &            status  = 'replace',         &
              &            form    = 'formatted'        )
            write(dunit, '(a)') 'Detailed timings for ' &
              &                 // trim(me%label%val(itimer))
            write(dunit, '(a8,1x,en16.6)') 'min:', minval(proctimer)
            write(dunit, '(a8,1x,en16.6)') 'max:', maxval(proctimer)
            write(dunit, '(a8,1x,en16.6)') 'sum:', sum(proctimer)
            write(dunit, *) ''
            do iProc=1,nProcs
              write(dunit, '(i7,a1,en16.6)') iProc, ':', proctimer(iProc)
            end do
            close(dunit)
          end if

        end select

      end do

      deallocate(proctimer)

      ! Write any user-defined, but unknown timers...
      if (myrank == 0) then
        do itimer=1,me%config%label%nVals
          if ( unmatched(iTimer)                       &
            &  .and. ( me%config%verbosity%val(iTimer) &
            &          /= tem_timer_ignored )          ) then
            write(funit,'(a1,a23,1x,a)') 'U',                               &
              &                          trim(me%config%label%val(itimer)), &
              &                          '    .:!UNKNOWN!:.'
          end if
        end do

        close(funit)
      end if

    end if

  end subroutine tem_timer_dumplabeled
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_timer_dump_glob(comm, myrank, nProcs)
    ! -------------------------------------------------------------------- !
    !> communicator handle
    integer, intent(in) :: comm
    !> MPI rank of the calling process.
    integer, intent(in) :: myrank
    !> Number of processes in the communicator.
    integer, intent(in) :: nProcs
    ! -------------------------------------------------------------------- !

    call tem_timer_dumplabeled(timer, comm, myrank, nProcs)

  end subroutine tem_timer_dump_glob
  ! ************************************************************************ !


  ! ***************************************************************************!
  !> This routine sets timer config required for apesmate
  subroutine tem_set_timerConfig(timerConfig)
    !---------------------------------------------------------------------------
    type(tem_timerconfig_type), intent(in)  :: timerConfig
    !---------------------------------------------------------------------------
    timer%config = timerConfig
  end subroutine tem_set_timerConfig
  ! ***************************************************************************!


  ! ***************************************************************************!
  !> This routine gets timer config required for apesmate
  function tem_get_timerConfig()result(timerConfig)
    !---------------------------------------------------------------------------
    type(tem_timerconfig_type) :: timerConfig
    !---------------------------------------------------------------------------
    timerConfig = timer%config
  end function tem_get_timerConfig
  ! ***************************************************************************!

end module tem_timer_module
! **************************************************************************** !
