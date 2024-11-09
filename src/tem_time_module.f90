! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2014, 2017, 2019-2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2013-2014 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014, 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
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
!> The module contains the time definition of treelm.
!!
!! It is useful to describe the progress in time of an explicit solver with
!! respect to different measurements.
!! Currently accounted are:
!!
!! - iterations (`iter`)
!! - simulation time (`sim`)
!! - running time (`clock`)
!!
!! Thus, the time definition takes the following general form:
!!
!!```lua
!!  {sim = 1.0, iter = 12, clock = 60.0}
!!```
!!
!! Any of these definitions may be left out. If the time definition is not
!! a table the given value is interpreted as simultion time, equivalent to
!!
!!```lua
!! {sim = 1.0}
!!```
!!
!! For the handling of time evolution in an explicit solver, a time definition
!! is needed to keep track of the current time.
!! This should be initialized with [[tem_time_reset]], and then advanced in each
!! iteration by tem_time_advance.
!! Time definitions in relation to this main time evolution object can then
!! be defined in the Lua configuration by [[tem_time_load]].
!! For an explicit time stepping solver this would at least require one
!! definition of time, when to stop the simulation.
!! A comparison of such points in time to the current time can be done by
!! [[tem_time_ge_trigger]].
!!
module tem_time_module
  use mpi

  use env_module, only: rk, labelLen

  use aotus_module, only: flu_State, &
    &                     aot_get_val
  use aot_out_module, only: aot_out_type,       &
    &                       aot_out_val,        &
    &                       aot_out_open_table, &
    &                       aot_out_close_table
  use aot_table_module, only: aot_table_open, &
    &                         aot_table_close

  implicit none

  private

  public :: tem_time_type

  public :: tem_time_load
  public :: tem_time_dump
  public :: tem_time_out
  public :: tem_time_ge_each
  public :: tem_time_ge_trigger
  public :: tem_time_gt_trigger
  public :: tem_time_reset
  public :: tem_time_set_clock
  public :: tem_time_advance
  public :: tem_time_never
  public :: tem_time_isDefined
  public :: tem_time_last_interval
  public :: tem_time_default_zero

  public :: operator(+)
  public :: operator(-)
  public :: max

  public :: tem_time_needs_reduce
  public :: tem_time_sim_stamp
  public :: tem_time_iter_stamp

  !> Number of available different time definitions.
  integer, parameter, public :: tem_time_n_ids = 3


  !> Index for the simulation time.
  integer, parameter, public :: tem_time_sim_id   = 1
  !> Index for the number of iterations.
  integer, parameter, public :: tem_time_iter_id  = 2
  !> Index for the clock time.
  integer, parameter, public :: tem_time_clock_id = 3


  !> Description of time
  !!
  !! This type provides the description of a point in
  !! time in different terms.
  !! Currently supported are:
  !! - number of iterations
  !! - simulation time
  !! - running time (clock)
  !!
  !! It can be used to describe certain points in time
  !! as well as periods of time.
  type tem_time_type
    !> Time in iterations.
    integer :: iter

    !> Time in terms of simulated time.
    real(kind=rk) :: sim

    !> Time passed in seconds of running time.
    real(kind=rk) :: clock

    !> Wtime of the last reset.
    real(kind=rk) :: clock_start
  end type tem_time_type


  interface operator(+)
    module procedure tem_time_add_time
  end interface

  interface operator(-)
    module procedure tem_time_subtract_time
  end interface

  interface max
    module procedure tem_time_max
  end interface


contains


  ! ************************************************************************ !
  !> Reading a time description from a Lua script given by conf.
  !!
  !! The time description has to be provided in the table given by thandle.
  !! The time can be given in terms of:
  !! - iter: number of iterations
  !! - sim: simulated time
  !! - clock: running time used to compute the simulation
  !! Thus, the configuration looks like:
  !! \verbatim
  !! time = { iter = 123, sim = 1.23, clock = 12.3 }
  !! \end verbatim
  !! Omitted time defintions are set to the maximal representable number,
  !! indicating something like never.
  !! If the time is not provided as a table, it is interpreted as a setting of
  !! the sim time, the other entries are set to never.
  !!
  !! The optional argument clock_start is useful to overwrite the current time
  !! by the settings in a loaded restart header.
  !! With clock_start the old reference wtime from tem_time_reset can be
  !! maintained, while iterations and simulation time can be inherited from
  !! the configuration.
  subroutine tem_time_load(me, conf, key, parent, clock_start)
    ! -------------------------------------------------------------------- !
    !> Time to be read from the Lua script
    type(tem_time_type), intent(out) :: me

    !> Handle to the Lua script.
    type(flu_state), intent(inout) :: conf

    !> Name of the table containing the time definition. Default: 'time'.
    character(len=*), intent(in), optional :: key

    !> Handle to the parent table.
    integer, intent(in), optional :: parent

    !> Reference MPI_Wtime mark to use for computation of wallclock
    !!
    !! This is basically only useful to update the time object for a read in
    !! restart without distorting the clock measurement.
    real(kind=rk), intent(in), optional :: clock_start
    ! -------------------------------------------------------------------- !
    integer :: iErr
    character(len=labelLen) :: loc_key
    integer :: thandle
    ! -------------------------------------------------------------------- !

    if ( present(key) ) then
      loc_key = trim(key)
    else
      loc_key = 'time'
    end if

    call aot_table_open(L       = conf,    &
      &                 parent  = parent,  &
      &                 thandle = thandle, &
      &                 key     = loc_key  )

    if (thandle /= 0) then
      ! The time is given as a table load its components accordingly.
      call aot_get_val(L       = conf,         &
        &              thandle = thandle,      &
        &              val     = me%sim,       &
        &              key     = 'sim',        &
        &              default = huge(me%sim), &
        &              ErrCode = iErr          )

      call aot_get_val(L       = conf,          &
        &              thandle = thandle,       &
        &              val     = me%iter,       &
        &              key     = 'iter',        &
        &              default = huge(me%iter), &
        &              ErrCode = iErr           )

      call aot_get_val(L       = conf,           &
        &              thandle = thandle,        &
        &              val     = me%clock,       &
        &              key     = 'clock',        &
        &              default = huge(me%clock), &
        &              ErrCode = iErr            )
    else
      ! The time is not given as a table, try to interpret it as a setting for
      ! the simtime.
      call aot_get_val(L       = conf,         &
        &              thandle = parent,       &
        &              key     = loc_key,      &
        &              val     = me%sim,       &
        &              default = huge(me%sim), &
        &              ErrCode = iErr          )
      me%iter  = huge(me%iter)
      me%clock = huge(me%clock)
    end if

    call aot_table_close(conf, thandle)

    !>\todo HK: maybe provide the possibility to add up run times from loaded
    !!          data instead of resetting it.
    if (present(clock_start)) then
      me%clock_start = clock_start
    else
      me%clock_start = MPI_Wtime()
    end if

  end subroutine tem_time_load
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Write the time into a Lua table
  !!
  !! Entries, which are set to the maximal represantable values (meaning never),
  !! are ignored and not written.
  subroutine tem_time_out(me, conf, key)
    ! -------------------------------------------------------------------- !
    !> The time descritpion to write out as Lua table.
    type(tem_time_type), intent(in) :: me

    !> Output handle for the script to write to.
    type(aot_out_type), intent(inout) :: conf

    !> Name for the tabel to write this time description to.
    character(len=*), intent(in), optional :: key
    ! -------------------------------------------------------------------- !
    character(len=labelLen) :: loc_key
    ! -------------------------------------------------------------------- !
    if (present(key)) then
      loc_key = trim(key)
    else
      loc_key = 'time'
    end if

    call aot_out_open_table( put_conf = conf, tname = loc_key )
    if (me%sim < huge(me%sim)) then
      call aot_out_val( put_conf = conf,   &
        &               val      = me%sim, &
        &               vname    = 'sim'   )
    end if
    if (me%iter < huge(me%iter)) then
      call aot_out_val( put_conf = conf,    &
        &               val      = me%iter, &
        &               vname    = 'iter'   )
    end if
    if (me%clock < huge(me%clock)) then
      call aot_out_val( put_conf = conf,     &
        &               val      = me%clock, &
        &               vname    = 'clock'   )
    end if
    call aot_out_close_table( put_conf = conf )

  end subroutine tem_time_out
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Dump the given time to outUnit.
  !!
  !! Values which are not set, are omitted.
  subroutine tem_time_dump(me, outUnit)
    ! -------------------------------------------------------------------- !
    !> The time that should be written to outunit.
    type(tem_time_type), intent(inout) :: me

    !> The unit to write to.
    integer, intent(in) :: outUnit
    ! -------------------------------------------------------------------- !

    if (me%iter  < huge(me%iter))  write(outUnit, *) ' iterations: ', me%iter
    if (me%sim   < huge(me%sim))   write(outUnit, *) ' simTime   : ', me%sim
    if (me%clock < huge(me%clock)) write(outUnit, *) ' wallClock : ', me%clock

  end subroutine tem_time_dump
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Function to find a default for time, if it is not defined already.
  !!
  !! Depending on dependency set the value in time to a zero time,
  !! if it is set to never.
  !! This is useful to find a suitable minimum, where an interval is
  !! defined in the timeControl, but no minimum.
  pure function tem_time_default_zero(time, dependency) result(zeroed)
    ! -------------------------------------------------------------------- !
    !> The time to find a default 0 for, if not define, but dependency is.
    type(tem_time_type), intent(in) :: time

    !> A time definition where we set a 0 default in time for, if the
    !! corresponding component is defined, but the time component is
    !! not.
    type(tem_time_type), intent(in) :: dependency

    !> Resulting time, with zeroed components where dependency is defined,
    !! but time not.
    type(tem_time_type) :: zeroed
    ! -------------------------------------------------------------------- !

    zeroed = time

    if (.not. time%sim < huge(time%sim)) then
      ! Time set to never.
      if (dependency%sim < huge(dependency%sim)) then
        ! But dependency not!
        zeroed%sim = 0.0_rk
      end if
    end if
    if (.not. time%iter < huge(time%iter)) then
      ! Time set to never.
      if (dependency%iter < huge(dependency%iter)) then
        ! But dependency not!
        zeroed%iter = 0
      end if
    end if
    if (.not. time%clock < huge(time%clock)) then
      ! Time set to never.
      if (dependency%clock < huge(dependency%clock)) then
        ! But dependency not!
        zeroed%clock = 0.0_rk
      end if
    end if

  end function tem_time_default_zero
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Return the (>=) comparison of each time measurement between the two
  !! arguments.
  !!
  !! This returns an array of logicals with the length tem_time_n_ids, with
  !! each one indicating the comparison result for the individual measures.
  pure function tem_time_ge_each(time, trigger) result(ge)
    ! -------------------------------------------------------------------- !
    !> A fully defined time definition (all components should have meaningful
    !! settings).
    type(tem_time_type), intent(in) :: time

    !> A comparison time definition, where some entries might be set to never.
    !! If any of the time components is larger, the result will be true.
    type(tem_time_type), intent(in) :: trigger

    !> Result of the comparison for each of the time specifications.
    logical :: ge(tem_time_n_ids)
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    ge(tem_time_sim_id)   = (time%sim   >= trigger%sim)
    ge(tem_time_iter_id)  = (time%iter  >= trigger%iter)
    ge(tem_time_clock_id) = (time%clock >= trigger%clock)

  end function tem_time_ge_each
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Compare (>=) a complete time defintion against a trigger.
  !!
  !! This returns true, if any of the time definitions in time is greater or
  !! equal to the corresponding time given in trigger.
  !! The time argument should be completely defined, while the trigger might
  !! have some of its definitions set to never.
  !! Due to this definition this comparison operator is not useful to define
  !! a unique ordering of time definitions!
  elemental function tem_time_ge_trigger(time, trigger) result(ge)
    ! -------------------------------------------------------------------- !
    !> A fully defined time definition (all components should have meaningful
    !! settings).
    type(tem_time_type), intent(in) :: time

    !> A comparison time definition, where some entries might be set to never.
    !! If any of the time components is larger, the result will be true.
    type(tem_time_type), intent(in) :: trigger

    !> Result of the comparison, true if any of the time specifications in time
    !! is larger or equal to the ones in trigger.
    logical :: ge
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    ge = (      (time%sim   >= trigger%sim)   &
      &    .or. (time%iter  >= trigger%iter)  &
      &    .or. (time%clock >= trigger%clock) )

  end function tem_time_ge_trigger
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Compare (>) a complete time defintion against a trigger.
  !!
  !! This returns true, if any of the time definitions in time is greater
  !! than the corresponding time given in trigger.
  !! The time argument should be completely defined, while the trigger might
  !! have some of its definitions set to never.
  !! Due to this definition this comparison operator is not useful to define
  !! a unique ordering of time definitions!
  elemental function tem_time_gt_trigger(time, trigger) result(gt)
    ! -------------------------------------------------------------------- !
    !> A fully defined time definition (all components should have meaningful
    !! settings).
    type(tem_time_type), intent(in) :: time

    !> A comparison time definition, where some entries might be set to never.
    !! If any of the time components is larger, the result will be true.
    type(tem_time_type), intent(in) :: trigger

    !> Result of the comparison, true if any of the time specifications in time
    !! is larger or equal to the ones in trigger.
    logical :: gt
    ! -------------------------------------------------------------------- !

    gt = (      (time%sim   > trigger%sim)   &
      &    .or. (time%iter  > trigger%iter)  &
      &    .or. (time%clock > trigger%clock) )

  end function tem_time_gt_trigger
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Add two time definitions to each other (provides the plus operator).
  !!
  !! Entries set to huge (never) will be kept as this, all other entries are
  !! simply added to each other.
  elemental function tem_time_add_time(left, right) result(res)
    ! -------------------------------------------------------------------- !
    !> First operant in the addition of times.
    type(tem_time_type), intent(in) :: left

    !> Second operant in the addition of times.
    type(tem_time_type), intent(in) :: right

    !> Resulting sum of left and right.
    type(tem_time_type) :: res
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: nvrReal
    integer :: nvrInt
    ! -------------------------------------------------------------------- !

    nvrReal = huge(left%sim)
    nvrInt = huge(left%iter)

    if (       (left%sim  < nvrReal)            &
      &  .and. (right%sim < (nvrReal-left%sim)) ) then
      res%sim = left%sim + right%sim
    else
      res%sim = nvrReal
    end if

    if (       (left%iter  < nvrInt)             &
      &  .and. (right%iter < (nvrInt-left%iter)) ) then
      res%iter = left%iter + right%iter
    else
      res%iter = nvrInt
    end if

    if (       (left%clock  < nvrReal)              &
      &  .and. (right%clock < (nvrReal-left%clock)) ) then
      res%clock = left%clock + right%clock
    else
      res%clock = nvrReal
    end if

  end function tem_time_add_time
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Subtract right time definition from left (provides the minus operator).
  !!
  !! Entries set to huge (never) will be kept as this, all other entries are
  !! simply subtracted.
  elemental function tem_time_subtract_time(left, right) result(res)
    ! -------------------------------------------------------------------- !
    !> First operant in the addition of times.
    type(tem_time_type), intent(in) :: left

    !> Second operant in the addition of times.
    type(tem_time_type), intent(in) :: right

    !> Resulting sum of left and right.
    type(tem_time_type) :: res
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: nvrReal
    integer :: nvrInt
    ! -------------------------------------------------------------------- !

    nvrReal = huge(left%sim)
    nvrInt = huge(left%iter)

    if ( (left%sim  < nvrReal) ) then
      if (right%sim < nvrReal) then
         res%sim = left%sim - right%sim
      else
         res%sim = left%sim
      end if
    else
      res%sim = nvrReal
    end if

    if ( (left%iter  < nvrInt) ) then
      if (right%iter < nvrInt) then
         res%iter = left%iter - right%iter
      else
         res%iter = left%iter
      end if
    else
      res%iter = nvrInt
    end if

    if ( (left%clock  < nvrReal) ) then
      if (right%clock < nvrReal) then
         res%clock = left%clock - right%clock
      else
         res%clock = left%clock
      end if
    else
      res%clock = nvrReal
    end if

  end function tem_time_subtract_time
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Returns the last full interval before now.
  !!
  !! If interval is 0 or smaller, the result is put to 0. This holds for
  !! time component individually.
  elemental function tem_time_last_interval(now, interval) result(last)
    ! -------------------------------------------------------------------- !
    type(tem_time_type), intent(in) :: now
    type(tem_time_type), intent(in) :: interval
    type(tem_time_type) :: last
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    last%sim   = 0.0_rk
    last%iter  = 0
    last%clock = 0.0_rk

    if (interval%sim > 0.0_rk) then
      last%sim = floor(now%sim / interval%sim) * interval%sim
    end if
    if (interval%iter > 0) then
      last%iter = (now%iter / interval%iter) * interval%iter
    end if
    if (interval%clock > 0.0_rk) then
      last%clock = floor(now%clock / interval%clock) * interval%clock
    end if
  end function tem_time_last_interval
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Get the maximum of two time definitions.
  !!
  !! This is simply the maximum in each component.
  elemental function tem_time_max(left, right) result(res)
    ! -------------------------------------------------------------------- !
    !> First time operant to compare.
    type(tem_time_type), intent(in) :: left

    !> Second time operant to compare.
    type(tem_time_type), intent(in) :: right

    !> Resulting maximum time of left and right.
    type(tem_time_type) :: res
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    res%sim   = max(left%sim, right%sim)
    res%iter  = max(left%iter, right%iter)
    res%clock = max(left%clock, right%clock)

  end function tem_time_max
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Reset the time definition.
  !!
  !! All counters are reset to 0, and the starting clock is set to the current
  !! time.
  subroutine tem_time_reset(me)
    ! -------------------------------------------------------------------- !
    !> Time type to reset
    type(tem_time_type), intent(out) :: me
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    me%sim   = 0.0_rk
    me%iter  = 0
    me%clock = 0.0_rk

    me%clock_start = MPI_Wtime()

  end subroutine tem_time_reset
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Update the clock measurement in the time description.
  !!
  !! The clock component in me is updated with the help of MPI_Wtime.
  subroutine tem_time_set_clock(me)
    ! -------------------------------------------------------------------- !
    !> Time setting to update.
    type(tem_time_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    me%clock = MPI_Wtime() - me%clock_start

  end subroutine tem_time_set_clock
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Advance the time definition.
  !!
  !! Advance the time object by the given simulation time difference.
  !! Optionally, iterations might be advanced by more than one. If iter is not
  !! provided, iterations will be advanced by one.
  !! The running time is always automatically updated using MPI_Wtime.
  !! \note This is only valid, if me has been properly set by tem_time_reset!
  subroutine tem_time_advance(me, sim_dt, iter)
    ! -------------------------------------------------------------------- !
    !> Time definition to advance
    type(tem_time_type), intent(inout) :: me

    !> Increment in simulation time to add
    real(kind=rk) :: sim_dt

    !> If the number of iterations should not be increased by one,
    !! this optional parameter can be used to define the increment for the
    !! iterations. Default: 1.
    integer, optional :: iter
    ! -------------------------------------------------------------------- !
    integer :: iter_inc
    ! -------------------------------------------------------------------- !

    iter_inc = 1
    if (present(iter)) iter_inc = iter

    me%sim   = me%sim + sim_dt
    me%iter  = me%iter + iter_inc
    me%clock = MPI_Wtime() - me%clock_start

  end subroutine tem_time_advance
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Representation of a point in time, that should never happen.
  function tem_time_never()
    ! -------------------------------------------------------------------- !
    !> A time that represents never in all definitions.
    type(tem_time_type) :: tem_time_never
    ! -------------------------------------------------------------------- !

    tem_time_never%sim   = huge(tem_time_never%sim)
    tem_time_never%iter  = huge(tem_time_never%iter)
    tem_time_never%clock = huge(tem_time_never%clock)

  end function tem_time_never
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function returns, if the given time definition requires a reduction
  !! operation if used as a trigger.
  !!
  !! If there is a clock setting involved, a reduction is needed, as the
  !! wtimes are not synchronous across all processes.
  elemental function tem_time_needs_reduce(me) result(needs_reduce)
    ! -------------------------------------------------------------------- !
    !> Time definition to check of its need of a reduction.
    type(tem_time_type), intent(in) :: me

    !> Flag indicating, if this time setting requires a reduction.
    logical :: needs_reduce
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    needs_reduce = me%clock < huge(me%clock)

  end function tem_time_needs_reduce
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function returns true if tem_time_type is defined either as iter,
  !! sim or clock.
  pure function tem_time_isDefined(me) result(isDefined)
    ! -------------------------------------------------------------------- !
    type(tem_time_type), intent(in) :: me
    logical :: isDefined
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    isDefined = (      (me%iter < huge(me%iter)) &
      &           .or. (me%sim < huge(me%sim))   &
      &           .or. (me%clock < huge(me%clock)) )

  end function tem_time_isDefined
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Generate a time stamp from the simulation time in the given time
  !! definition.
  !!
  !! This basically generates a string identifying the solution time and
  !! writing it in a meaningful format, so that it can be easily recognized.
  function tem_time_sim_stamp(time) result(timeStamp)
    ! -------------------------------------------------------------------- !
    !> Time definition to create the stamp off.
    type(tem_time_type), intent(in) :: time

    !> String representation of the given simulation time.
    character(len=12) :: timeStamp
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    write(timeStamp, '(EN12.3)') time%sim

    ! remove leading empty spaces in the timestamp
    timeStamp = adjustl(timeStamp)

  end function  tem_time_sim_stamp
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Generate a time stamp from the iteration in the given time
  !! definition.
  !!
  !! This basically generates a string identifying the solution time and
  !! writing it in a meaningful format, so that it can be easily recognized.
  function tem_time_iter_stamp(time) result(timeStamp)
    ! -------------------------------------------------------------------- !
    !> Time definition to create the stamp off.
    type(tem_time_type), intent(in) :: time

    !> String representation of the given simulation time.
    character(len=labelLen) :: timeStamp
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    write(timeStamp, '(I0)') time%iter

    ! remove leading empty spaces in the timestamp
    timeStamp = adjustl(timeStamp)

  end function  tem_time_iter_stamp
  ! ************************************************************************ !

end module tem_time_module
