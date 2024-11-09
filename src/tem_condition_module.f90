! Copyright (c) 2012 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013, 2019-2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2015 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
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
!
!! This module contains type definition and simple routine to load conditions
!! Condition type has a `threshold` and an `operator` against which a quantity
!! can be compared. This can be used for flow properties based geometry
!! changes as well as for assessing the convergence.
!!
!! There may be multiple conditions bundled together in one table.
!! A single condition takes the following form:
!!
!!```lua
!!  condition = {threshold = 2.0e-10, operator = '<='}
!!```
!!
!! Multiple conditions can be defined as shown in the following
!! example.
!!
!!```lua
!!  condition = {
!!    {threshold = 1.e-15, operator = '<=' },
!!    {threshold = 1.e-12, operator = '<=' }
!!  }
!!```
!!
!! `threshold` needs to be a number, while `operator` needs to be a string.
!! The following options are possible for `operator`:
!! * '<'
!! * '='
!! * '>'
!! * '<='
!! * '>='
!!
module tem_condition_module

  ! include treelm modules
  use env_module,         only: rk, labelLen
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit
  use tem_float_module,   only: operator(.feq.)

  ! include aotus modules
  use aotus_module, only: flu_State,        &
    &                     aot_get_val,      &
    &                     aoterr_Fatal,     &
    &                     aoterr_WrongType, &
    &                     aoterr_NonExistent
  use aot_table_module, only: aot_table_open,   &
    &                         aot_table_close,  &
    &                         aot_table_length, &
    &                         aot_get_val
  use aot_out_module, only: aot_out_type,       &
    &                       aot_out_open,       &
    &                       aot_out_close,      &
    &                       aot_out_val,        &
    &                       aot_out_open_table, &
    &                       aot_out_close_table

  implicit none

  private

  public :: tem_condition_type
  public :: tem_load_condition
  public :: tem_condition_out
  public :: tem_condition_dump
  public :: tem_comparator


  !> Datatype containing different conditions to be checked
  !! Currently only threshold and operator are defined
  type tem_condition_type
    !> Contains the threshold value defined in lua file
    real(kind=rk) :: threshold
    !> Contains the operator defined in lua file
    character(len=2) :: operation
  end type tem_condition_type

  interface tem_condition_dump
    module procedure tem_condition_dump_vector
    module procedure tem_condition_dump_single
  end interface tem_condition_dump

  interface tem_condition_out
    module procedure tem_condition_out_vector
    module procedure tem_condition_out_single
  end interface tem_condition_out


contains


  ! ************************************************************************ !
  !> Load the condition table in case of convergence
  !!
  !! Example:
  !!```lua
  !! condition = {threshold = 2.0e-10, operator = '<='}
  !!```
  !!
  subroutine tem_load_condition( me, conf, parent )
    ! -------------------------------------------------------------------- !
    !>
    type(tem_condition_type), allocatable, intent(inout) :: me(:)
    !>
    type(flu_state), intent(in)             :: conf
    !>
    integer, intent(in)                     :: parent
    ! -------------------------------------------------------------------- !
    integer :: cond_handle          ! handle for the condition table
    integer :: sub_cond_handle      ! handle for subtables inside condition
    integer :: nCond                ! number of conditions
    integer :: iCond                ! index for condition loop
    ! -------------------------------------------------------------------- !
    !! Open the condition table
    call aot_table_open( L       = conf,        &
      &                  parent  = parent,      &
      &                  thandle = cond_handle, &
      &                  key     = 'condition'  )
    !! The tables inside condition table should be equal to the nVars
    !! If thats not the case we return an error message
    nCond = aot_table_length( L=conf, thandle=cond_handle )
    !! check single or multiple table
    call aot_table_open( L       = conf,            &
      &                  parent  = cond_handle,     &
      &                  thandle = sub_cond_handle, &
      &                  pos     = 1                )

    if (sub_cond_handle == 0) then
      call aot_table_close( L = conf, thandle = sub_cond_handle )
      ! just one table within
      allocate ( me(1) )
      call tem_load_cond_single( me(1), conf, cond_handle )
    else
      ! IF there are more tables within condition
      call aot_table_close( L = conf, thandle = sub_cond_handle )
      allocate ( me(nCond) )
      do iCond = 1, nCond
        call aot_table_open( L       = conf,            &
          &                  parent  = cond_handle,     &
          &                  thandle = sub_cond_handle, &
          &                  pos     = iCond            )
        call tem_load_cond_single( me(iCond), conf, sub_cond_handle )
        call aot_table_close( L = conf, thandle = sub_cond_handle )
      end do
    end if ! sub condition check
    call aot_table_close(L=conf, thandle=cond_handle )

  end subroutine tem_load_condition
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Load the conditions for geomIncr and convergence check within convergence
  !! conditions mean the operator and threshold against which the macroscopic
  !! variable has to be compared
  !!
  subroutine tem_load_cond_single( cond, conf, thandle )
    ! -------------------------------------------------------------------- !
    !>
    type(tem_condition_type), intent(inout) :: cond
    !>
    type(flu_state), intent(in)             :: conf
    !>
    integer, intent(in)                     :: thandle
    ! -------------------------------------------------------------------- !
    integer :: iError
    ! -------------------------------------------------------------------- !
    call aot_get_val( L       = conf,           &
      &               thandle = thandle,        &
      &               val     = cond%threshold, &
      &               ErrCode = iError,         &
      &               key     = 'threshold'     )
    if ( btest(iError, aoterr_Fatal) ) then
      write(logUnit(0),*) "Fatal Error: In reading 'threshold' in condition"
      if ( btest(iError, aoterr_NonExistent) ) &
        & write(logUnit(0),*) 'NonExistent.'
      if ( btest(iError, aoterr_WrongType) ) write(logUnit(0),*) 'WrongType.'
      call tem_abort()
    end if

    call aot_get_val( L       = conf,           &
      &               thandle = thandle,        &
      &               val     = cond%operation, &
      &               ErrCode = iError,         &
      &               key     = 'operator'      )
    if ( btest(iError, aoterr_Fatal) ) then
      write(logUnit(0),*) "Fatal Error: In reading 'operator' for condition"
      if ( btest(iError, aoterr_NonExistent) ) &
        & write(logUnit(0),*) 'NonExistent.'
      if( btest(iError, aoterr_WrongType) ) write(logUnit(0),*) 'WrongType.'
      call tem_abort()
    end if

  end subroutine tem_load_cond_single
  ! ************************************************************************ !


  ! ************************************************************************ !
  !! Return a logical if the input relation holds where the relation is:
  !! val _operation_ threshold
  !!
  function tem_comparator( val, operation, threshold ) result(comp)
    ! -------------------------------------------------------------------- !
    !>
    real(kind=rk), intent(in)      :: val
    !>
    character(len=2), intent(in)   :: operation
    !>
    real(kind=rk), intent(in)      :: threshold
    !> return value
    logical                        :: comp
    ! -------------------------------------------------------------------- !

    comp = .false.

    select case( trim(operation))
    case ('<')
      comp = (val < threshold)

    case ('<=')
      comp = (val <= threshold)

    case ('>')
      comp = (val > threshold)

    case ('>=')
      comp = (val >= threshold)

    case ('=')
      comp = (val .feq. threshold)

    end select

  end function tem_comparator
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Dumps array of condition to given unit
  subroutine tem_condition_dump_vector(me, outUnit)
    ! -------------------------------------------------------------------- !
    !> condition to write into the lua file
    type(tem_condition_type), intent(in) :: me(:)
    !> unit to write to
    integer, intent(in) :: outUnit
    ! -------------------------------------------------------------------- !
    ! aotus type handling the output to the file in lua format
    type(aot_out_type) :: conf
    ! -------------------------------------------------------------------- !

    call aot_out_open( put_conf = conf, outUnit = outUnit )
    call tem_condition_out_vector( me, conf )
    call aot_out_close( put_conf = conf )

  end subroutine tem_condition_dump_vector
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Dump single condition to given unit
  subroutine tem_condition_dump_single(me, outUnit)
    ! -------------------------------------------------------------------- !
    !> condition to write into the lua file
    type(tem_condition_type), intent(in) :: me
    !> unit to write to
    integer, intent(in) :: outUnit
    ! -------------------------------------------------------------------- !
    ! aotus type handling the output to the file in lua format
    type(aot_out_type) :: conf
    ! -------------------------------------------------------------------- !

    call aot_out_open( put_conf = conf, outUnit = outUnit )
    call tem_condition_out_single( me, conf )
    call aot_out_close( put_conf = conf )

  end subroutine tem_condition_dump_single
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Allows the output of array of condition to lua out
  subroutine tem_condition_out_vector(me, conf)
    ! -------------------------------------------------------------------- !
    !> condition to write into the lua file
    type(tem_condition_type), intent(in) :: me(:)
    !> aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    ! -------------------------------------------------------------------- !
    integer :: iCond
    ! -------------------------------------------------------------------- !

    call aot_out_open_table( put_conf = conf, tname='condition' )
    do iCond = 1,size(me)
      call tem_condition_out_single( me(iCond), conf, level=1 )
    end do
    call aot_out_close_table( put_conf = conf )

  end subroutine tem_condition_out_vector
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Allows the output of the single condition to lua out.
  !!
  !! The data is written into the file, the lunit is connected to.
  !! It is formatted as a Lua table.
  !!
  subroutine tem_condition_out_single(me, conf, level)
    ! -------------------------------------------------------------------- !
    !> condition to write into the lua file
    type(tem_condition_type), intent(in) :: me
    !> aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    !> to dump variable with key or without key
    integer, optional, intent(in) :: level
    ! -------------------------------------------------------------------- !
    integer :: level_loc
    ! -------------------------------------------------------------------- !

    if (present(level)) then
      level_loc = level
    else
      level_loc = 0
    end if

    if( level_loc == 0) then
      call aot_out_open_table( put_conf = conf, tname = 'condition' )
    else
      call aot_out_open_table( put_conf = conf )
    end if

    call aot_out_val( put_conf = conf,         &
      &               val      = me%operation, &
      &               vname    = 'operator'    )

    call aot_out_val( put_conf = conf,         &
      &               val      = me%threshold, &
      &               vname    = 'threshold'   )

    call aot_out_close_table( put_conf = conf )

  end subroutine tem_condition_out_single
  ! ************************************************************************ !

end module tem_condition_module
