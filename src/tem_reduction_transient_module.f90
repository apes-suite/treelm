! Copyright (c) 2012 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012-2013, 2015, 2019, 2021 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2013 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2013-2015, 2017, 2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015 Julia Moos <julia.moos@student.uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> author: Kannan Masilamani
!! summary: This module reduces values in time
!!
!! Transient reduction enables the definition of quantities that are reduced
!! over a given interval of iterations.
!! It is mainly meant to provide temporally averaged values, but minimum,
!! maximum and sum are also available.
!! For this a double buffer is used, where the current inverval is first filled
!! with the intermediate computation, while the resulting reduced value resides
!! in the other storage place.
!! This double buffer is provided for all elements in the domain.
!! When the value of the transient reduction is requested (via the getElement
!! function) the previously completed interval is returned.
!! This means the data is only useful after the first completed interval and
!! remains constant until the next interval is completed.
!!
!! A variable providing the temporal average of a variable called pressure would
!! for example be defined as follows:
!!
!!```lua
!!  variable = {
!!    {
!!       name = 'press_timeavg',
!!       ncomponents = 1,
!!       vartype = "operation",
!!       operation = {
!!         kind='reduction_transient',
!!         input_varname={'pressure'},
!!         reduction_transient = {
!!           kind    = 'average',
!!           nrecord = 1000
!!         }
!!       }
!!    }
!!  }
!!```
!!
!! Note that each transient variable definition increases the memory consumption
!! by two reals for every degree of freedom of the input variable in the domain.

module tem_reduction_transient_module

  ! include treelm modules
  use env_module,         only: rk, labelLen
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit
  use tem_time_module,    only: tem_time_type

  ! include aotus modules
  use aotus_module,     only: aoterr_Fatal, aoterr_NonExistent, flu_State, &
    &                         aoterr_WrongType, aot_get_val
  use aot_table_module, only: aot_table_open, aot_table_close, aot_table_length

  implicit none

  private

  public :: tem_reduction_transient_load
  public :: tem_reduction_transient_config_type
  public :: tem_reduction_transient_type
  public :: tem_reduction_transient_init
  public :: tem_reduction_transient_reset
  public :: tem_reduction_transient_update
  public :: tem_reduction_transient_getElement

  !> Contains transient reduction info loaded from variable table for
  !! reduction_transient operation kind
  type tem_reduction_transient_config_type
    !> Type of time reduction operation
    character(len=labelLen) :: reduceType

    !> number of interations to record
    integer :: nRecord = 0
  end type tem_reduction_transient_config_type

  !> all data needed for a transient reduction,
  !! operation to perform and necessary data from previous timesteps
  type tem_reduction_transient_type
    !> reduction info loaded from config file
    type(tem_reduction_transient_config_type) :: config

    !> Number of "recorded" previous iterations
    integer :: nTimes = 0

    !> Number of values to store in the double buffer val
    integer :: nEntries = 0

    !> number of components
    integer :: nComponents = 0

    !> Number of degrees of freedom
    integer :: nDofs

    !> Index of the storage for the currently filling (incomplete) position
    !! in val
    integer :: curr

    !> Index of the storage location with the previously completed interval
    !! holding the reduced value over that interval.
    integer :: last

    !> Double buffer to store data from previous timesteps
    !! size (nComponents*nDofs*tree%nElems,2)
    !! 2nd index is used to maintain last valid reduced value when
    !! nRecord is reached. It will be swapped to avoid copy operations
    real(kind=rk), allocatable :: val(:,:)
  end type tem_reduction_transient_type


contains


  ! ------------------------------------------------------------------------ !
  !> Read time reduction info from reduction_transient operation variable
  !!
  subroutine tem_reduction_transient_load( me, conf, parent )
    ! -------------------------------------------------------------------- !
    !> time reduction
    type(tem_reduction_transient_config_type), intent(out) :: me
    !> handle for lua file
    type(flu_State), intent(inout)                    :: conf
    !> operation table handle
    integer, intent(in)                               :: parent
    ! -------------------------------------------------------------------- !
    integer :: iError, reduce_handle
    ! -------------------------------------------------------------------- !
    write(logUnit(10),*) 'Loading reduction transient'
    call aot_table_open( L       = conf,                 &
      &                  thandle = reduce_handle,        &
      &                  parent  = parent,               &
      &                  key     = 'reduction_transient' )

    call aot_get_val( L       = conf,          &
      &               thandle = reduce_handle, &
      &               val     = me%nRecord,    &
      &               ErrCode = iError,        &
      &               key     = "nrecord"      )

    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*)'FATAL Error occured, while retrieving nRecord:'
      if (btest(iError, aoterr_NonExistent))        &
        & write(logUnit(1),*)'Variable not existent!'
      if (btest(iError, aoterr_WrongType))            &
        & write(logUnit(1),*)'Variable has wrong type!'
      write(logUnit(1),*)'STOPPING'
      call tem_abort()
    end if

    call aot_get_val( L       = conf,          &
      &               thandle = reduce_handle, &
      &               val     = me%reduceType, &
      &               ErrCode = iError,        &
      &               default = 'average',     &
      &               key     = "kind"         )

    write(logUnit(10),*) 'Reduction type: '//trim(me%reduceType)
    write(logUnit(10),*) 'nIter to record: ', me%nRecord
    call aot_table_close(L=conf, thandle=reduce_handle)

  end subroutine tem_reduction_transient_load
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Initialize time reduction.
  subroutine tem_reduction_transient_init(me, nElems, nComponents, nDofs)
    ! -------------------------------------------------------------------- !
    !> current time reduction
    type(tem_reduction_transient_type), intent(inout) :: me
    !> Number of elements in tree
    integer, intent(in) :: nElems
    !> nComponents of operation variable
    integer, intent(in) :: nComponents
    !> Number of degrees of freedom
    integer, intent(in) :: nDofs
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    me%nTimes = 0
    me%last = 1
    me%curr = 2
    me%nDofs = nDofs
    me%nComponents = nComponents
    me%nEntries = nElems * nComponents * nDofs

    allocate(me%val(me%nEntries,2))

    select case(trim(me%config%reduceType))
      case('min')
        me%val(:,me%curr) = huge(1.0_rk)

      case('max')
        me%val(:,me%curr) = tiny(1.0_rk)

      case('sum', 'average')
        me%val(:,me%curr) = 0.0_rk

      case default
        write(logUnit(1),*)'Unknown time reduction operation: '//       &
          &            trim(me%config%reduceType)
        call tem_abort()
    end select

  end subroutine tem_reduction_transient_init
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Reset time reduction.
  subroutine tem_reduction_transient_reset(me)
    ! -------------------------------------------------------------------- !
    !> current time reduction
    type(tem_reduction_transient_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !
    me%nTimes = 0
    select case(trim(me%config%reduceType))
      case('min')
        me%val(:, me%curr) = huge(1.0_rk)
      case('max')
        me%val(:, me%curr) = tiny(1.0_rk)
      case('sum', 'average')
        me%val(:, me%curr) = 0.0_rk
    end select

  end subroutine tem_reduction_transient_reset
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Apply time reduction on given res
  subroutine tem_reduction_transient_update(me, res)
    ! -------------------------------------------------------------------- !
    !> Time reduction data to update
    type(tem_reduction_transient_type), intent(inout) :: me

    !> Current values of the variable to reduce.
    real(kind=rk), intent(in) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: idx
    ! -------------------------------------------------------------------- !
    me%nTimes = me%nTimes + 1

    select case( trim(me%config%reduceType) )
      case('min')
        do idx=1,me%nEntries
          me%val(idx, me%curr) = min( me%val(idx, me%curr),  res(idx) )
        end do

      case('max')
        do idx=1,me%nEntries
          me%val(idx, me%curr) = max( me%val(idx, me%curr),  res(idx) )
        end do

      case('sum', 'average')
        me%val(:, me%curr) = me%val(:, me%curr) + res(:)
    end select

    ! Check whether interval is completed
    if (me%nTimes == me%config%nRecord) then
      select case( trim(me%config%reduceType) )
        case('average')
          me%val(:, me%curr) = me%val(:, me%curr) / me%nTimes
      end select

      ! swap curr and last
      me%curr = mod(me%curr,2)+1
      me%last = mod(me%last,2)+1

      ! Reset current val reduction when nTimes has reached nRecord
      call tem_reduction_transient_reset(me)
    end if

  end subroutine tem_reduction_transient_update
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This routine returns the time reduced value for given elemPos
  subroutine tem_reduction_transient_getElement(me, elemPos, nElems, nDofs, res)
    ! -------------------------------------------------------------------- !
    !> current time reduction
    type(tem_reduction_transient_type), intent(in) :: me
    !> Position of elements in global tree is same as me%val
    integer, intent(in) :: elemPos(:)
    !> Number of elements to return
    integer, intent(in) :: nElems
    !> Number of degrees of freedom to return
    integer, intent(in) :: nDofs
    !> Result array
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: iElem, iDof, iComp
    integer :: eSize_val, eSize_res, offset
    ! -------------------------------------------------------------------- !
    eSize_val = me%nComponents*me%nDofs
    eSize_res = me%nComponents*nDofs
    do iElem = 1, nElems
      do iDof = 1, nDofs
        do iComp = 1, me%nComponents
          offset = (iDof-1)*me%nComponents + iComp
          res( (iElem-1)*eSize_res + offset )                          &
            & = me%val( (elemPos(iElem)-1)*eSize_val + offset, me%last )
        end do
      end do
    end do

  end subroutine tem_reduction_transient_getElement
  ! ------------------------------------------------------------------------ !

end module tem_reduction_transient_module
! **************************************************************************** !
