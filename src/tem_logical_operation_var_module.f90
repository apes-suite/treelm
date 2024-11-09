! Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2019 Harald Klimach <harald.klimach@uni-siegen.de>
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
! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015, 2018, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
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
!--------------------------------------------
!    A O S - Array of structures layout new
!-------------------------------------------
! Access to get_point value output
! Access to get_element value output
! ****************************************************************************** !
!> This module provides logical operations for variables in the variable
!! system.
!!
module tem_logical_operation_var_module

  use, intrinsic :: iso_c_binding,  only: c_f_pointer

  use env_module,                   only: rk, eps
  use tem_float_module,             only: operator(.feq.), &
    &                                     operator(.fne.), &
    &                                     operator(.fgt.), &
    &                                     operator(.fge.), &
    &                                     operator(.flt.), &
    &                                     operator(.fle.)
  use tem_operation_module,         only: tem_varSys_op_data_type
  use tem_varSys_module,            only: tem_varSys_type,   &
    &                                     tem_varSys_op_type
  use tem_time_module,              only: tem_time_type
  use treelmesh_module,             only: treelmesh_type
  use tem_logging_module,           only: logUnit, llwarning
  use tem_aux_module,               only: tem_abort
  use tem_grow_array_module,        only: grw_intArray_type

  implicit none

  private

  public :: evalLogicalAnd_forPoint
  public :: evalLogicalOr_forPoint
  public :: evalLogicalGreater_forPoint
  public :: evalLogicalGreaterOrEqual_forPoint
  public :: evalLogicalLess_forPoint
  public :: evalLogicalLessOrEqual_forPoint
  public :: evalLogicalEqual_forPoint
  public :: evalLogicalNotEqual_forPoint

  public :: evalLogicalAnd_forElement
  public :: evalLogicalOr_forElement
  public :: evalLogicalGreater_forElement
  public :: evalLogicalGreaterOrEqual_forElement
  public :: evalLogicalLess_forElement
  public :: evalLogicalLessOrEqual_forElement
  public :: evalLogicalEqual_forElement
  public :: evalLogicalNotEqual_forElement


  public :: evalLogicalAnd_fromIndex
  public :: evalLogicalOr_fromIndex
  public :: evalLogicalGreater_fromIndex
  public :: evalLogicalGreaterOrEqual_fromIndex
  public :: evalLogicalLess_fromIndex
  public :: evalLogicalLessOrEqual_fromIndex
  public :: evalLogicalEqual_fromIndex
  public :: evalLogicalNotEqual_fromIndex

  ! The following routines have been made public solely for unit tests
  public :: logicalToReal
  public :: logicalToRealArray
  public :: realToLogical
  public :: realToLogicalArray

  public :: numFalse, numTrue

  real(kind=rk), parameter :: numTrue = 1._rk
  real(kind=rk), parameter :: numFalse = 0._rk

contains


  ! ************************************************************************** !
  !> Converts a logical into a real.
  !!
  !! 0 equals to false, everything else equals to true.
  function logicalToReal(value) result(res)
    ! ------------------------------------------------------------------------ !
    !> The value to interpret as boolean
    logical, intent(in) :: value
    !> The value interpreted as boolean
    real(kind=rk) :: res
    ! ------------------------------------------------------------------------ !

    if (value) then
      res = numTrue
    else
      res = numFalse
    end if

  end function logicalToReal
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Converts an array of logicals into an array of reals.
  !!
  !! 0 equals to false, everything else equals to true.
  function logicalToRealArray(value,n) result(res)
    ! ------------------------------------------------------------------------ !
    !> The number of values in the input array
    integer, intent(in) :: n
    !> The value to interpret as boolean
    logical, intent(in) :: value(n)
    !> The values interpreted as booleans
    real(kind=rk) :: res(n)
    ! ------------------------------------------------------------------------ !

    res = numFalse
    where(value) res = numTrue

  end function logicalToRealArray
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Converts a real into a logical.
  !!
  !! 0 equals to false, everything else equals to true.
  function realToLogical(value) result(res)
    ! ------------------------------------------------------------------------ !
    !> The value to interpret as boolean
    real(kind=rk), intent(in) :: value
    !> The value interpreted as boolean
    logical :: res
    ! ------------------------------------------------------------------------ !

    res = .not. (value .feq. numFalse)

  end function realToLogical
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Converts an array of reals into an array of logicals.
  !!
  !! 0 equals to false, everything else equals to true.
  function realToLogicalArray(value,n) result(res)
    ! ------------------------------------------------------------------------ !
    !> The number of values in the input array
    integer, intent(in) :: n
    !> The value to interpret as boolean
    real(kind=rk), intent(in) :: value(n)
    !> The values interpreted as booleans
    logical :: res(n)
    ! ------------------------------------------------------------------------ !
    integer :: ii
    ! ------------------------------------------------------------------------ !

    do ii = 1, n
      res(ii) = .not. (value(ii) .feq. numFalse)
    end do

  end function realToLogicalArray
  ! ************************************************************************** !


  subroutine evalLogicalAnd_forPoint(fun, varsys, point, time, tree, nPnts, res)
    !---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !---------------------------------------------------------------------------
    !> first dimension is nElems, second dimension is the number of input
    !! variables.
    real(kind=rk), allocatable :: input_varRes(:,:)
    integer :: iDep, posDepVar
    !---------------------------------------------------------------------------

    allocate(input_varRes(nPnts*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_point( &
        & varSys = varSys,                         &
        & point  = point,                          &
        & time   = time,                           &
        & tree   = tree,                           &
        & nPnts  = nPnts,                          &
        & res    = input_varRes(:, iDep)           )

    end do

    res = logicalToRealArray(                                  &
      & value = realToLogicalArray(input_varRes(:, 1), nPnts)  &
      &   .and. realToLogicalArray(input_varRes(:, 2), nPnts), &
      & n     = nPnts                                          )

    deallocate(input_varRes)

  end subroutine evalLogicalAnd_forPoint


  subroutine evalLogicalAnd_forElement( fun, varsys, elempos, time, tree, &
      &                                 nElems, nDofs, res                )
    !--------------------------------------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Position of the TreeID of the element to get the variable for in the
    !! global treeID list.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !--------------------------------------------------------------------------!
    integer :: iDep, posDepVar, iComp, iElem, lIdx, uIdx
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    !--------------------------------------------------------------------------!


    ! nInputs must be two
    allocate(input_varRes(nElems*nDofs*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_element(                        &
        &                                   varSys  = varSys,               &
        &                                   elemPos = elemPos,              &
        &                                   time    = time,                 &
        &                                   tree    = tree,                 &
        &                                   nElems  = nElems,               &
        &                                   nDofs   = nDofs,                &
        &                                   res     = input_varRes(:, iDep) )

    end do

    do iComp = 1, fun%nComponents
      do iElem = 1, nElems
        ! We only take the first mode into consideration, but store the result
        ! value for all modes.
        lIdx = (( ielem-1)* fun%ncomponents* ndofs+( 1-1)* fun%ncomponents+icomp)
        uIdx = (( ielem-1)* fun%ncomponents* ndofs+( ndofs-1)* fun%ncomponents+icomp)

        res( lIdx : uIdx ) = logicalToReal(                  &
          & realToLogical( input_varRes( lIdx, 1 ) )         &
          &   .and. realToLogical( input_varRes( lIdx, 2 ) ) )
      end do
    end do

    deallocate(input_varRes)

  end subroutine evalLogicalAnd_forElement


  subroutine evalLogicalAnd_fromIndex( fun, varSys, time, iLevel, idx, &
      &                                idxLen, nVals, res              )
    !--------------------------------------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !--------------------------------------------------------------------------!
    integer :: iDep, posDepVar
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    type(tem_varSys_op_data_type), pointer :: fPtr
    !--------------------------------------------------------------------------!

    call C_F_POINTER( fun%method_Data, fPtr )


    ! nInputs must be two
    allocate(input_varRes(nVals*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable for indices
      !rent derive dependent variable
      call varSys%method%val(posDepVar)%get_valOfIndex( &
        & varSys  = varSys,                             &
        & time    = time,                               &
        & iLevel  = iLevel,                             &
        & idx     = fPtr%input_pntIndex(iDep)           &
        &           %indexLvl(iLevel)%val( idx(:) ),    &
        & idxLen  = idxLen,                             &
        & nVals   = nVals,                              &
        & res     = input_varRes(:,iDep)                )
    end do

    res = logicalToRealArray(                                  &
      & value = realToLogicalArray(input_varRes(:, 1), nVals)  &
      &   .and. realToLogicalArray(input_varRes(:, 2), nVals), &
      & n     = nVals                                          )

    deallocate(input_varRes)

  end subroutine evalLogicalAnd_fromIndex


  subroutine evalLogicalOr_forPoint(fun, varsys, point, time, tree, nPnts, res)
    !---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !---------------------------------------------------------------------------
    !> first dimension is nElems, second dimension is the number of input
    !! variables.
    real(kind=rk), allocatable :: input_varRes(:,:)
    integer :: iDep, posDepVar
    !---------------------------------------------------------------------------

    allocate(input_varRes(nPnts*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_point( &
        & varSys = varSys,                         &
        & point  = point,                          &
        & time   = time,                           &
        & tree   = tree,                           &
        & nPnts  = nPnts,                          &
        & res    = input_varRes(:, iDep)           )

    end do

    res = logicalToRealArray(                                 &
      & value = realToLogicalArray(input_varRes(:, 1), nPnts) &
      &   .or. realToLogicalArray(input_varRes(:, 2), nPnts), &
      & n     = nPnts                                         )

    deallocate(input_varRes)

  end subroutine evalLogicalOr_forPoint


  subroutine evalLogicalOr_forElement( fun, varsys, elempos, time, tree, &
      &                                nElems, nDofs, res                )
    !--------------------------------------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Position of the TreeID of the element to get the variable for in the
    !! global treeID list.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !--------------------------------------------------------------------------!
    integer :: iDep, posDepVar, iComp, iElem, lIdx, uIdx
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    !--------------------------------------------------------------------------!


    ! nInputs must be two
    allocate(input_varRes(nElems*nDofs*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_element(                        &
        &                                   varSys  = varSys,               &
        &                                   elemPos = elemPos,              &
        &                                   time    = time,                 &
        &                                   tree    = tree,                 &
        &                                   nElems  = nElems,               &
        &                                   nDofs   = nDofs,                &
        &                                   res     = input_varRes(:, iDep) )

    end do

    do iComp = 1, fun%nComponents
      do iElem = 1, nElems
        ! We only take the first mode into consideration, but store the result
        ! value for all modes.
        lIdx = (( ielem-1)* fun%ncomponents* ndofs+( 1-1)* fun%ncomponents+icomp)
        uIdx = (( ielem-1)* fun%ncomponents* ndofs+( ndofs-1)* fun%ncomponents+icomp)

        res( lIdx : uIdx ) = logicalToReal(                 &
          & realToLogical( input_varRes( lIdx, 1 ) )        &
          &   .or. realToLogical( input_varRes( lIdx, 2 ) ) )
      end do
    end do

    deallocate(input_varRes)

  end subroutine evalLogicalOr_forElement


  subroutine evalLogicalOr_fromIndex( fun, varSys, time, iLevel, idx, &
      &                               idxLen, nVals, res                    )
    !--------------------------------------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !--------------------------------------------------------------------------!
    integer :: iDep, posDepVar
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    type(tem_varSys_op_data_type), pointer :: fPtr
    !--------------------------------------------------------------------------!

    call C_F_POINTER( fun%method_Data, fPtr )


    ! nInputs must be two
    allocate(input_varRes(nVals*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable for indices
      !rent derive dependent variable
      call varSys%method%val(posDepVar)%get_valOfIndex( &
        & varSys  = varSys,                             &
        & time    = time,                               &
        & iLevel  = iLevel,                             &
        & idx     = fPtr%input_pntIndex(iDep)           &
        &           %indexLvl(iLevel)%val( idx(:) ),    &
        & idxLen  = idxLen,                             &
        & nVals   = nVals,                              &
        & res     = input_varRes(:,iDep)                )
    end do

    res = logicalToRealArray(                                 &
      & value = realToLogicalArray(input_varRes(:, 1), nVals) &
      &   .or. realToLogicalArray(input_varRes(:, 2), nVals), &
      & n     = nVals                                         )

    deallocate(input_varRes)

  end subroutine evalLogicalOr_fromIndex


  subroutine evalLogicalGreater_forPoint( fun, varsys, point, time, tree, &
      &                                   nPnts, res                      )
    !---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !---------------------------------------------------------------------------
    !> first dimension is nElems, second dimension is the number of input
    !! variables.
    real(kind=rk), allocatable :: input_varRes(:,:)
    integer :: iDep, posDepVar
    !---------------------------------------------------------------------------

    allocate(input_varRes(nPnts*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_point( &
        & varSys = varSys,                         &
        & point  = point,                          &
        & time   = time,                           &
        & tree   = tree,                           &
        & nPnts  = nPnts,                          &
        & res    = input_varRes(:, iDep)           )

    end do

    res = logicalToRealArray(                            &
      & value = input_varRes(:, 1) > input_varRes(:, 2), &
      & n     = nPnts                                    )

    deallocate(input_varRes)

  end subroutine evalLogicalGreater_forPoint


  subroutine evalLogicalGreater_forElement( fun, varsys, elempos, time, tree, &
      &                                     nElems, nDofs, res                )
    !--------------------------------------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Position of the TreeID of the element to get the variable for in the
    !! global treeID list.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !--------------------------------------------------------------------------!
    integer :: iDep, posDepVar, iComp, iElem, lIdx, uIdx
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    !--------------------------------------------------------------------------!


    ! nInputs must be two
    allocate(input_varRes(nElems*nDofs*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_element(                        &
        &                                   varSys  = varSys,               &
        &                                   elemPos = elemPos,              &
        &                                   time    = time,                 &
        &                                   tree    = tree,                 &
        &                                   nElems  = nElems,               &
        &                                   nDofs   = nDofs,                &
        &                                   res     = input_varRes(:, iDep) )

    end do

    do iComp = 1, fun%nComponents
      do iElem = 1, nElems
        ! We only take the first mode into consideration, but store the result
        ! value for all modes.
        lIdx = (( ielem-1)* fun%ncomponents* ndofs+( 1-1)* fun%ncomponents+icomp)
        uIdx = (( ielem-1)* fun%ncomponents* ndofs+( ndofs-1)* fun%ncomponents+icomp)

        res( lIdx : uIdx ) = logicalToReal(                       &
          & input_varRes( lIdx, 1 ) .fgt. input_varRes( lIdx, 2 ) )
      end do
    end do

    deallocate(input_varRes)


  end subroutine evalLogicalGreater_forElement


  subroutine evalLogicalGreater_fromIndex( fun, varSys, time, iLevel, &
      &                                    idx, idxLen, nVals, res    )
    !--------------------------------------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !--------------------------------------------------------------------------!
    integer :: iDep, posDepVar
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    type(tem_varSys_op_data_type), pointer :: fPtr
    !--------------------------------------------------------------------------!

    call C_F_POINTER( fun%method_Data, fPtr )


    ! nInputs must be two
    allocate(input_varRes(nVals*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable for indices
      !rent derive dependent variable
      call varSys%method%val(posDepVar)%get_valOfIndex( &
        & varSys  = varSys,                             &
        & time    = time,                               &
        & iLevel  = iLevel,                             &
        & idx     = fPtr%input_pntIndex(iDep)           &
        &           %indexLvl(iLevel)%val( idx(:) ),    &
        & idxLen  = idxLen,                             &
        & nVals   = nVals,                              &
        & res     = input_varRes(:,iDep)                )
    end do

    res = logicalToRealArray(                            &
      & value = input_varRes(:, 1) > input_varRes(:, 2), &
      & n     = nVals                                    )

    deallocate(input_varRes)

  end subroutine evalLogicalGreater_fromIndex


  subroutine evalLogicalGreaterOrEqual_forPoint( fun, varsys, point, time, &
      &                                          tree, nPnts, res          )
    !---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !---------------------------------------------------------------------------
    !> first dimension is nElems, second dimension is the number of input
    !! variables.
    real(kind=rk), allocatable :: input_varRes(:,:)
    integer :: iDep, posDepVar
    !---------------------------------------------------------------------------

    allocate(input_varRes(nPnts*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_point( &
        & varSys = varSys,                         &
        & point  = point,                          &
        & time   = time,                           &
        & tree   = tree,                           &
        & nPnts  = nPnts,                          &
        & res    = input_varRes(:, iDep)           )

    end do

    res = logicalToRealArray(                             &
      & value = input_varRes(:, 1) >= input_varRes(:, 2), &
      & n     = nPnts                                     )

    deallocate(input_varRes)

  end subroutine evalLogicalGreaterOrEqual_forPoint


  subroutine evalLogicalGreaterOrEqual_forElement( fun, varsys, elempos, time, &
      &                                            tree, nElems, nDofs, res    )
    !--------------------------------------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Position of the TreeID of the element to get the variable for in the
    !! global treeID list.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !--------------------------------------------------------------------------!
    integer :: iDep, posDepVar, iComp, iElem, lIdx, uIdx
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    !--------------------------------------------------------------------------!


    ! nInputs must be two
    allocate(input_varRes(nElems*nDofs*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_element(                        &
        &                                   varSys  = varSys,               &
        &                                   elemPos = elemPos,              &
        &                                   time    = time,                 &
        &                                   tree    = tree,                 &
        &                                   nElems  = nElems,               &
        &                                   nDofs   = nDofs,                &
        &                                   res     = input_varRes(:, iDep) )

    end do

    do iComp = 1, fun%nComponents
      do iElem = 1, nElems
        ! We only take the first mode into consideration, but store the result
        ! value for all modes.
        lIdx = (( ielem-1)* fun%ncomponents* ndofs+( 1-1)* fun%ncomponents+icomp)
        uIdx = (( ielem-1)* fun%ncomponents* ndofs+( ndofs-1)* fun%ncomponents+icomp)

        res( lIdx : uIdx ) = logicalToReal(                       &
          & input_varRes( lIdx, 1 ) .fge. input_varRes( lIdx, 2 ) )
      end do
    end do

    deallocate(input_varRes)

  end subroutine evalLogicalGreaterOrEqual_forElement


  subroutine evalLogicalGreaterOrEqual_fromIndex( fun, varSys, time, iLevel, &
      &                                           idx, idxLen, nVals, res    )
    !--------------------------------------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !--------------------------------------------------------------------------!
    integer :: iDep, posDepVar
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    type(tem_varSys_op_data_type), pointer :: fPtr
    !--------------------------------------------------------------------------!

    call C_F_POINTER( fun%method_Data, fPtr )


    ! nInputs must be two
    allocate(input_varRes(nVals*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable for indices
      !rent derive dependent variable
      call varSys%method%val(posDepVar)%get_valOfIndex( &
        & varSys  = varSys,                             &
        & time    = time,                               &
        & iLevel  = iLevel,                             &
        & idx     = fPtr%input_pntIndex(iDep)           &
        &           %indexLvl(iLevel)%val( idx(:) ),    &
        & idxLen  = idxLen,                             &
        & nVals   = nVals,                              &
        & res     = input_varRes(:,iDep)                )
    end do

    res = logicalToRealArray(                             &
      & value = input_varRes(:, 1) >= input_varRes(:, 2), &
      & n     = nVals                                     )

    deallocate(input_varRes)

  end subroutine evalLogicalGreaterOrEqual_fromIndex


  subroutine evalLogicalLess_forPoint( fun, varsys, point, time, tree, nPnts, &
      &                                res                                    )
    !---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !---------------------------------------------------------------------------
    !> first dimension is nElems, second dimension is the number of input
    !! variables.
    real(kind=rk), allocatable :: input_varRes(:,:)
    integer :: iDep, posDepVar
    !---------------------------------------------------------------------------

    allocate(input_varRes(nPnts*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_point( &
        & varSys = varSys,                         &
        & point  = point,                          &
        & time   = time,                           &
        & tree   = tree,                           &
        & nPnts  = nPnts,                          &
        & res    = input_varRes(:, iDep)           )

    end do

    res = logicalToRealArray(                            &
      & value = input_varRes(:, 1) < input_varRes(:, 2), &
      & n     = nPnts                                    )

    deallocate(input_varRes)

  end subroutine evalLogicalLess_forPoint


  subroutine evalLogicalLess_forElement( fun, varsys, elempos, time, tree, &
      &                                  nElems, nDofs, res                )
    !--------------------------------------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Position of the TreeID of the element to get the variable for in the
    !! global treeID list.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !--------------------------------------------------------------------------!
    integer :: iDep, posDepVar, iComp, iElem, lIdx, uIdx
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    !--------------------------------------------------------------------------!


    ! nInputs must be two
    allocate(input_varRes(nElems*nDofs*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_element(                        &
        &                                   varSys  = varSys,               &
        &                                   elemPos = elemPos,              &
        &                                   time    = time,                 &
        &                                   tree    = tree,                 &
        &                                   nElems  = nElems,               &
        &                                   nDofs   = nDofs,                &
        &                                   res     = input_varRes(:, iDep) )

    end do

    do iComp = 1, fun%nComponents
      do iElem = 1, nElems
        ! We only take the first mode into consideration, but store the result
        ! value for all modes.
        lIdx = (( ielem-1)* fun%ncomponents* ndofs+( 1-1)* fun%ncomponents+icomp)
        uIdx = (( ielem-1)* fun%ncomponents* ndofs+( ndofs-1)* fun%ncomponents+icomp)

        res( lIdx : uIdx ) = logicalToReal(                       &
          & input_varRes( lIdx, 1 ) .flt. input_varRes( lIdx, 2 ) )
      end do
    end do

    deallocate(input_varRes)

  end subroutine evalLogicalLess_forElement


  subroutine evalLogicalLess_fromIndex( fun, varSys, time, iLevel, idx, &
      &                                 idxLen, nVals, res                    )
    !--------------------------------------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !--------------------------------------------------------------------------!
    integer :: iDep, posDepVar
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    type(tem_varSys_op_data_type), pointer :: fPtr
    !--------------------------------------------------------------------------!

    call C_F_POINTER( fun%method_Data, fPtr )


    ! nInputs must be two
    allocate(input_varRes(nVals*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable for indices
      !rent derive dependent variable
      call varSys%method%val(posDepVar)%get_valOfIndex( &
        & varSys  = varSys,                             &
        & time    = time,                               &
        & iLevel  = iLevel,                             &
        & idx     = fPtr%input_pntIndex(iDep)           &
        &           %indexLvl(iLevel)%val( idx(:) ),    &
        & idxLen  = idxLen,                             &
        & nVals   = nVals,                              &
        & res     = input_varRes(:,iDep)                )
    end do

    res = logicalToRealArray(                            &
      & value = input_varRes(:, 1) < input_varRes(:, 2), &
      & n     = nVals                                    )

    deallocate(input_varRes)

  end subroutine evalLogicalLess_fromIndex


  subroutine evalLogicalLessOrEqual_forPoint( fun, varsys, point, time, tree, &
      &                                       nPnts, res                      )
    !---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !---------------------------------------------------------------------------
    !> first dimension is nElems, second dimension is the number of input
    !! variables.
    real(kind=rk), allocatable :: input_varRes(:,:)
    integer :: iDep, posDepVar
    !---------------------------------------------------------------------------

    allocate(input_varRes(nPnts*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_point( &
        & varSys = varSys,                         &
        & point  = point,                          &
        & time   = time,                           &
        & tree   = tree,                           &
        & nPnts  = nPnts,                          &
        & res    = input_varRes(:, iDep)           )

    end do

    res = logicalToRealArray(                             &
      & value = input_varRes(:, 1) <= input_varRes(:, 2), &
      & n     = nPnts                                     )

    deallocate(input_varRes)

  end subroutine evalLogicalLessOrEqual_forPoint


  subroutine evalLogicalLessOrEqual_forElement( fun, varsys, elempos, time, &
      &                                         tree, nElems, nDofs, res    )
    !--------------------------------------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Position of the TreeID of the element to get the variable for in the
    !! global treeID list.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !--------------------------------------------------------------------------!
    integer :: iDep, posDepVar, iComp, iElem, lIdx, uIdx
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    !--------------------------------------------------------------------------!


    ! nInputs must be two
    allocate(input_varRes(nElems*nDofs*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_element(                        &
        &                                   varSys  = varSys,               &
        &                                   elemPos = elemPos,              &
        &                                   time    = time,                 &
        &                                   tree    = tree,                 &
        &                                   nElems  = nElems,               &
        &                                   nDofs   = nDofs,                &
        &                                   res     = input_varRes(:, iDep) )

    end do

    do iComp = 1, fun%nComponents
      do iElem = 1, nElems
        ! We only take the first mode into consideration, but store the result
        ! value for all modes.
        lIdx = (( ielem-1)* fun%ncomponents* ndofs+( 1-1)* fun%ncomponents+icomp)
        uIdx = (( ielem-1)* fun%ncomponents* ndofs+( ndofs-1)* fun%ncomponents+icomp)

        res( lIdx : uIdx ) = logicalToReal(                       &
          & input_varRes( lIdx, 1 ) .fle. input_varRes( lIdx, 2 ) )
      end do
    end do

    deallocate(input_varRes)

  end subroutine evalLogicalLessOrEqual_forElement


  subroutine evalLogicalLessOrEqual_fromIndex( fun, varSys, time, iLevel,    &
      &                                        idx, idxLen, nVals, res )
    !--------------------------------------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !--------------------------------------------------------------------------!
    integer :: iDep, posDepVar
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    type(tem_varSys_op_data_type), pointer :: fPtr
    !--------------------------------------------------------------------------!

    call C_F_POINTER( fun%method_Data, fPtr )


    ! nInputs must be two
    allocate(input_varRes(nVals*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable for indices
      !rent derive dependent variable
      call varSys%method%val(posDepVar)%get_valOfIndex( &
        & varSys  = varSys,                             &
        & time    = time,                               &
        & iLevel  = iLevel,                             &
        & idx     = fPtr%input_pntIndex(iDep)           &
        &           %indexLvl(iLevel)%val( idx(:) ),    &
        & idxLen  = idxLen,                             &
        & nVals   = nVals,                              &
        & res     = input_varRes(:,iDep)                )
    end do

    res = logicalToRealArray(                             &
      & value = input_varRes(:, 1) <= input_varRes(:, 2), &
      & n     = nVals                                     )

    deallocate(input_varRes)

  end subroutine evalLogicalLessOrEqual_fromIndex


  subroutine evalLogicalEqual_forPoint( fun, varsys, point, time, tree, nPnts, &
      &                                 res                                    )
    !---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !---------------------------------------------------------------------------
    !> first dimension is nElems, second dimension is the number of input
    !! variables.
    real(kind=rk), allocatable :: input_varRes(:,:)
    integer :: iDep, posDepVar
    !---------------------------------------------------------------------------

    allocate(input_varRes(nPnts*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_point( &
        & varSys = varSys,                         &
        & point  = point,                          &
        & time   = time,                           &
        & tree   = tree,                           &
        & nPnts  = nPnts,                          &
        & res    = input_varRes(:, iDep)           )

    end do

    res = 1.0
    where((input_varRes(:, 1) - input_varRes(:, 2)) > eps ) res = 0.0

    deallocate(input_varRes)

  end subroutine evalLogicalEqual_forPoint


  subroutine evalLogicalEqual_forElement( fun, varsys, elempos, time, tree, &
      &                                   nElems, nDofs, res                )
    !--------------------------------------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Position of the TreeID of the element to get the variable for in the
    !! global treeID list.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !--------------------------------------------------------------------------!
    integer :: iDep, posDepVar, iComp, iElem, lIdx, uIdx
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    !--------------------------------------------------------------------------!


    ! nInputs must be two
    allocate(input_varRes(nElems*nDofs*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_element(                        &
        &                                   varSys  = varSys,               &
        &                                   elemPos = elemPos,              &
        &                                   time    = time,                 &
        &                                   tree    = tree,                 &
        &                                   nElems  = nElems,               &
        &                                   nDofs   = nDofs,                &
        &                                   res     = input_varRes(:, iDep) )

    end do

    do iComp = 1, fun%nComponents
      do iElem = 1, nElems
        ! We only take the first mode into consideration, but store the result
        ! value for all modes.
        lIdx = (( ielem-1)* fun%ncomponents* ndofs+( 1-1)* fun%ncomponents+icomp)
        uIdx = (( ielem-1)* fun%ncomponents* ndofs+( ndofs-1)* fun%ncomponents+icomp)

        res( lIdx : uIdx ) = logicalToReal(                       &
          & input_varRes( lIdx, 1 ) .feq. input_varRes( lIdx, 2 ) )
      end do
    end do

    deallocate(input_varRes)

  end subroutine evalLogicalEqual_forElement


  subroutine evalLogicalEqual_fromIndex( fun, varSys, time, iLevel, idx, &
      &                                  idxLen, nVals, res                    )
    !--------------------------------------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !--------------------------------------------------------------------------!
    integer :: iDep, posDepVar
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    type(tem_varSys_op_data_type), pointer :: fPtr
    !--------------------------------------------------------------------------!

    call C_F_POINTER( fun%method_Data, fPtr )


    ! nInputs must be two
    allocate(input_varRes(nVals*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable for indices
      !rent derive dependent variable
      call varSys%method%val(posDepVar)%get_valOfIndex( &
        & varSys  = varSys,                             &
        & time    = time,                               &
        & iLevel  = iLevel,                             &
        & idx     = fPtr%input_pntIndex(iDep)           &
        &           %indexLvl(iLevel)%val( idx(:) ),    &
        & idxLen  = idxLen,                             &
        & nVals   = nVals,                              &
        & res     = input_varRes(:,iDep)                )
    end do

    res = numTrue
    where((input_varRes(:, 1) - input_varRes(:, 2)) > eps ) res = numFalse

    deallocate(input_varRes)

  end subroutine evalLogicalEqual_fromIndex


  subroutine evalLogicalNotEqual_forPoint( fun, varsys, point, time, tree, &
      &                                    nPnts, res                      )
    !---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !---------------------------------------------------------------------------
    !> first dimension is nElems, second dimension is the number of input
    !! variables.
    real(kind=rk), allocatable :: input_varRes(:,:)
    integer :: iDep, posDepVar
    !---------------------------------------------------------------------------

    allocate(input_varRes(nPnts*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_point( &
        & varSys = varSys,                         &
        & point  = point,                          &
        & time   = time,                           &
        & tree   = tree,                           &
        & nPnts  = nPnts,                          &
        & res    = input_varRes(:, iDep)           )

    end do

    res = numFalse
    where((input_varRes(:, 1) - input_varRes(:, 2)) > eps ) res = numTrue

    deallocate(input_varRes)

  end subroutine evalLogicalNotEqual_forPoint


  subroutine evalLogicalNotEqual_forElement( fun, varsys, elempos, time, tree, &
      &                                      nElems, nDofs, res                )
    !--------------------------------------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Position of the TreeID of the element to get the variable for in the
    !! global treeID list.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !--------------------------------------------------------------------------!
    integer :: iDep, posDepVar, iComp, iElem, lIdx, uIdx
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    !--------------------------------------------------------------------------!


    ! nInputs must be two
    allocate(input_varRes(nElems*nDofs*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_element(                        &
        &                                   varSys  = varSys,               &
        &                                   elemPos = elemPos,              &
        &                                   time    = time,                 &
        &                                   tree    = tree,                 &
        &                                   nElems  = nElems,               &
        &                                   nDofs   = nDofs,                &
        &                                   res     = input_varRes(:, iDep) )

    end do

    do iComp = 1, fun%nComponents
      do iElem = 1, nElems
        ! We only take the first mode into consideration, but store the result
        ! value for all modes.
        lIdx = (( ielem-1)* fun%ncomponents* ndofs+( 1-1)* fun%ncomponents+icomp)
        uIdx = (( ielem-1)* fun%ncomponents* ndofs+( ndofs-1)* fun%ncomponents+icomp)

        res( lIdx : uIdx ) = logicalToReal(                       &
          & input_varRes( lIdx, 1 ) .fne. input_varRes( lIdx, 2 ) )
      end do
    end do

    deallocate(input_varRes)


  end subroutine evalLogicalNotEqual_forElement


  subroutine evalLogicalNotEqual_fromIndex( fun, varSys, time, iLevel, &
      &                                     idx, idxLen, nVals, res    )
    !--------------------------------------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !--------------------------------------------------------------------------!
    integer :: iDep, posDepVar
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    type(tem_varSys_op_data_type), pointer :: fPtr
    !--------------------------------------------------------------------------!

    call C_F_POINTER( fun%method_Data, fPtr )


    ! nInputs must be two
    allocate(input_varRes(nVals*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable for indices
      !rent derive dependent variable
      call varSys%method%val(posDepVar)%get_valOfIndex( &
        & varSys  = varSys,                             &
        & time    = time,                               &
        & iLevel  = iLevel,                             &
        & idx     = fPtr%input_pntIndex(iDep)           &
        &           %indexLvl(iLevel)%val( idx(:) ),    &
        & idxLen  = idxLen,                             &
        & nVals   = nVals,                              &
        & res     = input_varRes(:,iDep)                )
    end do

    res = numFalse
    where((input_varRes(:, 1) - input_varRes(:, 2)) > eps ) res = numTrue

    deallocate(input_varRes)

  end subroutine evalLogicalNotEqual_fromIndex


end module tem_logical_operation_var_module
