! Copyright (c) 2016-2017, 2019, 2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016-2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016-2017 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016-2017, 2019, 2021 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Jana Gericke <jana.gericke@student.uni-siegen.de>
! Copyright (c) 2017 Michael Gaida  <michael.gaida@student.uni-siegen.de>
! Copyright (c) 2017 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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

! *************************************************************************** !
!> This module provides a mechanism to define new variables in the simulation
!! by applying some operator on existing variables.
!!
module tem_operation_var_module
  use, intrinsic :: iso_c_binding,  only: c_ptr, c_loc, c_f_pointer, &
    &                                     c_null_ptr

  ! include treelm modules
  use env_module,                   only: rk, long_k, labelLen
  use tem_varSys_module,            only: tem_varSys_type,               &
    &                                     tem_varSys_check_inArgs,       &
    &                                     tem_get_element_chunk,         &
    &                                     tem_get_point_chunk,           &
    &                                     tem_varSys_op_type,            &
    &                                     tem_varSys_append_derVar,      &
    &                                     tem_varSys_proc_point,         &
    &                                     tem_varSys_proc_element,       &
    &                                     tem_varSys_proc_setParams,     &
    &                                     tem_varSys_proc_getParams,     &
    &                                     tem_varSys_proc_setupIndices,  &
    &                                     tem_varSys_proc_getValOfIndex, &
    &                                     tem_varSys_solverData_evalElem_type
  use tem_variable_module,          only: tem_variable_type
  use tem_varMap_module,            only: tem_varMap_type, tem_create_varMap
  use tem_time_module,              only: tem_time_type
  use treelmesh_module,             only: treelmesh_type
  use tem_logging_module,           only: logUnit, llerror
  use tem_aux_module,               only: tem_abort
  use tem_dyn_array_module,         only: PositionOfVal
  use tem_grow_array_module,        only: append, init, &
    &                                     truncate, grw_intArray_type, &
    &                                     grw_labelArray_type
  use tem_operation_module,         only: tem_varSys_op_data_type, &
    &                                     tem_indexLvl_type
  use tem_reduction_transient_module, only: tem_reduction_transient_type,    &
    &                                       tem_reduction_transient_init,    &
    &                                       tem_reduction_transient_reset,   &
    &                                       tem_reduction_transient_update,  &
    &                                       tem_reduction_transient_getElement
  use tem_logical_operation_var_module,                                        &
    &                              only: evalLogicalAnd_forPoint,              &
    &                                    evalLogicalOr_forPoint,               &
    &                                    evalLogicalGreater_forPoint,          &
    &                                    evalLogicalGreaterOrEqual_forPoint,   &
    &                                    evalLogicalLess_forPoint,             &
    &                                    evalLogicalLessOrEqual_forPoint,      &
    &                                    evalLogicalEqual_forPoint,            &
    &                                    evalLogicalNotEqual_forPoint,         &
    &                                    evalLogicalAnd_forElement,            &
    &                                    evalLogicalOr_forElement,             &
    &                                    evalLogicalGreater_forElement,        &
    &                                    evalLogicalGreaterOrEqual_forElement, &
    &                                    evalLogicalLess_forElement,           &
    &                                    evalLogicalLessOrEqual_forElement,    &
    &                                    evalLogicalEqual_forElement,          &
    &                                    evalLogicalNotEqual_forElement,       &
    &                                    evalLogicalAnd_fromIndex,             &
    &                                    evalLogicalOr_fromIndex,              &
    &                                    evalLogicalGreater_fromIndex,         &
    &                                    evalLogicalGreaterOrEqual_fromIndex,  &
    &                                    evalLogicalLess_fromIndex,            &
    &                                    evalLogicalLessOrEqual_fromIndex,     &
    &                                    evalLogicalEqual_fromIndex,           &
    &                                    evalLogicalNotEqual_fromIndex

  implicit none

  private

  public :: tem_divideVecByScal_forPoint
  public :: tem_divideVecByScal_fromIndex
  public :: tem_varSys_append_operVar
  public :: tem_opVar_setupIndices
  public :: tem_opVar_fill_inputIndex
  public :: tem_varSys_op_data_type
  public :: tem_get_new_varSys_data_ptr
  public :: tem_free_varSys_data_ptr
  public :: tem_evalMag_forPoint, tem_evalMag_forElement, tem_evalMag_fromIndex
  public :: tem_evalAdd_forPoint, tem_evalAdd_forElement, tem_evalAdd_fromIndex
  public :: tem_evalDiff_forPoint, tem_evalDiff_forElement, &
    &       tem_evalDiff_fromIndex
  public :: tem_evalMultiply_forPoint, tem_evalMultiply_forElement, &
    &       tem_evalMultiply_fromIndex
  public :: tem_division_forPoint, tem_division_forElement, &
    &       tem_division_fromIndex
  public :: tem_opVar_getParams, tem_opVar_setParams

  public :: tem_opVar_reduction_transient_init
  public :: tem_opVar_reduction_transient_update


contains


  ! ************************************************************************** !
  !> Routine to get a pointer to a new instance of method_data for an operation
  !! variable
  function tem_get_new_varSys_data_ptr(solver_bundle) result(resPtr)
    ! ---------------------------------------------------------------------- !
    !> Pointer to the newly created instance.
    type(c_ptr) ::  resPtr
    !> Optional solver data to store in tem_varSys_op_data_type
    type(c_ptr), optional, intent(in) :: solver_bundle
    ! ---------------------------------------------------------------------- !
    !> Local variable to allocate a new instance.
    type(tem_varSys_op_data_type), pointer :: res
    ! ---------------------------------------------------------------------- !

    allocate(res)
    if (present(solver_bundle)) res%solver_bundle = solver_bundle
    resPtr = c_loc(res)

  end function tem_get_new_varSys_data_ptr
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Free a method data structure again.
  subroutine tem_free_varSys_data_ptr(vardat_ptr)
    ! ---------------------------------------------------------------------- !
    !> Data pointer to free
    type(c_ptr), intent(inout) :: vardat_ptr
    ! ---------------------------------------------------------------------- !
    type(tem_varSys_op_data_type), pointer :: vardat
    ! ---------------------------------------------------------------------- !

    call c_f_pointer(vardat_ptr, vardat)
    deallocate(vardat)
    vardat_ptr = c_null_ptr

  end subroutine tem_free_varSys_data_ptr
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> subroutine to add the variables from the input lua script to the varsys
  subroutine tem_varSys_append_operVar( operVar, varSys, pos,     &
    &                                   solverData_evalElem       )
    !--------------------------------------------------------------------------
    !> variables defined in the lua file
    type(tem_variable_type), intent(in)           :: operVar

    !> global variable system to which operVar to be appended
    type(tem_varSys_type), intent(inout)            :: varSys

    !> Position of the variable in the system.
    integer, optional, intent(out) :: pos

    !> A setter routine that allows the caller to define routine for the
    !! construction of an element representation.
    type(tem_varSys_solverData_evalElem_type), &
      &  optional, intent(in) :: solverData_evalElem
    !--------------------------------------------------------------------------
    integer :: addedPos
    integer :: nComps, nInputs
    integer, allocatable :: inPos(:)
    logical :: wasAdded
    character(len=labelLen), allocatable :: input_varname(:)
    integer, allocatable :: input_varIndex(:)
    logical :: isSatisfied
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setParams), pointer :: set_params => NULL()
    procedure(tem_varSys_proc_getParams), pointer :: get_params => NULL()
    procedure(tem_varSys_proc_setupIndices), pointer :: &
      &                                      setup_indices => NULL()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: &
      &                                       get_valOfIndex => NULL()
    type(c_ptr) ::  method_data
    type(tem_varSys_op_data_type), pointer :: opData
    !--------------------------------------------------------------------------
    nullify(get_point, get_element, set_params, get_params, setup_indices, &
      &     get_valOfIndex)

    nComps = operVar%nComponents

    nInputs = size(operVar%input_varname)
    allocate(input_varname(nInputs))
    input_varname = operVar%input_varname

    ! check if input_varnames satisfy requirements for operType
    ! and correct user defined nComps if it does not match operType
    call check_opVar_prerequisites( operType      = operVar%operType, &
      &                             nInputs       = nInputs,          &
      &                             input_varname = input_varname,    &
      &                             varSys        = varSys,           &
      &                             nComps        = nComps,           &
      &                             isSatisfied   = isSatisfied       )

    !> If not satisfied then it is not possible to append current variable
    !! to varSys
    if (.not. isSatisfied) then
      write(logUnit(1),*) 'WARNING: input varnames does not satisfy'
      write(logUnit(1),*) '         requirements for operType '&
        &               //trim(operVar%operType)
      write(logUnit(1),*) 'Variable: "'//trim(operVar%label) &
        &              //'" is not appended.'
      return
    end if

    ! for operation variables, set_params, get_params and setup_indices are same
    ! since they send information to depend variables
    set_params => tem_opVar_setParams
    get_params => tem_opVar_getParams
    setup_indices => tem_opVar_setupIndices
    ! Get method data container to store indices for getValOfIndex
    ! Overwrite this method data with solver method data if operation
    ! is solver-specific
    method_data = tem_get_new_varSys_data_ptr()
    call C_F_POINTER(method_data, opData)

    select case(trim(operVar%operType))
    ! magnitude, division, multiplication are solver specific
    ! so set using getEvalFuncionsCallback
    case( 'difference' )
      get_element => tem_evalDiff_forElement
      get_point => tem_evalDiff_forPoint
      get_valOfIndex => tem_evalDiff_fromIndex

    case( 'rel_difference' )
      get_element => evalRelDiff_forElement
      get_point => evalRelDiff_forPoint
      get_valOfIndex => evalRelDiff_fromIndex

    case( 'addition' )
      get_element => tem_evalAdd_forElement
      get_point => tem_evalAdd_forPoint
      get_valOfIndex => tem_evalAdd_fromIndex

    case ('multiplication')
      get_point => tem_evalMultiply_forPoint
      get_element => tem_evalMultiply_forElement
      get_valOfindex => tem_evalMultiply_fromIndex

    case( 'multiply_scalar_times_vector' )
      get_point => tem_multiplyScalTimesVec_forPoint
      get_element => multiplyScalTimesVec_forElement
      get_valOfindex => tem_multiplyScalTimesVec_fromIndex

    case( 'division', 'div' )
      get_point => tem_division_forPoint
      get_element => tem_division_forElement
      get_valOfindex => tem_division_fromIndex

    case( 'divide_vector_by_scalar' )
      get_point => tem_divideVecByScal_forPoint
      get_element => divideVecByScal_forElement
      get_valOfindex => tem_divideVecByScal_fromIndex

    case( 'gradient', 'grad', 'gradientX', 'gradX','gradientY','gradY', &
      &   'gradientZ', 'gradZ'                                          )
    ! Pointers set by the solvers using opVar_setter callback, see below

    case( 'magnitude' )
      get_point => tem_evalMag_forPoint
      get_element => tem_evalMag_forElement
      get_valOfindex => tem_evalMag_fromIndex

    case( 'extract' )
      allocate( input_varIndex(size(operVar%input_varIndex)) )
      input_varIndex = operVar%input_varIndex

      get_element => extract_forElement
      get_point => extract_forPoint
      get_valOfIndex => extract_fromIndex

    case( 'combine' )
      get_element => combine_forElement
      get_point => combine_forPoint
      get_valOfIndex => combine_fromIndex

    case( 'greater_than', 'gt', '>' )
      get_point => evalLogicalGreater_forPoint
      get_element => evalLogicalGreater_forElement
      get_valOfIndex => evalLogicalGreater_fromIndex

    case( 'greater_than_or_equal', 'ge', '>=' )
      get_point => evalLogicalGreaterOrEqual_forPoint
      get_element => evalLogicalGreaterOrEqual_forElement
      get_valOfIndex => evalLogicalGreaterOrEqual_fromIndex

    case( 'less_than', 'lt', '<' )
      get_point => evalLogicalLess_forPoint
      get_element => evalLogicalLess_forElement
      get_valOfIndex => evalLogicalLess_fromIndex

    case( 'less_than_or_equal', 'le', '<=' )
      get_point => evalLogicalLessOrEqual_forPoint
      get_Element => evalLogicalLessOrEqual_forElement
      get_valOfIndex => evalLogicalLessOrEqual_fromIndex

    case( 'equal', 'eq', '=' )
      get_point => evalLogicalEqual_forPoint
      get_Element => evalLogicalEqual_forElement
      get_valOfIndex => evalLogicalEqual_fromIndex

    case( 'not_equal', 'ne', '/=' )
      get_point => evalLogicalNotEqual_forPoint
      get_Element => evalLogicalNotEqual_forElement
      get_valOfIndex => evalLogicalNotEqual_fromIndex

    case( 'and' )
      get_point => evalLogicalAnd_forPoint
      get_Element => evalLogicalAnd_forElement
      get_valOfIndex => evalLogicalAnd_fromIndex

    case( 'or' )
      get_point => evalLogicalOr_forPoint
      get_Element => evalLogicalOr_forElement
      get_valOfIndex => evalLogicalOr_fromIndex

    case('reduction_transient')
      opData%redTrans%config = operVar%redTransConfig
      get_point => reductionTransient_forPoint
      get_Element => reductionTransient_forElement
      get_valOfIndex => reductionTransient_fromIndex

    case default
      if (.not. (associated(get_point)      &
        &   .or. associated(get_element)    &
        &   .or. associated(set_params)     &
        &   .or. associated(get_params)     &
        &   .or. associated(setup_indices)  &
        &   .or. associated(get_valOfIndex))) then
        write(logUnit(3),*) 'operType: '                 &
          & // trim(operVar%operType)                    &
          & // ' not supported. Variable is not appended.'
        return ! go to next variable
      end if
    end select

    ! Workaround for Intel 15 compiler
    if ( .not. allocated(input_varname) ) then
      allocate( input_varname(0) )
    end if

    if ( .not. allocated(input_varIndex) ) then
      allocate( input_varIndex(0) )
    end if

    ! append variable to varSys
    call tem_varSys_append_derVar(          &
      &  me             = varSys,           &
      &  varName        = operVar%label,    &
      &  operType       = operVar%operType, &
      &  nComponents    = nComps,           &
      &  input_varname  = input_varname,    &
      &  input_varIndex = input_varIndex,   &
      &  method_data    = method_data,      &
      &  get_point      = get_point,        &
      &  get_element    = get_element,      &
      &  set_params     = set_params,       &
      &  get_params     = get_params,       &
      &  setup_indices  = setup_indices,    &
      &  get_valOfIndex = get_valOfIndex,   &
      &  pos            = addedPos,         &
      &  wasAdded       = wasAdded          )

    if (wasAdded) then
      if (present(solverData_evalElem)) then
        ! If an solverData_evalElem operation is provided, override the
        ! get_element pointer and use the provided setter solverData_evalElem
        ! instead to define the get_element routine.
        call solverData_evalElem%opVar_setter(varSys%method%val(addedPos))
      end if
      write(logUnit(9),*) 'Successfully appended variable "' &
        & // trim(operVar%label)// '" to the variable system'
    else if (addedpos < 1) then
      write(logUnit(7),*) 'WARNING: variable '//trim(operVar%label)// &
        &                 ' is not added to variable system'
    end if

    if (present(pos)) pos = addedPos

    ! deallocate here to be allocated for next variable
    deallocate(input_varname)
    if (allocated(input_varIndex)) deallocate(input_varIndex)
    if (allocated(inPos)) deallocate(inPos)

  end subroutine tem_varSys_append_operVar
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This subroutine checks whether input variables satisfy requirements for
  !! opertype.
  !!
  !! For example: nComponents, number of inputs etc.
  !! If input_varname not found in varSys then this function returns false.
  !! If user defined nComps = -1 then nComps is set according to operType.
  !!
  !! This subroutine checks operations used in both treelm and solvers
  subroutine check_opVar_prerequisites( operType, nInputs, input_varname, &
    &                                   varSys, nComps, isSatisfied )
    !--------------------------------------------------------------------------!
    !> Operation type
    character(len=*), intent(in) :: operType
    !> Number of inputs
    integer, intent(in) :: nInputs
    !> Input varnames for current operation
    character(len=*), intent(in) :: input_varname(nInputs)
    !> Variable system to look for input_varname
    type(tem_varSys_type), intent(in) :: varSys
    !> Number of components defined for operation variable.
    !! If nComps == -1 then current nComps is set here
    integer, intent(inout) :: nComps
    !> true if all requirements for opertype are satisfied
    logical, intent(out) :: isSatisfied
    !--------------------------------------------------------------------------!
    integer :: iIn, total_input_nComps
    integer :: inpos(nInputs), input_nComps(nInputs)
    !--------------------------------------------------------------------------!
    write(logUnit(7),*) 'Checking prerequisites for opertype: '//trim(operType)
    isSatisfied = .true.

    ! Position of input variables in varSys
    do iIn = 1, nInputs
      inPos(iIn) = PositionofVal(varSys%varname, input_varname(iIn))
    end do

    if ( all(inPos > 0) ) then
      input_nComps(:) = varSys%method%val(inPos(:))%nComponents
    else
      isSatisfied = .false.
      write(logUnit(1),*) 'Error: input varnames not found in varsys'
      return
    end if

    select case (trim(operType))
    case ('addition', 'difference', 'rel_difference', 'multiplication', &
      &   'division')
      ! operations which require two inputs and both to have same number of
      ! components
      if (size(input_varname) /= 2) then
        write(logUnit(1),*) 'Error: In operation type: '//trim(operType)
        write(logUnit(1),*) 'Number of input_varname /= 2'
        call tem_abort()
      end if

      if (nComps == -1) then
        nComps = input_nComps(1)
        write(logUnit(7),*) 'INFO: nComponents is not defined by user.'
        write(logUnit(7),*)  '     nComps set to nComps of 1st input variable'
      end if

      if (.not. all(input_nComps == nComps)) then
        write(logUnit(1),*) 'Error: nComps of operation variable does not ' &
          &               //'match with input variables nComps'
        write(logUnit(1),*) 'Input nComps: ', input_nComps
        call tem_abort()
      end if

    case ('multiply_scalar_times_vector')
      ! operations which require two inputs and 1st input must be a scalar
      if (size(input_varname) /= 2) then
        write(logUnit(1),*) 'Error: In operation type: '//trim(operType)
        write(logUnit(1),*) 'Number of input_varname /= 2'
        call tem_abort()
      end if

      ! 1st input variable must be scalar
      if (input_nComps(1) /= 1) then
        write(logUnit(1),*) 'Error: In operation type: '//trim(operType)
        write(logUnit(1),*) '1st input variable is not a scalar'
        call tem_abort()
      end if

      if (nComps == -1 .or. nComps /= input_nComps(2)) then
        write(logUnit(1),*) 'Warning: nComps for operVar is wrong.'
        write(logUnit(1),*) '         Setting nComps =', input_nComps(2)
        nComps = input_nComps(2)
      end if

    case ('divide_vector_by_scalar')
      ! operations which require two inputs and 1st input must be a scalar
      if (size(input_varname) /= 2) then
        write(logUnit(1),*) 'Error: In operation type: '//trim(operType)
        write(logUnit(1),*) 'Number of input_varname /= 2'
        call tem_abort()
      end if

      ! 2nd input variable must be scalar
      if (input_nComps(2) /= 1) then
        write(logUnit(1),*) 'Error: In operation type: '//trim(operType)
        write(logUnit(1),*) '1st input variable is not a scalar'
        call tem_abort()
      end if

      if (nComps == -1 .or. nComps /= input_nComps(1)) then
        write(logUnit(1),*) 'Warning: nComps for operVar is wrong.'
        write(logUnit(1),*) '         Setting nComps =', input_nComps(1)
        nComps = input_nComps(1)
      end if

    case ('magnitude')
      ! operations which require only one input and ncomponents must be 1
      if (nInputs /= 1 ) then
        write(logUnit(1),*) 'Error: In operation type: '//trim(operType)
        write(logUnit(1),*) 'Number of input_varname /= 1'
        call tem_abort()
      end if

      if ( nComps /= 1 ) then
        write(logUnit(1),*) 'Warning: nComps /= 1. Setting nComps = 1'
        nComps = 1
      end if

    case ('meansquare')
      ! operations which require only one input and ncomps == input_ncomps
      if (nInputs /= 1 ) then
        write(logUnit(1),*) 'Error: In operation type: '//trim(operType)
        write(logUnit(1),*) 'Number of input_varname /= 1'
        call tem_abort()
      end if

      if ( nComps /= input_nComps(1) ) then
        write(logUnit(1),*) 'Warning: nComps /= input, should be:', &
          &                 input_nComps(1)
        write(logUnit(1),*) 'but is:', nComps
        write(logUnit(1),*) 'Setting it to ', input_nComps(1)
        nComps = input_nComps(1)
      end if

    case ('locall2mean','deviation')
      ! operations which require only one input and ncomps == input_ncomps
      if (nInputs /= 1 ) then
        write(logUnit(1),*) 'Error: In operation type: '//trim(operType)
        write(logUnit(1),*) 'Number of input_varname /= 1'
        call tem_abort()
      end if

      if ( nComps /= input_nComps(1) ) then
        write(logUnit(1),*) 'Warning: nComps /= input, should be:', &
          &                 input_nComps(1)
        write(logUnit(1),*) 'but is:', nComps
        write(logUnit(1),*) 'Setting it to ', input_nComps(1)
        nComps = input_nComps(1)
      end if


    case ('extract', 'gradient', 'gradientX', 'gradientY', 'gradientZ', &
      &   'reduction_transient'                                         )
      ! operations which require only one input and ncomponents can be >= 1
      if (nInputs /= 1 ) then
        write(logUnit(1),*) 'Error: In operation type: '//trim(operType)
        write(logUnit(1),*) 'Number of input_varname /= 1'
        call tem_abort()
      end if

      if (nComps == -1) then
        write(logUnit(1),*) 'Error: nComps for operVar is not defined.'
        call tem_abort()
      end if

    case ('combine')
      ! operations which require more than two inputs and
      ! nComps is sum of nComps of input vars
      total_input_nComps = sum(input_nComps)
      if (nComps /= total_input_nComps) then
        write(logUnit(llerror),*) 'Warning: In variable operation combine'
        write(logUnit(llerror),*) '  User defined nComps ' &
          &                    // '/= sum(input_nComponents).'
        write(logUnit(llerror),*) &
          &  '  So, setting nComps = sum(input_nComponents)'
        write(logUnit(llerror),*) '  nComps: ', total_input_nComps
        nComps = total_input_nComps
      end if

    case( 'greater_than', 'gt', '>',         &
      & 'greater_than_or_equal', 'ge', '>=', &
      & 'less_than', 'lt', '<',              &
      & 'less_than_or_equal', 'le', '<=',    &
      & 'equal', 'eq', '=',                  &
      & 'not_equal', 'ne', '/=',             &
      & 'and',                               &
      & 'or'                                 )
      if (size(input_varname) /= 2) then
        write(logUnit(1),*) 'Error: In operation type: '//trim(operType)
        write(logUnit(1),*) 'Number of input_varname /= 2'
        call tem_abort()
      end if

    case default
      write(logUnit(1),*) 'ERROR: operType: ' // trim(operType) &
        & // ' not supported. Variable is not appended.'
      call tem_abort()
    end select

  end subroutine check_opVar_prerequisites
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Initialize time reduction operation variable
  !! Loop over all variable in varSys and allocate redTrans%val for
  !! reduction_transient operation variable with nElems
  subroutine tem_opVar_reduction_transient_init(varSys, tree, redTransVarMap,&
    &                                           nDofs, time)
    ! -------------------------------------------------------------------------
    !> Global variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> treelmesh_type
    type(treelmesh_type), intent(in) :: tree
    !> position of time reduction variable in varSys
    type(tem_varMap_type), intent(out) :: redTransVarMap
    !> Solver nDegrees of freedom
    integer, intent(in), optional :: nDofs
    !> Current time
    type(tem_time_type), intent(in) :: time
    ! -------------------------------------------------------------------------
    type(tem_varSys_op_data_type), pointer :: opData
    integer :: iVar, iElem, varPos, nDofs_loc, posDepVar, nCompMax, idxMax
    integer :: nRedVars
    type(grw_labelArray_type) :: redTransVarName
    integer :: elemPos(tree%nElems)
    real(kind=rk), allocatable :: input_varRes(:)
    ! -------------------------------------------------------------------------
    if (present(nDofs)) then
      nDofs_loc = nDofs
    else
      nDofs_loc = 1
    end if

    ! Gather list of variable names which has reduction_transient operation
    call init(redTransVarName)
    do iVar = 1, varSys%varName%nVals
      if (trim(varSys%method%val(iVar)%operType) == 'reduction_transient') then
        call append(me=redTransVarName, val=varSys%varName%val(iVar))
      end if
    end do

    ! create varMap to store position of reduction_transient variable in varSys
    call tem_create_varMap(                                                 &
      &             varName = redTransVarName%val(1:redTransVarName%nVals), &
      &             varSys  = varSys,                                       &
      &             varMap  = redTransVarMap                                )

    elemPos(1:tree%nElems) = (/ (iElem, iElem=1, tree%nElems) /)
    nRedVars =  redTransVarMap%varPos%nVals
    nCompMax = maxval(varSys%method%val(redTransVarMap%varPos%val(1:nRedVars)) &
      &                            %nComponents)
    allocate(input_varRes(tree%nElems*nCompMax*nDofs_loc))

    ! Initialize time reduction
    do iVar = 1, redTransVarMap%varPos%nVals
      varPos = redTransVarmap%varPos%val(iVar)
      call C_F_POINTER(varSys%method%val(varPos)%method_data, opData)

      call tem_reduction_transient_init(                            &
        &      me          = opData%redTrans,                       &
        &      nElems      = tree%nElems,                           &
        &      nComponents = varSys%method%val(varPos)%nComponents, &
        &      nDofs       = nDofs_loc                              )
      ! Fill last
      posDepVar = varSys%method%val(varPos)%input_varPos(1)
      call varSys%method%val(posDepVar)%get_element(                        &
        &                                   varSys  = varSys,               &
        &                                   elemPos = elemPos,              &
        &                                   time    = time,                 &
        &                                   tree    = tree,                 &
        &                                   nElems  = tree%nElems,          &
        &                                   nDofs   = nDofs_loc,            &
        &                                   res     = input_varRes(:)       )

      idxMax = opData%redTrans%nEntries
      opData%redTrans%val(:, opData%redTrans%last) = input_varRes(1:idxMax)
    end do

  end subroutine tem_opVar_reduction_transient_init
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Update all time reduction operation variables for entire domain
  subroutine tem_opVar_reduction_transient_update(redTransVarPos, varSys, tree,&
    &                                             time)
    ! ---------------------------------------------------------------------- !
    !> Position of time reduction variables in varSys
    integer, intent(in) :: redTransVarPos(:)
    !> Global Variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> treelmesh_type
    type(treelmesh_type), intent(in) :: tree
    !> Current time
    type(tem_time_type), intent(in) :: time
    ! ---------------------------------------------------------------------- !
    integer :: elemPos(tree%nElems)
    real(kind=rk), allocatable :: input_varRes(:)
    integer :: iVar, iElem, posDepVar, nCompMax, idxMax, nDofs
    type(tem_varSys_op_data_type), pointer :: opData
    ! ---------------------------------------------------------------------- !

    ! Only need to do anything here if at least one variable with reduction
    ! is to be computed.
    if (size(redTransVarPos) < 1) RETURN

    ! nDofs of solver are stored in opData%redTrans which is same for all
    ! variables so get from any reduction_transient variable
    call C_F_Pointer(varSys%method%val(redTransVarPos(1))%method_data, &
      &              opData                                            )
    nDofs = opData%redTrans%nDofs

    elemPos(1:tree%nElems) = (/ (iElem, iElem=1, tree%nElems) /)

    nCompMax = maxval(varSys%method%val(redTransVarPos(:))%nComponents)
    allocate(input_varRes(tree%nElems*nCompMax*nDofs))

    do iVar = 1, size(redTransVarPos)
      call C_F_Pointer(varSys%method%val(redTransVarPos(iVar))%method_data, &
        &              opData                                               )

      posDepVar = varSys%method%val(redTransVarPos(iVar))%input_varPos(1)
      call varSys%method%val(posDepVar)%get_element(                        &
        &                                   varSys  = varSys,               &
        &                                   elemPos = elemPos,              &
        &                                   time    = time,                 &
        &                                   tree    = tree,                 &
        &                                   nElems  = tree%nElems,          &
        &                                   nDofs   = nDofs,                &
        &                                   res     = input_varRes(:)       )

      idxMax = opData%redTrans%nEntries

      call tem_reduction_transient_update(me     = opData%redTrans,       &
        &                                 res    = input_varRes(1:idxMax) )
    end do
    deallocate(input_varRes)
  end subroutine tem_opVar_reduction_transient_update
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Evaluate the function pointers of the dependent variables,
  !! and then calculate the difference between these two. ( scalar or vector )
  !! In lua file, first define new variable with varType operation kind as
  !! "difference" and provide two dependent variable via input_varname.
  !! If input_varname variable is not part of predefined solver variables then
  !! add also that variable via variable table for example spacetime function
  !! variable.
  !! For example: Define a variable called difference, which depend on density
  !! and spacetime. one can get an error between simulation
  !! results and analytical solution.
  !!
  !! \verbatim
  !! -- in lua file, one can define as following:
  !! variable = {{
  !!   name = 'dens_reference',
  !!   ncomponents = 1,
  !!   vartype = "st_fun",
  !!   st_fun =  luaFun
  !!   },
  !!   {
  !!   name = 'dens_difference',
  !!   ncomponents = 1,
  !!   vartype = "operation",
  !!   operation = {kind='difference', input_varname={density, dens_reference}}
  !!   }
  !!   ...
  !!  }
  !! tracking = {
  !!   variable = {'dens_difference'},
  !!   folder = 'tracking/',
  !!   shape = {kind = 'canoND', object = {origin = {3.0,3.1,3.0} } },
  !!   format = 'ascii',
  !!   time = {min = 0, max = tmax, interval = 1},
  !!  }
  !! \endverbatim
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_element.
  !!
  recursive subroutine tem_evalDiff_forElement( fun, varsys, elempos, time,   &
    &                                           tree, nElems, nDofs, res                )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    integer :: iDep, posDepVar
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    ! ---------------------------------------------------------------------- !
    integer :: nTotal
    ! ---------------------------------------------------------------------- !

    nTotal = nElems*nDofs*fun%nComponents

    ! nInputs must be two
    allocate(input_varRes(nTotal, fun%nInputs))

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

    res(:nTotal) = input_varRes(:, 1) - input_varRes(:, 2)

    deallocate(input_varRes)

  end subroutine tem_evalDiff_forElement
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Same as evalDiff_forElement except it evaluate diff from points
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_point
  recursive subroutine tem_evalDiff_forPoint( fun, varsys, point, time, tree, &
    &                                         nPnts, res                      )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    integer :: iDep, posDepVar
    integer :: nTotal
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    ! ---------------------------------------------------------------------- !

    nTotal = nPnts*fun%nComponents

    ! nInputs must be two
    allocate(input_varRes(nPnts*fun%nComponents, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_point(                         &
        &                                   varSys = varSys,               &
        &                                   point  = point,                &
        &                                   time   = time,                 &
        &                                   tree   = tree,                 &
        &                                   nPnts  = nPnts,                &
        &                                   res    = input_varRes(:, iDep) )

    end do

    res(:nTotal) = input_varRes(:, 1) - input_varRes(:, 2)

    deallocate(input_varRes)

  end subroutine tem_evalDiff_forPoint
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Same as evalDiff_forPoint except it evaluate diff from points via index
  !! which are setup before hand
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_point
  recursive subroutine tem_evalDiff_fromIndex( fun, varSys, time, iLevel, idx, &
    &                                          idxLen, nVals, res              )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    integer :: iDep, posDepVar
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    ! ---------------------------------------------------------------------- !
    type(tem_varSys_op_data_type), pointer :: opData
    integer :: nTotal
    ! ---------------------------------------------------------------------- !

    call C_F_POINTER( fun%method_Data, opData )

    ! check if number of index are the same as number of values asked for
    call tem_varSys_check_inArgs( fun, varSys, time, iLevel, idx, idxLen, &
      &    nVals, label = 'tem_evalDiff_fromIndex'                        )

    nTotal = nVals*fun%nComponents
    ! nInputs must be two
    allocate(input_varRes(nTotal, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable for indices
      !rent derive dependent variable
      call varSys%method%val(posDepVar)%get_valOfIndex( &
          & varSys  = varSys,                           &
          & time    = time,                             &
          & iLevel  = iLevel,                           &
          & idx     = opData%input_pntIndex(iDep)       &
          &           %indexLvl(iLevel)%val( idx(:) ),  &
          & nVals   = nVals,                            &
          & res     = input_varRes(:,iDep)              )
    end do

    res(:nTotal) = input_varRes(:, 1) - input_varRes(:, 2)

    deallocate(input_varRes)

  end subroutine tem_evalDiff_fromIndex
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Same as evalDiff but output of 2nd dependent variable is used to compute
  !! relative difference.
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_element.
  !!
  recursive subroutine evalRelDiff_forElement( fun, varsys, elempos, time, &
    &                                          tree, nElems, nDofs, res    )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    integer :: iDep, posDepVar
    integer :: nTotal
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    ! ---------------------------------------------------------------------- !

    nTotal = nElems*nDofs*fun%nComponents
    ! nInputs must be two
    allocate(input_varRes(nTotal, fun%nInputs))

    if (nDofs > 1) then
      write(*,*) 'TODO: evalreldiff does not work for polynomial data yet'
      write(*,*) '      It makes use of a division, which can not directly'
      write(*,*) '      be done in modal space!'
      write(*,*) ''
      write(*,*) 'Need to replace this routine in Ateles!'
      write(*,*) 'Stopping now.'
      call tem_abort()
    end if

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

    res(:nTotal) = ( input_varRes(:, 1) - input_varRes(:, 2) ) &
      &          / input_varRes(:, 2)

    deallocate(input_varRes)

  end subroutine evalRelDiff_forElement
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Same as evalRelDiff_forElement except it evaluate relatice diff
  !! from points
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_point.
  recursive subroutine evalRelDiff_forPoint( fun, varsys, point, time, tree, &
    &                                        nPnts, res                      )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    integer :: iDep, posDepVar
    integer :: nTotal
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    ! ---------------------------------------------------------------------- !

    nTotal = nPnts*fun%nComponents
    ! nInputs must be two
    allocate(input_varRes(nTotal, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_point(                         &
        &                                   varSys = varSys,               &
        &                                   point  = point,                &
        &                                   time   = time,                 &
        &                                   tree   = tree,                 &
        &                                   nPnts  = nPnts,                &
        &                                   res    = input_varRes(:, iDep) )

    end do

    res(:nTotal) = ( input_varRes(:, 1) - input_varRes(:, 2) ) &
      &            / input_varRes(:, 2)

    deallocate(input_varRes)

  end subroutine evalRelDiff_forPoint
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Same as evalRelDiff_forElement except it evaluate relative diff
  !! from points via index which are setup before
  !!
  !!
  recursive subroutine evalRelDiff_fromIndex( fun, varSys, time, iLevel, idx, &
    &                                         idxLen, nVals, res              )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    integer :: iDep, posDepVar
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk), allocatable :: input_varRes(:,:)
    ! ---------------------------------------------------------------------- !
    type(tem_varSys_op_data_type), pointer :: opData
    integer :: nTotal
    ! ---------------------------------------------------------------------- !

    call C_F_POINTER( fun%method_Data, opData )

    nTotal = nVals*fun%nComponents

    call tem_varSys_check_inArgs( fun, varSys, time, iLevel, idx, idxLen, &
      &    nVals, label = 'evalRelDiff_fromIndex'                         )

    ! nInputs must be two
    allocate(input_varRes(nTotal, fun%nInputs))

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable for indices
      call varSys%method%val(posDepVar)%get_valOfIndex( &
          & varSys  = varSys,                           &
          & time    = time,                             &
          & iLevel  = iLevel,                           &
          & idx     = opData%input_pntIndex(iDep)       &
          &           %indexLvl(iLevel)%val( idx(:) ),  &
          & nVals   = nVals,                            &
          & res     = input_varRes(:,iDep)              )
    end do

    res(:nTotal) = ( input_varRes(:, 1) - input_varRes(:, 2) ) &
      &           / input_varRes(:, 2)

    deallocate(input_varRes)

  end subroutine evalRelDiff_fromIndex
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Evaluate the function pointers of the dependent variables,
  !! and then calculate the addition or arbinary number of input variables
  !! with all input variables have same nComponents. ( scalar or vector )
  !! In lua file, first define new variable with varType operation kind as
  !! "addition" and provide two dependent variable via input_varname.
  !! If input_varname variable is not part of predefined solver variables then
  !! add also that variable via variable table for example spacetime function
  !! variable.
  !!
  !! \verbatim
  !! -- in lua file, one can define as following:
  !! variable = {{
  !!   name = 'dens_reference',
  !!   ncomponents = 1,
  !!   vartype = "st_fun",
  !!   st_fun =  luaFun
  !!   },
  !!   {
  !!   name = 'dens_difference',
  !!   ncomponents = 1,
  !!   vartype = "operation",
  !!   operation = {kind='addition', input_varname={density, dens_reference}}
  !!   }
  !!   ...
  !!  }
  !! \endverbatim
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_element.
  !!
  recursive subroutine tem_evalAdd_forElement( fun, varsys, elempos, time, &
    &                                          tree, nElems, nDofs, res    )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    integer :: iDep, posDepVar
    real(kind=rk), allocatable :: input_varRes(:)
    integer :: nTotal
    ! ---------------------------------------------------------------------- !

    nTotal = nElems*nDofs*fun%nComponents
    allocate(input_varRes(nTotal))

    res = 0.0_rk
    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_element(               &
        &                                   varSys  = varSys,      &
        &                                   elemPos = elemPos,     &
        &                                   time    = time,        &
        &                                   tree    = tree,        &
        &                                   nElems  = nElems,      &
        &                                   nDofs   = nDofs,       &
        &                                   res     = input_varRes )

      res(:nTotal) = res(:nTotal) + input_varRes
    end do

    deallocate(input_varRes)

  end subroutine tem_evalAdd_forElement
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Same as tem_evalAdd_forElement except it evaluate addition from points
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_point
  recursive subroutine tem_evalAdd_forPoint( fun, varsys, point, time, tree, &
    &                                        nPnts, res                      )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    integer :: iDep, posDepVar
    real(kind=rk) :: input_varRes(nPnts*fun%nComponents)
    ! ---------------------------------------------------------------------- !
    res = 0.0_rk
    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_point( varSys = varSys,      &
        &                                          point  = point,       &
        &                                          time   = time,        &
        &                                          tree   = tree,        &
        &                                          nPnts  = nPnts,       &
        &                                          res    = input_varRes )

      res(:nPnts*fun%nComponents) = res(:nPnts*fun%nComponents) &
        &                         + input_varRes
    end do

  end subroutine tem_evalAdd_forPoint
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Same as tem_evalAdd_forPoint except it evaluate addition from points via
  !! index which are setup before hand
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_point
  recursive subroutine tem_evalAdd_fromIndex( fun, varSys, time, iLevel, idx, &
    &                                         idxLen, nVals, res              )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    integer :: iDep, posDepVar
    ! first dimension:nElems,
    ! second dimension:size of difference variable in input_varname
    real(kind=rk) :: input_varRes(nVals * fun%nComponents)
    ! ---------------------------------------------------------------------- !
    type(tem_varSys_op_data_type), pointer :: opData
    ! ---------------------------------------------------------------------- !

    call C_F_POINTER( fun%method_Data, opData )

    call tem_varSys_check_inArgs( fun, varSys, time, iLevel, idx, idxLen, &
      &    nVals, label = 'evalAdd_fromIndex'                             )

    res = 0.0_rk
    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable for indices
      !rent derive dependent variable
      call varSys%method%val(posDepVar)%get_valOfIndex( &
          & varSys  = varSys,                           &
          & time    = time,                             &
          & iLevel  = iLevel,                           &
          & idx     = opData%input_pntIndex(iDep)       &
          &           %indexLvl(iLevel)%val( idx(:) ),  &
          & nVals   = nVals,                            &
          & res     = input_varRes                      )
      res(:nVals*fun%nComponents) = res(:nVals*fun%nComponents) + input_varRes
    end do

  end subroutine tem_evalAdd_fromIndex
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Evaluate magnitude of any vectorial variable.
  !! In lua file, first define new variable with varType operation kind as
  !! "magnitude" and provide name of the variable from which magnitude
  !! to be derived in input_varname.
  !! If input_varname variable is not part of predefined solver variables then
  !! add also that variable via variable table.
  !!
  !! \verbatim
  !! -- in lua file, one can define as following:
  !! variable = {{
  !!   name = 'velMag',
  !!   ncomponents = 1,
  !!   vartype = "operation",
  !!   operation = {kind='magnitude',input_varname={'velocity'}}
  !!   },
  !!  }
  !! tracking = {
  !!   variable = {'velMag'},
  !!   folder = 'tracking/',
  !!   shape = {kind = 'canoND', object = {origin = {3.0,3.1,3.0} } },
  !!   format = 'ascii',
  !!   time = {min = 0, max = tmax, interval = 1},
  !!  }
  !! \endverbatim
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_element.
  !!
  recursive subroutine tem_evalMag_forElement( fun, varsys, elempos, time, &
    &                                          tree, nElems, nDofs, res    )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    integer :: iElem, depVar_pos, depVar_nComps, offset, iDof
    real(kind=rk), allocatable :: input_varRes(:)
    ! ---------------------------------------------------------------------- !

    ! get the position of dependent variable, assuming only one
    depVar_pos = fun%input_varPos(1)

    ! get the nComp of dependent variable
    depVar_nComps = varSys%method%val( depVar_pos )%nComponents

    ! allocate array to save results
    allocate(input_varRes( nElems * nDofs * depVar_nComps ))

    ! derive dependent variable
    call varSys%method%val(depVar_pos)%get_element(                 &
      &                                   varSys  = varSys,         &
      &                                   elemPos = elemPos,        &
      &                                   time    = time,           &
      &                                   tree    = tree,           &
      &                                   nElems  = nElems,         &
      &                                   nDofs   = nDofs,          &
      &                                   res     = input_varRes(:) )

    ! compute magnitude
    ! assuming fun%nComponents = 1
    input_varRes = input_varRes * input_varRes
    do iElem = 1, nElems
      do iDof = 1, nDofs
        offset = (( ielem-1)* depvar_ncomps* ndofs+( idof-1)* depvar_ncomps+0)
        res((( ielem-1)* 1* ndofs+( idof-1)* 1+1) )                    &
          & = sqrt(sum(input_varRes(offset+1 : offset+depVar_nComps)))
      end do
    end do

    deallocate(input_varRes)

  end subroutine tem_evalMag_forElement
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Same as evalMag_forElement except it evaluate magnitude from points
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_element.
  recursive subroutine tem_evalMag_forPoint( fun, varsys, point, time, tree, &
    &                                        nPnts, res                      )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    integer :: iPnt, depVar_pos, depVar_nComps, offset
    real(kind=rk), allocatable :: input_varRes(:)
    ! ---------------------------------------------------------------------- !

    ! get the position of dependent variable, assuming only one
    depVar_pos = fun%input_varPos(1)

    ! get the nComp of dependent variable
    depVar_nComps = varSys%method%val( depVar_pos )%nComponents

    ! allocate array to save results
    allocate(input_varRes( nPnts * depVar_nComps ))

    ! derive dependent variable
    call varSys%method%val(depVar_pos)%get_point(                  &
      &                                   varSys = varSys,         &
      &                                   point  = point,          &
      &                                   time   = time,           &
      &                                   tree   = tree,           &
      &                                   nPnts  = nPnts,          &
      &                                   res    = input_varRes(:) )

    ! compute magnitude
    ! assuming fun%nComponents = 1
    input_varRes = input_varRes * input_varRes
    do iPnt = 1, nPnts
      offset = (( ipnt-1)* depvar_ncomps+0)
      res(iPnt) = sqrt(sum(input_varRes(offset+1 : offset+depVar_nComps)))
    end do

    deallocate(input_varRes)

  end subroutine tem_evalMag_forPoint
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Same as evalMag_forPoint except it evaluate magnitude from points via
  !! indices which need to be setup before
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_element.
  recursive subroutine tem_evalMag_fromIndex( fun, varSys, time, iLevel, idx, &
    &                                         idxLen, nVals, res              )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    integer :: iVal, depVar_pos, depVar_nComps, offset
    real(kind=rk), allocatable :: input_varRes(:)
    ! ---------------------------------------------------------------------- !
    type(tem_varSys_op_data_type), pointer :: opData
    ! ---------------------------------------------------------------------- !

    call C_F_POINTER( fun%method_Data, opData )

    call tem_varSys_check_inArgs( fun, varSys, time, iLevel, idx, idxLen, &
      &    nVals, label = 'evalMag_fromIndex'                             )

    ! get the position of dependent variable, assuming only one
    depVar_pos = fun%input_varPos(1)

    ! get the nComp of dependent variable
    depVar_nComps = varSys%method%val( depVar_pos )%nComponents

    ! allocate array to save results
    allocate(input_varRes( nVals * depVar_nComps ))

    ! derive dependent variable
    call varSys%method%val(depVar_pos)%get_valOfIndex(            &
      &                varSys  = varSys,                          &
      &                time    = time,                            &
      &                iLevel  = iLevel,                          &
      &                idx     = opData%input_pntIndex(1)         &
      &                          %indexLvl(iLevel)%val( idx(:) ), &
      &                nVals   = nVals,                           &
      &                res     = input_varRes(:)                  )

    ! compute magnitude
    ! assuming fun%nComponents = 1
    input_varRes = input_varRes * input_varRes
    do iVal = 1, nVals
      offset = (( ival-1)* depvar_ncomps+0)
      res(iVal) = sqrt(sum(input_varRes(offset+1 : offset+depVar_nComps)))
    end do

    deallocate(input_varRes)

  end subroutine tem_evalMag_fromIndex
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Routine to multiply variables if all variables have same number of
  !! components.
  !!
  !! \verbatim
  !! -- in lua file, one can define as following:
  !! variable = {{
  !!   name = 'coeff',
  !!   ncomponents = 1,
  !!   vartype = "st_fun",
  !!   st_fun =  0.25
  !!   },
  !!   {
  !!   name = 'newVel',
  !!   ncomponents = 1,
  !!   vartype = "operation",
  !!   operation = {kind='multiplication',
  !!                input_varname={coeff, vel_mag}}
  !!   }
  !!   ...
  !!  }
  !! \endverbatim
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_element.
  !!
  recursive subroutine tem_evalMultiply_forElement( fun, varsys, elempos,    &
    &                                               time, tree, nElems,      &
    &                                               nDofs, res               )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    integer :: iDep, posDepVar
    integer :: nTotal
    real(kind=rk), allocatable :: input_varRes(:)
    ! ---------------------------------------------------------------------- !

    if (nDofs > 1) then
      write(*,*) 'TODO: evalmultiply does not work for polynomial data yet'
      write(*,*) '      It makes use of a multiplication, which can not'
      write(*,*) '      directly be done that simple in modal space!'
      write(*,*) ''
      write(*,*) 'Need to replace this routine in Ateles!'
      write(*,*) 'Stopping now.'
      call tem_abort()
    end if

    nTotal = nElems*nDofs*fun%nComponents
    allocate( input_varRes(nTotal) )

    res = 1.0_rk

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_element(               &
        &                                   varSys  = varSys,      &
        &                                   elemPos = elemPos,     &
        &                                   time    = time,        &
        &                                   tree    = tree,        &
        &                                   nElems  = nElems,      &
        &                                   nDofs   = nDofs,       &
        &                                   res     = input_varRes )

      res(:nTotal) = res(:nTotal) * input_varRes

    end do

    deallocate(input_varRes)

  end subroutine tem_evalMultiply_forElement
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Same as tem_evalMultiply_forElement except it evaluate it multiply values
  !! from points
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_point.
  recursive subroutine tem_evalMultiply_forPoint( fun, varsys, point, time,   &
    &                                             tree, nPnts, res            )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    integer :: iDep, posDepVar
    integer :: nTotal
    real(kind=rk), allocatable :: input_varRes(:)
    ! ---------------------------------------------------------------------- !

    nTotal = nPnts*fun%nComponents
    allocate( input_varRes(nTotal) )

    res = 1.0_rk

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_point(                &
        &                                   varSys = varSys,      &
        &                                   point  = point,       &
        &                                   time   = time,        &
        &                                   tree   = tree,        &
        &                                   nPnts  = nPnts,       &
        &                                   res    = input_varRes )

      res(:nTotal) = res(:nTotal) * input_varRes

    end do

    deallocate(input_varRes)

  end subroutine tem_evalMultiply_forPoint
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Same as tem_evalMultiply_forPoint except it multiply values from points via
  !! indices which are setup before
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_getvalofindex.
  recursive subroutine tem_evalMultiply_fromIndex( fun, varSys, time, iLevel,  &
    &                                              idx, idxLen, nVals, res     )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    integer :: iDep, posDepVar
    real(kind=rk), allocatable :: input_varRes(:)
    ! ---------------------------------------------------------------------- !
    type(tem_varSys_op_data_type), pointer :: opData
    integer :: nTotal
    ! ---------------------------------------------------------------------- !

    call C_F_POINTER( fun%method_Data, opData )

    call tem_varSys_check_inArgs( fun, varSys, time, iLevel, idx, idxLen, &
      &    nVals, label = 'tem_evalMultiply_fromIndex'                    )

    nTotal = nVals*fun%nComponents
    allocate( input_varRes(nTotal) )

    res = 1.0_rk

    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! derive dependent variable
      call varSys%method%val(posDepVar)%get_valOfIndex( &
          & varSys  = varSys,                           &
          & time    = time,                             &
          & iLevel  = iLevel,                           &
          & idx     = opData%input_pntIndex(iDep)       &
          &           %indexLvl(iLevel)%val( idx(:) ),  &
          & nVals   = nVals,                            &
          & res     = input_varRes                      )

      res(:nTotal) = res(:nTotal) * input_varRes

    end do

    deallocate(input_varRes)

  end subroutine tem_evalMultiply_fromIndex
  ! ************************************************************************** !


  ! ***************************************************************************!
  !> Routine to divide variables if all variables have same number of
  !! components.
  !!
  !! \verbatim
  !! -- in lua file, one can define as following:
  !! variable = {{
  !!   name = 'coeff',
  !!   ncomponents = 3,
  !!   vartype = "st_fun",
  !!   st_fun =  0.25
  !!   },
  !!   {
  !!   name = 'newVel',
  !!   ncomponents = 1,
  !!   vartype = "operation",
  !!   operation = {kind='division',
  !!                input_varname={velocity, coeff}} -- numerator, denominator
  !!   }
  !!   ...
  !!  }
  !! \endverbatim
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_element.
  !!
  recursive subroutine tem_division_forElement( fun, varsys, elempos, time,   &
    &                                           tree, nElems, nDofs, res      )
    ! -------------------------------------------------------------------- !
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
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: divisor(nElems*nDofs*fun%nComponents)
    real(kind=rk) :: dividend(nElems*nDofs*fun%nComponents)
    integer :: iComp, iDof, idx, iElem
    ! -------------------------------------------------------------------- !

    call varSys%method%val(fun%input_varPos(1))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = dividend                                   )

    call varSys%method%val(fun%input_varPos(2))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = divisor                                    )

    do iElem = 1, nElems
      do iDof = 1, nDofs
        do iComp = 1, fun%nComponents
          idx = (( ielem-1)* fun%ncomponents* ndofs+( idof-1)* fun%ncomponents+icomp)
          res(idx) = dividend(idx) / divisor(idx)
        end do
      end do
    end do

  end subroutine tem_division_forElement
  ! ***************************************************************************!


  ! ************************************************************************** !
  !> Same as division_forElement except it evaluate it division values
  !! from points
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_point.
  recursive subroutine tem_division_forPoint( fun, varsys, point, time,    &
    &                                         tree, nPnts, res             )
    ! -------------------------------------------------------------------- !
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
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: divisor(nPnts*fun%nComponents)
    real(kind=rk) :: dividend(nPnts*fun%nComponents)
    integer :: iComp, iPnt, idx
    ! -------------------------------------------------------------------- !

    call varSys%method%val(fun%input_varPos(1))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = dividend                                 )

    call varSys%method%val(fun%input_varPos(2))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = divisor                                  )

    do iPnt = 1, nPnts
      do iComp = 1, fun%nComponents
        idx = (( ipnt-1)* fun%ncomponents+icomp)
        res(idx) = dividend(idx) / divisor(idx)
      end do
    end do

  end subroutine tem_division_forPoint
  ! ***************************************************************************!


  ! ************************************************************************** !
  !> Same as division_forElement except it division values from points via
  !! indices which are setup before
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_getvalofindex.
  recursive subroutine tem_division_fromIndex( fun, varSys, time, iLevel,  &
    &                                          idx, idxLen, nVals, res     )
    ! -------------------------------------------------------------------- !
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
    ! -------------------------------------------------------------------- !
    type(tem_varSys_op_data_type), pointer :: opData
    real(kind=rk) :: divisor(nVals*fun%nComponents)
    real(kind=rk) :: dividend(nVals*fun%nComponents)
    integer :: iComp, iVal
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, opData )

    ! get the nodal values for the 2 inputs
    !>TODO make it working for idxLen and contiguous access of index array
    call tem_varSys_check_inArgs( fun, varSys, time, iLevel, idx, idxLen, &
      &    nVals, label = 'tem_division_fromIndex'                         )

    call varSys%method%val(fun%input_varPos(1))%get_valOfIndex( &
      & varSys  = varSys,                                       &
      & time    = time,                                         &
      & iLevel  = iLevel,                                       &
      & idx     = opData%input_pntIndex(1)                      &
      &           %indexLvl(iLevel)%val( idx(:) ),              &
      & nVals   = nVals,                                        &
      & res     = dividend                                      )

    call varSys%method%val(fun%input_varPos(2))%get_valOfIndex( &
      & varSys  = varSys,                                       &
      & time    = time,                                         &
      & iLevel  = iLevel,                                       &
      & idx     = opData%input_pntIndex(2)                      &
      &           %indexLvl(iLevel)%val( idx(:) ),              &
      & nVals   = nVals,                                        &
      & res     = divisor                                       )

    do iVal = 1, nVals
      do iComp = 1, fun%nComponents
        res((iVal-1)*fun%nComponents + iComp) =         &
          & dividend((iVal-1)*fun%nComponents + iComp)  &
          & / divisor((iVal-1)*fun%nComponents + iComp)
      end do
    end do

  end subroutine tem_division_fromIndex
  ! ***************************************************************************!


  ! ***************************************************************************!
  !> Routine to multiply scalat times vector
  !!
  !! \verbatim
  !! -- in lua file, one can define as following:
  !! variable = {{
  !!   name = 'coeff',
  !!   ncomponents = 1,
  !!   vartype = "st_fun",
  !!   st_fun =  0.25
  !!   },
  !!   {
  !!   name = 'newVel',
  !!   ncomponents = 3,
  !!   vartype = "operation",
  !!   operation = {kind='multiply_scalar_times_vector',
  !!                input_varname={coeff, velocity }} -- scalar, vector
  !!   }
  !!   ...
  !!  }
  !! \endverbatim
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_element.
  !!
  recursive subroutine multiplyScalTimesVec_forElement( fun, varsys, elempos, &
    &                                                   time, tree, nElems,   &
    &                                                   nDofs, res            )
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
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: scalar(nElems*nDofs)
    real(kind=rk) :: vector(nElems*nDofs*fun%nComponents)
    integer :: iComp, iElem, iDof
    ! -------------------------------------------------------------------- !

    call varSys%method%val(fun%input_varPos(1))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = scalar                                     )

    call varSys%method%val(fun%input_varPos(2))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = vector                                     )

    do iElem = 1, nElems
      do iDof = 1, nDofs
        do iComp = 1, fun%nComponents
          res((( ielem-1)* fun%ncomponents* ndofs+( idof-1)* fun%ncomponents+icomp)) =        &
            & scalar((iElem-1)*nDofs + iDof)                                  &
            & * vector((( ielem-1)* fun%ncomponents* ndofs+( idof-1)* fun%ncomponents+icomp))
        end do
      end do
    end do

  end subroutine multiplyScalTimesVec_forElement
  ! ***************************************************************************!


  ! ***************************************************************************!
  !> Same as multiplyScalTimesVec_forElement except it multiply values for
  !! given corrdinate points
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_getvalofindex.
  recursive subroutine tem_multiplyScalTimesVec_forPoint( fun, varsys, point, &
    &                                                     time, tree, nPnts,  &
    &                                                     res                 )
    !--------------------------------------------------------------------------!
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
    !--------------------------------------------------------------------------!
    real(kind=rk) :: scalar(nPnts)
    real(kind=rk) :: vector(nPnts*fun%nComponents)
    integer :: iComp, iPnt
    !--------------------------------------------------------------------------!

    call varSys%method%val(fun%input_varPos(1))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = scalar                                   )

    call varSys%method%val(fun%input_varPos(2))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = vector                                   )

    do iPnt = 1, nPnts
      do iComp = 1, fun%nComponents
        res((( ipnt-1)* fun%ncomponents+icomp)) =                   &
          & scalar(iPnt) * vector((( ipnt-1)* fun%ncomponents+icomp))
      end do
    end do

  end subroutine tem_multiplyScalTimesVec_forPoint
  ! ***************************************************************************!


  ! ***************************************************************************!
  !> Same as multiplyScalTimesVec_fromElement except it multiply values
  !! from points via indices which are setup before
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_getvalofindex.
  recursive subroutine tem_multiplyScalTimesVec_fromIndex( fun, varSys, time, &
    &                                                      iLevel, idx,       &
    &                                                      idxLen, nVals, res )
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
    type(tem_varSys_op_data_type), pointer :: opData
    real(kind=rk) :: scalar(nVals)
    real(kind=rk) :: vector(nVals*fun%nComponents)
    integer :: iComp, iVal
    !--------------------------------------------------------------------------!
    write(logUnit(10),*) 'Get the values of indices for derived variable',&
      &                  'by multiplyScalTimesVec ',                      &
      &                  trim(varSys%varname%val(fun%myPos))
    call C_F_POINTER( fun%method_Data, opData )

    ! get the nodal values for the 2 inputs
    !>TODO make it working for idxLen and contiguous access of index array
    call tem_varSys_check_inArgs( fun, varSys, time, iLevel, idx, idxLen, &
      &    nVals, label = 'tem_multiplyScalTimesVec_fromIndex'            )

    call varSys%method%val(fun%input_varPos(1))%get_valOfIndex( &
      & varSys  = varSys,                                       &
      & time    = time,                                         &
      & iLevel  = iLevel,                                       &
      & idx     = opData%input_pntIndex(1)                      &
      &           %indexLvl(iLevel)%val( idx(:) ),              &
      & nVals   = nVals,                                        &
      & res     = scalar                                        )

    call varSys%method%val(fun%input_varPos(2))%get_valOfIndex( &
      & varSys  = varSys,                                       &
      & time    = time,                                         &
      & iLevel  = iLevel,                                       &
      & idx     = opData%input_pntIndex(2)                      &
      &           %indexLvl(iLevel)%val( idx(:) ),              &
      & nVals   = nVals,                                        &
      & res     = vector                                        )

    do iVal = 1, nVals
      do iComp = 1, fun%nComponents
        res((iVal-1)*fun%nComponents + iComp) =                   &
          & scalar(iVal) * vector((iVal-1)*fun%nComponents + iComp)
      end do
    end do
  end subroutine tem_multiplyScalTimesVec_fromIndex
  ! ***************************************************************************!


  ! ***************************************************************************!
  !> Routine to divide vector by scalar
  !!
  !! \verbatim
  !! -- in lua file, one can define as following:
  !! variable = {{
  !!   name = 'coeff',
  !!   ncomponents = 1,
  !!   vartype = "st_fun",
  !!   st_fun =  0.25
  !!   },
  !!   {
  !!   name = 'newVel',
  !!   ncomponents = 3,
  !!   vartype = "operation",
  !!   operation = {kind='division',
  !!                input_varname={velocity, coeff}} -- numerator, denominator
  !!   }
  !!   ...
  !!  }
  !! \endverbatim
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_element.
  !!
  recursive subroutine divideVecByScal_forElement( fun, varsys, elempos, time, &
    &                                              tree, nElems, nDofs, res    )
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
    real(kind=rk) :: divisor(nElems*nDofs)
    real(kind=rk) :: dividend(nElems*nDofs*fun%nComponents)
    integer :: iComp, iElem, iDof
    !--------------------------------------------------------------------------!

    call varSys%method%val(fun%input_varPos(1))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = dividend                                   )

    call varSys%method%val(fun%input_varPos(2))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = divisor                                    )

    do iElem = 1, nElems
      do iDof = 1, nDofs
        do iComp = 1, fun%nComponents
          res((( ielem-1)* fun%ncomponents* ndofs+( idof-1)* fun%ncomponents+icomp)) =        &
            & dividend((( ielem-1)* fun%ncomponents* ndofs+( idof-1)* fun%ncomponents+icomp)) &
            & / divisor((iElem-1)*nDofs + iDof)
        end do
      end do
    end do

  end subroutine divideVecByScal_forElement
  ! ***************************************************************************!


  ! ***************************************************************************!
  !> Same as divideVecByScal_forElement except it divides values from points
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_getvalofindex.
  recursive subroutine tem_divideVecByScal_forPoint( fun, varsys, point, time, &
    &                                                tree, nPnts, res          )
    !--------------------------------------------------------------------------!
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
    !--------------------------------------------------------------------------!
    real(kind=rk) :: divisor(nPnts)
    real(kind=rk) :: dividend(nPnts*fun%nComponents)
    integer :: iComp, iPnt
    !--------------------------------------------------------------------------!

    call varSys%method%val(fun%input_varPos(1))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = dividend                                 )

    call varSys%method%val(fun%input_varPos(2))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = divisor                                  )

    do iPnt = 1, nPnts
      do iComp = 1, fun%nComponents
        res((( ipnt-1)* fun%ncomponents+icomp)) =                      &
          & dividend((( ipnt-1)* fun%ncomponents+icomp)) / divisor(iPnt)
      end do
    end do

  end subroutine tem_divideVecByScal_forPoint
  ! ***************************************************************************!


  ! ***************************************************************************!
  !> Same as divideVecByScal_fromElement except it divide values from points
  !! via indices which are setup before
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_getvalofindex.
  recursive subroutine tem_divideVecByScal_fromIndex( fun, varSys, time,   &
    &                                                 iLevel, idx, idxLen, &
    &                                                 nVals, res           )
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
    type(tem_varSys_op_data_type), pointer :: opData
    real(kind=rk) :: divisor(nVals)
    real(kind=rk) :: dividend(nVals*fun%nComponents)
    integer :: iComp, iVal
    !--------------------------------------------------------------------------!
    write(logUnit(4),*) 'Get the values of indices for derived variable',&
      &                  'by divideVecByScal ',                          &
      &                  trim(varSys%varname%val(fun%myPos))
    call C_F_POINTER( fun%method_Data, opData )

    ! get the nodal values for the 2 inputs
    !>TODO make it working for idxLen and contiguous access of index array
    call tem_varSys_check_inArgs( fun, varSys, time, iLevel, idx, idxLen, &
      &    nVals, label = 'tem_divideVecByScal_fromIndex'            )

    call varSys%method%val(fun%input_varPos(1))%get_valOfIndex( &
      & varSys  = varSys,                                       &
      & time    = time,                                         &
      & iLevel  = iLevel,                                       &
      & idx     = opData%input_pntIndex(1)                      &
      &           %indexLvl(iLevel)%val( idx(:) ),              &
      & nVals   = nVals,                                        &
      & res     = dividend                                      )

    call varSys%method%val(fun%input_varPos(2))%get_valOfIndex( &
      & varSys  = varSys,                                       &
      & time    = time,                                         &
      & iLevel  = iLevel,                                       &
      & idx     = opData%input_pntIndex(2)                      &
      &           %indexLvl(iLevel)%val( idx(:) ),              &
      & nVals   = nVals,                                        &
      & res     = divisor                                       )

    do iVal = 1, nVals
      do iComp = 1, fun%nComponents
        res((iVal-1)*fun%nComponents + iComp) =                      &
          & dividend((iVal-1)*fun%nComponents + iComp) / divisor(iVal)
      end do
    end do
  end subroutine tem_divideVecByScal_fromIndex
  ! ***************************************************************************!


  ! ************************************************************************** !
  !> Extract component index of any vectorial variable.
  !! In lua file, first define new variable with varType operation kind as
  !! "extract" and provide name of the variable from which to extract
  !! component index via input_varname (it must be single variable) and
  !! index to extract via input_varIndex.
  !! If input_varname variable is not part of predefined solver variables then
  !! add also that variable via variable table.
  !!
  !! \verbatim
  !! -- in lua file, one can define as following:
  !! variable = {{
  !!   name = 'vel_y',
  !!   ncomponents = 1,
  !!   vartype = "operation",
  !!   operation = {kind='extract',
  !!                input_varname={'velocity'},
  !!                input_varindex = {2}
  !!               }
  !!   },
  !!  }
  !! \endverbatim
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_element.
  !!
  recursive subroutine extract_forElement( fun, varsys, elempos, time, tree, &
    &                                      nElems, nDofs, res                )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    integer :: iElem, depVar_pos, depVar_nComps, iDof, comp_index, iComp
    real(kind=rk), allocatable :: input_varRes(:)
    ! ---------------------------------------------------------------------- !

    ! get the position of dependent variable, assuming only one
    depVar_pos = fun%input_varPos(1)

    ! get the nComp of dependent variable
    depVar_nComps = varSys%method%val( depVar_pos )%nComponents

    ! allocate array to save results
    allocate(input_varRes( nElems * nDofs * depVar_nComps ))

    ! derive dependent variable
    call varSys%method%val(depVar_pos)%get_element(              &
      &                                   varSys  = varSys,      &
      &                                   elemPos = elemPos,     &
      &                                   time    = time,        &
      &                                   tree    = tree,        &
      &                                   nElems  = nElems,      &
      &                                   nDofs   = nDofs,       &
      &                                   res     = input_varRes )

    ! Extract component index from input_varRes
    do iElem = 1, nElems
      do iDof = 1, nDofs
        do iComp = 1, fun%nComponents
          comp_index = fun%input_varIndex(iComp)
          res((( ielem-1)* fun%ncomponents* ndofs+( idof-1)* fun%ncomponents+icomp) )  &
            & = input_varRes(                                          &
            & (( ielem-1)* depvar_ncomps* ndofs+( idof-1)* depvar_ncomps+comp_index))
        end do
      end do
    end do

    deallocate(input_varRes)

  end subroutine extract_forElement
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Same as extract_forElement except it extract from points
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_element.
  !!
  recursive subroutine extract_forPoint( fun, varsys, point, time, tree, &
    &                                    nPnts, res                      )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    integer :: iPnt, depVar_pos, depVar_nComps, iComp, comp_index
    real(kind=rk), allocatable :: input_varRes(:)
    ! ---------------------------------------------------------------------- !

    ! get the position of dependent variable, assuming only one
    depVar_pos = fun%input_varPos(1)

    ! get the nComp of dependent variable
    depVar_nComps = varSys%method%val( depVar_pos )%nComponents

    ! allocate array to save results
    allocate(input_varRes( nPnts * depVar_nComps ))

    ! derive dependent variable
    call varSys%method%val(depVar_pos)%get_point(                &
      &                                   varSys  = varSys,      &
      &                                   point   = point,       &
      &                                   time    = time,        &
      &                                   tree    = tree,        &
      &                                   nPnts   = nPnts,       &
      &                                   res     = input_varRes )

    ! Extract component index from input_varRes
    do iPnt = 1, nPnts
      do iComp = 1, fun%nComponents
        comp_index = fun%input_varIndex(iComp)
        res((( ipnt-1)* fun%ncomponents+icomp) )                  &
          & = input_varRes( (( ipnt-1)* depvar_ncomps+comp_index) )
      end do
    end do

    deallocate(input_varRes)

  end subroutine extract_forPoint
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Same as extract_from Point except it extract from points via indices
  !! which are setup before
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_getvalofIndex.
  !!
  recursive subroutine extract_fromIndex( fun, varSys, time, iLevel, idx, &
    &                                     idxLen, nVals, res              )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    integer :: iVal, depVar_pos, depVar_nComps, iComp, comp_index
    real(kind=rk), allocatable :: input_varRes(:)
    ! ---------------------------------------------------------------------- !
    type(tem_varSys_op_data_type), pointer :: opData
    ! ---------------------------------------------------------------------- !

    call C_F_POINTER( fun%method_Data, opData )

    call tem_varSys_check_inArgs( fun, varSys, time, iLevel, idx, idxLen, &
      &    nVals, label = 'extract_fromIndex'                             )

    ! get the position of dependent variable, assuming only one
    depVar_pos = fun%input_varPos(1)

    ! get the nComp of dependent variable
    depVar_nComps = varSys%method%val( depVar_pos )%nComponents

    ! allocate array to save results
    allocate(input_varRes( nVals * depVar_nComps ))
    ! derive dependent variable
    call varSys%method%val(depVar_pos)%get_valOfIndex( &
        & varSys  = varSys,                            &
        & time    = time,                              &
        & iLevel  = iLevel,                            &
        & idx     = opData%input_pntIndex(1)           &
        &           %indexLvl(iLevel)%val( idx(:) ),   &
        & nVals   = nVals,                             &
        & res     = input_varRes                       )

    ! Extract component index from input_varRes
    do iVal = 1, nVals
      do iComp = 1, fun%nComponents
        comp_index = fun%input_varIndex(iComp)
        res((( ival-1)* fun%ncomponents+icomp) )                  &
          & = input_varRes( (( ival-1)* depvar_ncomps+comp_index) )
      end do
    end do

    deallocate(input_varRes)

  end subroutine extract_fromIndex
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Combine multiple variables into single variable with nComponent of
  !! output variable as sum of all input variables nComponents.
  !! In lua file, first define new variable with varType operation kind as
  !! "combine" and provide name of the variable from which to extract
  !! component index via input_varname (it must be single variable) and
  !! index to combine via input_varIndex.
  !! If input_varname variable is not part of predefined solver variables then
  !! add also that variable via variable table.
  !!
  !! \verbatim
  !! -- in lua file, one can define as following:
  !! variable = {{
  !!   name = 'dens_and_vel',
  !!   ncomponents = 4,
  !!   vartype = "operation",
  !!   operation = {kind='combine',
  !!                input_varname={'density','velocity'}
  !!               }
  !!   },
  !!  }
  !! \endverbatim
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_element.
  !!
  recursive subroutine combine_forElement( fun, varsys, elempos, time, tree, &
    &                                      nElems, nDofs, res                )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !

    ! This routine combines multiple variable together
    call tem_get_element_chunk( varSys  = varSys,           &
      &                         varPos  = fun%input_varPos, &
      &                         elemPos = elemPos,          &
      &                         time    = time,             &
      &                         tree    = tree,             &
      &                         nElems  = nElems,           &
      &                         nDofs   = nDofs,            &
      &                         res     = res               )

  end subroutine combine_forElement
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> same as combine_fromelement except it extract from points
  !!
  !! the interface has to comply to the abstract interface
  !! tem_varsys_module#tem_varsys_proc_element.
  !!
  recursive subroutine combine_forPoint( fun, varsys, point, time, tree, &
    &                                    nPnts, res                      )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    ! This routine combines multiple variable together
    call tem_get_point_chunk( varSys = varSys,           &
      &                       varPos = fun%input_varPos, &
      &                       point  = point,            &
      &                       time   = time,             &
      &                       tree   = tree,             &
      &                       nPnts  = nPnts,            &
      &                       res    = res               )

  end subroutine combine_forPoint
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Same as combine_from Point except it combine from points via indices
  !! which are setup before
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_getvalofIndex.
  !!
  recursive subroutine combine_fromIndex( fun, varSys, time, iLevel, idx, &
    &                                     idxLen, nVals, res              )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    integer :: iVal, iDep
    integer :: maxComponents, compOff, dep_nComps, depVar_pos
    real(kind=rk), allocatable :: input_varRes(:)
    integer :: e_start, t_start, res_size
    type(tem_varSys_op_data_type), pointer :: opData
    ! ---------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, opData )

    call tem_varSys_check_inArgs( fun, varSys, time, iLevel, idx, idxLen, &
      &    nVals, label = 'combine_fromIndex'                             )

    ! Need to obtain the data variable for variable, and store it in an
    ! intermediate array, because all components should be put together in the
    ! res array.
    ! The temporary array therefore needs to be sufficiently large to store the
    ! maximal number of components.
    maxComponents = maxval(varSys%method%val(fun%input_varPos(:))%nComponents)

    ! Using a temporary array to store the variables and transfer them to res
    ! in the correct ordering afterwards.
    allocate(input_varRes(nVals*maxComponents))

    compOff = 0
    do iDep = 1, fun%nInputs

      ! get the position of dependent variable
      depVar_pos = fun%input_varPos(iDep)

      ! get the number of components for variable iVar
      dep_nComps = varSys%method%val(depVar_pos)%nComponents

      ! get the size of the needed part of the res array
      res_size = nVals * dep_nComps

      ! derive dependent variable
      call varSys%method%val(depVar_pos)%get_valOfIndex( &
          & varSys  = varSys,                            &
          & time    = time,                              &
          & iLevel  = iLevel,                            &
          & idx     = opData%input_pntIndex(iDep)        &
          &           %indexLvl(iLevel)%val( idx(:) ),   &
          & nVals   = nVals,                             &
          & res     = input_varRes(:res_size)            )

      ! copy the information to the right positions in the result array
      ! res contains results for all variables,
      ! input_varRes is only for one variable
      do iVal = 1, nVals
        e_start = (iVal-1)*fun%nComponents + compOff
        t_start = (iVal-1)*dep_nComps
        res( (e_start+1) : (e_start+dep_nComps) )               &
          &  = input_varRes( t_start + 1 : t_start + dep_nComps )
      end do
      ! Increase the component offset for the next variables.
      compOff = compOff + dep_nComps
    end do !iDep

    deallocate(input_varRes)

  end subroutine combine_fromIndex
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Routine to return time reduction variables
  !!
  !! \verbatim
  !! -- in lua file, one can define as following:
  !! variable = {
  !!   {
  !!     name = 'press_timeavg',
  !!     ncomponents = 1,
  !!     vartype = "operation",
  !!     operation = {
  !!       kind='reduction_transient',
  !!       input_varname={'pressure'},
  !!       reduction_transient = {kind = 'average', nrecord = 1000}
  !!     }
  !!   }
  !!   ...
  !!  }
  !! \endverbatim
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_element.
  !!
  recursive subroutine reductionTransient_forElement( fun, varsys, elempos, &
    &                                                 time, tree, nElems,   &
    &                                                 nDofs, res            )
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
    ! ------------------------------------------------------------------------ !
    type(tem_varSys_op_data_type), pointer :: opData
    ! ------------------------------------------------------------------------ !
    !call C_F_POINTER( fun%method_Data, opData )
    ! Avoid warning about unused varsys dummy argument
    call C_F_POINTER( varsys%method%val(fun%mypos)%method_Data, opData )

    if (time%sim < 0.0_rk) then
      write(logunit(10),*) 'Avoid unused argument warning for time'
      write(logunit(10),*) 'tree%nElems:', tree%nElems
    end if

    res = 0.0_rk
    call tem_reduction_transient_getElement( me      = opData%redTrans, &
      &                                      elemPos = elemPos,         &
      &                                      nElems  = nElems,          &
      &                                      nDofs   = nDofs,           &
      &                                      res     = res              )

  end subroutine reductionTransient_forElement
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Same as reductionTransient_forElement except it evaluate it multiply values
  !! from points
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_point.
  recursive subroutine reductionTransient_forPoint( fun, varsys, point, time, &
    &                                               tree, nPnts, res          )
    ! ---------------------------------------------------------------------- !
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
    ! ---------------------------------------------------------------------- !
    res = 0.0_rk
    write(logunit(1),*) 'Variable: ', trim(varsys%varname%val(fun%mypos))
    write(logunit(1),*) 'nPnts=', nPnts
    write(logunit(1),*) 'tree%nElems=', tree%nElems
    write(logunit(1),*) 'time%sim=', time%sim
    write(logunit(1),*) 'size(point)=', size(point)
    flush(logunit(1))
    call tem_abort('Reduction_transient for Point is not implemented yet')

  end subroutine reductionTransient_forPoint
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Same as reductionTransient_forPoint except it multiply values from points
  !! via indices which are setup before
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_getvalofindex.
  recursive subroutine reductionTransient_fromIndex( fun, varSys, time,   &
    &                                                iLevel, idx, idxLen, &
    &                                                nVals, res           )
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
    ! ---------------------------------------------------------------------- !
    res = 0.0_rk
    call tem_varSys_check_inArgs( fun, varSys, time, iLevel, idx, idxLen, &
      &    nVals, label = 'reductionTransient_fromIndex'                  )
    call tem_abort('Reduction_transient from Index is not implemented yet')
  end subroutine reductionTransient_fromIndex
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This subroutine call set_params of input_variable
  !!
  !! the interface has to comply to the abstract interface
  !! tem_varsys_module#tem_varsys_proc_setParams.
  !!
  recursive subroutine tem_opVar_setParams(fun, varSys, instring)
    ! ---------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Input string with parameter to set in method_data
    character(len=*), intent(in) :: instring
    ! ---------------------------------------------------------------------- !
    integer :: iDep, posDepVar
    ! ---------------------------------------------------------------------- !
    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! set params in dependent variable
      call varSys%method%val(posDepVar)%set_params( varSys   = varSys,  &
        &                                           inString = inString )
    end do

  end subroutine tem_opVar_setParams
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This subroutine call get_params of input_variable
  !!
  !! the interface has to comply to the abstract interface
  !! tem_varsys_module#tem_varsys_proc_setParams.
  !!
  recursive subroutine tem_opVar_getParams( fun, varSys, instring, outstring )
    ! ---------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Input string with parameter to set in method_data
    character(len=*), intent(in) :: instring

    !> Output string with requested parameter value from method_data
    character(len=*), intent(out) :: outstring
    ! ---------------------------------------------------------------------- !
    integer :: iDep, posDepVar
    ! ---------------------------------------------------------------------- !
    do iDep = 1, fun%nInputs
      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! get params from dependent variable
      call varSys%method%val(posDepVar)%get_params( varSys    = varSys,   &
        &                                           inString  = inString, &
        &                                           outString = outString )

      ! if outString is filled by any dependent variable then exit loop
      if (len(trim(outString)) > 0) then
        outString = trim(outString)//'_oper'
        exit
      end if
    end do

  end subroutine tem_opVar_getParams
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This subroutine call setup indices of input_variable
  !!
  !! the interface has to comply to the abstract interface
  !! tem_varsys_module#tem_varsys_proc_setupIndices.
  !!
  recursive subroutine tem_opVar_setupIndices( fun, varSys, point, offset_bit, &
    &                                          iLevel, tree, nPnts, idx        )
    ! ---------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> List of space coordinate points to store as growing array in
    !! method_data
    real(kind=rk), intent(in) :: point(:,:)

    !> Offset bit encoded as character for every point.
    !!
    !! Offset integer coord(3) is converted into a character with
    !! offset_bit = achar( (coord(1)+1) + (coord(2)+1)*4 + (coord(3)+1)*16 )
    !! Backward transformation form character to 3 integer:
    !! coord(1) = mod(ichar(offset_bit),4) - 1
    !! coord(2) = mod(ichar(offset_bit),16)/4 - 1
    !! coord(3) = ichar(offset_bit)/16 - 1
    !!
    !! If not present default is to center i.e offset_bit = achar(1+4+16)
    character, optional, intent(in) :: offset_bit(:)

    !> Level to which input points belong to
    integer, intent(in) :: iLevel

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of points to add in method_data of this variable
    integer, intent(in) :: nPnts

    !> Index of points in the growing array and variable val array.
    !! Size: nPoints
    !!
    !! This must be stored in boundary or source depends on who
    !! calls this routine.
    !! This index is required to return a value using getValOfIndex.
    integer, intent(out) :: idx(:)
    ! ---------------------------------------------------------------------- !
    type(tem_varSys_op_data_type), pointer :: opData
    integer :: iPnt, iDep
    type(grw_intArray_type), allocatable :: inputIndex_loc(:)
    integer, allocatable :: idxPerPnt(:)
    ! ---------------------------------------------------------------------- !

    call C_F_POINTER( fun%method_Data, opData )

    ! allcoate the index array for all inpits
    if (.not. allocated(opData%input_pntIndex)) then
      allocate( opData%input_pntIndex(fun%nInputs) )
    end if
    ! allocate temporary inputIndex with size of nInputs and initialize
    ! growing array with length nPnts
    allocate(inputIndex_loc(fun%nInputs))

    ! Now fill in the index arrays for the inputs
    call tem_opVar_fill_inputIndex( fun        = fun,                &
      &                             varSys     = varSys,             &
      &                             point      = point,              &
      &                             offset_bit = offset_bit,         &
      &                             iLevel     = iLevel,             &
      &                             tree       = tree,               &
      &                             nPnts      = nPnts,              &
      &                             inputIndex = inputIndex_loc      )

    ! KM: Workaround for intel compiler in SuperMUC, the pointer
    ! gets corrupted after recursive call to depend variable
    call C_F_POINTER( fun%method_Data, opData )
    ! initialize index with zero to identify points which does not
    ! belong to subTree
    allocate(idxPerPnt(fun%nInputs))
    idx = 0
    do iPnt = 1, nPnts
      do iDep = 1, fun%nInputs
        idxPerPnt(iDep) = inputIndex_loc(iDep)%val(iPnt)
      end do
      ! set index only when any of dependent variable has valid index
      if (any(idxPerPnt > 0)) then
        do iDep = 1, fun%nInputs
          call append(me  = opData%input_pntIndex(iDep)%indexLvl(iLevel), &
            &         val = inputIndex_loc(iDep)%val(iPnt)              )
        end do
        ! set index to last position in input_pntIndex of dep var 1 of
        ! indexLvl of iLevel
        idx(iPnt) = opData%input_pntIndex(1)%indexLvl(iLevel)%nVals
      end if
    end do

    do iDep = 1, fun%nInputs
      call truncate (opData%input_pntIndex(iDep)%indexLvl(iLevel) )
    end do

  end subroutine tem_opVar_setupIndices
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> subroutine to fill index for the setuo Index routine called for operation
  !! variables, it is also used by the solver
  recursive subroutine tem_opVar_fill_inputIndex( fun, varSys, point,       &
    &                                             offset_bit, iLevel, tree, &
    &                                             nPnts, inputIndex         )
    !---------------------------`----------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> List of space coordinate points to store as growing array in
    !! method_data
    real(kind=rk), intent(in) :: point(:,:)

    !> Offset bit encoded as character for every point.
    !!
    !! Offset integer coord(3) is converted into a character with
    !! offset_bit = achar( (coord(1)+1) + (coord(2)+1)*4 + (coord(3)+1)*16 )
    !! Backward transformation form character to 3 integer:
    !! coord(1) = mod(ichar(offset_bit),4) - 1
    !! coord(2) = mod(ichar(offset_bit),16)/4 - 1
    !! coord(3) = ichar(offset_bit)/16 - 1
    !!
    !! If not present default is to center i.e offset_bit = achar(1+4+16)
    character, optional, intent(in) :: offset_bit(:)

    !> Level to which input points belong to
    integer, intent(in) :: iLevel

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of points to add in method_data of this variable
    integer, intent(in) :: nPnts

    !> input index for dependent variables
    !! size: fun%nInputs
    type(grw_intArray_type), intent(out) :: inputIndex(:)
    !--------------------------------------------------------------------------!
    integer, allocatable :: idx_loc(:)
    integer :: iDep, posDepVar
    !--------------------------------------------------------------------------!
    ! allocate local index array, needed for setup_indices call
    allocate( idx_loc(nPnts) )

    do iDep = 1, fun%nInputs
      idx_loc = 0

      ! get position of dependent var in varSys
      posDepVar = fun%input_varPos(iDep)

      ! setup indices in dependent variable.
      ! idx output from any variables will be same so it just overwrites and
      ! returns the idx of last dependent variable
      call varSys%method%val(posDepVar)%setup_indices( &
        & varSys     = varSys,                         &
        & point      = point,                          &
        & offset_bit = offset_bit,                     &
        & iLevel     = iLevel,                         &
        & tree       = tree,                           &
        & nPnts      = nPnts,                          &
        & idx        = idx_loc                         )

      call append( me  = inputIndex(iDep), &
        &          val = idx_loc           )

      call truncate (inputIndex(iDep))
      write(logUnit(9),*) 'nIndex on level for input variable ',  &
        &                 trim(varSys%varname%val(posDepVar)),    &
        &                 ' on Lvl ', iLevel,                     &
        &                 ' are = ', inputIndex(iDep)%nVals
    end do
    deallocate(idx_loc)

  end subroutine tem_opVar_fill_inputIndex

end module tem_operation_var_module
! **************************************************************************** !
