! Copyright (c) 2016-2017, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016-2017 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
! ***************************************************************************** !
!> This module provides variable definitions for space time functions.
!!
module tem_spacetime_var_module
  use, intrinsic :: iso_c_binding,  only: c_ptr,       &
    &                                     c_loc,       &
    &                                     c_f_pointer

  ! include treelm modules
  use env_module,               only: rk, long_k, labelLen
  use tem_param_module,         only: qOffset_inChar, q000
  use tem_varSys_module,        only: tem_varSys_type,               &
    &                                 tem_varSys_op_type,            &
    &                                 tem_varSys_check_inArgs,       &
    &                                 tem_varSys_append_derVar,      &
    &                                 tem_varSys_proc_point,         &
    &                                 tem_varSys_proc_element,       &
    &                                 tem_varSys_proc_setParams,     &
    &                                 tem_varSys_proc_getParams,     &
    &                                 tem_varSys_proc_setupIndices,  &
    &                                 tem_varSys_proc_getValOfIndex, &
    &                                 tem_varSys_solverData_evalElem_type
  use tem_tools_module,         only: tem_PositionInSorted
  use tem_geometry_module,      only: tem_CoordOfReal, &
    &                                 tem_posOfId
  use tem_topology_module,      only: tem_IdOfCoord
  use tem_variable_module,      only: tem_variable_type
  use tem_time_module,          only: tem_time_type
  use tem_spacetime_fun_module, only: tem_spacetime_for,          &
    &                                 tem_st_fun_linkedList_type, &
    &                                 append,                     &
    &                                 tem_st_fun_listElem_type,   &
    &                                 tem_spacetime_fun_type
  use tem_spatial_module,       only: tem_spatial_storeVal
  use treelmesh_module,         only: treelmesh_type
  use tem_logging_module,       only: logUnit, llerror
  use tem_aux_module,           only: tem_abort
  use tem_dyn_array_module,     only: truncate, PositionOfVal
  use tem_grow_array_module,    only: init, append, truncate
  use tem_pointData_module,     only: init, append, truncate

  ! aotus modules
  use aotus_module,     only: aot_get_val, aot_top_get_val,                    &
    &                         aoterr_Fatal, aoterr_WrongType,                  &
    &                         aoterr_NonExistent, flu_State, open_config_chunk

  implicit none

  private

  public :: tem_varSys_append_stFun
  public :: evaluate_add_spacetime_scalarByTreeID
  public :: evaluate_add_spacetime_vectorByTreeID
  public :: evaluate_add_spacetime_scalarByCoordinate
  public :: evaluate_add_spacetime_vectorByCoordinate
  public :: evaluate_first_spacetime_scalarByCoordinate
  public :: evaluate_first_spacetime_scalarByTreeID
  public :: evaluate_first_spacetime_vectorByCoordinate
  public :: evaluate_first_spacetime_vectorByTreeID
  public :: get_valOfIndex_add_scalar_spacetime
  public :: get_valOfIndex_add_vector_spacetime
  public :: get_valOfIndex_first_scalar_spacetime
  public :: get_valOfIndex_first_vector_spacetime
  public :: set_params_spacetime
  public :: get_params_spacetime
  public :: setup_indices_spacetime
  public :: evaluate_FOAG_spacetime_scalarByTreeID
  public :: evaluate_FOAG_spacetime_vectorByTreeID
  public :: evaluate_FOAG_spacetime_scalarByCoordinate
  public :: evaluate_FOAG_spacetime_vectorByCoordinate
  public :: get_valOfIndex_FOAG_scalar_spacetime
  public :: get_valOfIndex_FOAG_vector_spacetime

  interface tem_varSys_append_stfun
    module procedure tem_varSys_append_stFunVar
    module procedure tem_varSys_append_stFun_raw
  end interface tem_varSys_append_stfun


contains


  ! *************************************************************************** !
  !> Returns the get_point and get_element pointer according to the requested
  !! evaluation type.
  subroutine tem_varSys_assignEvalType(evaltype, nComp, get_point, &
    &                                  get_element, get_valOfIndex )
    ! -------------------------------------------------------------------------
    character(len=*), intent(in) :: evaltype
    integer, intent(in) :: nComp
    !> The function pointer to the get_point subroutine for the given
    !! operation.
    procedure(tem_varSys_proc_point), pointer, intent(out) :: get_point
    !> The function pointer to the get_element subroutine for the given
    !! operation.
    procedure(tem_varSys_proc_element), pointer, intent(out) :: get_element
    !> The function pointer to the get_valOfIndex subroutine for the given
    !! operation.
    procedure(tem_varSys_proc_getValOfIndex), pointer, intent(out) :: &
      &                                       get_valOfIndex
    ! -------------------------------------------------------------------------

    select case(evalType)

    case ('add')
      ! if nComp = 1, use scalar function else use vector function
      if ( nComp == 1 ) then
        get_point => evaluate_add_spacetime_scalarByCoordinate
        get_element => evaluate_add_spacetime_scalarByTreeID
        get_valOfIndex => get_valOfIndex_add_scalar_spacetime
      else
        get_point => evaluate_add_spacetime_vectorByCoordinate
        get_element => evaluate_add_spacetime_vectorByTreeID
        get_valOfIndex => get_valOfIndex_add_vector_spacetime
      end if

    case ('first')
      ! if nComp = 1, use scalar function else use vector function
      if ( nComp == 1 ) then
        get_point => evaluate_first_spacetime_scalarByCoordinate
        get_element => evaluate_first_spacetime_scalarByTreeID
        get_valOfIndex => get_valOfIndex_first_scalar_spacetime
      else
        get_point => evaluate_first_spacetime_vectorByCoordinate
        get_element => evaluate_first_spacetime_vectorByTreeID
        get_valOfIndex => get_valOfIndex_first_vector_spacetime
      end if
    case ('firstonly_asglobal')
      ! Anonymous spacetime function variable created for boundary
      ! or source variable assumes global shape and single st-fun
      !
      ! if nComp = 1, use scalar function else use vector function
      if ( nComp == 1 ) then
        get_point => evaluate_FOAG_spacetime_scalarByCoordinate
        get_element => evaluate_FOAG_spacetime_scalarByTreeID
        get_valOfIndex => get_valOfIndex_FOAG_scalar_spacetime
      else
        get_point => evaluate_FOAG_spacetime_vectorByCoordinate
        get_element => evaluate_FOAG_spacetime_vectorByTreeID
        get_valOfIndex => get_valOfIndex_FOAG_vector_spacetime
      end if

    end select

  end subroutine tem_varSys_assignEvalType
  ! *************************************************************************** !


  ! *************************************************************************** !
  !> subroutine to add the variables from the input lua script to the varsys
  subroutine tem_varSys_append_stFunVar( stFunVar, varSys, st_funList, &
    &                                    solverData_evalElem           )
    ! -------------------------------------------------------------------------
    !> variables defined in the lua file
    type(tem_variable_type), intent(in)           :: stFunVar

    !> global variable system to which stFunVar to be appended
    type(tem_varSys_type), intent(inout)            :: varSys

    !> contains spacetime functions of all variables
    type(tem_st_fun_linkedList_type), intent(inout) :: st_funList

    !> A setter routine that allows the caller to define routine for the
    !! construction of an element representation.
    type(tem_varSys_solverData_evalElem_type), &
      &  optional, intent(in) :: solverData_evalElem
    ! -------------------------------------------------------------------------
    integer :: addedPos
    integer :: nComp
    logical :: wasAdded
    type(c_ptr) :: method_data
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setParams), pointer :: set_params => null()
    procedure(tem_varSys_proc_getParams), pointer :: get_params => null()
    procedure(tem_varSys_proc_setupIndices), pointer :: &
      &                                      setup_indices => null()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: &
      &                                       get_valOfIndex => null()
    type(tem_st_fun_listElem_type), pointer  :: newElem
    ! -------------------------------------------------------------------------
    nullify(get_point, get_element, set_params, get_params, setup_indices, &
      &     get_valOfIndex)
    nComp = stFunVar%nComponents

    ! append space time function to linked list of spacetime functions
    call append( st_funList, stFunVar%st_fun, newElem )

    ! c pointer to list of spacetime functions of current variable
    !method_data = c_loc(st_funList%current)
    method_data = c_loc(newElem)

    ! assign function pointer depends on evaluation type
    call tem_varSys_assignEvalType( evaltype       = stfunVar%evaltype,    &
      &                             nComp          = stfunVar%nComponents, &
      &                             get_point      = get_point,            &
      &                             get_element    = get_element,          &
      &                             get_valOfIndex = get_valOfIndex        )

    set_params => set_params_spacetime
    get_params => get_params_spacetime
    setup_indices => setup_indices_spacetime

    if (.not. associated(get_point)) then
      call tem_abort( 'Error: No evaluation is defined for variable ' &
        & // trim(stfunvar%label) )
    end if

    ! append variable to varSys
    call tem_varSys_append_derVar( me             = varSys,         &
      &                            varName        = stFunVar%label, &
      &                            operType       = 'st_fun',       &
      &                            nComponents    = nComp,          &
      &                            method_data    = method_data,    &
      &                            get_point      = get_point,      &
      &                            get_element    = get_element,    &
      &                            set_params     = set_params,     &
      &                            get_params     = get_params,     &
      &                            setup_indices  = setup_indices,  &
      &                            get_valOfIndex = get_valOfIndex, &
      &                            pos            = addedPos,       &
      &                            wasAdded       = wasAdded        )

    if (wasAdded) then
      if (present(solverData_evalElem)) then
        ! If an solverData_evalElem function is provided,
        ! override the get_element pointer and use the provided setter
        ! solverData_evalElem instead to define the get_element routine.
        call solverData_evalElem%stFun_setter(varSys%method%val(addedPos))
      end if
      write(logUnit(9),*) 'Successfully appended variable "' &
        & // trim(stFunVar%label) // '" to the variable system'
    else if (addedpos < 1) then
      write(logUnit(1),*) 'WARNING: variable '//trim(stFunVar%label)// &
        &                 ' is not added to variable system'
    end if

  end subroutine tem_varSys_append_stFunVar
  ! *************************************************************************** !


  ! *************************************************************************** !
  !> subroutine to add the variables from the input lua script to the varsys
  subroutine tem_varSys_append_stFun_raw( varSys, stFun, varname, nComp, &
    &                                     evaltype, st_funList,          &
    &                                     solverData_evalElem            )
    ! -------------------------------------------------------------------------
    !> global variable system to which stFunVar to be appended
    type(tem_varSys_type), intent(inout) :: varSys

    !> variables defined in the lua file
    type(tem_spacetime_fun_type), pointer, intent(in) :: stFun(:)

    character(len=*), intent(in) :: varname

    integer, intent(in), optional :: nComp

    character(len=*), intent(in), optional :: evaltype

    !> contains spacetime functions of all variables
    type(tem_st_fun_linkedList_type), intent(inout), optional :: st_funList

    !> A setter routine that allows the caller to define routine for the
    !! construction of an element representation.
    type(tem_varSys_solverData_evalElem_type), &
      &  optional, intent(in) :: solverData_evalElem
    ! -------------------------------------------------------------------------
    type(tem_st_fun_listElem_type), pointer :: stfun_listelem
    integer :: addedPos
    logical :: wasAdded
    type(c_ptr) :: method_data
    integer :: ncomp_loc
    character(len=labelLen) :: evaltype_loc
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setParams), pointer :: set_params => null()
    procedure(tem_varSys_proc_getParams), pointer :: get_params => null()
    procedure(tem_varSys_proc_setupIndices), pointer :: &
      &                                      setup_indices => null()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: &
      &                                       get_valOfIndex => null()
    type(tem_st_fun_listElem_type),pointer :: newElem
    ! -------------------------------------------------------------------------
    nullify(get_point, get_element, set_params, get_params, setup_indices, &
      &     get_valOfIndex)

    nComp_loc = 1
    if (present(nComp)) nComp_loc = nComp

    evaltype_loc = 'add'
    if (present(evaltype)) evaltype_loc = evaltype

    if (present(st_funList)) then
      ! append space time function to linked list of spacetime functions
      call append( st_funList, stfun, newElem )

      ! c pointer to list of spacetime functions of current variable
      method_data = c_loc(newElem)
    else
      allocate(stfun_listelem)
      stfun_listelem%val => stfun
      stfun_listelem%nvals = size(stfun)
      method_data = c_loc(stfun_Listelem)
    end if


    ! assign function pointer depends on evaluation type
    call tem_varSys_assignEvalType( evaltype       = evaltype_loc,  &
      &                             nComp          = nComp_loc,     &
      &                             get_point      = get_point,     &
      &                             get_element    = get_element,   &
      &                             get_valOfIndex = get_valOfIndex )

    set_params => set_params_spacetime
    get_params => get_params_spacetime
    setup_indices => setup_indices_spacetime

    if (.not. associated(get_point)) then
      write(logUnit(1),*) 'Error: No evaluation is defined for variable '//&
        &                 trim(varname)
      call tem_abort()
    end if

    ! append variable to varSys
    call tem_varSys_append_derVar( me             = varSys,         &
      &                            varName        = varName,        &
      &                            operType       = 'st_fun',       &
      &                            nComponents    = nComp_loc,      &
      &                            method_data    = method_data,    &
      &                            get_point      = get_point,      &
      &                            get_element    = get_element,    &
      &                            set_params     = set_params,     &
      &                            get_params     = get_params,     &
      &                            setup_indices  = setup_indices,  &
      &                            get_valOfIndex = get_valOfIndex, &
      &                            pos            = addedPos,       &
      &                            wasAdded       = wasAdded        )

    if (wasAdded) then
      if (present(solverData_evalElem)) then
        ! If an solverData_evalElem function is provided, override
        ! the get_element pointer and use the provided setter
        ! solverData_evalElem instead to define the get_element routine.
        call solverData_evalElem%stFun_setter(varSys%method%val(addedPos))
      end if
      write(logUnit(9),*) 'Successfully appended variable "' &
        & // trim(varname) // '" to the variable system'

    else if (addedpos < 1) then
      write(logUnit(1),*) 'WARNING: variable '//trim(varname)// &
        &                 ' is not added to variable system'
    end if

  end subroutine tem_varSys_append_stFun_raw
  ! *************************************************************************** !


  ! *************************************************************************** !
  !> This subroutine directly evaluate spacetime function for given coordinate
  !! points on 1st st_fun assuming global shape.
  recursive subroutine evaluate_FOAG_spacetime_scalarByCoordinate( fun, &
    & varsys, point, time, tree, nPnts, res                             )
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
    ! --------------------------------------------------------------------------
    type(tem_st_fun_listElem_type), pointer :: fPtr
    ! --------------------------------------------------------------------------
    !call C_F_POINTER( fun%method_Data, fPtr )
    ! Avoid compiler warning about unused varsys
    call C_F_POINTER( varsys%method%val(fun%mypos)%method_Data, fPtr )
    res = tem_spacetime_for( me    = fPtr%val(1), &
      &                      coord = point,       &
      &                      time  = time,        &
      &                      n     = nPnts        )

    if (tree%nElems < 0) write(logunit(10),*) 'Avoid unused dummy argument'

  end subroutine evaluate_FOAG_spacetime_scalarByCoordinate
  ! *************************************************************************** !


  ! *************************************************************************** !
  !> This subroutine add multiple spacetime function for given coordinate
  !! points. Spacetime function are evaluated only on the coordinate which
  !! belong to st_fun shape.
  recursive subroutine evaluate_add_spacetime_scalarByCoordinate( fun, varsys, &
    & point, time, tree, nPnts, res                                            )
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
    ! --------------------------------------------------------------------------
    integer :: iStFun, iPoint, elemPos, pos
    type(tem_st_fun_listElem_type), pointer :: fPtr
    integer :: coord(4)
    integer(kind=long_k):: treeID
    real(kind=rk) :: st_res(1)
    ! --------------------------------------------------------------------------
    !call C_F_POINTER( fun%method_Data, fPtr )
    ! Avoid compiler warning about unused varsys
    call C_F_POINTER( varsys%method%val(fun%mypos)%method_Data, fPtr )
    res = 0.0_rk

    ! The space time function has a custom shape, thus we have to check
    ! whether the coordinates are covered by the space time function.
    ! As we are dealing with vectorial access here, we have to loop over all
    ! requested points.
    stLoop: do iStFun = 1, fPtr%nVals
      if(fPtr%val(iStFun)%subTree%useGlobalMesh) then
        res = res + tem_spacetime_for( me      = fPtr%val(iStFun), &
          &                            coord   = point,            &
          &                            time    = time,             &
          &                            n       = nPnts             )
      else
        do iPoint = 1, nPnts
          ! The space time function has a custom shape, thus we have to check
          ! whether the coordinates are covered by the space time function.
          coord =  tem_CoordOfReal(tree, point(iPoint,:), &
            &                      tree%global%maxLevel)
          treeId = tem_IdOfCoord(coord)
          elemPos = tem_PosofId(treeId, tree%treeID)
          pos = tem_PositionInSorted(                      &
            & me    = fPtr%val(iStFun)%subTree%map2global, &
            & val   = elemPos                              )
          ! When the elemenet is covered by the space time function's shape,
          ! we evaluate it at the given coordinates.
          if (pos > 0) then
            st_res = tem_spacetime_for(                            &
              & me      = fPtr%val(iStFun),                        &
              & coord   = reshape( source = (/ point(iPoint,1),    &
              &                                point(iPoint,2),    &
              &                                point(iPoint,3) /), &
              &                    shape  = (/ 1, 3 /) ),          &
              & time    = time,                                    &
              & n       = 1                                        )

            res(iPoint) = res(iPoint) + st_res(1)
          end if ! point part of subtree
        end do !iPoint
      end if ! global mesh
    end do stLoop

  end subroutine evaluate_add_spacetime_scalarByCoordinate
  ! *************************************************************************** !


  ! *************************************************************************** !
  !> This subroutine returns first spacetime function on overlapping
  !! regions fo given coordinate points.
  !! Spacetime function are evaluated only on the coordinate which
  !! belong to st_fun shape.
  recursive subroutine evaluate_first_spacetime_scalarByCoordinate( fun, &
    & varsys, point, time, tree, nPnts, res                              )
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
    ! --------------------------------------------------------------------------
    integer :: iStFun, iPnt, elemPos, pos
    logical :: found
    type(tem_st_fun_listElem_type), pointer :: fPtr
    integer :: coord(4)
    integer(kind=long_k):: treeID
    real(kind=rk) :: st_res(1)
    ! --------------------------------------------------------------------------
    !call C_F_POINTER( fun%method_Data, fPtr )
    ! avoid compiler warning about unused varsys
    call C_F_POINTER( varsys%method%val(fun%mypos)%method_Data, fPtr )
    res = 0.0_rk

    do iPnt = 1, nPnts
      found = .false.
      ! Loop over all the space-time functions, if value is set by one then go
      ! next point
      stLoop: do iStFun = 1, fPtr%nVals
        if(fPtr%val(iStFun)%subTree%useGlobalMesh) then
          found = .true.
        else
          ! The space time function has a custom shape, thus we have to check
          ! whether the coordinates are covered by the space time function.
          coord =  tem_CoordOfReal(tree, point(iPnt,:), tree%global%maxLevel )
          treeId = tem_IdOfCoord(coord)
          elemPos = tem_PosofId(treeId, tree%treeID)
          pos = tem_PositionInSorted(                      &
            & me    = fPtr%val(iStFun)%subTree%map2global, &
            & val   = elemPos                              )
          ! When the elemenet is covered by the space time function's shape, we
          ! evaluate it at the given coordinates.
          if (pos > 0) then
            found = .true.
          end if ! point part of subtree
        end if ! global mesh

        ! When the first space time function was able to fulfill the request, we
        ! stop looping over them, because we only want one value.
        if (found) then
          st_res = tem_spacetime_for(                          &
            & me      = fPtr%val(iStFun),                      &
            & coord   = reshape( source = (/ point(iPnt,1),    &
            &                                point(iPnt,2),    &
            &                                point(iPnt,3) /), &
            &                    shape  = (/ 1, 3 /) ),        &
            & time    = time,                                  &
            & n       = 1                                      )

          res(iPnt) = st_res(1)
          exit stLoop
        endif
      end do stLoop
    end do !iPnt

  end subroutine evaluate_first_spacetime_scalarByCoordinate
  ! *************************************************************************** !


  ! *************************************************************************** !
  !> This subroutine directly evaluate spacetime function for given coordinate
  !! points on 1st st_fun assuming global shape for vectorial variable
  recursive subroutine evaluate_FOAG_spacetime_vectorByCoordinate( fun, &
    & varsys, point, time, tree, nPnts, res                             )
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
    ! --------------------------------------------------------------------------
    type(tem_st_fun_listElem_type), pointer :: fPtr
    integer :: iPoint
    real(kind=rk) :: st_res(nPnts, fun%nComponents)
    ! --------------------------------------------------------------------------
    !call C_F_POINTER( fun%method_Data, fPtr )
    ! avoid compiler warning about unused varsys
    call C_F_POINTER( varsys%method%val(fun%mypos)%method_Data, fPtr )

    st_res = tem_spacetime_for( me      = fPtr%val(1),    &
      &                         coord   = point,          &
      &                         time    = time,           &
      &                         n       = nPnts,          &
      &                         nComp   = fun%nComponents )

    do iPoint = 1, nPnts
      res( (iPoint-1)*fun%nComponents + 1 : iPoint*fun%nComponents ) &
        & = st_res(iPoint, :)
    end do
    if (tree%nElems < 0) write(logunit(10),*) 'Avoid unused dummy argument'

  end subroutine evaluate_FOAG_spacetime_vectorByCoordinate
  ! *************************************************************************** !


  ! *************************************************************************** !
  !> This subroutine add multiple spacetime function for given coordinate
  !! points. Spacetime function are evaluated only on the coordinate which
  !! belong to st_fun shape for vectorial variable.
  recursive subroutine evaluate_add_spacetime_vectorByCoordinate( fun, varsys, &
    & point, time, tree, nPnts, res                                            )
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
    ! --------------------------------------------------------------------------
    integer :: iStFun, iPoint, elemPos, pos
    type(tem_st_fun_listElem_type), pointer :: fPtr
    integer :: coord(4)
    integer(kind=long_k):: treeID
    real(kind=rk) :: st_res_pnt(1, fun%nComponents)
    real(kind=rk) :: st_res(nPnts, fun%nComponents)
    ! --------------------------------------------------------------------------
    !call C_F_POINTER( fun%method_Data, fPtr )
    ! avoid compiler warning about unused varsys
    call C_F_POINTER( varsys%method%val(fun%mypos)%method_Data, fPtr )

    st_res = 0.0_rk

    ! The space time function has a custom shape, thus we have to check
    ! whether the coordinates are covered by the space time function.
    ! As we are dealing with vectorial access here, we have to loop over all
    ! requested points.
    stLoop: do iStFun = 1, fPtr%nVals
      if (fPtr%val(iStFun)%subTree%useGlobalMesh) then
        st_res = st_res + tem_spacetime_for( me      = fPtr%val(iStFun), &
          &                                  coord   = point,            &
          &                                  time    = time,             &
          &                                  n       = nPnts,            &
          &                                  nComp   = fun%nComponents   )
      else
        ! The space time function has a custom shape, thus we have to check
        ! whether the coordinates are covered by the space time function.
        do iPoint = 1, nPnts
          coord =  tem_CoordOfReal(tree, point(iPoint,:), tree%global%maxLevel)
          treeId = tem_IdOfCoord(coord)
          elemPos = tem_PosofId(treeId, tree%treeID)
          pos = tem_PositionInSorted(                      &
            & me    = fPtr%val(iStFun)%subTree%map2global, &
            & val   = elemPos                              )
          ! When the element is covered by the space time function's shape, we
          ! evaluate it at the given coordinates.
          if (pos > 0) then
            st_res_pnt = tem_spacetime_for(                        &
              & me      = fPtr%val(iStFun),                        &
              & coord   = reshape( source = (/ point(iPoint,1),    &
              &                                point(iPoint,2),    &
              &                                point(iPoint,3) /), &
              &                    shape  = (/ 1, 3 /) ),          &
              & time    = time,                                    &
              & n       = 1,                                       &
              & nComp   = fun%nComponents                          )

            st_res(iPoint,:) = st_res(iPoint,:) +  st_res_pnt(1,:)
          end if ! point part of subtree
        end do !iPoint
      end if ! global mesh
    end do stLoop

    do iPoint = 1, nPnts
      res( (iPoint-1)*fun%nComponents + 1 : iPoint*fun%nComponents ) &
        & = st_res(iPoint, :)
    end do

  end subroutine evaluate_add_spacetime_vectorByCoordinate
  ! *************************************************************************** !


  ! *************************************************************************** !
  !> This subroutine returns first spacetime function on overlapping
  !! regions fo given coordinate points for vectorial variable.
  !! Spacetime function are evaluated only on the coordinate which
  !! belong to st_fun shape.
  recursive subroutine evaluate_first_spacetime_vectorByCoordinate( fun, &
    & varsys, point, time, tree, nPnts, res                              )
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
    ! --------------------------------------------------------------------------
    integer :: iStFun, iPoint, elemPos, pos
    logical :: found
    type(tem_st_fun_listElem_type), pointer :: fPtr
    integer :: coord(4)
    integer(kind=long_k):: treeID
    real(kind=rk) :: st_res(1, fun%nComponents)
    ! --------------------------------------------------------------------------

    !call C_F_POINTER( fun%method_Data, fPtr )
    ! avoid compiler warning about unused varsys
    call C_F_POINTER( varsys%method%val(fun%mypos)%method_Data, fPtr )

    res = 0.0_rk

    ! We have to loop over all requested points. For each point we have to
    ! figure out which is the first space time function that can fulfill the
    ! request. Thus a vectorized access is not possible here.
    do iPoint = 1, nPnts
      found = .false.
      ! The space time function has a custom shape, thus we have to check
      ! whether the coordinates are covered by the space time function.
      ! As we are dealing with vectorial access here, we have to loop over all
      ! requested points.
      stLoop: do iStFun = 1, fPtr%nVals
        if(fPtr%val(iStFun)%subTree%useGlobalMesh) then
          found = .true.
        else
          ! The space time function has a custom shape, thus we have to check
          ! whether the coordinates are covered by the space time function.
          coord =  tem_CoordOfReal(tree, point(iPoint,:), tree%global%maxLevel)
          treeId = tem_IdOfCoord(coord)
          elemPos = tem_PosofId(treeId, tree%treeID)
          pos = tem_PositionInSorted(                      &
            & me    = fPtr%val(iStFun)%subTree%map2global, &
            & val   = elemPos                              )
          ! When the elemenet is covered by the space time function's shape, we
          ! evaluate it at the given coordinates.
          if (pos > 0) then
            found = .true.
          end if ! point part of subtree
        end if ! global mesh

        ! When the first space time function was able to fulfill the request, we
        ! stop looping over them, because we only want one value.
        if (found) then
          st_res = tem_spacetime_for(                            &
            & me      = fPtr%val(iStFun),                        &
            & coord   = reshape( source = (/ point(iPoint,1),    &
            &                                point(iPoint,2),    &
            &                                point(iPoint,3) /), &
            &                    shape  = (/ 1, 3 /) ),          &
            & time    = time,                                    &
            & n       = 1,                                       &
            & nComp   = fun%nComponents                          )

          res( (iPoint-1)*fun%nComponents + 1 : iPoint*fun%nComponents ) &
            & = st_res(1, :)
          exit stLoop
        end if ! found
      end do stLoop
    end do !iVal

  end subroutine evaluate_first_spacetime_vectorByCoordinate
  ! *************************************************************************** !


  ! *************************************************************************** !
  !> This subroutine directly evaluate spacetime function for given element
  !! position in global TreeID list on 1st st_fun assuming global shape.
  recursive subroutine evaluate_FOAG_spacetime_scalarByTreeID( fun, varsys, &
    & elempos, time, tree, nElems, nDofs, res                               )
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
    ! --------------------------------------------------------------------------
    type(tem_st_fun_listElem_type), pointer :: fPtr
    ! -------------------------------------------------------------------------!
    !call C_F_POINTER( fun%method_Data, fPtr )
    ! avoid compiler warning about unused varsys
    call C_F_POINTER( varsys%method%val(fun%mypos)%method_Data, fPtr )

    ! assuming nDofs = 1
    if (nDofs > 1) then
      write(logUnit(1),*) 'Spacetime function does not support nDofs>1'
      call tem_abort()
    end if

    res(:) = tem_spacetime_for( me      = fPtr%val(1),             &
         &                      treeIDs = tree%treeID(elemPos(:)), &
         &                      tree    = tree,                    &
         &                      time    = time,                    &
         &                      n       = nElems                   )

  end subroutine evaluate_FOAG_spacetime_scalarByTreeID
  ! *************************************************************************** !


  ! *************************************************************************** !
  !> Call the spacetime function for scalar variable,
  !! which are stored in method_data which intern
  !! points to tem_st_fun_listelem_type.
  !! this spacetime function can be used as a analytical solution for reference.
  !! note: ndofs is not used by this routine. so, to evaluate spacetime function
  !! for ndofs in an element, use tem_varsys_proc_point interface.
  !!
  !! \verbatim
  !! -- in lua file, one can define it as following:
  !! variable = {{
  !!   name = 'reference',
  !!   ncomponents = 1,
  !!   vartype = st_fun,
  !!   st_fun =  luafun
  !!   },
  !!   ...
  !!  }
  !! \endverbatim
  !!
  !! the interface has to comply to the abstract interface
  !! tem_varsys_module#tem_varsys_proc_element.
  !!
  recursive subroutine evaluate_add_spacetime_scalarByTreeID( fun, varsys, &
    & elempos, time, tree, nElems, nDofs, res                              )
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
    ! --------------------------------------------------------------------------
    integer :: iStFun, iElem, pos, oldPos
    type(tem_st_fun_listElem_type), pointer :: fPtr
    ! element wise output
    real(kind=rk) :: st_res(1)
    integer(kind=long_k) :: treeID(1)
    ! -------------------------------------------------------------------------!
    !call C_F_POINTER( fun%method_Data, fPtr )

    ! avoid compiler warning about unused varsys
    call C_F_POINTER( varsys%method%val(fun%mypos)%method_Data, fPtr )

    res = 0.0_rk
    ! assuming nDofs = 1
    if (nDofs > 1) then
      write(logUnit(1),*) 'Spacetime function does not support nDofs>1'
      call tem_abort()
    end if

    do iStFun = 1, fPtr%nVals
      if (fPtr%val(iStFun)%subTree%useGlobalMesh) then
        res(:) = res(:)                                                &
         &     + tem_spacetime_for( me      = fPtr%val(iStFun),        &
         &                          treeIDs = tree%treeID(elemPos(:)), &
         &                          tree    = tree,                    &
         &                          time    = time,                    &
         &                          n       = nElems                   )
      else
        oldPos = 1
        do iElem = 1, nElems
          ! position of element matches with any of map2global then
          ! this element is part of this subTree
          pos = tem_PositionInSorted(  &
            &                 me    = fPtr%val(iStFun)%subTree%map2global, &
            &                 val   = elemPos(iElem),                      &
            &                 lower = oldPos                               )
          if (pos > 0) then
            oldPos = pos
            treeID = tree%treeID(elemPos(iElem))
            st_res = tem_spacetime_for( me      = fPtr%val(iStFun),     &
              &                         treeIDs = treeID,               &
              &                         tree    = tree,                 &
              &                         time    = time,                 &
              &                         n       = 1                     )
            res(iElem) = res(iElem) + st_res(1)
          end if ! elemPos(iElem) part of subtree
        end do ! iElem
      end if ! global mesh
    end do ! nr. st_funs

  end subroutine evaluate_add_spacetime_scalarByTreeID
  ! *************************************************************************** !

  ! *************************************************************************** !
  !> Call the spacetime function for scalar variable,
  !! which are stored in method_data which intern
  !! points to tem_st_fun_listelem_type.
  !! this spacetime function can be used as a analytical solution for reference.
  !! note: ndofs is not used by this routine. so, to evaluate spacetime function
  !! for ndofs in an element, use tem_varsys_proc_point interface.
  !!
  !! \verbatim
  !! -- in lua file, one can define it as following:
  !! variable = {{
  !!   name = 'reference',
  !!   ncomponents = 1,
  !!   vartype = st_fun,
  !!   st_fun =  luafun
  !!   },
  !!   ...
  !!  }
  !! \endverbatim
  !!
  !! the interface has to comply to the abstract interface
  !! tem_varsys_module#tem_varsys_proc_element.
  !!
  recursive subroutine evaluate_first_spacetime_scalarByTreeID( fun, varsys, &
    & elempos, time, tree, nElems, nDofs, res                                )
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
    ! --------------------------------------------------------------------------
    integer :: iStFun, iElem, pos, oldPos
    logical :: found = .false.
    type(tem_st_fun_listElem_type), pointer :: fPtr
    ! element wise output
    real(kind=rk) :: st_res(1)
    integer(kind=long_k) :: treeID(1)
    ! -------------------------------------------------------------------------!
    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk
    ! assuming nDofs = 1
    if (nDofs > 1) then

    ! avoid compiler warning about unused varsys
      write(logunit(1),*) 'Variable:', trim(varsys%varname%val(fun%mypos))
      write(logunit(1),*) 'tree%nElems:', tree%nElems
      write(logUnit(1),*) 'Spacetime function does not support nDofs>1'
      call tem_abort()
    end if
    do iElem = 1, nElems
      found = .false.
      stLoop: do iStFun = 1, fPtr%nVals
        ! If the space time function has a global shape, we can leave out the
        ! search for the elemPos.
        if(fPtr%val(iStFun)%subTree%useGlobalMesh) then
          found = .true.
        else
          oldPos = 1
          ! position of element matches with any of map2global then
          ! this element is part of this subTree
          pos = tem_PositionInSorted(                                      &
            &                 me    = fPtr%val(iStFun)%subTree%map2global, &
            &                 val   = elemPos(iElem),                      &
            &                 lower = oldPos                               )
          if (pos > 0) then
            oldPos = pos
            found = .true.
          end if ! elemPos(iElem) part of subtree
        end if ! global mesh

        ! We found the one value we were looking for, so we can stop now
        if (found) then
          treeID = tree%treeID(elemPos(iElem))
          st_res = tem_spacetime_for( me      = fPtr%val(iStFun),     &
            &                         treeIDs = treeID,               &
            &                         tree    = tree,                 &
            &                         time    = time,                 &
            &                         n       = 1                     )
          res(iElem) = st_res(1)
          exit stLoop
        end if !found
      end do stLoop
    end do ! iElem
  end subroutine evaluate_first_spacetime_scalarByTreeID
  ! *************************************************************************** !


  ! *************************************************************************** !
  !> This subroutine directly evaluate spacetime function for given element
  !! position in global TreeID list on 1st st_fun assuming global shape for
  !! vectorial variable.
  recursive subroutine evaluate_FOAG_spacetime_vectorByTreeID( fun, varsys, &
    & elempos, time, tree, nElems, nDofs, res                               )
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
    ! -------------------------------------------------------------------------!
    type(tem_st_fun_listElem_type), pointer :: fPtr
    integer :: iElem
    real(kind=rk) :: st_res(nElems, fun%nComponents)
    ! -------------------------------------------------------------------------!
    call C_F_POINTER( fun%method_Data, fPtr )
    ! assuming nDofs = 1
    if (nDofs > 1) then
      write(logunit(1),*) 'Variable:', trim(varsys%varname%val(fun%mypos))
      write(logunit(1),*) 'tree%nElems:', tree%nElems
      write(logUnit(1),*) 'Spacetime function does not support nDofs>1'
      call tem_abort()
    end if

    st_res = tem_spacetime_for( me      = fPtr%val(1),             &
      &                         treeIDs = tree%treeID(elemPos(:)), &
      &                         tree    = tree,                    &
      &                         time    = time,                    &
      &                         nComp   = fun%nComponents,         &
      &                         n       = nElems                   )


    do iElem = 1, nElems
      res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents )             &
        & = st_res(iElem, :)
    end do

  end subroutine evaluate_FOAG_spacetime_vectorByTreeID
  ! *************************************************************************** !


  ! *************************************************************************** !
  !> Call the spacetime function for vector variable,
  !! which are stored in method_data which intern
  !! points to tem_st_fun_listElem_type.
  !! This spacetime function can be used as a analytical solution for reference.
  !! NOTE: nDofs is not used by this routine. So, to evaluate spacetime function
  !! for nDofs in an element, use tem_varSys_proc_point interface.
  !!
  !! \verbatim
  !! -- in lua file, one can define it as following:
  !! variable = {{
  !!   name = 'reference',
  !!   ncomponents = 3,
  !!   vartype = st_fun,
  !!   st_fun =  luaFun
  !!   },
  !!   ...
  !!  }
  !! \endverbatim
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_element.
  !!
  recursive subroutine evaluate_add_spacetime_vectorByTreeID( fun, varsys, &
    & elempos, time, tree, nElems, nDofs, res                              )
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
    ! -------------------------------------------------------------------------!
    integer :: iStFun, iElem, pos, oldPos
    type(tem_st_fun_listElem_type), pointer :: fPtr
    ! element wise output
    real(kind=rk), allocatable :: st_res(:,:)
    real(kind=rk) :: st_res_elem(1, fun%nComponents)
    integer(kind=long_k) :: treeID(1)
    ! -------------------------------------------------------------------------!
    !write(*,*) 'derive variable label: ', trim(varSys%varname%val(fun%myPos))
    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk
    ! assuming nDofs = 1
    if (nDofs > 1) then
      write(logunit(1),*) 'Variable:', trim(varsys%varname%val(fun%mypos))
      write(logUnit(1),*) 'Spacetime function does not support nDofs>1'
      call tem_abort()
    end if

    allocate(st_res(nElems, fun%nComponents))
    st_res = 0.0_rk
    do iStFun = 1, fPtr%nVals
      if (fPtr%val(iStFun)%subTree%useGlobalMesh ) then
        !write(*,*) 'global mesh'
        st_res = st_res +                                         &
         &  tem_spacetime_for( me      = fPtr%val(iStFun),        &
         &                     treeIDs = tree%treeID(elemPos(:)), &
         &                     tree    = tree,                    &
         &                     time    = time,                    &
         &                     nComp   = fun%nComponents,         &
         &                     n       = nElems                   )
      else
        !write(*,*) 'subTree Mesh'
        oldPos = 1
        do iElem = 1, nElems
          ! position of element matches with any of map2global then
          ! this element is part of this subTree
          pos = tem_PositionInSorted(  &
            &                 me    = fPtr%val(iStFun)%subTree%map2global, &
            &                 val   = elemPos(iElem),                      &
            &                 lower = oldPos                               )
          ! pos>0 if elemPos(iElem) is part of subTree
          if (pos > 0) then
            oldPos = pos
            treeID = tree%treeID(elemPos(iElem))
            st_res_elem =                                           &
              &  tem_spacetime_for( me      = fPtr%val(iStFun),     &
              &                     treeIDs = treeID,               &
              &                     tree    = tree,                 &
              &                     time    = time,                 &
              &                     nComp   = fun%nComponents,      &
              &                     n       = 1                     )
            st_res(iElem,:) = st_res(iElem,:) +  st_res_elem(1,:)
          end if
        end do
      end if ! global mesh
    end do ! iStFun

    do iElem = 1, nElems
      res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents )             &
        & = st_res(iElem, :)
    end do

  end subroutine evaluate_add_spacetime_vectorByTreeID
  ! *************************************************************************** !


  ! *************************************************************************** !
  !> Call the spacetime function for vector variable,
  !! which are stored in method_data which intern
  !! points to tem_st_fun_listElem_type.
  !! This spacetime function can be used as a analytical solution for reference.
  !! NOTE: nDofs is not used by this routine. So, to evaluate spacetime function
  !! for nDofs in an element, use tem_varSys_proc_point interface.
  !!
  !! \verbatim
  !! -- in lua file, one can define it as following:
  !! variable = {{
  !!   name = 'reference',
  !!   ncomponents = 3,
  !!   vartype = st_fun,
  !!   st_fun =  luaFun
  !!   },
  !!   ...
  !!  }
  !! \endverbatim
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_element.
  !!
  recursive subroutine evaluate_first_spacetime_vectorByTreeID( fun, varsys, &
    & elempos, time, tree, nElems, nDofs, res                                )
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
    ! -------------------------------------------------------------------------!
    integer :: iStFun, iElem, pos, oldPos
    logical :: found
    type(tem_st_fun_listElem_type), pointer :: fPtr
    ! element wise output
    real(kind=rk) :: st_res(1, fun%nComponents)
    integer(kind=long_k) :: treeID(1)
    ! -------------------------------------------------------------------------!
    !write(*,*) 'derive variable label: ', trim(varSys%varname%val(fun%myPos))
    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk
    ! assuming nDofs = 1
    if (nDofs > 1) then
      write(logunit(1),*) 'Variable:', trim(varsys%varname%val(fun%mypos))
      write(logUnit(1),*) 'Spacetime function does not support nDofs>1'
      call tem_abort()
    end if

    do iElem = 1, nElems
      found = .false.
      stLoop: do iStFun = 1, fPtr%nVals
        if (fPtr%val(iStFun)%subTree%useGlobalMesh) then
          found = .true.
        else
          oldPos = 1
          ! position of element matches with any of map2global then
          ! this element is part of this subTree
          pos = tem_PositionInSorted(                      &
            & me    = fPtr%val(iStFun)%subTree%map2global, &
            & val   = elemPos(iElem),                      &
            & lower = oldPos                               )
          ! pos>0 if elemPos(iElem) is part of subTree
          if (pos > 0) then
            oldPos = pos
            found = .true.
          end if
        end if ! global mesh

        ! We found the one value we were looking for, so we can stop now
        if (found) then
          treeID = tree%treeID(elemPos(iElem))
          st_res = tem_spacetime_for( me      = fPtr%val(iStFun), &
            &                         treeIDs = treeID,           &
            &                         tree    = tree,             &
            &                         time    = time,             &
            &                         nComp   = fun%nComponents,  &
            &                         n       = 1                 )

          res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents ) &
            & = st_res(1,:)

          exit stLoop
        end if ! found

      end do stLoop
    end do !iElem

  end subroutine evaluate_first_spacetime_vectorByTreeID
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine process instring and store information in spacetime function.
  !!
  !! Abort if same spacetime function is used as boundaries and sources.
  recursive subroutine set_params_spacetime( fun, varSys, instring )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Input string with parameter to set in method_data
    character(len=*), intent(in) :: instring
    ! -------------------------------------------------------------------- !
    type(flu_state) :: conf
    type(tem_st_fun_listElem_type), pointer :: fPtr
    logical :: isSurface_loc
    integer :: iStFun
    integer :: iError
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    call open_config_chunk(L = conf, chunk = instring)
    call aot_get_val( L       = conf,          &
      &               key     = 'isSurface',   &
      &               val     = isSurface_loc, &
      &               ErrCode = iError         )

    if (btest(iError,  aoterr_Fatal)) then
      write(logUnit(1),*) 'Error: In set_params to load "isSurface" '//&
        &                 'for spacetime variable ', &
        &                 trim(varSys%varname%val(fun%mypos))
      if (btest(iError, aoterr_WrongType)) then
        write(logUnit(1),*) '       "isSurface" is wrong type'
      end if
      if (btest(iError, aoterr_NonExistent)) then
        write(logUnit(1),*) '       "isSurface" does not exist'
      end if
      call tem_abort()
    end if

    ! KM Store bc_id and field id only for apesmate coupling
    do iStFun = 1, fPtr%nVals
      select case (trim(fPtr%val(iStFun)%fun_kind))
      case ('apesmate')
        if (fPtr%isSurface == -1) then
          if (isSurface_loc) then
            fPtr%isSurface = 0
          else
            fPtr%isSurface = 1
          end if
        else
          write(logUnit(10),*) 'Warning: ' &
            &                  // 'Set params stfun: isSurface is already set'
          if (isSurface_loc) then
            if (fPtr%isSurface == 1) then
              write(logUnit(1),*) 'Error: In set params for spacetime function'
              write(logUnit(1),*) '       This st_fun is already used for' &
                &                 // ' volume'
              write(logUnit(1),*) '       so cannot be used for surface.'
              call tem_abort()
            end if
          else
            if (fPtr%isSurface == 0) then
              write(logUnit(1),*) 'Error: In set params for spacetime function'
              write(logUnit(1),*) '       This st_fun is already used for' &
                &                 // ' surface'
              write(logUnit(1),*) '       so cannot be used for volume.'
              call tem_abort()
            end if
          end if
        end if

        ! store surface or volume information in coupling type
        fPtr%val(iStFun)%aps_coupling%isSurface = fPtr%isSurface

      end select
    end do

  end subroutine set_params_spacetime
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine process instring and return string with requested info
  !! from spacetime function
  recursive subroutine get_params_spacetime(fun, varSys, instring, outstring)
    ! -------------------------------------------------------------------- !
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
    ! -------------------------------------------------------------------- !
    select case (trim(instring))
    case ('vartype')
      ! return vartype
      outstring = 'st_fun'
    case default
      outstring = trim(varsys%varname%val(fun%mypos)) // '_UNKNOWN'
    end select
  end subroutine get_params_spacetime
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine stores provided points in method_data of spacetime_listElem
  !! and return the indices of points or evaluated value in the growing array.
  !! If spacetime function is time-independent then pre-compute values
  !! and store in growing array of evalVal in tem_pointData_type.
  recursive subroutine setup_indices_spacetime( fun, varSys, point, &
    & offset_bit, iLevel, tree, nPnts, idx                          )
    !---------------------------`-----------------------------------------------!
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
    !--------------------------------------------------------------------------!
    ! -------------------------------------------------------------------------!
    type(tem_st_fun_listElem_type), pointer :: fPtr
    integer :: iStFun, iPnt, posInTree, nUniquePnts, iVar
    character :: offset_bit_local
    logical, allocatable :: storePnt(:), storeOffsetBit(:), storeVal(:)
    integer :: elemPos
    integer(kind=long_k) :: treeID
    logical :: addPoint, wasAdded
    real(kind=rk) :: uniquePnts(nPnts,3)
    ! -------------------------------------------------------------------------!


    call C_F_POINTER( fun%method_Data, fPtr )

    allocate(storePnt(fPtr%nVals))
    allocate(storeVal(fPtr%nVals))
    allocate(storeOffsetBit(fPtr%nVals))

    ! Store points only for time dependent spacetime functions
    do iStFun = 1, fPtr%nVals
      select case (trim(fPtr%val(iStFun)%fun_kind))
      case ('none', 'const')
        ! time independent and no need to store points
        storePnt(iStFun) = .false.
        storeOffsetBit(iStFun) = .false.
        storeVal(iStFun) = .false.
      case ('combined')
        ! spatial is time independent so compute spatial value and store it.
        storePnt(iStFun) = .false.
        storeOffsetBit(iStFun) = .false.
        storeVal(iStFun) = .true.

      case ('lua_fun', 'miescatter_displacementfieldz',               &
        &   'miescatter_magneticfieldx', 'miescatter_magneticfieldy', &
        &   'cylindrical_wave')
        storePnt(iStFun) = .true.
        storeOffsetBit(iStFun) = .false.
        storeVal(iStFun) = .false.
      case ('apesmate')
        storePnt(iStFun) = .true.
        storeOffsetBit(iStFun) = (fPtr%val(iStFun)%aps_coupling%isSurface==0)
        storeVal(iStFun) = .false.
      case ('precice')
        storePnt(iStFun) = .true.
        storeOffsetBit(iStFun) = .true.
        storeVal(iStFun) = .true.
        ! for writing to precice, we need the position of the variables in
        ! the variable system
        do iVar = 1, fPtr%val(iStFun)%precice_coupling%writeVar%nVars
          fPtr%val(iStFun)%precice_coupling%writeVar%varPos(iVar) =     &
             PositionOfVal( me = varSys%varName,                        &
               &            val = fPtr%val(iStFun)%precice_coupling     &
               &                                  %writeVar%names(iVar) )
             if (fPtr%val(iStFun)%precice_coupling%writeVar%varPos(iVar) == 0) &
               & then
               write(*,*) 'position in Varsys for writing variable ',          &
                 & trim(fPtr%val(iStFun)%precice_coupling%writeVar%names(iVar)),&
                 & ' not found'
                 call tem_abort()
             end if
        end do

        ! same for reading from precice
        do iVar = 1, fPtr%val(iStFun)%precice_coupling%readVar%nVars
          fPtr%val(iStFun)%precice_coupling%readVar%varPos(iVar) =     &
             PositionOfVal( me = varSys%varName,                       &
               &            val = fPtr%val(iStFun)%precice_coupling    &
               &                                  %readVar%names(iVar) )
             if (fPtr%val(iStFun)%precice_coupling%readVar%varPos(iVar) == 0) &
               & then
               write(*,*) 'position in Varsys for reading ',                   &
                 & trim(fPtr%val(iStFun)%precice_coupling%readVar%names(iVar)),&
                 & 'not found'
                 call tem_abort()
             end if
        end do

      case default
        write(logUnit(1),*)'ERROR: Unknown spatial function in '// &
          &                'setup_indices.'
        call tem_abort()
      end select
    end do !iStFun

    ! initialize index with zero to identify points which does not
    ! belong to subTree
    idx = 0

    ! number of unique points added
    nUniquePnts = 0
    do iPnt = 1, nPnts
      addPoint = .false.

      ! get treeID from globalmaxLevel since points from ghost elements are
      ! also passed to this routine
      treeID = tem_IdOfCoord( tem_CoordOfReal(tree, point(iPnt, :)) )
      elemPos = tem_PosOfId( treeID, tree%treeID )
      !if (elemPos == 0)  then
       ! call tem_abort('Error: treeID not found in st-fun setup_indices')
     ! end if

      ! if any spacetime function has useGlobalMesh then store all points
      ! else store only points which belong to subTree of spacetime variable.
      if ( any(fPtr%val(:)%subTree%useGlobalMesh) ) then
        addPoint = .true.
      else
        stFunLoop: do iStFun = 1, fPtr%nVals
          posInTree = tem_PositionInSorted(                                &
            &                 me    = fPtr%val(iStFun)%subTree%map2global, &
            &                 val   = elemPos                              )
          if (posInTree > 0) then
            addPoint = .true.
            exit stFunLoop
          end if
        end do stFunLoop
      end if !use global mesh

      if (addPoint) then
        ! use center offset bit as default
        if (present(offset_bit)) then
          offset_bit_local = offset_bit(iPnt)
        else
          offset_bit_local = qOffset_inChar(q000)
        end if

        ! append point, offset_bit and elemPos to pointData type
        call append(me             = fPtr%pntData%pntLvl(iLevel), &
          &         point          = point(iPnt,:),               &
          &         storePnt       = any(storePnt),               &
          &         offset_bit     = offset_bit_local,            &
          &         storeOffsetBit = any(storeOffsetBit),         &
          &         elemPos        = elemPos,                     &
          &         tree           = tree,                        &
          &         pos            = idx(iPnt),                   &
          &         wasAdded       = wasAdded                     )

        if (wasAdded) then
          nUniquePnts = nUniquePnts + 1
          uniquePnts(nUniquePnts,:) = point(iPnt,:)
        end if

      end if ! add point
    end do !iPnt

    if (any(storePnt)) call truncate(fPtr%pntData%pntLvl(iLevel)%grwPnt)
    if (any(storeOffsetBit)) &
      & call truncate(fPtr%pntData%pntLvl(iLevel)%offset_bit)

    deallocate(storePnt)
    deallocate(storeOffsetBit)
    call truncate(fPtr%pntData%pntLvl(iLevel)%treeID)
    call truncate(fPtr%pntData%pntLvl(iLevel)%elemPos)

    !!! Store spatial value for unique points depends on stFun type
    if ( any(storeVal) .and. nUniquePnts > 0 ) then
      do iStFun = 1, fPtr%nVals
        select case (trim(fPtr%val(iStFun)%fun_kind))
        case ('combined')
          if (fun%nComponents == 1) then
            call tem_spatial_storeVal( me     = fPtr%val(iStFun)%spatial,    &
              &                        coord  = uniquePnts(1:nUniquePnts,:), &
              &                        nVals  = nUniquePnts,                 &
              &                        iLevel = iLevel                       )
          else
            call tem_spatial_storeVal( me     = fPtr%val(iStFun)%spatial,    &
              &                        coord  = uniquePnts(1:nUniquePnts,:), &
              &                        nVals  = nUniquePnts,                 &
              &                        iLevel = iLevel,                      &
              &                        nComps = fun%nComponents              )
          end if
        end select
      end do !iStFun
    end if !store spatial value
    deallocate(storeVal)

  end subroutine setup_indices_spacetime
  ! *************************************************************************** !


  ! ************************************************************************** !
  !> This routine returns value at the given index. If value is pre-computed
  !! and value at given index is returned else value is computing on the
  !! points for given index.
  !!
  !! If a variable has more than one spacetime function, then result of 1st
  !! spacetime function is returned with assumed global shape
  recursive subroutine get_valOfIndex_FOAG_scalar_spacetime( fun, varSys,  &
    &                                                        time, iLevel, &
    &                                                        idx, idxLen,  &
    &                                                        nVals, res    )
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
    ! -------------------------------------------------------------------------!
    type(tem_st_fun_listElem_type), pointer :: fPtr
    ! -------------------------------------------------------------------------!
    call C_F_POINTER( fun%method_Data, fPtr )

    ! check if number of index are the same as number of values asked for
    call tem_varsys_check_inArgs( fun, varSys, time, iLevel, idx, idxLen, &
      &    nVals, label = 'get_valOfIndex_FOAG_scalar_spacetime'          )

    res(:) = tem_spacetime_for( me     = fPtr%val(1),                &
      &                         grwPnt = fPtr%pntData%pntLvl(iLevel) &
      &                                      %grwPnt,                &
      &                         idx    = idx,                        &
      &                         nVals  = nVals,                      &
      &                         iLevel = ilevel,                     &
      &                         time   = time                        )
  end subroutine get_valOfIndex_FOAG_scalar_spacetime
  ! ************************************************************************ !


  ! *************************************************************************** !
  !> This routine returns value at the given index. If value is pre-computed
  !! and value at given index is returned else value is computing on the
  !! points for given index.
  !!
  !! If a variable has more than one spacetime function, then result of 1st
  !! spacetime function is returned with assumed global shape
  !!
  !! For vectorial variable.
  recursive subroutine get_valOfIndex_FOAG_vector_spacetime( fun, varSys,  &
    &                                                        time, iLevel, &
    &                                                        idx, idxLen,  &
    &                                                        nVals, res    )
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
    ! -------------------------------------------------------------------------!
    type(tem_st_fun_listElem_type), pointer :: fPtr
    integer :: iVal
    real(kind=rk) :: st_res(nVals, fun%nComponents)
    ! -------------------------------------------------------------------------!
    call C_F_POINTER( fun%method_Data, fPtr )

    ! check if number of index are the same as number of values asked for
    call tem_varSys_check_inArgs( fun, varSys, time, iLevel, idx, idxLen, &
      &    nVals, label = 'get_valOfIndex_FOAG_vector_spacetime'          )

    st_res = tem_spacetime_for( me     = fPtr%val(1),                &
      &                         grwPnt = fPtr%pntData%pntLvl(iLevel) &
      &                                      %grwPnt,                &
      &                         idx    = idx,                        &
      &                         nVals  = nVals,                      &
      &                         iLevel = ilevel,                     &
      &                         time   = time,                       &
      &                         nComps = fun%nComponents             )

    do iVal = 1, nVals
      res( (iVal-1)*fun%nComponents + 1 : iVal*fun%nComponents ) &
        & = st_res(iVal, :)
    end do

  end subroutine get_valOfIndex_FOAG_vector_spacetime
  ! ************************************************************************ !


  ! *************************************************************************** !
  !> This routine returns value at the given index. If value is pre-computed
  !! and value at given index is returned else value is computing on the
  !! points for given index.
  !!
  !! If a variable has more than one spacetime function, then result is sum
  !! of all spacetime functions
  recursive subroutine get_valOfIndex_add_scalar_spacetime( fun, varSys, time, &
    &                                                       iLevel, idx,       &
    &                                                       idxLen, nVals, res )
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
    ! -------------------------------------------------------------------------!
    integer :: iStFun, iVal, pos, oldPos
    type(tem_st_fun_listElem_type), pointer :: fPtr
    ! element wise output
    real(kind=rk) :: st_res(1)
    ! -------------------------------------------------------------------------!
    call C_F_POINTER( fun%method_Data, fPtr )

    ! check if number of index are the same as number of values asked for
    call tem_varsys_check_inArgs( fun, varSys, time, iLevel, idx, idxLen, &
      &    nVals, label = 'get_valOfIndex_add_scalar_spacetime'           )

    res = 0.0_rk
    do iStFun = 1, fPtr%nVals
      if (fPtr%val(iStFun)%subTree%useGlobalMesh) then
        res(:) = res(:)                                                  &
          &    + tem_spacetime_for( me     = fPtr%val(iStFun),           &
          &                         grwPnt = fPtr%pntData%pntLvl(iLevel) &
          &                                      %grwPnt,                &
          &                         idx    = idx,                        &
          &                         nVals  = nVals,                      &
          &                         iLevel = ilevel,                     &
          &                         time   = time                        )
      else
        oldPos = 1
        do iVal = 1, nVals

          if (idx(iVal) > 0) then
            ! position of element matches with any of map2global then
            ! this element is part of this subTree
            pos = tem_PositionInSorted(  &
              &                 me    = fPtr%val(iStFun)%subTree%map2global, &
              &                 val   = fPtr%pntData%pntLvl(iLevel)          &
              &                             %elemPos%val(idx(iVal)),         &
              &                 lower = oldPos                               )
            if (pos > 0) then
              oldPos = pos
              st_res = tem_spacetime_for( me     = fPtr%val(iStFun),           &
                &                         grwPnt = fPtr%pntData%pntLvl(iLevel) &
                &                                      %grwPnt,                &
                &                         idx    = (/ idx(iVal) /),            &
                &                         nVals  = 1,                          &
                &                         iLevel = ilevel,                     &
                &                         time   = time                        )
              res(iVal) = res(iVal) + st_res(1)
            end if ! elemPos(iElem) part of subtree
          end if !idx > 0
        end do !iVal
      end if ! global mesh
    end do ! nr. st_funs

  end subroutine get_valOfIndex_add_scalar_spacetime
  ! ************************************************************************ !


  ! *************************************************************************** !
  !> This routine returns value at the given index. If value is pre-computed
  !! and value at given index is returned else value is computing on the
  !! points for given index.
  !!
  !! If a variable has more than one spacetime function, then result is sum
  !! of all spacetime functions
  !!
  !! For vectorial variable
  recursive subroutine get_valOfIndex_add_vector_spacetime( fun, varSys, time, &
    &                                                       iLevel, idx,       &
    &                                                       idxLen, nVals, res )
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
    integer :: iStFun, iVal, pos, oldPos
    type(tem_st_fun_listElem_type), pointer :: fPtr
    real(kind=rk) :: st_res(nVals, fun%nComponents)
    real(kind=rk) :: st_res_elem(1, fun%nComponents)
    ! -------------------------------------------------------------------------!
    call C_F_POINTER( fun%method_Data, fPtr )

    ! check if number of index are the same as number of values asked for
    call tem_varSys_check_inArgs( fun, varSys, time, iLevel, idx, idxLen, &
      &    nVals, label = 'get_valOfIndex_add_vector_spacetime'           )

    st_res = 0.0_rk
    do iStFun = 1, fPtr%nVals
      if (fPtr%val(iStFun)%subTree%useGlobalMesh) then
        st_res = st_res                                                  &
          &    + tem_spacetime_for( me     = fPtr%val(iStFun),           &
          &                         grwPnt = fPtr%pntData%pntLvl(iLevel) &
          &                                      %grwPnt,                &
          &                         idx    = idx,                        &
          &                         nVals  = nVals,                      &
          &                         iLevel = ilevel,                     &
          &                         time   = time,                       &
          &                         nComps = fun%nComponents             )
      else
        oldPos = 1
        do iVal = 1, nVals
          if (idx(iVal) > 0) then
            ! position of element matches with any of map2global then
            ! this element is part of this subTree
            pos = tem_PositionInSorted(  &
              &                 me    = fPtr%val(iStFun)%subTree%map2global, &
              &                 val   = fPtr%pntData%pntLvl(iLevel)          &
              &                             %elemPos%val(idx(iVal)),         &
              &                 lower = oldPos                               )
            if (pos > 0) then
              oldPos = pos
              st_res_elem =                                               &
                & tem_spacetime_for( me     = fPtr%val(iStFun),           &
                &                    grwPnt = fPtr%pntData%pntLvl(iLevel) &
                &                                 %grwPnt,                &
                &                    idx    = (/ idx(iVal) /),            &
                &                    nVals  = 1,                          &
                &                    iLevel = ilevel,                     &
                &                    time   = time,                       &
                &                    nComps = fun%nComponents             )
              st_res(iVal,:) = st_res(iVal,:) + st_res_elem(1,:)
            end if ! elemPos(iElem) part of subtree
          end if
        end do !iVal
      end if ! global mesh
    end do ! nr. st_funs

    do iVal = 1, nVals
      res( (iVal-1)*fun%nComponents + 1 : iVal*fun%nComponents ) &
        & = st_res(iVal, :)
    end do

  end subroutine get_valOfIndex_add_vector_spacetime
  ! *************************************************************************** !


  ! *************************************************************************** !
  !> This routine returns value at the given index. If value is pre-computed
  !! and value at given index is returned else value is computing on the
  !! points for given index.
  !!
  !! If a variable has more than one spacetime function, then result is
  !! first spacetime function in the overlapping region
  recursive subroutine get_valOfIndex_first_scalar_spacetime( fun, varSys,  &
    &                                                         time, iLevel, &
    &                                                         idx, idxLen,  &
    &                                                         nVals, res    )
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
    ! -------------------------------------------------------------------------!
    integer :: iStFun, iVal, pos
    type(tem_st_fun_listElem_type), pointer :: fPtr
    ! element wise output
    real(kind=rk) :: st_res(1)
    logical :: found
    ! -------------------------------------------------------------------------!
    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk

    ! check if number of index are the same as number of values asked for
    call tem_varSys_check_inArgs( fun, varSys, time, iLevel, idx, idxLen, &
      &    nVals, label = 'get_valOfIndex_first_scalar_spacetime'         )

    do iVal = 1, nVals
      found = .false.
      stLoop: do iStFun = 1, fPtr%nVals
        if (fPtr%val(iStFun)%subTree%useGlobalMesh) then
          found = .true.
        else
          if (idx(iVal) > 0) then
            ! position of element matches with any of map2global then
            ! this element is part of this subTree
            pos = tem_PositionInSorted(  &
              &                 me    = fPtr%val(iStFun)%subTree%map2global, &
              &                 val   = fPtr%pntData%pntLvl(iLevel)          &
              &                             %elemPos%val(idx(iVal))          )

            found = (pos > 0)
          end if !idx>0
        end if ! global mesh
        ! We found the one value we were looking for, so we can stop now
        if (found) then
          st_res = tem_spacetime_for( me     = fPtr%val(iStFun),           &
            &                         grwPnt = fPtr%pntData%pntLvl(iLevel) &
            &                                      %grwPnt,                &
            &                         idx    = (/ idx(iVal) /),            &
            &                         nVals  = 1,                          &
            &                         iLevel = ilevel,                     &
            &                         time   = time                        )
          res(iVal) = st_res(1)
          exit stLoop
        end if
      end do stLoop
    end do !iVal

  end subroutine get_valOfIndex_first_scalar_spacetime
  ! *************************************************************************** !


  ! *************************************************************************** !
  !> This routine returns value at the given index. If value is pre-computed
  !! and value at given index is returned else value is computing on the
  !! points for given index.
  !!
  !! If a variable has more than one spacetime function, then result is
  !! first spacetime function in the overlapping region.
  !!
  !! For vectorial variable
  recursive subroutine get_valOfIndex_first_vector_spacetime( fun, varSys,  &
    &                                                         time, iLevel, &
    &                                                         idx, idxLen,  &
    &                                                         nVals, res    )
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
    ! -------------------------------------------------------------------------!
    integer :: iStFun, iVal, pos
    type(tem_st_fun_listElem_type), pointer :: fPtr
    ! element wise output
    real(kind=rk) :: st_res(1, fun%nComponents)
    logical :: found
    ! -------------------------------------------------------------------------!
    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk

    ! check if number of index are the same as number of values asked for
    call tem_varsys_check_inArgs( fun, varSys, time, iLevel, idx, idxLen, &
      &    nVals, label = 'get_valOfIndex_first_vector_spacetime'         )

    do iVal = 1, nVals
      found = .false.
      stLoop: do iStFun = 1, fPtr%nVals
        if (fPtr%val(iStFun)%subTree%useGlobalMesh) then
          found = .true.
        else
          if (idx(iVal) > 0) then
            ! position of element matches with any of map2global then
            ! this element is part of this subTree
            pos = tem_PositionInSorted(  &
              &                 me    = fPtr%val(iStFun)%subTree%map2global, &
              &                 val   = fPtr%pntData%pntLvl(iLevel)          &
              &                             %elemPos%val(idx(iVal))          )

            found = (pos > 0)
          end if
        end if ! global mesh

        ! We found the one value we were looking for, so we can stop now
        if (found) then
          st_res = tem_spacetime_for( me     = fPtr%val(iStFun),           &
            &                         grwPnt = fPtr%pntData%pntLvl(iLevel) &
            &                                      %grwPnt,                &
            &                         idx    = (/ idx(iVal) /),            &
            &                         nVals  = 1,                          &
            &                         iLevel = ilevel,                     &
            &                         time   = time,                       &
            &                         nComps = fun%nComponents             )

          res( (iVal-1)*fun%nComponents + 1 : iVal*fun%nComponents ) &
            & = st_res(1, :)

          exit stLoop
        end if ! found
      end do stLoop
    end do !iVal

  end subroutine get_valOfIndex_first_vector_spacetime
  ! *************************************************************************** !


end module tem_spacetime_var_module
! ****************************************************************************** !
