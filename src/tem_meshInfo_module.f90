! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013, 2015, 2017-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016-2017 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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
! ---------------------------------------------------------------------------- !
!> author: Simon Zimny
!! This module provides the functionality to derive mesh dependent information
!! for tracking output.
module tem_meshInfo_module
  use, intrinsic :: iso_c_binding, only: c_null_ptr

  ! include treelm modules
  use env_module,              only: rk, long_k, labelLen
  use tem_aux_module,          only: tem_abort
  use tem_varSys_module,       only: tem_varSys_type, tem_varSys_op_type, &
    &                                tem_varSys_append_derVar,            &
    &                                tem_varSys_proc_point,               &
    &                                tem_varSys_proc_element,             &
    &                                tem_varSys_proc_setParams,           &
    &                                tem_varSys_proc_getParams,           &
    &                                tem_varSys_proc_setupIndices,        &
    &                                tem_varSys_proc_getValOfIndex
  use tem_time_module,         only: tem_time_type
  use treelmesh_module,        only: treelmesh_type, tem_load_weights
  use tem_geometry_module,     only: tem_ElemSize
  use tem_property_module,     only: prp_hasIBM, prp_solid, prp_sendHalo
  use tem_logging_module,      only: logUnit
  use tem_topology_module,     only: tem_levelOf

  ! include aotus modules
  use aotus_module, only: flu_state

  implicit none

  private

  public :: tem_varSys_append_meshInfoVar

  contains


  ! ------------------------------------------------------------------------ !
  !> This subroutine appends the list of meshInfo variables (e.g. element
  !! volume, element volume fraction, treeID, ... )
  !!
  subroutine tem_varSys_append_meshInfoVar( varSys )
    ! -------------------------------------------------------------------- !
    !> global variable system to which meshInfoVar to be appended
    type(tem_varSys_type), intent(inout)        :: varSys
    ! -------------------------------------------------------------------- !
    integer :: iVar, addedPos, nComponents
    logical :: wasAdded
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setParams), pointer :: set_params => null()
    procedure(tem_varSys_proc_getParams), pointer :: get_params => null()
    procedure(tem_varSys_proc_setupIndices), pointer :: &
      &                                      setup_indices => null()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: &
      &                                       get_valOfIndex => null()
    ! Mesh info variables
    integer :: nMeshVars
    character(len=labelLen), allocatable ::  meshInfoVar(:)
    ! -------------------------------------------------------------------- !
    write(logUnit(5),*) 'Append meshInfo variables to varSys'

    nMeshVars = 9
    allocate(meshInfoVar(nMeshVars))
    meshInfoVar = [  'treeid       ', 'process      ', &
      &              'weight       ',                  &
      &              'elem_vol     ', 'vol_frac     ', &
      &              'level        ', 'solidified   ', &
      &              'has_ibm      ', 'has_sendhalos'  ]

    do iVar = 1, nMeshVars
      ! all meshInfo variables are scalar.
      ! If not change nComp for that variable
      nComponents = 1
      select case( trim(adjustl(meshInfoVar(iVar))) )
      case ('treeid')
        get_element => getTreeID
      case ('weight')
        get_element => getElemWeight
      case ('process')
        get_element => getMPIproc
      case ('elem_vol')
        get_element => deriveElemVol
      case ('vol_frac')
        get_element => deriveVolFrac
      case ('solidified', 'has_ibm', 'has_sendhalos')
        get_element => deriveProperty
      case ('level')
        get_element => deriveLevel
      case default
        write(logUnit(1),*) 'WARNING: Cannot append meshInfo variable: ' &
          &                 // trim(meshInfoVar(iVar))
        write(logUnit(1),*) 'without variable operation routine'
        cycle !go to next variable
      end select

      ! append variable to varSys
      call tem_varSys_append_derVar( me             = varSys,                  &
        &                            varName        = trim(meshInfoVar(iVar)), &
        &                            operType       = 'st_fun',                &
        &                            nComponents    = nComponents,             &
        &                            method_data    = c_null_ptr,              &
        &                            get_point      = get_point,               &
        &                            get_element    = get_element,             &
        &                            set_params     = set_params,              &
        &                            get_params     = get_params,              &
        &                            setup_indices  = setup_indices,           &
        &                            get_valOfIndex = get_valOfIndex,          &
        &                            pos            = addedPos,                &
        &                            wasAdded       = wasAdded                 )

      if (wasAdded) then
        write(logUnit(10),*) ' Appended variable:'//trim(meshInfoVar(iVar))
      else if (addedpos < 1) then
        write(logUnit(1),*) 'Error: variable '//trim(meshInfoVar(iVar)) &
          &                 // ' is not added to variable system'
      end if

    end do ! iVar

  end subroutine tem_varSys_append_meshInfoVar
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Get the treeID of an element
  recursive subroutine getTreeID( fun, varsys, elempos, time, tree, nElems, &
    &                             nDofs, res                                )
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
    ! local variable
    integer :: iElem
    ! -------------------------------------------------------------------- !

    res = 0.0_rk
    do iElem = 1, nElems
      res((iElem-1)*nDofs+1) = real( tree%treeID( elemPos(iElem) ), kind=rk )
    end do  ! iElem

    if (time%sim < 0.0_rk) then
      write(logUnit(10),*) 'Negative time for variable ', &
        &                  trim(varSys%varName%val(fun%myPos))
    end if

  end subroutine getTreeID
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Get the computational weight of elements
  recursive subroutine getElemWeight( fun, varsys, elempos, time, tree, &
    &                                 nElems, nDofs, res                )
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
    ! local variable
    integer :: iElem
    real(kind=rk), allocatable :: weights(:)
    logical :: has_weights
    ! -------------------------------------------------------------------- !

    allocate(weights(tree%nElems))
    call tem_load_weights( me      = tree,       &
      &                    weights = weights,    &
      &                    success = has_weights )
    res = 0.0_rk
    if (has_weights) then
      do iElem = 1, nElems
        res((iElem-1)*nDofs+1) = weights( elemPos(iElem) )
      end do  ! iElem
    end if

    if (time%sim < 0.0_rk) then
      write(logUnit(10),*) 'Negative time for variable ', &
        &                  trim(varSys%varName%val(fun%myPos))
    end if

  end subroutine getElemWeight
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Get the MPI rank to which the corresponding treeID belongs.
  !! This routine is used to visualize the rank affinity of the mesh elements
  recursive subroutine getMPIproc( fun, varsys, elempos, time, tree, nElems, &
    &                              nDofs, res                                )
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
    ! local variable
    integer :: iElem
    ! -------------------------------------------------------------------- !

    res = 0.0_rk
    do iElem = 1, nElems
      res((iElem-1)*nDofs+1) = real( tree%global%myPart, kind=rk )
    end do  ! iElem

    if (time%sim < 0.0_rk) then
      write(logUnit(10),*) 'Negative time for variable ', &
        &                  trim(varSys%varName%val(fun%myPos))
      write(logUnit(10),*) 'elempos(1):', elempos(1)
    end if

  end subroutine getMPIproc
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Get the level of an element
  recursive subroutine deriveLevel( fun, varsys, elempos, time, tree, nElems, &
    &                               nDofs, res                                )
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
    ! local variable
    integer :: iElem
    ! -------------------------------------------------------------------- !

    res = 0.0_rk
    do iElem = 1, nElems
      res((iElem-1)*nDofs+1) = tem_levelOf( tree%treeID(elempos(iElem)) )
    end do  ! iElem

    if (time%sim < 0.0_rk) then
      write(logUnit(10),*) 'Negative time for variable ', &
        &                  trim(varSys%varName%val(fun%myPos))
    end if

  end subroutine deriveLevel
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Get the volume of an element
  recursive subroutine deriveElemVol( fun, varsys, elempos, time, tree, &
    &                                 nElems, nDofs, res                )
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
    integer :: iElem
    real(kind=rk) :: dx
    ! -------------------------------------------------------------------- !

    res = 0.0_rk
    do iElem = 1, nElems
      if ( .not. btest( tree%ElemPropertyBits( elemPos(iElem) ), &
        &               prp_solid )                              ) then
        ! calculate without qVals
        dx = tem_ElemSize( tree, tree%treeID( elemPos(iElem) ) )
        res( (iElem-1)*nDofs+1 ) = dx*dx*dx
      end if
    end do  ! iElem

    if (time%sim < 0.0_rk) then
      write(logUnit(10),*) 'Negative time for variable ', &
        &                  trim(varSys%varName%val(fun%myPos))
    end if

  end subroutine deriveElemVol
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> compute the volume of the current element as a fraction of a reference
  !! bounding hexahedral volume (stored in tree%global%effLength).
  !! The sum of all element fractions then gives the fractional fluid volume
  !! occupation inside the reference volume, i.e. the porosity
  recursive subroutine deriveVolFrac( fun, varsys, elempos, time, tree, &
    &                                 nElems, nDofs, res                )
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
    integer :: iElem
    real(kind=rk) :: dx
    real(kind=rk) :: boundVol
    ! -------------------------------------------------------------------- !

    ! calculate the volume of the shape defined by the tree (only rectangulars)
    boundVol = tree%global%effLength(1) * tree%global%effLength(2) &
      &                                 * tree%global%effLength(3)

    res = 0.0_rk
    if (boundvol > 0.0_rk) then
      do iElem = 1, nElems
        if ( .not. btest( tree%ElemPropertyBits( elemPos(iElem) ), &
          &               prp_solid ) ) then
          dx = tem_ElemSize( tree, tree%treeID( elemPos(iElem) ))
          res( (iElem-1)*nDofs+1 ) = dx*dx*dx/boundVol
        end if
      end do  ! iElem
    end if

    ! Avoid warnings about unused arguments
    if (time%sim < 0.0_rk) then
      write(logUnit(10),*) 'Negative time for variable ', &
        &                  trim(varSys%varName%val(fun%myPos))
    end if

  end subroutine deriveVolFrac
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !>  Evaluate if the element is fluidifiable or not
  !!
  recursive subroutine deriveProperty( fun, varsys, elempos, time, tree, &
    &                                  nElems, nDofs, res                )
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
    integer :: iElem, check_prp
    ! -------------------------------------------------------------------- !

    check_prp = -999
    select case(trim(varSys%varname%val(fun%myPos)))
    case ('solidified')
      check_prp = prp_solid
    case ('has_ibm')
      check_prp = prp_hasIBM
    case ('has_sendhalos')
      check_prp = prp_sendHalo
    case default
      write(logUnit(1),*) 'Property unknown, check_prp != (prp_solid,' &
        &                 // ' prp_hasIBM, prp_sendHalo)'
      call tem_abort()
      write(logUnit(1),*) time%sim
    end select

    res = 0.0_rk
    do iElem = 1, nElems
      if( btest( tree%ElemPropertyBits( elemPos(iElem) ), check_prp) ) then
        res( (iElem-1)*nDofs+1 ) = 1.0_rk
      end if
    end do  ! iElem

  end subroutine deriveProperty
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


end module tem_meshInfo_module
! ---------------------------------------------------------------------------- !
! ---------------------------------------------------------------------------- !
