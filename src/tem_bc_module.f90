! Copyright (c) 2011-2013, 2016, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2014, 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2015-2017 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> Definition of the boundary property for elements
!!
!! The has-boundary [[tem_property_module]] highlights an element
!! that has at least one side with a boundary condition.
!! Therefore the neighbor on this side is not part of the computational domain,
!! and needs to be treated accordingly.
!! For all elements with the has-boundary property, the boundary conditions are
!! stored for all sides.
!! Yet this is done in a elementwise manner to allow arbitrary distribution of
!! the elments to parallel processing units. The appropiate handling of the
!! boundaries is solver specific and left to them.
!!
!! Each boundary is identified by a unique number (ID), which starts from 1 and
!! ranges up to the number of existing boundaries (nBCtypes).
!! In the solvers, the labels of the boundaries attached to the boundary IDs in
!! the array of strings <em>BC_label</em> are analyzed in order to treat the
!! boundaries adequately.
!!
!! For each element with the has-boundary property, a boundary ID is assigned
!! for each of the 26 direct neighbor directions.
!! This is done in the two-dimensional array <em>boundary_ID</em>,
!! where the first dimension has nSides entries and the second dimension has as
!! many entries as there elements with the has-boundary property on the local
!! partition.
!! For regular fluid neighbors, the boundary ID is set to zero in the respective
!! direction.
!!
!! ![boundary data structure](|media|/boundary_data_structure.png)
!!
!! To attach boundary conditions a bit is set for an element to indicate that
!! any of its sides has a boundary condition attached to it.
!! For all elements with this property and only those elements, boundary
!! conditions on every side of the element are stored.
!! This data is put in a separate file, containing the boundary conditions per
!! side for all elements with the boundary property.
!! The binary file has the same number of integer entries for all elements and
!! can therefore be easily accessed on an elemental basis.
!! In this list of elements the ordering from the complete mesh is preserved.
!! To describe the meaning of each boundary condition, indicated by an integer
!! number, there is an additional header file with labels for each type of
!! boundary condition in the mesh.
!! Only positive integers are considered as boundary entries, zeros indicate no
!! boundary in the corresponding direction, and negative numbers are reserved
!! for references on <em>treeIDs</em> as neighboring element.
!! Directly referencing a certain <em>treeID</em> instead of a boundary condition
!! on one of the sides of an element is used to implement periodic interfaces
!! within the universe cube.
!! Due to this layout of the boundary description, the boundaries can be read in
!! parallel with the following procedure:
!!
!! - Count the local number of elements with boundary conditions (given by the
!!   property bits defined for all elements).
!! - Build a prefix summation across all partitions to find the offset of the
!!   local set of boundary conditions in the sparse set of global boundary
!!   conditions.
!! - Read this set of elemental boundary conditions from the file accordingly
!!   from the offset up to the number of local elements with boundary
!!   conditions.
!!
!! With the strict ordering and uniform elemental treatment of data, the sparse
!! property information can be accessed efficiently.
!! Just a single collective prefix operation is required in the parallel system,
!! while the consumed disk storage is kept at a minimum.
!!
module tem_bc_module

  ! include treelm modules
  use env_module,               only: labelLen, rk
  use tem_spacetime_fun_module, only: tem_spacetime_fun_type
  use tem_logging_module,       only: logUnit
  use tem_aux_module,           only: tem_abort
  use tem_varSys_module,        only: tem_varSys_type, &
    &                                 tem_varSys_solverData_evalElem_type
  use tem_stringKeyValuePair_module, only: tem_stringKeyValuePair_type,      &
    &                                      grw_stringKeyValuePairArray_type, &
    &                                      append
  use tem_varMap_module,         only: tem_variable_loadMapping
  use tem_operation_module,      only: tem_indexLvl_type

  ! include aotus modules
  use aot_table_module, only: aot_table_open, aot_table_close, aot_get_val, &
    &                         aot_exists
  use aotus_module,     only: flu_state, aoterr_Fatal, aoterr_WrongType

  implicit none

  private

  public :: tem_bc_state_type
  public :: tem_load_bc_state

  !> boundary state type definition for boundary state variable
  !!
  !! State variables like  density, velocityX can be
  !! constant, temporal, spatial, dirichlet/neumann.
  !! res = f(x,y,z)*g(t) or res = f(x,y,z,t)
  type tem_bc_state_type
    !> Name of the state, this boundary condition applies to
    character(len=LabelLen) :: state_name
    !> Style of this boundary condition
    !! dirichlet = set value itself
    !! neumann = set derivative of value
    character(len=LabelLen) :: style
    !> Number of Components in this boundary variable.
    integer :: nComponents
    !> A flag to indicate that the state is properly defined
    logical :: isDefined
    !> Position of variable defined for the state_name in the varSys
    integer :: varPos
    !> Indices for points on the boundary, required for setup_index,
    !! getvalof_Index
    type(tem_indexLvl_type) :: pntIndex
  end type tem_bc_state_type

contains

! ****************************************************************************** !
  !> Function to load spatial, constant or temporal boundary conditions to
  !! the boundary state type (or combinations of them).
  !! \(res = f(x,y,z)*g(t)\) or \(res = f(x,y,z,t)\)
  !! Valid definitions:
  !! Example given for state variable velocityX in variable table
  !!```lua
  !! variable = {
  !!   { name = 'bc_pressure',
  !!     ncomponents = 1,
  !!     var_type = 'st_fun',
  !!     st_fun = 0.01
  !!   },
  !!   { name = 'bc_velocity',
  !!     ncomponents = 3,
  !!     var_type = 'st_fun',
  !!     st_fun = {0.01,0.0,0.0}
  !!   },
  !!   { name = 'bc_velocity_2',
  !!     ncomponents = 3,
  !!     var_type = 'st_fun',
  !!     st_fun = {
  !!                predefined = 'combined',
  !!                temporal = {predefined="linear", min_factor = 0.0,
  !!                             max_factor = 1.0, from_time = 0, to_time = 1000},
  !!                spatial = {predefined='parabol', shape = {kind = 'canoND',
  !!                           object = { origin={-2.0,0.0,0.0},vec={0.0,1.0,0.0}}},
  !!                           amplitude = {0.0,1.0,0.0}
  !!                          }
  !!                 }
  !!       }
  !! }
  !! boundary_condition = {
  !!  {
  !!    label = 'from_seeder',
  !!    kind = 'bc_kind',
  !!    style = 'dirichlet',
  !!    pressure = 'bc_pressure',
  !!    velocity = 'bc_velocity' --'bc_velocity_2'
  !!  }
  !! }
  !!```
  !!
  subroutine tem_load_bc_state( bc, state_name, nComp, style, conf, bc_handle, &
    &                           varDict, varSys, solverData_evalElem, ErrCode        )
    !---------------------------------------------------------------------------
    !> The boundary to fill
    type(tem_bc_state_type), intent(inout) :: bc
    !> The state variable to set with this boundary condition
    character(len=*), intent(in) :: state_name
    !> Number of Components in this boundary variable.
    integer, intent(in), optional :: nComp
    !> Style of this boundary condition
    !! dirichlet = set value itself
    !! neumann = set derivative of value
    character(len=*), optional, intent(in) :: style
    type(flu_State),intent(in) :: conf !< Lua state
    !> Handle to the table describing the boundary
    integer, intent(in) :: bc_handle
    !> The dictionary that contains the mapping between expected variables
    !! and the actual variables defined by the user.
    type(grw_stringKeyValuePairArray_type), intent(inout) :: varDict
    type(tem_varSys_type), intent(inout) :: varSys
    !> A routine that allows the construction of an element representation
    !! from a point values.
    type(tem_varSys_solverData_evalElem_type), &
      &  optional, intent(in) :: solverData_evalElem
    !> Error code
    integer, optional, intent(out) :: ErrCode
    ! ---------------------------------------------------------------------------
    integer :: iError
    integer :: state_handle
    character(len=labelLen) :: def_style, varname
    logical :: varExist
    type(tem_stringKeyValuePair_type) :: kvp
    real(kind=rk) :: numtest
    ! ---------------------------------------------------------------------------
    ! Assume undefined boundary condition for this state.
    bc%isDefined = .false.
    bc%state_name = trim(state_name)

    if (present(style)) then
      def_style = trim(style)
    else
      def_style = 'dirichlet'
    end if

    if (present(nComp)) then
      bc%nComponents = nComp
    else
      bc%nComponents = 1
    end if

    write(logUnit(1),*)' Loading bc state for '// trim(state_name)

    ! If table is defined load spacetime function variable name from table
    ! else load from bc_handle
    call aot_table_open( L       = conf,         &
      &                  parent  = bc_handle,    &
      &                  thandle = state_handle, &
      &                  key     = state_name    )

    varExist = .false.
    ! If boundary variable is defined as a table, then load boundary
    ! style and check if varname exist.
    if (state_handle /= 0) then
      ! Found a table, now check for the style of this boundary condition
      ! Default to dirichlet.
      call aot_get_val( L       = conf,         &
        &               thandle = state_handle, &
        &               val     = bc%style,     &
        &               ErrCode = iError,       &
        &               key     = 'style',      &
        &               default = def_style     )

      ! try to load variable as string
      varExist = aot_exists( L       = conf,         &
        &                    thandle = state_handle, &
        &                    key     = 'varname'     )
    end if ! variable defined as table

    !If varname exist load refered variable name from 'varname'
    ! and append the value to dictionary
    if (varExist) then
      bc%isDefined = .true.

      ! First check, whether this variable definition is a number.
      ! (They also satisfy strings).
      ! We do not accept numbers as variable names, instead this
      ! will be read as constant stfun.
      call aot_get_val( L       = conf,         &
        &               thandle = state_handle, &
        &               key     = 'varname',    &
        &               val     = numtest,      &
        &               ErrCode = iError        )
      if (btest(iError, aoterr_WrongType)) then
        ! Not a number, try to interpret it as a string.
        call aot_get_val( L       = conf,         &
           &              thandle = state_handle, &
           &              key     = 'varname',    &
           &              val     = varname,      &
           &              ErrCode = iError        )
        if( iError == 0 ) then
          ! Found a string, use it to refer to a variable.
          write(logUnit(3),*) 'Corresponding variable for ' &
            & // trim(state_name) // ' found: ' // trim(varname)
          kvp%key = state_name
          kvp%value = varname
          call append( me = varDict, val = kvp )
        else
          write(logUnit(1),*) 'Error: Unable to load state name "'// &
            & trim(state_name)//'" for boundary "'
          write(logUnit(1),*) '"'//trim(state_name)//'" is defined as table with'
          write(logUnit(1),*) 'varname, but failed to load variable name '// &
            &                 'as string'
          call tem_abort()
        end if
      end if
    else
      bc%isDefined = .true.
      bc%style = def_style

      ! Load spacetime function variable name
      call tem_variable_loadMapping( expectedName        = trim(state_name),   &
        &                            conf                = conf,               &
        &                            thandle             = bc_handle,          &
        &                            varDict             = varDict,            &
        &                            varSys              = varSys,             &
        &                            nComp               = bc%nComponents,     &
        &                            solverData_evalElem = solverData_evalElem,&
        &                            ErrCode             = iError              )

      if (iError /= 0) then
        write(logUnit(1),*) 'Error: Unable to load state name "'// &
          & trim(state_name)//'" for boundary '
        call tem_abort()
      end if

    end if

    call aot_table_close(L=conf, thandle=state_handle)

    if (present(ErrCode)) ErrCode = iError

  end subroutine tem_load_bc_state
! ****************************************************************************** !


end module tem_bc_module
! ****************************************************************************** !
