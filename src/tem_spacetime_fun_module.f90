! Copyright (c) 2012-2016,2019-2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2016, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012 Daniel Harlacher <d.harlacher@grs-sim.de>
! Copyright (c) 2012 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2017, 2019 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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
!> A module to deal with generic space-time functions.
!! ---
!!
!! It allows the definition of arbitrary space-time
!! dependent functions, and the multiplication of
!! a spatial function with a temporal function as
!! a often required special case.
!!
!! Space-Time functions might be predefined as Fortran functions, simple
!! constants or Lua functions.
!! They might return scalars or one-dimensional arrays.
!! If only a single value is expected from a Lua function, the Lua function
!! is supposed to return a scalar value. Irregardless if the function is
!! invoked by an array valued routine or a plain scalar one.
!! Otherwise, the Lua function has to return a table with the correct number
!! of entries identified by position.
!!
!! Note, that this makes the interface in the Lua script the same:
!!
!! - whenever a single return value is expected, a plain scalar should be
!!   returned.
!! - otherwise a table should be returned
!!
module tem_spacetime_fun_module
  use iso_c_binding,              only: c_ptr

  ! include treelm modules
  use env_module,                 only: LabelLen, rk, long_k, globalMaxLevels
  use treelmesh_module,           only: treelmesh_type
  use tem_bc_prop_module,         only: tem_bc_prop_type
  use tem_logging_module,         only: tem_toStr, logUnit
  use tem_subTree_type_module,    only: tem_subTree_type, tem_destroy_subTree
  use tem_subTree_module,         only: tem_create_subTree_of
  use tem_aux_module,             only: tem_abort
  use tem_geometry_module,        only: tem_baryOfId
  use tem_spatial_module,         only: tem_spatial_type, tem_load_spatial, &
    &                                   tem_spatial_for
  use tem_tools_module,           only: upper_to_lower, tem_horizontalSpacer
  use tem_temporal_module,        only: tem_temporal_type, tem_load_temporal, &
    &                                   tem_temporal_for
  use tem_shape_module,           only: tem_shape_type, tem_load_shape, &
    &                                   tem_global_shape
  use tem_time_module,            only: tem_time_type
  use tem_grow_array_module,      only: grw_intArray_type
  use tem_pointData_module,       only: tem_grwPoints_type,   &
    &                                   tem_pointData_list_type
  use tem_coupling_module,        only: tem_aps_coupling_type, &
    &                                   tem_aps_load_coupling
  use tem_polygon_material_module,only: tem_polygon_material_movement_single, &
   &                                    tem_polygon_material_movement_multi,     &
   &                                    tem_polygon_material_single_load,     &
   &                                    tem_polygon_material_type, &
   &                                    tem_polygon_material_load, &
   &                                    tem_polygon_material_multi_load


 use tem_miescatter_module,                         &
    & only: tem_miescatter_field_type,              &
    &       tem_load_miescatter_displacementfieldz, &
    &       tem_load_miescatter_magneticfieldx,     &
    &       tem_load_miescatter_magneticfieldy,     &
    &       tem_eval_miescatter_magnx,              &
    &       tem_eval_miescatter_magny,              &
    &       tem_eval_miescatter_displz
  use tem_cylindricalWave_module, only: tem_cylindricalWave_type, &
    &                                   tem_load_cylindricalWave, &
    &                                   tem_eval_cylindricalWave
  use tem_acoustic_pulse_module,  only: tem_acoustic_pulse_type, &
    &                                   tem_load_acoustic_pulse, &
    &                                   tem_eval_acoustic_pulse
  use tem_stencil_module,         only: tem_stencilHeader_type
  use tem_precice_module,         only: tem_precice_coupling_type, &
    &                                   tem_precice_load_coupling, &
    &                                   tem_precice_read,          &
    &                                   precice_available

  ! include aotus modules
  use aotus_module,               only: aot_get_val, aot_top_get_val,   &
    &                                   aoterr_Fatal, aoterr_WrongType, &
    &                                   aoterr_NonExistent, flu_State
  use aot_fun_module,             only: aot_fun_type, aot_fun_open, &
    &                                   aot_fun_close, aot_fun_put, &
    &                                   aot_fun_do, aot_fun_id
  use aot_table_module,           only: aot_table_open, aot_table_close, &
    &                                   aot_table_length, aot_exists,    &
    &                                   aot_type_of
  use aot_references_module,      only: aot_reference_for, aot_reference_to_top
  use flu_binding,                only: FLU_TNUMBER, FLU_TFUNCTION, &
    &                                   FLU_TTABLE, FLU_TSTRING

  implicit none

  private

  public :: tem_spacetime_fun_type
  public :: tem_load_spacetime
  public :: tem_spacetime_for
  public :: tem_st_fun_linkedList_type
  public :: tem_st_fun_listElem_type
  public :: append
  public :: tem_create_subTree_of_st_funList
  public :: tem_destroy_subTree_of_st_funList
  public :: tem_spacetime_hash_id

  !> Contains space time function definition
  type tem_spacetime_fun_type
    !> The function kind
    !!
    !! Should be either:
    !!
    !! - 'none': Not defined at all
    !! - 'const': Constant for all (x,y,z,t)
    !! - 'combined': This returns spatial(x,y,z)*temporal(t)
    !! - 'lua_fun': Function defined in the Lua script
    !! - Add predefined functions here
    character(len=labelLen) :: fun_kind

    !> spatial restrictions
    type(tem_shape_type), allocatable :: geom(:)

    !> subTree build from the shapes
    type(tem_subTree_type) :: subTree

    !> Number of components
    integer :: nComps

    !> constant value for nComponents
    real(kind=rk), allocatable :: const(:)

    !> Description of composite spatial fun
    type(tem_spatial_type) :: spatial

    !> Composite temporal fun
    type(tem_temporal_type) :: temporal

    !> Lua state handle to evaluate space time Lua function
    type(flu_State) :: conf

    !> Reference to the Lua function, if the st_fun is a Lua function.
    integer :: lua_fun_ref = 0

    !> type for the movement of the polygon
    type(tem_polygon_material_type) :: polygon_material
    !> Space-time function for Mie-series solution
    !! of electrodynamic scattering at dielectric sphere.
    type(tem_miescatter_field_type) :: mie_fun

    !> type for a scalar cylindrical wave.
    type(tem_cylindricalWave_type) :: cylindricalWave

    !> Description of an acoustic pulse.
    type(tem_acoustic_pulse_type) :: acoustic_pulse

    !> Apesmate coupling description
    type(tem_aps_coupling_type) :: aps_coupling

    !> preCICE coupling description
    type(tem_precice_coupling_type) :: precice_coupling
  end type tem_spacetime_fun_type


  !> An element for a spacetime function within a linked list.
  !!
  !! Besides the actual list of spacetime definitions that are provided, there
  !! is a pointer to the next element in the list.
  type tem_st_fun_listElem_type

    !> Number of space time functions
    integer :: nVals

    !> Space time function target which C_ptr will point to
    !
    !! We maintain a list of spacetime functions here, as each one might be
    !! restricted to a subtree, and multiple of those locally different function
    !! definitions might be used to define a single variable like a source term.
    type(tem_spacetime_fun_type), dimension(:), pointer :: val => NULL()

    !> Points data containing space coordinates or evaluated
    !! values for time-indepentent functions
    type(tem_pointData_list_type) :: pntData

    !> A pointer to possibly additional solver data.
    !!
    !! This is for example used to keep a link to the projection data
    !! in Ateles to enable the construction of element data from the
    !! point data provided by the space-time function.
    type(c_ptr) :: solver_bundle

    !> Used to decided whether this spacetime functions are used
    !! for surface or volume i.e boundary or source.
    !! Boundary is treated as surface and source as volume
    !! coupling type can be rather surface or volume.
    !! For boundary. isSurface = 0
    !! For volume, isSurface = 1
    integer :: isSurface = -1

    !> Pointer to next space time function
    type(tem_st_fun_listElem_type), pointer :: next => NULL()

  end type tem_st_fun_listElem_type


  !> Type used to create linked list of space time function
  !! (tem_st_fun_listElem_type)
  !!
  !! We need to point to the spacetime functions in variable declarations and
  !! therefore need to have dynamic data structure to keep all space time
  !! functions. Most likely this data structure is only filled once and never
  !! iterated through afterwards. Direct pointers to the entries of the list
  !! are maintained wherever necessary. Therefore, it should be fine to use
  !! a linked list here.
  type tem_st_fun_linkedList_type
    !> Pointer to the first entry in the linked list
    type(tem_st_fun_listElem_type), pointer :: head => NULL()
  end type tem_st_fun_linkedList_type

  !> Maximum length for constant vectors to read from the configuration
  integer, parameter :: maxveclen = 10

  interface append
    module procedure append_stFunSingle_ToLinkList
    module procedure append_stFunArray_ToLinkList
  end interface append

  interface tem_spacetime_lua_for
    module procedure tem_spacetime_lua_for_treeIds
    module procedure tem_spacetime_lua_vector_for_treeIds
    module procedure tem_spacetime_lua_for_coord
    module procedure tem_spacetime_lua_vector_for_coord
  end interface tem_spacetime_lua_for

  interface tem_spacetime_for
    module procedure tem_spacetime_for_treeIds
    module procedure tem_spacetime_vector_for_treeIds
    module procedure tem_spacetime_for_coord
    module procedure tem_spacetime_vector_for_coord
    module procedure tem_spacetime_for_stcoord
    module procedure tem_spacetime_scalar_for_index
    module procedure tem_spacetime_vector_for_index
  end interface tem_spacetime_for

  interface tem_load_spacetime
    module procedure tem_load_spacetime_single
    module procedure tem_load_spacetime_table
  end interface tem_load_spacetime


contains


  ! ************************************************************************ !
  !> This routine appends a new array of space time functions st_fun to the
  !! linked list me.
  !!
  !! HK: It might be useful to return a pointer to the appended new stfun entry
  !!     from this routine.
  subroutine append_stFunArray_ToLinkList(me, st_fun, new)
    ! -------------------------------------------------------------------- !
    !> Linked list to append the array of spacetime functions to.
    type(tem_st_fun_linkedList_type), intent(inout) :: me

    !> Spacetime fun information to add to the list.
    type(tem_spacetime_fun_type), intent(in) :: st_fun(:)

    type(tem_st_fun_listElem_type), optional, pointer, intent(out) :: new
    ! -------------------------------------------------------------------- !
    type(tem_st_fun_listElem_type), pointer :: current
    type(tem_st_fun_listElem_type), pointer :: lnew
    ! -------------------------------------------------------------------- !

    current => NULL()
    allocate(lnew)
    lnew%nVals = size(st_fun)
    allocate(lnew%val(lnew%nVals))
    lnew%val = st_fun

    !look for the last element in the linked list
    if (associated(me%head)) then
      current => me%head
      do while (associated(current%next))
        current => current%next
      enddo

      allocate(current%next)
      current%next => lnew
    else
      ! If current does not point anywhere yet, this has to be the first entry,
      ! allocate that and store it as the head element.
      allocate(me%head)
      me%head => lnew
    endif

    if (present(new)) new => lnew

  end subroutine append_stFunArray_ToLinkList
  ! ************************************************************************ !



  ! ************************************************************************ !
  !> This routine appends new space time function to linked list of
  !! tem_st_fun_linkedList
  subroutine append_stFunSingle_ToLinkList(me, st_fun, new)
    ! -------------------------------------------------------------------- !
    !> Linked list to append the spacetime function to.
    type(tem_st_fun_linkedList_type), intent(inout) :: me

    !> Spacetime fun information to add to the list.
    type(tem_spacetime_fun_type), intent(in) :: st_fun
    type(tem_st_fun_listElem_type), pointer, optional, intent(out) :: new
    ! -------------------------------------------------------------------- !
    type(tem_st_fun_listElem_type), pointer :: current
    type(tem_st_fun_listElem_type), pointer :: lnew
    ! -------------------------------------------------------------------- !

    allocate(lnew)
    lnew%nVals = 1
    allocate(lnew%val(lnew%nVals))
    lnew%val = st_fun
    lnew%next => null()

    ! Look for the last element in linked list
    if (associated(me%head)) then
      current => me%head
      do
       if (.not. associated(current%next)) EXIT
         current => current%next
      end do
      current%next => lnew
    else
      allocate(me%head)
      me%head => lnew
    end if

    if (present(new)) new => lnew

  end subroutine append_stFunSingle_ToLinkList
  ! ************************************************************************ !
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> create subtree for shapes defined in each spacetime functions
  subroutine tem_create_subTree_of_st_funList( me, tree, bc_prop, stencil )
    ! -------------------------------------------------------------------- !
    !> Linked list to append the spacetime function to.
    type( tem_st_fun_linkedList_type ), intent(inout) :: me
    !> Global treelmesh
    type( treelmesh_type ),     intent(in) :: tree
    !> bc property
    type( tem_bc_prop_type ), intent(in) :: bc_prop
    !> stencil
    type( tem_stencilHeader_type ), optional, intent(in) :: stencil

    ! -------------------------------------------------------------------- !
    type(tem_st_fun_listElem_type), pointer :: st_fun
    integer :: iSt, iList
    ! -------------------------------------------------------------------- !
    call tem_horizontalSpacer( fUnit = logUnit(3))
    write(logUnit(3),*) 'Create subtree for all space time functions stored '
    write(logUnit(3),*) 'in linked list of spacetime function'

    st_fun => me%head
    iList = 0
    do
      if (.not. associated(st_fun)) EXIT
      iList = iList + 1
      do iSt = 1, st_fun%nVals
        call tem_create_subTree_of( inTree  = tree,                    &
          &                         subTree = st_fun%val(iSt)%subTree, &
          &                         bc_prop = bc_prop,                 &
          &                         stencil = stencil,                 &
          &                         inShape = st_fun%val(iSt)%geom     )
      end do
      st_fun => st_fun%next
    end do
    write(logUnit(3),'(a,i3,a)') ' Done creating subtree for ', iList, &
      &                          ' spacetime functions'
    call tem_horizontalSpacer( fUnit = logUnit(3))

  end subroutine tem_create_subTree_of_st_funList
  ! ************************************************************************ !
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> destroy subtree for shapes defined in each spacetime functions for
  !! dynamic load balancing
  subroutine tem_destroy_subTree_of_st_funList( me )
    ! -------------------------------------------------------------------- !
    !> Linked list to append the spacetime function to.
    type( tem_st_fun_linkedList_type ), intent(inout) :: me
    ! -------------------------------------------------------------------- !
    type(tem_st_fun_listElem_type), pointer :: st_fun
    integer :: iSt, iList
    ! -------------------------------------------------------------------- !
    call tem_horizontalSpacer( fUnit = logUnit(3))
    write(logUnit(3),*) 'Create subtree for all space time functions stored '
    write(logUnit(3),*) 'in linked list of spacetime function'

    st_fun => me%head
    iList = 0
    do
      if (.not. associated(st_fun)) EXIT
      iList = iList + 1
      do iSt = 1, st_fun%nVals
        call tem_destroy_subTree(me = st_fun%val(iSt)%subTree)
      end do
      st_fun => st_fun%next
    end do
    write(logUnit(3),'(a,i3,a)') ' Done creating subtree for ', iList, &
      &                          ' spacetime functions'
    call tem_horizontalSpacer( fUnit = logUnit(3))

  end subroutine tem_destroy_subTree_of_st_funList
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine loads table of spacetime functions from the given key or pos
  !!
  !! NOTE: If any of the entries in the table can not be interpreted as a
  !!       space-time function, none will be returne at all, and an error code
  !!       of -1 will be set. "Me" will be deallocated in this case.
  !!       The routine first attempts to read the given key as a single
  !!       space-time function definition, only if that fails, it tries to read
  !!       it as a table of functions.
  subroutine tem_load_spacetime_table( me, conf, parent, key, nComp, &
    &                                  errCode )
    ! -------------------------------------------------------------------- !
    !> spacetime fun information
    type(tem_spacetime_fun_type), allocatable, intent(out) :: me(:)
    !> lua state handle
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(inout), optional :: parent
    !> name of the variable which is defined as spacetime function
    character(len=*), intent(in) :: key
    !> number of components of the variable
    integer, intent(in), optional :: nComp
    !> errCode /=0, space time function fails
    !! use errCode to abort code outside this routine call
    integer, optional, intent(out) :: errCode
    ! -------------------------------------------------------------------- !
    !> aotus table handle
    integer :: thandle
    ! counter variables
    integer :: nSt, iSt
    integer :: errCode_loc
    character(len=labelLen) :: buffer
    type(tem_st_fun_listElem_type), pointer :: current
    ! -------------------------------------------------------------------- !
    current => NULL()
    nSt = 1
    allocate(me(nSt))
    ! read in a single spacetime function with given key
    call tem_load_spacetime_single( me      = me(1),      &
      &                             conf    = conf,       &
      &                             parent  = parent,     &
      &                             key     = key,        &
      &                             nComp   = nComp,      &
      &                             errCode = errCode_loc )

    if (errCode_loc /= 0) then
      write(logUnit(3),*) 'Error loading spacetime function from key ' &
       &                  // trim(key)
      write(logUnit(3),*) 'Try to load it from table'

      ! Error loading spacetime function directly via key
      ! Try opening as table of space time functions
      call aot_table_open( L       = conf,     &
        &                  thandle = thandle,  &
        &                  parent  = parent,   &
        &                  key     = trim(key) )

      nSt = aot_table_length( L=conf, thandle=thandle )
      write(logUnit(3),*) 'Multiple spacetime functions are defined'
      write(buffer,'(i3)') nSt
      write(logUnit(3),*) 'Number of spacetime fun tables '//trim(buffer)
      deallocate(me)
      allocate(me(nSt))
      do iSt = 1, nSt
        write(buffer,'(i3)') iSt
        write(logUnit(3),*)
        write(logUnit(3),*) 'loading space time function at pos: '//trim(buffer)
        call tem_load_spacetime_single( me      = me(iSt),    &
          &                             conf    = conf,       &
          &                             parent  = thandle,    &
          &                             pos     = iSt,        &
          &                             nComp   = nComp,      &
          &                             errCode = errCode_loc )

        if (errCode_loc /= 0) then
          write(logUnit(3),*) 'Error loading spacetime function at pos ' &
            &                 //trim(buffer)
          write(logUnit(3),*) 'Aborting the attempt to more functions '
          EXIT
        else
        end if
      end do
      call aot_table_close( L = conf, thandle = thandle )
    end if

    if (errCode_loc /= 0) then
      deallocate(me)
    else
      write(logUnit(1),*) 'Space-time function for key '//trim(key)//':'
      do iSt=1,nSt
          write(logUnit(1),*) '      pos ', iSt, &
            &                 ' is defined as ' // trim(me(iSt)%fun_kind)
      end do
    end if

    if (present(errCode)) errCode = errCode_loc

  end subroutine tem_load_spacetime_table
  ! ************************************************************************ !
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine loads the single spacetime function from the
  !! given key or position
  !!
  !! If spacetime is defined as block than read block for key word
  !! predefined/fun/const and load shape inside a block else
  !! define directly as lua function or constant.
  !! If predefined is defined inside a block, define other neccessary
  !! parameters for predefined.
  !! If shape table is not defined, shape is set to "all"
  !!
  !! Valid definitions:
  !! - Constant
  !!~~~~~~~~~~~~~~~~~~~~~{.lua}
  !!st_fun = 1.0
  !!or
  !!st_fun = {const = 1.0, shape = {..}}
  !!~~~~~~~~~~~~~~~~~~~~~
  !! - lua_function
  !!~~~~~~~~~~~~~~~~~~~~~{.lua}
  !!st_fun = lua_fun_name
  !! --or
  !!st_fun = {fun=lua_fun_name, shape={..}}
  !!~~~~~~~~~~~~~~~~~~~~~
  !! Note. Lua function take 4 input arguments (x,y,z,t) i.e barycentric
  !! coordinates of an element and time
  !! - Predefined Fortran function
  !!~~~~~~~~~~~~~~~~~~~~~{.lua}
  !! st_fun  = {predefined = "fun_name", fun_parameters}
  !!~~~~~~~~~~~~~~~~~~~~~
  ! For more examples, look at [[tem_spacetime_fun_test]]
  !! This definition can itself to be part of tables to define multiple
  !! space time functions.
  recursive subroutine tem_load_spacetime_single( me, conf, parent, key, pos, &
    &                                             nComp, errCode, recurred )
    ! -------------------------------------------------------------------- !
    !> spacetime fun information
    type(tem_spacetime_fun_type), intent(out) :: me
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in), optional :: parent
    !> name of the variable which is defined as spacetime function
    character(len=*), intent(in), optional :: key
    !> position of spacetime fun in a table
    integer, intent(in), optional :: pos
    !> number of components of the variable
    integer, intent(in), optional :: nComp
    !> errCode /=0, space time function fails
    !! use errCode to abort code outside this routine call
    integer, optional, intent(out) :: errCode
    !> Number of recursion steps done so far (defaults to 0)
    integer, optional, intent(in) :: recurred
    ! -------------------------------------------------------------------- !
    type(aot_fun_type) :: fun
    ! aotus handle
    integer :: thandle
    ! error variables
    integer :: iError, iError_shape
    ! local ncomp
    logical :: stFunNotATable
    integer :: ltype
    logical :: has_key(3) ! There are three different possible keys we need
                          ! to check for.
    character(len=labelLen) :: fun_key
    integer :: loc_recurred
    ! -------------------------------------------------------------------- !
    loc_recurred = 0
    if (present(recurred)) loc_recurred = recurred

    iError = huge(iError)
    iError_shape = huge(iError_shape)

    if (present(ErrCode)) ErrCode = iError

    ! Do not allow more than 1 recursion step
    if (loc_recurred > 1) RETURN

    if (present(key)) then
      write(logUnit(3),*) 'loading space time function from key: ', trim(key)
    end if

    ! store conf to load lua space time function
    me%conf = conf

    ! default values
    stFunNotATable = .true.
    me%fun_kind = 'none'

    if (present(nComp)) then
      me%nComps = nComp
    else
      me%nComps = 1
    end if

    ltype = aot_type_of( L       = conf,   &
      &                  thandle = parent, &
      &                  key     = key,    &
      &                  pos     = pos     )

    select case(ltype)
    case(FLU_TNUMBER)
      write(logunit(9),*) 'Trying to load ST-Fun as a scalar constant...'
      ! Try to load the top of the stack as a constant value.
      call load_spacetime_asConst( me      = me,      &
        &                          conf    = conf,    &
        &                          errCode = iError,  &
        &                          nComp   = nComp    )

    case(FLU_TFUNCTION)
      ! Try to interpret the top of the stack as a Lua function.
      write(logunit(9),*) 'Trying to load ST-Fun as Lua function...'
      call aot_fun_open( L      = conf, &
        &                fun    = fun   )
      if (fun%handle /= 0) then
        write(logunit(9),*) '... ST-Fun is a Lua function!'
        ! There is a function defined in Lua.
        me%fun_kind = 'lua_fun'
        ! Store a reference to this function.
        me%lua_fun_ref = aot_reference_for(conf)
        call aot_fun_close( L=conf, fun=fun )
        iError = 0
      else
        iError = -1
      end if

    case(FLU_TSTRING)
      if (loc_recurred == 1) then
        write(logunit(9),*) 'Trying to load ST-Fun as predefined function...'
        call aot_get_val( L       = conf,        &
          &               val     = me%fun_kind, &
          &               default = 'none',      &
          &               ErrCode = iError       )
        if (iError == 0) then
          call load_spacetime_predefined( me      = me,     &
            &                             conf    = conf,   &
            &                             thandle = parent, &
            &                             nComp   = nComp   )
        end if
      else
        ! A predefined spacetime function is not possible without an embedding
        ! table, return an error if we are not inside a table!
        iError = -1
      end if

    case(FLU_TTABLE)
      ! First, try to interpret the table as a vectorial constant.
      write(logunit(9),*) 'Trying to load ST-Fun as a vectorial constant...'
      ! Try to load the top of the stack as a constant value.
      call load_spacetime_asConst( me      = me,      &
        &                          conf    = conf,    &
        &                          errCode = iError,  &
        &                          nComp   = nComp    )

      if (iError < 0) then
        write(logunit(9),*) '... not a vectorial constant.'
        call aot_table_open( L       = conf,    &
          &                  thandle = thandle, &
          &                  parent  = parent,  &
          &                  key     = key,     &
          &                  pos     = pos      )

        recursion: if (loc_recurred == 0) then

          write(logunit(9),*) 'Trying to obtain spacetime function definition' &
            & // ' within the provided table.'
          stFunNotATable = .false.

          ! For backwards compatibility we have several options to use as
          ! keywords for the function definition.
          ! Exactly one of them has to be defined.
          has_key(1) = aot_exists( L       = conf,    &
            &                      thandle = thandle, &
            &                      key     = 'const'  )
          if (has_key(1)) fun_key = 'const'

          has_key(2) = aot_exists( L       = conf,    &
            &                      thandle = thandle, &
            &                      key     = 'fun'    )
          if (has_key(2)) fun_key = 'fun'

          has_key(3) = aot_exists( L       = conf,        &
            &                      thandle = thandle,     &
            &                      key     = 'predefined' )
          if (has_key(3)) fun_key = 'predefined'

          ! Only if exactly one key is defined, we proceed and try to load
          ! that as a space-time function itself.
          if ( count(has_key) == 1 ) then
            call tem_load_spacetime_single( me       = me,              &
              &                             conf     = conf,            &
              &                             parent   = thandle,         &
              &                             key      = trim(fun_key),   &
              &                             nComp    = nComp,           &
              &                             errCode  = iError,          &
              &                             recurred = loc_recurred + 1 )
          end if

          ! Only during first call try to load the shape for the function, and
          ! identify function itself by one of the keywords.
          ! As the definition is a table, there might be a
          ! shape defined to restrict the area of the function.
          ! Shape either has to be given via the keyword 'shape'.
          write(logunit(9),*) 'Trying to obtain the shape...'
          call tem_load_shape( me      = me%geom,     &
            &                  conf    = conf,        &
            &                  parent  = thandle,     &
            &                  key     = 'shape',     &
            &                  iError  = iError_shape )
        else recursion

          ! Loading predefined space time function from a subtable.
          write(logunit(9),*) '... failed loading vectorial constant'
          write(logunit(9),*) 'Attempting to load a predefined function in' &
            & // ' a subtable.'
          call aot_table_open( L       = conf,    &
            &                  thandle = thandle, &
            &                  parent  = parent,  &
            &                  key     = key,     &
            &                  pos     = pos      )
          call aot_get_val( L       = conf,        &
            &               val     = me%fun_kind, &
            &               thandle = thandle,     &
            &               pos     = 1,           &
            &               default = 'none',      &
            &               ErrCode = iError       )
          if (iError == 0) then
            call load_spacetime_predefined( me      = me,      &
              &                             conf    = conf,    &
              &                             thandle = thandle, &
              &                             nComp   = nComp    )
          end if

        end if recursion

      end if

    end select


    if ( trim(me%fun_kind) == 'const') then
      if ( me%nComps /= size(me%const) ) then
        write(logUnit(1),*) 'WARNING: In loading stFun, nComps of const ' &
          &                 //'loaded:', size(me%const)
        write(logUnit(1),*) '         does not match argumental nComps: ', &
          &                 me%nComps
        write(logUnit(1),*) '         Setting nComps to size(const).'
        me%nComps = size(me%const)
      end if
    end if

    if (loc_recurred == 0) then
      if (iError /= 0) then
        me%fun_kind = 'none'
        write(logunit(3), *) 'Could not load spacetime function!'
      end if

      if (trim(me%fun_kind) == 'const') &
        &  write(logUnit(3),*)          &
        &           'Spacetime function is a const value: ', me%const

      ! if shape is defined inside a table but function type is not
      ! defined with key word "fun", "predefined", "const" then
      ! terminate code with error message
      if (iError_shape == 0 .and. iError == -1) then
        write(logUnit(1),*) 'ERROR: Shape is defined inside a table but' &
          & // ' spacetime function is unidentified.'
        write(logUnit(1),*) 'Provide spacetime function kind via key word:' &
          & // ' "fun" / "predefined" / "const"'
        call tem_abort()
      end if
      ! if shape table is not defined
      if (stFunNotATable .and. iError /= -1) then
        write(logUnit(1),*) 'St-fun is not a table, thus setting global shape.'
        if (allocated(me%geom)) deallocate(me%geom)
        allocate(me%geom(1))
        me%geom(1)%kind = 'all'
        me%geom(1)%shapeID = tem_global_shape
      end if
    end if

    if (present(errCode)) errCode = iError

  end subroutine tem_load_spacetime_single
  ! ************************************************************************ !
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Load predefined space time function
  !!
  !! @todo Provide an error code to return instead of aborting.
  subroutine load_spacetime_predefined( me, conf, thandle, nComp )
    ! --------------------------------------------------------------------------
    !> spacetime fun information
    type(tem_spacetime_fun_type), intent(inout) :: me
    !> lua state type
    type(flu_State) :: conf
    !> spacetime function handle
    integer, intent(in) :: thandle
    !> number of components of the variable
    integer, optional, intent(in) :: nComp
    ! --------------------------------------------------------------------------
    character(len=labelLen) :: funkind
    integer :: iError
    ! --------------------------------------------------------------------------

    funkind = adjustl(me%fun_kind)
    funkind = upper_to_lower(funkind)
    select case(trim(funkind))
    case('combined')
      ! A typical case of a space-time function that can be represented
      ! a multiplicative combination of a temporal and a spatial part.
      call tem_load_spatial( me     = me%spatial, &
        &                    conf   = conf,       &
        &                    parent = thandle,    &
        &                    errCode= iError,     &
        &                    nComp  = nComp       )
      if (me%spatial%kind == 'none') then
        call tem_abort('Error in loading combined space-time function:' &
          &            //' no spatial function defined!'                )
      end if
      call tem_load_temporal( me     = me%temporal, &
        &                     conf   = conf,        &
        &                     parent = thandle      )
      if (me%temporal%kind == 'none') then
        call tem_abort('Error in loading combined space-time function:' &
          &            //' no temporal function defined!'               )
      end if

    case('miescatter_displacementfieldz')
      call tem_load_miescatter_displacementfieldz( conf    = conf,      &
        &                                          thandle = thandle,   &
        &                                          me      = me%mie_fun )

    case('miescatter_magneticfieldx')
      call tem_load_miescatter_magneticfieldx( conf    = conf,      &
        &                                      thandle = thandle,   &
        &                                      me      = me%mie_fun )

    case('miescatter_magneticfieldy')
      call tem_load_miescatter_magneticfieldy( conf    = conf,      &
        &                                      thandle = thandle,   &
        &                                      me      = me%mie_fun )

    case('cylindricalwave')
      call tem_load_cylindricalWave( conf    = conf,              &
        &                            thandle = thandle,           &
        &                            me      = me%cylindricalWave )

    case('acoustic_pulse')
      call tem_load_acoustic_pulse( conf    = conf,             &
        &                           thandle = thandle,          &
        &                           me      = me%acoustic_pulse )

    case('polygon_body_2d','polygon_body_3d')
      call tem_polygon_material_single_load( &
        & conf    = conf,                    &
        & thandle = thandle,                 &
        & me      = me%polygon_material      )

    case('polygon_multi_body_2d', 'polygon_multi_body_3d')
      call tem_polygon_material_multi_load( &
        & conf    = conf,                   &
        & thandle = thandle,                &
        & me      = me%polygon_material     )
    case('apesmate')
      call tem_aps_load_coupling( me      = me%aps_coupling, &
        &                         thandle = thandle,         &
        &                         conf    = conf             )

    case('precice')
      if (.not. precice_available) then
        call tem_abort(' Error: Spacetime function predefinded = precice not &
          &  available if not compiled with preCICE support, stopping... '  )
      end if
      call tem_precice_load_coupling( me      = me%precice_coupling, &
        &                             conf    = conf,                &
        &                             thandle = thandle              )

    case default
      ! If you introduce new predefined functions, add their loading
      ! routine in a seperate case branch here.
      write(logUnit(1),*) 'ERROR in definition of a space-time function:'
      write(logUnit(1),*) 'Unknown "predefined" space-time function: '// &
        &                 trim( me%fun_kind )
      call tem_abort()
    end select

  end subroutine load_spacetime_predefined
  ! ************************************************************************ !
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Load space time function as constant
  subroutine load_spacetime_asConst( me, conf, errCode, parent, key, pos, &
    &                                nComp )
    ! -------------------------------------------------------------------- !
    !> spacetime fun information
    type(tem_spacetime_fun_type), intent(inout) :: me
    !> lua state type
    type(flu_State) :: conf
    !> errCode = -1, if space time function is not defined as constant
    integer, intent(out) :: errCode
    !> aotus parent handle
    integer, intent(in), optional :: parent
    !> name of the variable which is defined as spacetime function
    character(len=*), intent(in), optional :: key
    !> position of spacetime fun in a table
    integer, intent(in), optional :: pos
    !> number of components of the variable
    integer, optional, intent(in) :: nComp
    ! -------------------------------------------------------------------- !
    integer, allocatable :: vError(:)
    logical :: check_varlen
    ! -------------------------------------------------------------------- !

    errCode = 0
    check_varlen = .true.

    if (allocated(me%const)) deallocate(me%const)

    if (present(nComp)) then
      ! There is a given number of components expected to be filled.

      if (nComp > 1) then
        ! Only load the fixed sized vector if more than a single component is
        ! expected. Otherwise, we use the loading for arrays of unknown length
        ! below, as that als takes care of scalar number definitions.

        check_varlen = .false. ! Check for fixed sized array.
        allocate(me%const(nComp), vError(nComp))
        write(logUnit(9),"(A,I0,A)") 'Trying to read constant as a vector with ', &
          &                          nComp, ' components.'
        if (present(key)) write(logUnit(9),*) ' key ', trim(key)
        call aot_get_val( L         = conf,     &
          &               thandle   = parent,   &
          &               val       = me%const, &
          &               key       = key,      &
          &               pos       = pos,      &
          &               ErrCode   = vError    )
        if ( any(btest(vError, aoterr_Fatal)) ) then
          write(logUnit(6),*) 'Attempted interpretation of spacetime function'
          write(logUnit(6),*) 'as a vectorial constant failed.'
          errCode = -1
        else
          me%fun_kind = 'const'
        end if

      end if
    end if


    if (check_varlen) then
      ! If ncomp is not provided, or ncomp == 1, we proceed and try to
      ! load the constant as a vector of unknown length.

      write(logUnit(9),*) 'Trying to read constant as a vector'
      write(logUnit(9),*) 'with unknown or at most 1 components.'
      call aot_get_val( L         = conf,     &
        &               thandle   = parent,   &
        &               maxLength = 10000,    &
        &               val       = me%const, &
        &               key       = key,      &
        &               pos       = pos,      &
        &               ErrCode   = vError    )
      if ( any(btest(vError, aoterr_Fatal)) ) then
        write(logUnit(6),*) 'Attempted interpretation of spacetime function'
        write(logUnit(6),*) 'as a vectorial constant failed.'
        errCode = -1
      else
        me%fun_kind = 'const'
      end if

    end if

    ! If a certain number of components has been requested, ensure that we
    ! got exactly that number of components.
    if (present(nComp)) then
      if (.not. nComp == size(me%const)) then
        me%fun_kind = 'none'
        errCode = -1
      end if
    end if

  end subroutine load_spacetime_asConst
  ! ************************************************************************ !
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function compute space time function for given list of treeIDs
  !! @todo pass subtree with treeIDs lists instead of treeIds and tree
  !!
  function tem_spacetime_for_treeIDs( me, treeIDs, time, tree, n ) result(res)
    ! -------------------------------------------------------------------- !
    !> Spacetime function to evaluate
    type(tem_spacetime_fun_type) :: me
    !> Global treelmesh to look for positions in
    type(treelmesh_type), intent(in) ::tree
    !> Number of values to return
    integer, intent(in) :: n
    !> TreeIDs where to evaluate the function
    integer(kind=long_k), intent( in ) :: treeIDs(n)
    !> timer object incl. the current time information
    type(tem_time_type), intent( in )  :: time
    !> return value of the function
    real(kind=rk) :: res(n)
    ! -------------------------------------------------------------------- !

    select case(trim(adjustl(me%fun_kind)))
    case('none')
      res = 0.0_rk
    case('const')
      res = me%const(1)
    case('lua_fun')
      res = tem_spacetime_lua_for( me%lua_fun_ref, treeIDs, time, tree, n, &
        &                          me%conf)
    case('combined')
      res = tem_spatial_for( me      = me%spatial,    &
        &                    treeIDs = treeIDs,       &
        &                    tree    = tree,          &
        &                    n       = n           )  &
        & * tem_temporal_for( temporal = me%temporal, &
        &                     time     = time         )
    case default
      write(logUnit(1),*)'ERROR: Unknown spatial function in '//               &
        &                'tem_spacetime_for_treeIDs.'
      write(logUnit(1),*) 'tem_spacetime_fun_type%fun_kind = ', &
        &                 trim(me%fun_kind)
      call tem_abort()
    end select
  end function tem_spacetime_for_treeIDs
  ! ************************************************************************ !
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function compute space time function that give bach a table of results
  !! for given list of treeIDs
  !! @todo pass subtree with treeIDs lists instead of treeIds and tree
  !!
  function tem_spacetime_vector_for_treeIDs( me, treeIDs, time, &
    &                                        tree, n, ncomp) result(res)
    ! -------------------------------------------------------------------- !
    !> Spacetime function to evaluate
    type(tem_spacetime_fun_type) :: me
    !> Global treelmesh to look for positions in
    type(treelmesh_type), intent(in) ::tree
    !> Number of tables to return
    integer, intent(in) :: n
    !> Number of values in a Table
    integer, intent(in) :: ncomp
    !> TreeIDs where to evaluate the function
    integer(kind=long_k), intent(in) :: treeIDs(n)
    !> timer object incl. the current time information
    type(tem_time_type), intent(in)  :: time
    !> return value of the function
    real(kind=rk) :: res(n,ncomp)
    ! -------------------------------------------------------------------- !
    ! counter
    integer :: i
    ! -------------------------------------------------------------------- !


    select case(trim(adjustl(me%fun_kind)))
    case('none')
      res = 0.0_rk
    case('const')
      do i = 1, nComp
        res(:,i) = me%const(i)
      end do
    case('lua_fun')
      res = tem_spacetime_lua_for( fun_ref = me%lua_fun_ref, &
        &                          treeIDs = treeIDs,        &
        &                          time    = time,           &
        &                          tree    = tree,           &
        &                          n       = n,              &
        &                          nComp   = nComp,          &
        &                          conf    = me%conf         )
    case('combined')
      res = tem_spatial_for( me      = me%spatial,    &
        &                    treeIDs = treeIDs,       &
        &                    tree    = tree,          &
        &                    n       = n,             &
        &                    ncomp   = nComp     )    &
        & * tem_temporal_for( temporal = me%temporal, &
        &                     time     = time         )
    case default
      call tem_abort('ERROR: Unknown spatial function in' &
        & // 'tem_spacetime_vector_for_treeIDs.'          )
    end select
  end function tem_spacetime_vector_for_treeIDs
  ! ************************************************************************ !
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function computes the space time function for given list of space-time
  !! coordinates.
  !!
  function tem_spacetime_for_stcoord( me, stCoord, n  ) result(res)
    ! -------------------------------------------------------------------- !
    !> Spacetime function to evaluate
    type(tem_spacetime_fun_type) :: me
    !> Number of values to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z,t coordinates
    real(kind=rk), intent(in) :: stCoord(n,4)
    !> return value of the function
    real(kind=rk) :: res(n)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: spaceCoord(1,3)
    integer :: iPoint
    type(tem_time_type) :: time
    ! -------------------------------------------------------------------- !

    !> @todo JZ: i added the subroutine for the space-time formulation.
    !! Below you see, that we create a tem_time_type variable in each loop.
    !! The code can be faster if we implement this routine in a nicer way.


    do iPoint = 1, n
      time%sim = stcoord(iPoint, 4)
      spaceCoord(1,1:3) = stcoord(iPoint, 1:3)
      res(iPoint:iPoint)                                                       &
        & = tem_spacetime_for_coord( me    = me,         &
        &                            coord = spaceCoord, &
        &                            time  = time,       &
        &                            n     = 1           )
    end do

  end function tem_spacetime_for_stcoord
  ! ************************************************************************ !
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function computes the space time function for a given list of
  !! coordinates
  !!
  function tem_spacetime_for_coord(me, coord, time, n) result(res)
    ! -------------------------------------------------------------------- !
    !> Spacetime function to evaluate
    type(tem_spacetime_fun_type) :: me
    !> Number of values to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent(in) :: coord(n,3)
    !> timer object incl. the current time information
    type(tem_time_type), intent(in)  :: time
    !> return value of the function
    real(kind=rk) :: res(n)
    ! -------------------------------------------------------------------- !

    select case(trim(adjustl(me%fun_kind)))
    case('none')
      res = 0.0_rk
    case('const')
      res = me%const(1)
    case('lua_fun')
      res = tem_spacetime_lua_for( fun_ref = me%lua_fun_ref, &
        &                          coord   = coord,          &
        &                          time    = time,           &
        &                          n       = n,              &
        &                          conf    = me%conf         )
    case('combined')
      res = tem_spatial_for( me    = me%spatial,  &
        &                    coord = coord,       &
        &                    n     = n          ) &
        & * tem_temporal_for( temporal   = me%temporal, &
        &                     time       = time         )
    case('miescatter_displacementfieldz')
      res = tem_eval_miescatter_displz( me    = me%mie_fun, &
        &                               coord = coord,      &
        &                               time  = time%sim,   &
        &                               n     = n           )
    case('miescatter_magneticfieldx')
      res = tem_eval_miescatter_magnx( me    = me%mie_fun, &
        &                              coord = coord,      &
        &                              time  = time%sim,   &
        &                              n     = n           )
    case('miescatter_magneticfieldy')
      res = tem_eval_miescatter_magny( me    = me%mie_fun, &
        &                              coord = coord,      &
        &                              time  = time%sim,   &
        &                              n     = n           )
    case('cylindricalwave')
      res = tem_eval_cylindricalWave( me    = me%cylindricalWave, &
        &                             coord = coord,              &
        &                             time  = time%sim,           &
        &                             n     = n                   )
    case('acoustic_pulse')
      res = tem_eval_acoustic_pulse( me    = me%acoustic_pulse, &
        &                            coord = coord,             &
        &                            time  = time%sim,          &
        &                            n     = n                  )
    case('polygon_body_2d', 'polygon_body_3d')
      res = tem_polygon_material_movement_single( &
        & me     = me%polygon_material,           &
        & coord  = coord,                         &
        & time   = time%sim,                      &
        & nPoint = n                              )
    case('polygon_multi_body_2d', 'polygon_multi_body_3d')
      res = tem_polygon_material_movement_multi( &
        & me     = me%polygon_material,           &
        & coord  = coord,                         &
        & time   = time%sim,                      &
        & nPoint = n                              )
    case default
      call tem_abort('ERROR: Unknown spatial function in' &
        & // ' tem_spacetime_for_coord.'                  )
    end select

  end function tem_spacetime_for_coord
  ! ************************************************************************ !
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function computes the space time function that gives back an array
  !! for a given list of coordinates
  !!
  function tem_spacetime_vector_for_coord( me, coord, time, n, ncomp)   &
    &                             result(res)
    ! -------------------------------------------------------------------- !
    !> Spacetime function to evaluate
    type(tem_spacetime_fun_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> Number of entrys in each array
    integer, intent(in) :: ncomp
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent(in) :: coord(n,3)
    !> timer object incl. the current time information
    type(tem_time_type), intent(in)  :: time
    !> return value of the function
    real(kind=rk) :: res(n,ncomp)
    ! -------------------------------------------------------------------- !
    ! counter
    integer :: i
    real(kind=rk) :: trans
    ! -------------------------------------------------------------------- !

    select case(trim(adjustl(me%fun_kind)))
    case('none')
      res = 0.0_rk
    case('const')
      do i = 1, ncomp
        res(:,i) = me%const(i)
      end do
    case('lua_fun')
      res = tem_spacetime_lua_for( fun_ref = me%lua_fun_ref, &
        &                          coord   = coord,          &
        &                          time    = time,           &
        &                          n       = n,              &
        &                          ncomp   = ncomp,          &
        &                          conf    = me%conf         )
    case('combined')
      trans = tem_temporal_for( temporal   = me%temporal, &
        &                       time       = time         )
      res = tem_spatial_for( me    = me%spatial, &
        &                    coord = coord,      &
        &                    n     = n,          &
        &                    ncomp = ncomp       )
      res = trans*res
    case default
      call tem_abort('ERROR: Unknown spatial function in' &
        & // ' tem_spacetime_vector_for_coord.'           )
    end select

  end function tem_spacetime_vector_for_coord
  ! ************************************************************************ !
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> \brief This function invokes the Lua function for barycentric coordinates
  !! of an element specified by treeIds
  !!
  function tem_spacetime_lua_for_treeIds(fun_ref, treeIds, time, tree, n, &
    &                                    conf  )      result(res)
    ! -------------------------------------------------------------------- !
    !> Reference of the function to open
    integer, intent(in) :: fun_ref
    !> global treelm mesh
    type(treelmesh_type), intent(in) ::tree
    !> number of return values
    integer, intent(in) :: n
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> timer object incl. the current time information
    type(tem_time_type), intent(in)  :: time
    !> return value
    real(kind=rk) :: res(n)
    !> lua state
    type(flu_State), intent(in) :: conf
    ! -------------------------------------------------------------------- !
    type(aot_fun_type) :: fun
    integer :: iError
    integer :: i, j
    real(kind=rk) :: coord(3)
    ! -------------------------------------------------------------------- !

    call aot_fun_open(L = conf, fun = fun, ref = fun_ref)

    do i=1,n
      coord = tem_BaryOfId( tree, treeIds(i) )
      do j=1,3
        call aot_fun_put(L=conf, fun=fun, arg=coord(j))
      end do
      call aot_fun_put(L=conf, fun=fun, arg=time%sim)
      call aot_fun_do(L=conf, fun=fun, nresults=1)
      call aot_top_get_val(L=conf, val=res(i), ErrCode=iError)
      if ( btest(iError,aoterr_Fatal) ) then
        write(logunit(0),*) "ERROR Obtaining a space time function"
        write(logunit(0),*) "Probably wrong number of components returned"
        write(logunit(0),*) "or a scalar was return as a lua table"
        write(logunit(0),*) 'Expected nComp: 1'
        write(logunit(0),*) 'ErrorCodes: ', iError
        write(logunit(0),*) "Check return values of your Lua functions!"
        call tem_abort()
      end if
    end do

    call aot_fun_close(L = conf, fun = fun)

  end function tem_spacetime_lua_for_treeIds
  ! ************************************************************************ !
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function invokes the Lua function for barycentric coordinates
  !! of an element specified by treeIds and returns an array with the given
  !! number of components.
  !!
  !! Note, that the returned object by the Lua function has to be a table,
  !! except if there is only one component. For arrays of length 1 the Lua
  !! return value has to be a simple scalar, not a table!
  function tem_spacetime_lua_vector_for_treeIds( fun_ref, treeIds, time, tree, &
    &                                            n, ncomp, conf ) result(res)
    ! -------------------------------------------------------------------- !
    !> Reference of the function to open
    integer, intent(in) :: fun_ref
    !> global treelm mesh
    type(treelmesh_type), intent(in) ::tree
    !> number of return values
    integer, intent(in) :: n
    !> number of components per returned value
    integer, intent(in) :: ncomp
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> timer object incl. the current time information
    type(tem_time_type), intent(in)  :: time
    !> return value
    real(kind=rk) :: res(n,ncomp)
    !> lua state
    type(flu_State), intent(in) :: conf
    ! -------------------------------------------------------------------- !
    type(aot_fun_type) :: fun
    integer :: iError(ncomp)
    integer :: i, j
    real(kind=rk) :: coord(3)
    ! -------------------------------------------------------------------- !

    call aot_fun_open(L = conf, fun = fun, ref = fun_ref)

    do i=1,n
      coord = tem_BaryOfId( tree, treeIds(i) )
      do j=1,3
        call aot_fun_put(L=conf, fun=fun, arg=coord(j))
      end do
      call aot_fun_put(L=conf, fun=fun, arg=time%sim)
      call aot_fun_do(L=conf, fun=fun, nresults=1)
      if (nComp == 1) then
        call aot_top_get_val(L=conf, val=res(i,1), ErrCode=iError(1))
        if ( any(btest(iError,aoterr_Fatal)) ) then
          write(logunit(0),*) "ERROR Obtaining a space time function"
          write(logunit(0),*) "Probably wrong number of components returned"
          write(logunit(0),*) "or a scalar was return as a lua table"
          write(logunit(0),*) 'Expected nComp: ', nComp
          write(logunit(0),*) 'ErrorCodes: ', iError
          write(logunit(0),*) "Check return values of your Lua functions!"
          call tem_abort()
        end if
      else
        call aot_top_get_val(L=conf, val=res(i,:), ErrCode=iError)
        if ( any(btest(iError,aoterr_Fatal)) ) then
          write(logunit(0),*) "ERROR Obtaining a space time function"
          write(logunit(0),*) "Probably wrong number of components returned"
          write(logunit(0),*) "or a scalar was return as a lua table"
          write(logunit(0),*) 'Expected nComp: ', nComp
          write(logunit(0),*) 'ErrorCodes: ', iError
          write(logunit(0),*) "Check return values of your Lua functions!"
          call tem_abort()
        end if
      end if
    end do

    call aot_fun_close(L = conf, fun = fun)

  end function tem_spacetime_lua_vector_for_treeIds
  ! ************************************************************************ !
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> \brief This function invokes the Lua function for a given coordinate.
  !!
  function tem_spacetime_lua_for_coord( fun_ref, coord, time, n, conf ) &
    &                                 result(res)
    ! -------------------------------------------------------------------- !
    !> Reference of the function to open
    integer, intent(in) :: fun_ref
    !> number of return values
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent(in) :: coord(n,3)
    !> timer object incl. the current time information
    type(tem_time_type), intent(in)  :: time
    !> return value
    real(kind=rk) :: res(n)
    !> lua state
    type(flu_State), intent(in) :: conf
    ! -------------------------------------------------------------------- !
    type(aot_fun_type) :: fun
    integer :: iError
    integer :: i, j
    ! -------------------------------------------------------------------- !

    call aot_fun_open(L = conf, fun = fun, ref = fun_ref)

    do i=1,n
      do j=1,3
        call aot_fun_put(L=conf, fun=fun, arg=coord(i,j))
      end do
      call aot_fun_put(L=conf, fun=fun, arg=time%sim)
      call aot_fun_do(L=conf, fun=fun, nresults=1)
      call aot_top_get_val(L=conf, val=res(i), ErrCode=iError)
      if ( btest(iError,aoterr_Fatal) ) then
        write(logunit(0),*) "ERROR Obtaining a space time function"
        write(logunit(0),*) "Probably wrong number of components returned"
        write(logunit(0),*) "or a scalar was return as a lua table"
        write(logunit(0),*) 'Expected nComp: 1'
        write(logunit(0),*) 'ErrorCodes: ', iError
        write(logunit(0),*) "Check return values of your Lua functions!"
        call tem_abort()
      end if
    end do

    call aot_fun_close(L = conf, fun = fun)

  end function tem_spacetime_lua_for_coord
  ! ************************************************************************ !
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function invokes the Lua function for a given coordinate and returns
  !! an array valued result.
  !!
  !! Note, that the returned object by the Lua function has to be a table,
  !! except if there is only one component. For arrays of length 1 the Lua
  !! return value has to be a simple scalar, not a table!
  function tem_spacetime_lua_vector_for_coord( fun_ref, coord, time, n, ncomp, &
    &                                          conf ) result(res)
    ! -------------------------------------------------------------------- !
    !> Reference of the function to open
    integer, intent(in) :: fun_ref
    !> number of return values
    integer, intent(in) :: n
    !> number of components returned for each value
    integer, intent(in) :: ncomp
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent(in) :: coord(n,3)
    !> timer object incl. the current time information
    type(tem_time_type), intent(in)  :: time
    !> return value
    real(kind=rk) :: res(n,ncomp)
    !> lua state
    type(flu_State), intent(in), optional :: conf
    ! -------------------------------------------------------------------- !
    type(aot_fun_type) :: fun
    integer :: iError(ncomp)
    integer :: i, j
    ! -------------------------------------------------------------------- !
    call aot_fun_open(L = conf, fun = fun, ref = fun_ref)

    do i=1,n
      do j=1,3
        call aot_fun_put(L=conf, fun=fun, arg=coord(i,j))
      end do
      call aot_fun_put(L=conf, fun=fun, arg=time%sim)
      call aot_fun_do(L=conf, fun=fun, nresults=1)
      call aot_top_get_val(L=conf, val=res(i,:), ErrCode=iError)
      if ( any(btest(iError,aoterr_Fatal)) ) then
        write(logunit(0),*) "ERROR Obtaining a space time function"
        write(logunit(0),*) "Probably wrong number of components returned"
        write(logunit(0),*) "or a scalar was return as a lua table"
        write(logunit(0),*) 'Expected nComp: ', nComp
        write(logunit(0),*) 'ErrorCodes: ', iError
        write(logunit(0),*) "Check return values of your Lua functions!"
        call tem_abort()
      end if
    end do

    call aot_fun_close(L = conf, fun = fun)

  end function tem_spacetime_lua_vector_for_coord
  ! ************************************************************************ !
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> This function returns pre-stored value at given idx if spacetime function
  !! is predefined apesmate else evaluate a spacetime
  !! function for a point at given idx in growing array of points.
  !! Return value is a scalar.
  function tem_spacetime_scalar_for_index( me, grwPnt, idx, nVals, iLevel, &
    &                                      time ) result (res)
    ! -------------------------------------------------------------------- !
    !> spacetime type
    type(tem_spacetime_fun_type), intent(in) :: me
    !> number of return values
    integer, intent(in) :: nVals
    !> growing array of all spacetime point of a variable
    type(tem_grwPoints_type), intent(in) :: grwPnt
    !> Index position to return a pre-store value or to compute
    integer, intent(in) :: idx(nVals)
    !> return value of a function
    real( kind=rk ) :: res(nVals)
    !> Level to access stored value in aps_coupling
    integer, intent(in) :: iLevel
    !> timer object incl. the current time information
    type(tem_time_type), intent(in)  :: time
    ! -------------------------------------------------------------------- !
    integer :: iVal
    real(kind=rk) :: coord(1,3), res_tmp(1), trans
    ! -------------------------------------------------------------------- !
    select case (trim(me%fun_kind))
    case ('none')
      res = 0.0
    case ('const')
      res = me%const(1)
    case ('lua_fun')
      do iVal = 1, nVals
        coord(1,:) =  (/ grwPnt%coordX%val( idx(iVal) ), &
          &              grwPnt%coordY%val( idx(iVal) ), &
          &              grwPnt%coordZ%val( idx(iVal) ) /)

        res_tmp = tem_spacetime_lua_for( fun_ref = me%lua_fun_ref, &
          &                              coord   = coord,          &
          &                              time    = time,           &
          &                              n       = 1,              &
          &                              conf    = me%conf         )

        res(iVal) = res_tmp(1)
      end do
    case ('combined')
      trans = tem_temporal_for( temporal   = me%temporal, &
        &                       time       = time         )
      res = tem_spatial_for( me     = me%spatial, &
        &                    grwPnt = grwPnt,     &
        &                    idx    = idx,        &
        &                    nVals  = nVals,      &
        &                    iLevel = iLevel      )
      res = trans*res
    case ('apesmate')
      res(1:nVals) = me%aps_coupling%valOnLvl(iLevel)      &
        &                           %evalVal( idx(1:nVals) )
    case ('precice')
      res(1:nVals) = tem_precice_read(                            &
        & dataID   = me%precice_coupling%readVar%IDs(1),          &
        & npoints  = nVals,                                       &
        & posIDs   = me%precice_coupling%readVar%posIDLvl(iLevel) &
        &                              %posIDs(idx(:))            )
    case default
      do iVal = 1, nVals
        coord(1,:) =  (/ grwPnt%coordX%val( idx(iVal) ), &
          &              grwPnt%coordY%val( idx(iVal) ), &
          &              grwPnt%coordZ%val( idx(iVal) ) /)

        res_tmp = tem_spacetime_for_coord( me    = me,    &
          &                                coord = coord, &
          &                                time  = time,  &
          &                                n     = 1      )
        res(iVal) = res_tmp(1)
      end do
    end select

  end function tem_spacetime_scalar_for_index
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function returns pre-stored value at given idx if spacetime function
  !! is predefined apesmate else evaluate a spacetime
  !! function for a point at given idx in growing array of points.
  !! Return value is a vector.
  function tem_spacetime_vector_for_index( me, grwPnt, idx, nVals, iLevel, &
    &                                      time, nComps ) result (res)
    ! -------------------------------------------------------------------- !
    !> spacetime type
    type(tem_spacetime_fun_type), intent(in) :: me
    !> number of return values
    integer, intent(in) :: nVals
    !> number of components per returned value
    integer, intent(in) :: nComps
    !> growing array of all spacetime point of a variable
    type(tem_grwPoints_type), intent(in) :: grwPnt
    !> Index position to return a pre-store value or to compute
    integer, intent(in) :: idx(nVals)
    !> return value of a function
    real( kind=rk ) :: res(nVals, nComps)
    !> Level to which the evaluated values to be returned
    integer, intent(in) :: iLevel
    !> timer object incl. the current time information
    type(tem_time_type), intent(in)  :: time
    ! -------------------------------------------------------------------- !
    integer :: iVal, iComp, iVar,offset
    real(kind=rk) :: coord(1,3), res_tmp(1, nComps), trans
    real(kind=rk), allocatable :: temp(:)
    ! -------------------------------------------------------------------- !
    select case (trim(me%fun_kind))
    case ('none')
      res = 0.0_rk
    case ('const')
      do iComp = 1, nComps
        res(:, iComp) = me%const(iComp)
      end do
    case ('lua_fun')
      do iVal = 1, nVals
        coord(1,:) =  (/ grwPnt%coordX%val( idx(iVal) ), &
          &              grwPnt%coordY%val( idx(iVal) ), &
          &              grwPnt%coordZ%val( idx(iVal) ) /)

        res_tmp = tem_spacetime_lua_for( fun_ref = me%lua_fun_ref, &
          &                              coord   = coord,          &
          &                              time    = time,           &
          &                              n       = 1,              &
          &                              nComp   = nComps,         &
          &                              conf    = me%conf         )

        res(iVal,:) = res_tmp(1,:)
      end do
    case ('combined')
      trans = tem_temporal_for( temporal   = me%temporal, &
        &                       time       = time         )
      res = tem_spatial_for( me     = me%spatial, &
        &                    grwPnt = grwPnt,     &
        &                    idx    = idx,        &
        &                    nVals  = nVals,      &
        &                    iLevel = iLevel,     &
        &                    nComps = nComps      )
      res = trans*res
    case ('apesmate')
      do iVal = 1, nVals
        offset = (idx(iVal)-1)*nComps
        res(iVal, :) = me%aps_coupling%valOnLvl(iLevel)                  &
          &                           %evalVal( offset+1 : offset+nComps )
      end do
    case('precice')
      allocate(temp(nVals))
      do iVar = 1, me%precice_coupling%readVar%nVars
        temp(1:nVals) = tem_precice_read(                          &
          & dataID  = me%precice_coupling%readVar%IDs(iVar),       &
          & posIDs  = me%precice_coupling%readVar%posIDLvl(iLevel) &
          &                                      %posIDs(idx(:)),  &
          & npoints = nVals                                        )
        res(1:nVals, iVar) = temp
      end do
      deallocate(temp)
    case default
      do iVal = 1, nVals
        coord(1,:) =  (/ grwPnt%coordX%val( idx(iVal) ), &
          &              grwPnt%coordY%val( idx(iVal) ), &
          &              grwPnt%coordZ%val( idx(iVal) ) /)

        res_tmp = tem_spacetime_vector_for_coord( me    = me,    &
          &                                       coord = coord, &
          &                                       time  = time,  &
          &                                       n     = 1,     &
          &                                       nComp = nComps )
        res(iVal,:) = res_tmp(1,:)
      end do
    end select

  end function tem_spacetime_vector_for_index
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function create unique id to create anonymous variable in
  !! tem_variable_loadMapping
  function tem_spacetime_hash_id(me, conf) result(id)
    ! -------------------------------------------------------------------- !
    type(tem_spacetime_fun_type), intent(in) :: me
    type(flu_State), intent(in) :: conf
    character(len=labelLen) :: id
    ! -------------------------------------------------------------------- !
    type(aot_fun_type) :: fun
    integer :: i
    character(len=labelLen) :: tmp
    real :: rnd
    ! -------------------------------------------------------------------- !

    id = me%fun_kind
    select case(trim(me%fun_kind))
    case('none')
      id = trim(id) // 'NONE'

    case('const')
      id = trim(id) // 'const:'
      do i=1,ubound(me%const,1)
        write(tmp,'(en17.8)') me%const(i)
        id = trim(id) // trim(adjustl(tmp))
      end do

    case('lua_fun')
      id = trim(id) // 'lua_fun:'
      call aot_fun_open(L = conf, fun = fun, ref = me%lua_fun_ref)
      tmp = trim(aot_fun_id(fun))
      call aot_fun_close(L = conf, fun = fun)
      id = trim(id) // trim(tmp)

    case('combined')
      id = trim(id) // 'combined-TODO:'
      call random_number(rnd)
      write(tmp, '(en17.8)') rnd
      id = trim(id) // trim(tmp)
    case('miescatter_displacementfieldz')
      id = trim(id) // 'mieZ-TODO:'
      call random_number(rnd)
      write(tmp, '(en17.8)') rnd
      id = trim(id) // trim(tmp)
    case('miescatter_magneticfieldx')
      id = trim(id) // 'mieX-TODO:'
      call random_number(rnd)
      write(tmp, '(en17.8)') rnd
      id = trim(id) // trim(tmp)
    case('miescatter_magneticfieldy')
      id = trim(id) // 'mieY-TODO:'
      call random_number(rnd)
      write(tmp, '(en17.8)') rnd
      id = trim(id) // trim(tmp)
    case('cylindricalWave')
      id = trim(id) // 'cylwav-TODO:'
      call random_number(rnd)
      write(tmp, '(en17.8)') rnd
      id = trim(id) // trim(tmp)
    case('apesmate')
      id = trim(id) // 'apesmate-TODO:'
      call random_number(rnd)
      write(tmp, '(en17.8)') rnd
      id = trim(id) // trim(tmp)
    case('precice')
      id = trim(id) // 'precice-TODO:'
      call random_number(rnd)
      write(tmp, '(en17.8)') rnd
      id = trim(id) // trim(tmp)
    case default
      id = trim(id) // 'UNKNOWN:'
      call random_number(rnd)
      write(tmp, '(en17.8)') rnd
      id = trim(id) // trim(tmp)
    end select

  end function tem_spacetime_hash_id
  ! ************************************************************************ !
  ! ************************************************************************ !


end module tem_spacetime_fun_module
! **************************************************************************** !
