! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2016, 2019, 2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013, 2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2015-2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2017 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
! Copyright (c) 2018 Raphael Haupt <Raphael.Haupt@student.uni-siegen.de>
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
!> This module contains subroutines needed to create and update subTrees from
!! corresponding global trees.
!!
!! The subroutine tem_create_subTree_of creates a subTree from a global tree
!! using an array of shapes. This functionality is used when certain operations
!! (e.g. tracking, source terms) shall be executed only on a subset of the fluid
!! tree (or another tree).
!!
module tem_subTree_module

  ! include treelm modules
  use mpi
  use env_module,              only: long_k, globalMaxLevels, pathLen, &
    &                                long_k_mpi, labelLen
  use treelmesh_module,        only: treelmesh_type
  use tem_bc_prop_module,      only: tem_bc_prop_type
  use tem_aux_module,          only: tem_abort
  use tem_tools_module,        only: tem_horizontalSpacer
  use tem_debug_module,        only: tem_debug_type, main_debug
  use tem_dyn_array_module,    only: dyn_intArray_type, init, append, destroy
  use tem_grow_array_module,   only: grw_intArray_type, init, append, destroy
  use tem_property_module,     only: gather_property, prp_hasQVal
  use tem_prophead_module,     only: tem_prop_countelems
  use tem_shape_module,        only: tem_shape_type,   &
    &                                tem_global_shape, &
    &                                tem_local_shape,  &
    &                                tem_shape2subTree
  use tem_geometry_module,     only: tem_setEffBoundingBox, tem_baryOfId
  use tem_construction_module, only: tem_levelDesc_type
  use tem_subTree_type_module, only: tem_subTree_type, &
    &                                tem_dump_subTree
  use tem_stencil_module,      only: tem_stencilHeader_type
  use tem_logging_module,      only: logUnit, tem_log, tem_toStr
  use tem_pointData_module,    only: tem_grwPoints_type, destroy

  implicit none

  private

  public :: tem_updatePropertyBits
  public :: tem_create_subTree_of
  public :: tem_write_debugMesh
  public :: tem_create_tree_from_sub
  public :: tem_subTree_from

  !> interface for copying the property bits from the global tree or the
  !! level descriptor
  interface tem_copyPropertyBits
    module procedure tem_copyPropertyBitsFromLevelDesc
    module procedure tem_copyPropertyBitsFromTree
  end interface

  contains

! ****************************************************************************** !
  !> Update the property of subTree with the ones from inTree, if something
  !! changed update the logical meshChange.
  !!
  subroutine tem_updatePropertyBits( inTree, subTree )
    ! ---------------------------------------------------------------------------
    !> tree to get information from
    type(treelmesh_type), intent(in) :: inTree
    !> tree to pass information to
    type(tem_subTree_type), intent(inout) :: subTree
    ! ---------------------------------------------------------------------------
    integer :: iElem
    integer :: elemPos
    ! ---------------------------------------------------------------------------

    do iElem = 1, subTree%nElems
      ! get the position in the global tree
      elemPos = subTree%map2global(iElem)
      ! test the individual property bits
      if ( inTree%ElemPropertyBits(elempos)   &
        &  /= subTree%ElemPropertyBits(iElem) ) then
        ! update the propertyBits in the subTree
        subTree%ElemPropertyBits(iElem) = inTree%ElemPropertyBits(elemPos)
        subTree%global%meshChange = .true.
      end if
    end do

  end subroutine tem_updatePropertyBits
! ****************************************************************************** !


! ****************************************************************************** !
  !> Copy the properties of the level decriptor to the ones in subTree.
  !!
  subroutine tem_copyPropertyBitsFromLevelDesc( levelDesc, subTree )
    ! ---------------------------------------------------------------------------
    !> level descriptor including all elements (fluid, ghost, halo)
    type(tem_levelDesc_type ),intent(in) :: levelDesc(:)
    !> tree to pass information to
    type(tem_subTree_type), intent(inout) :: subTree
    ! ---------------------------------------------------------------------------
    integer :: iLevel
    integer :: nElems
    integer :: elemCounter
    ! ---------------------------------------------------------------------------
    allocate( subTree%ElemPropertyBits( subTree%nElems ))

    elemCounter = 0
    do iLevel = 1, size(levelDesc)
      nElems = size(levelDesc(iLevel)%property(:))
      subTree%ElemPropertyBits(elemCounter+1 : elemCounter+nElems) =           &
        &                                         levelDesc(iLevel)%property(:)
      elemCounter = elemCounter + nElems
    end do
  end subroutine tem_copyPropertyBitsFromLevelDesc
! ****************************************************************************** !


! ****************************************************************************** !
  !> Copy the properties of inTree to the ones in subTree.
  !!
  subroutine tem_copyPropertyBitsFromTree( inTree, subTree )
    ! ---------------------------------------------------------------------------
    !> tree to get information from
    type(treelmesh_type), intent(in) :: inTree
    !> tree to pass information to
    type(tem_subTree_type), intent(inout) :: subTree
    ! ---------------------------------------------------------------------------
    integer :: iElem
    integer :: elemPos
    integer :: iProp
    ! ---------------------------------------------------------------------------
    ! allocate the bitfield
    allocate(subTree%elemPropertyBits(subTree%nElems))

    ! initialize the bitfield with 0
    subTree%elemPropertyBits = 0

    do iElem = 1, subTree%nElems
      ! get the position in the global tree
      elemPos = subTree%map2global(iElem)
      ! copy the propertyBits in the subTree
      subTree%ElemPropertyBits(iElem) = inTree%ElemPropertyBits(elemPos)
    end do

    ! update the property information
    allocate(subTree%property(subTree%global%nProperties))

    do iProp=1,subTree%global%nProperties
      call gather_Property( Property = subTree%Property(iProp),                &
        &                   Header   = subTree%global%Property(iProp),         &
        &                   BitField = subTree%ElemPropertyBits,               &
        &                   comm     = subTree%global%comm )
    end do

  end subroutine tem_copyPropertyBitsFromTree
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine creates a subtree based on a provided map or list of
  !! treeIDs (in case a local shape is used) to the corresponding tree.
  !! Only processes in comm will be involved.
  subroutine tem_subTree_from( me, map2global, treeID, comm, dirname, grwPnts )
    ! ---------------------------------------------------------------------------
    !> subTree to be created from list of elements (map2global)
    type(tem_subTree_type), intent(inout) :: me
    !> position of the treeID in the global treeID list
    integer, optional, intent(in) :: map2global(:)
    !> list of treeIDs only use this in case a local shape is set
    integer(kind=long_k), optional, intent(in) :: treeID(:)
    !> mpi communicator to use, defaults to the one in me%global%comm if not
    !! specified
    integer, intent(in), optional :: comm
    !> directory to store the mesh in. is taken to be me%global%dirname if not
    !! specified
    character(len=*), intent(in), optional :: dirname
    !> array of point vaues that neeeds to be stored in the subtree
    type(tem_grwPoints_type), intent(in), optional :: grwPnts
    ! ---------------------------------------------------------------------------
    integer :: commloc
    integer(kind=long_k) :: offset, nElems, nPoints
    integer :: ierror
    integer :: nElemsList
    ! mpi variables
    integer :: commsize, rank, iPnt
    ! ---------------------------------------------------------------------------

    if (present(dirname)) then
      me%global%dirname = dirname
    end if

    if (present(comm)) then
      commloc = comm
    else
      commloc = me%global%comm
    end if

    ! allocate and store the points for get point tracking if
    ! applicable
    if (present(grwPnts)) then
      me%nPoints = grwPnts%coordX%nVals
      allocate(me%points(me%nPoints,3))
      do iPnt = 1, me%nPoints
        me%points(:, 1) = grwPnts%coordX%val(1:me%nPoints)
        me%points(:, 2) = grwPnts%coordY%val(1:me%nPoints)
        me%points(:, 3) = grwPnts%coordZ%val(1:me%nPoints)
      end do
    end if

    ! determine number of elements in the provided map2global
    nElemsList = 0
    if(present(map2global) .and. .not. present(treeID))then
      nElemsList  = size(map2global)
      allocate(me%map2global( nElemsList ))
      ! copy the map2global
      me%map2global = map2global
    else if(present(treeID) .and. .not. present(map2global))then
      nElemsList  = size(treeID)
      allocate(me%treeID( nElemsList ))
      ! copy the list of treeID
      me%treeID = treeID
      ! in case a treeID is set we have to set useLocalMesh to true
      me%useLocalMesh = .true.
    else
      write(logUnit(1),*)'error: either none or both of treeID or map2global '
      write(logUnit(1),*)'       are passed to tem_subTree_from. '//           &
        &                '       only one of them is allowed'
      call tem_abort()
    end if

    ! get size of communicator and my rank in it
    call mpi_comm_size(commloc, commsize, ierror)
    call mpi_comm_rank(commloc, rank, ierror)
    me%global%nparts = commsize
    me%global%mypart = rank

    ! assign treeID list to the new mesh. sort first
    offset = 0
    nElems  = int(nElemsList, kind=long_k)
    ! determine offset for each process (get nElems of ranks before)
    call MPI_Exscan( nElems, offset, 1, long_k_mpi, mpi_sum, commloc, ierror )
    me%elemOffset = offset
    me%nElems = nElemsList

    ! last process offset + its nElems is total number of elems
    nElems = offset + nElems

    ! distribute among new communicator
    call MPI_Bcast( nElems, 1, long_k_mpi, commsize-1, commloc, ierror )
    ! store in global number of elements
    me%global%nElems = nElems

    if (present(grwPnts)) then
      nPoints = me%nPoints
      call MPI_Reduce(nPoints, me%glob_nPoints, 1, long_k_mpi, mpi_sum, &
        &             0, commloc, ierror)
    end if

  end subroutine tem_subTree_from
! ****************************************************************************** !


! ****************************************************************************** !
  !> Create a subTree from a given inTree and an array of shapes.
  !!
  !! All elements to the inTree%treeID are collected into the subTree%map2global
  !! array and a sub-communicator with participating processes is created.
  !! Additionally how many and which elements exist on my local process and are
  !! requested from the shapes is identified.
  !!
  subroutine tem_create_subTree_of( inTree, subTree, inShape, levelDesc,       &
    &                               bc_prop, storePnts, stencil, prefix )
    ! ---------------------------------------------------------------------------
    !> Global mesh from which the elements are identified and then stored to
    type(treelmesh_type), intent(in)                     :: inTree
    !> new mesh
    type(tem_subTree_type), intent(out)                  :: subTree
    !> shape objects on which to work
    type(tem_shape_type),intent(in)                      :: inShape(:)
    !> optional level descriptor needed for local shape
    type(tem_levelDesc_type ), optional, intent(in)      :: levelDesc(:)
    !> bc property which is used to identify elements belong to certain BCs
    type( tem_bc_prop_type ), optional, intent(in)       :: bc_prop
    !> To store space points in subTree
    logical, optional, intent(in)                        :: storePnts
    !> stencil used to find bcID on certain links
    type( tem_stencilHeader_type ), optional, intent(in) :: stencil
    !> prefix for the subTree
    character(len=*), optional, intent(in)               :: prefix
    ! ---------------------------------------------------------------------------
    ! local growing array for the map to the global tree
    type(dyn_intArray_type) :: local_map2global
    type(tem_grwPoints_type) :: local_grwPnts
    integer, allocatable :: map2globalTemp(:)
    integer(kind=long_k), allocatable :: tmp_treeID(:)
    integer :: iShape, iElem
    ! mpi variables for parallel writing of the log-wise trees
    integer :: color, iError
    ! if the local rank has part of the tracking object
    logical :: participateInMesh, globalParticipateInMesh
    ! if any element be found for this iShape
    integer :: local_countElems( globalMaxLevels )
    integer :: local_countPnts
    logical :: local_storePnts
    ! --------------------------------------------------------------------------

    if (present(storePnts)) then
      local_storePnts = storePnts
    else
      local_storePnts = .false.
    end if

    subTree%useGlobalMesh = any(inShape(:)%shapeID == tem_global_shape)
    subTree%useLocalMesh = any(inShape(:)%shapeID == tem_local_shape)
    subTree%created_new_comm = .false.

    if (subTree%useGlobalMesh) then
      write(logUnit(5),*) 'Subtree encompasses all fluid elements'
    else
      if (subTree%useLocalMesh) then
        if (present(levelDesc)) then
          write(logUnit(5),*) 'Subtree encompasses all elements of local domain'
        else
          write(logUnit(1),*) 'ERROR:levelDesc is not passed to ' &
            &                 // 'tem_create_subTree_of'
          write(logUnit(1),*) '      levelDesc is required to create ' &
            &                 // 'subTree for tem_local_shape!'
          call tem_abort()
        end if
      end if
    end if

    ! Prepare output structure
    ! init here the sub-tree mesh
    subTree%global = inTree%global
    ! initialize the number of elements on this partition with 0
    subTree%nElems = 0
    subTree%global%nElems = 0
    if ( present(prefix) ) subTree%global%dirname = trim(prefix)

    check_shapes: if (.not. (subTree%useGlobalMesh .or. subTree%useLocalMesh) ) then
      local_countElems = 0
      local_countPnts  = 0
      call init( local_map2global, length = 128)

      do iShape = 1, size( inShape )

        call tem_shape2subTree(me         = inShape(iShape),  &
          &                    iShape     = iShape,           &
          &                    inTree     = inTree,           &
          &                    storePnts  = local_storePnts,  &
          &                    map2global = local_map2global, &
          &                    grwPnts    = local_grwPnts,    &
          &                    countElems = local_countElems, &
          &                    countPnts  = local_countPnts,  &
          &                    bcIDs      = subTree%bc_ID,    &
          &                    bc_prop    = bc_prop,          &
          &                    stencil    = stencil           )

      end do ! iShape

      participateInMesh = ( sum(local_countElems) > 0 .or. local_countPnts > 0 )

      ! abort if no elements are found for the given shape
      call MPI_ALLREDUCE( participateInMesh, globalParticipateInMesh, 1,   &
        &                 MPI_LOGICAL, MPI_LOR, inTree%global%comm, iError )
      if (.not. globalParticipateInMesh) then
        write(logUnit(1),*) 'Error: No elements found for ' &
          &                 // trim(subTree%global%dirname) // ' in the mesh'
        call tem_abort()
      end if

    else check_shapes
      ! For local or global shapes we always participate
      participateInMesh = .true.
    end if check_shapes

    if ( .not. subTree%useGlobalMesh ) then
      subTree%global%predefined = ''
    end if

    ! Exchange information about the mesh
    ! Now create a new communicator with ranks that own parts of the sub-mesh
    ! This is done with the color
    if (participateInMesh) then
      color = 1
    else
      ! Not participating, set color to undefined
      color = MPI_UNDEFINED
    end if

    ! Input rank as a key to keep sorting as in global communicator
    ! reorder process ranks into tracking communicator
    if (subTree%useLocalMesh) then
      subTree%global%comm = MPI_COMM_SELF
      call tem_shape_initLocal(levelDesc, tmp_treeID)
    else
      call mpi_comm_split( inTree%global%comm, color, inTree%global%myPart, &
        &                  subTree%global%comm, iError                      )
      if (iError .ne. mpi_success) then
        write(logUnit(1),*)' There was an error splitting the communicator'
        write(logUnit(1),*)' for the tracking subset. Aborting...'
        call tem_abort()
      end if
      subTree%created_new_comm = .true.
    end if

    if (participateInMesh) then
      call mpi_comm_size( subTree%global%comm, subTree%global%nParts, iError )
      call mpi_comm_rank( subTree%global%comm, subTree%global%myPart, iError )

      if (.not. subTree%useGlobalmesh) then

        if (subTree%useLocalMesh) then
          call tem_subTree_from( me     = subTree,   &
            &                    treeID = tmp_treeID )
          ! update the property list of the newly defined tree according to the
          ! the local shape (simply copy the bits)
          call tem_copyPropertyBits( levelDesc = levelDesc, subTree = subTree )
        else ! NOT useLocalMesh
          ! Copy the created dynamic local_tIDlist to a long integer array
          allocate( map2globalTemp(local_map2global%nVals) )
          ! Assign list of links to the temporary list, first sorting
          do iElem = 1, local_map2global%nVals
            map2globalTemp(iElem) = local_map2global                       &
              &                       %val( local_map2global%sorted(iElem) )
          end do

          call tem_subTree_from( me         = subTree,        &
            &                    map2global = map2globalTemp, &
            &                    grwPnts    = local_grwPnts   )
          call destroy( me = local_grwPnts )

          ! copy property bits from the global tree
          call tem_copyPropertyBits( inTree = inTree, subTree = subTree )
         deallocate( map2globalTemp )
        end if ! useLocalMesh with level desc

        ! Calculate the effective extents of the tracking sub-tree bounding box
        call tem_setEffBoundingBox( subTree, inTree )

      else ! useGlobalMesh

        ! assign the global properties to the subTree
        subTree%nElems = inTree%nElems
        subTree%global%nElems = inTree%global%nElems

      end if ! not global mesh?

    else ! NOT participating in subtree
      subTree%nElems = 0
      subTree%global%nElems = 0
      allocate( subTree%map2global(0) )
      subTree%nPoints = 0
      subTree%glob_nPoints = 0_long_k
      allocate( subTree%points(0,0) )
    end if ! Participate in subtree?

    write(logUnit(5),"(A,I0)") 'Found corresponding nElems: ', &
      &                        subTree%global%nElems

    call destroy( me = local_map2global )
    call tem_horizontalSpacer( fUnit= logUnit(5) )

  end subroutine tem_create_subTree_of
! ****************************************************************************** !


! ****************************************************************************** !
  !> Create newtree out of intree by restricting to the elements of subtree.
  !!
  !! The new mesh will have no properties
  subroutine tem_create_tree_from_sub(intree, subtree, newtree, keep_props)
    !> The tree on which the subtree is defined.
    type(treelmesh_type), intent(in) :: intree

    !> Subtree describing the part of the mesh to create a new mesh from.
    type(tem_subtree_type), intent(in) :: subtree

    !> Resulting new tree with the elements selected by subtree from newtree.
    type(treelmesh_type), intent(out) :: newtree

    !> Flag to indicate whether to keep properties from intree also in newtree.
    !!
    !! If this is true, the properties will be copied from the intree to the
    !! newtree. An actual copy is done, as we can not rely on the pointer
    !! targets in intree to exist further on.
    !! Default is `.false.`, which means all properties will be dropped and
    !! newtree will have no properties at all.
    logical, optional, intent(in) :: keep_props

    logical :: withprop

    integer(kind=long_k) :: nNewElems
    integer :: iProp
    integer :: iError

    withprop = .false.
    if (present(keep_props)) withprop = keep_props

    newtree%global%maxlevel = intree%global%maxlevel
    newtree%global%minlevel = intree%global%minlevel
    newtree%global%origin = intree%global%origin
    newtree%global%BoundingCubeLength = intree%global%BoundingCubeLength

    if (subtree%useGlobalMesh) then

      ! Copy complete tree, but ignore properties.
      newtree = intree

      nullify(newtree%global%property)
      nullify(newtree%property)

    else

      newtree%nelems = subtree%nElems
      nullify(newtree%global%property)
      nullify(newtree%property)
      newtree%global%comm = subtree%global%comm
      newtree%global%nparts = subtree%global%nparts
      newtree%global%myPart = subtree%global%myPart
      allocate(newtree%treeID(newtree%nelems))
      allocate(newtree%ElemPropertyBits(newtree%nelems))
      allocate(newtree%Part_First(newtree%global%nparts))
      allocate(newtree%Part_Last(newtree%global%nparts))
      newtree%treeID = intree%treeID(subtree%map2global)
      if (withprop) then
        newtree%ElemPropertyBits = intree%ElemPropertyBits(subtree%map2global)
      else
        newtree%ElemPropertyBits = 0_long_k
      end if
      nNewElems = int(newtree%nElems, kind=long_k)
      ! Overall number of elements in the new mesh and offsets.
      call MPI_Exscan(nNewelems, newtree%ElemOffset, 1, long_k_mpi, &
        &             MPI_SUM, newtree%global%comm, iError          )
      newtree%global%nElems = newtree%ElemOffset+nNewElems
      call MPI_Bcast( newtree%global%nElems, 1, long_k_mpi, &
        &             newtree%global%nParts-1,              &
        &             newtree%global%comm, iError           )

      call MPI_Allgather( newtree%treeID(1), 1, long_k_mpi,  &
        &                 newtree%Part_First, 1, long_k_mpi, &
        &                 newtree%global%comm, iError        )

      call MPI_Allgather( newtree%treeID(newtree%nElems), 1, long_k_mpi, &
        &                 newtree%Part_Last, 1, long_k_mpi,              &
        &                 newtree%global%comm, iError                    )

    end if

    if (withprop) then
      newtree%global%nProperties = intree%global%nProperties
      allocate(newtree%global%property(newtree%global%nProperties))
      allocate(newtree%property(newtree%global%nProperties))
      newtree%global%property = intree%global%property
      do iProp=1,newtree%global%nProperties
        ! In the new mesh there may be a different number of elements with
        ! this property, recount them and update the header information
        ! accordingly.
        call tem_prop_countelems( me               = newtree%global            &
          &                                                 %Property(iProp),  &
          &                       elempropertybits = newtree%Elempropertybits, &
          &                       comm             = newtree%global%comm       )
        ! Now create the process local information on the property.
        call gather_Property( Property = newtree%Property(iProp),        &
          &                   Header   = newtree%global%Property(iProp), &
          &                   BitField = newtree%ElemPropertyBits,       &
          &                   comm     = newtree%global%comm             )
      end do
    else
      newtree%global%nProperties = 0
      allocate(newtree%global%property(0))
      allocate(newtree%property(0))
    end if

  end subroutine tem_create_tree_from_sub
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine collects the treeIDs of all elements in the level
  !! descriptor on the local partition (e.g. used in debug mesh, no general
  !! option for now).
  !!
  subroutine tem_shape_initLocal( levelDesc, treeID )
    ! ---------------------------------------------------------------------------
    !> level descriptor including all elements (fluid, ghost, halo)
    type(tem_levelDesc_type ),intent(in) :: levelDesc(:)
    !> temporary array of treeIDs
    integer(kind=long_k), allocatable, intent(out) :: treeID(:)
    ! ---------------------------------------------------------------------------
    integer :: offset
    integer :: nTotalElems
    integer :: iLevel
    integer :: iElem
    ! ---------------------------------------------------------------------------

    write(logUnit(1),*) 'Initializing shape with all elements on this partition'

    ! get the total number of elements in the level descriptor
    nTotalElems = sum(levelDesc(:)%nElems)

    ! allocate the array of treeIDs
    allocate(treeID(nTotalElems))
    offset = 0

    ! copy the treeIDs to the array
    do iLevel = 1,size(levelDesc)
      do iElem = 1, size( levelDesc( iLevel )%total )
        treeID( offset + iElem ) = levelDesc(iLevel)%total(iElem)
      end do
      offset = offset + size(levelDesc( iLevel )%total)
    end do
  end subroutine tem_shape_initLocal
! ****************************************************************************** !


! ****************************************************************************** !
  !> Write the complete mesh including fluid, ghosts and halo elements to disk
  !!
  !! Define in the [[tem_debug_module:tem_load_debug]] table as
  !!```lua
  !!   debug = { debugMode = true,
  !!             debugMesh = true }
  !!```
  !!
  subroutine tem_write_debugMesh( globtree, levelDesc, debug, myPart )
    ! ---------------------------------------------------------------------------
    !> mesh to locate the point in
    type(treelmesh_type), intent(in) :: globtree
    !> current level descriptor
    type(tem_levelDesc_type ),intent(in) :: levelDesc(:)
    !> debug info
    type(tem_debug_type ), optional, intent(inout) :: debug
    !> Partition to use on the calling process (= MPI Rank in comm)
    integer, intent(in) :: myPart
    ! ---------------------------------------------------------------------------
    type( tem_shape_type ) :: inShape(1)
    character(len=pathLen) :: dirname
    type( tem_debug_type ) :: local_debug
    ! local subTree that is dumped to disc
    type( tem_subTree_type ) :: subTree
    ! dummy variable which is required by tem_create_subTree_of but indeed not used
    type( tem_bc_prop_type ) :: bc_prop
    ! dummy variable which is required by tem_create_subTree_of but indeed not used
    type( tem_stencilHeader_type ) :: stencil
    ! ---------------------------------------------------------------------------

    if( present( debug ))then
      local_debug = debug
    else
      local_debug = main_debug
    end if

    if( trim( local_debug%debugMesh ) .ne. '')then
      write(logUnit(1),*)' Writing debug mesh to disk '
      inShape(1)%shapeID = tem_local_shape
      write(dirname ,'(a,i6.6,a)')trim( local_debug%debugMesh ),myPart,'_'

      call tem_create_subTree_of( inTree    = globtree,                        &
        &                         subTree   = subTree,                         &
        &                         inShape   = inShape,                         &
        &                         levelDesc = levelDesc,                       &
        &                         bc_prop   = bc_prop,                         &
        &                         stencil   = stencil,                         &
        &                         prefix    = trim( dirname ))

      call tem_dump_subTree( subTree, globtree, root_only = .false. )

      write(logUnit(1),*) 'Done writing debugMesh'

    endif ! output debugmesh

  end subroutine tem_write_debugMesh
! ****************************************************************************** !

end module tem_subTree_module
! ****************************************************************************** !
