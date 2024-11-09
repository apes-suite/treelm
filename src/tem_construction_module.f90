! Copyright (c) 2011-2012 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011-2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2011 Laura Didinger <l.didinger@grs-sim.de>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2012, 2014-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014-2016, 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2017-2018 Raphael Haupt <Raphael.Haupt@student.uni-siegen.de>
! Copyright (c) 2021 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
!> author: Harald Klimach
!! author: Manuel Hasert
!! author: Jens Zudrop
!! This module creates the required treeID lists, from which then the state
!! vector and the neighbor information can be constructed.
!!
!! Also, the structures for the ghost and halo elements are created.
!! For the ghost cells, this means establishing the dependencies. The source
!! element positions for each ghost is stored. For the halo elements, the
!! information from where to get the element information is stored in the mpi
!! buffers.
!! See the example for stencil construction in the Documentation.
!!
module tem_construction_module

  ! include treelm modules
  use mpi
  use env_module,            only: long_k, rk, globalMaxLevels, &
    &                              io_buffer_size, pathLen
  use tem_aux_module,        only: tem_abort
  use tem_param_module,      only: childPosition
  use tem_dyn_array_module,  only: dyn_longArray_type, dyn_intArray_type, &
    &                              init, append, expand, destroy,         &
    &                              PositionOfVal
  use tem_grow_array_module, only: grw_longArray_type, grw_intArray_type, &
    &                              init, append, expand, destroy, empty
  use treelmesh_module,      only: treelmesh_type
  use tem_geometry_module,   only: tem_tIdInfo, tem_PosOfPath, tem_baryOfID
  use tem_tools_module,      only: tem_PositionInSorted, tem_horizontalSpacer, &
    &                              tem_printArray
  use tem_logging_module,    only: logUnit, tem_logging_isActive, tem_toStr, &
    &                              tem_log
  use tem_debug_module,      only: tem_debug, main_debug, dbgUnit
  use tem_comm_env_module,   only: tem_comm_env_type
  use tem_comm_module,       only: tem_communication_type, tem_comm_dumpType, &
    &                              tem_commPattern_type, tem_comm_init,       &
    &                              tem_comm_count
  use tem_topology_module,   only: tem_LevelOf, tem_FirstIdAtLevel,          &
    &                              tem_PathComparison, tem_ParentOf,         &
    &                              tem_directChildren, tem_childNumber,      &
    &                              tem_PathOf, tem_path_type, tem_coordOfId, &
    &                              tem_IDofCoord
  use tem_property_module,   only: prp_hasBnd, prp_solid, prp_fluid
  use tem_bc_prop_module,    only: tem_bc_prop_type
  use tem_stencil_module,    only: tem_stencilHeader_type,      &
    &                              tem_stencilElement_type,     &
    &                              tem_stencil_map_toTreelmDef, &
    &                              tem_stencil_getHeaderPos
  use tem_element_module,    only: tem_element_type, tem_element_dump,    &
    &                              tem_eTypeOfID,                         &
    &                              print_nElems,                          &
    &                              init, append, PositionOfVal, truncate, &
    &                              changeType,                            &
    &                              eT_labels,                             &
    &                              eT_fluid, eT_nonExisting,              &
    &                              eT_halo, destroy,                      &
    &                              eT_ghostFromCoarser, eT_minNumber,     &
    &                              eT_ghostFromFiner, eT_maxNumber,       &
    &                              eT_minRelevant, eT_maxRelevant,        &
    &                              eT_undefined,                          &
    &                              eT_distributedGhostFromFiner
  use tem_halo_module,       only: tem_haloList_type, tem_halo_append, &
    &                              tem_halo_init, tem_halo_destroy
  use tem_sparse_comm_module, only: tem_sparse_alltoall_int, use_sparse_alltoall

  implicit none

  private

  public :: tem_levelDesc_type
  public :: tem_levelNeighbor_type
  public :: depSource_type
  public :: tem_updateTree_properties

  public :: tem_find_allElements
  public :: tem_build_horizontalDependencies
  public :: tem_build_verticalDependencies
  public :: tem_cleanupDependencyArrays
  public :: tem_cleanup_arrays
  public :: tem_init_elemLevels
  public :: tem_treeIDinTotal
  public :: tem_find_depProc
  public :: tem_find_depProc_globSearch
  public :: tem_elemList_dump
  public :: tem_debug_HorizontalDependencies
  public :: tem_dumpTreeIDlists
  public :: identify_stencilNeigh
  public :: depSource_append

  !> limit for searching neighbors
  !! * For acoustic scaling (2 coarser elements required) set to 1
  !! * For diffusive scaling (4 coarser elements required) set to 3
  integer, save :: nestingLimit

  !> Hash to quickly identify if an element was reconstructed before.
  !! If so, it is part of the hash
  !! this array contains treeIDs of recently accessed elements
  !! allocated and used in tem_init_elemLevels
  !! also used in identify_elements
  integer(kind=long_k), allocatable, save :: hash(:)
  !! this array contains positions in element_type of recently accessed elements
  integer, allocatable, save :: hash_elemPos(:)
  !> Entries in the hash
  integer(kind=long_k) :: nHashes

  !> identification parameters for different lists
  ! integer, parameter :: prp_nLists    =  4

  !@todo replace by element%stencil%tID
  !> includes the direct neighbors of each tree ID
  type tem_levelNeighbor_type
    !> array of the neighbors in the resulting totalList.
    !! Use this one in the solver!
    !! size: stencil%QQN, nElems(to treat with stencil)
    integer,allocatable :: nghElems(:,:)
  end type tem_levelNeighbor_type

  !> Type to specify the dependencies of ghost and halo cells.
  !! E.g.: used to specify which cells have to be known to be able
  !! to interpolate a ghost/halo cell
  !! @todo: incorperate into element_type?
  type depSource_type
    !> the source level, from where the current ghost element
    !! gets the source elements for the interpolation
    integer                  ::  dependencyLevel = -1
    !> position of the source elements in the totalList
    type( grw_intArray_type ) :: elem
    !> position of the source elements in the all source elements list
    !! i.e. levelDesc( targetLevel )%sourceFromCoarser
    type( grw_intArray_type ) :: elemBuffer
    !> Interpolation weight for each source element specified above
    real(kind=rk), allocatable :: weight(:)
    real(kind=rk) :: coord(3)
    integer :: childNum
    !> Pointer to array of interpolation matrix calculated from available
    !! sources
    integer :: posInIntpMatLSF
    
    !logical, allocatable :: bitmask(:)
  end type depSource_type

  !> detailed information of a complete level of elements
  !! including all treeIDs, properties and neighbors as well as informations
  !! about ghost/halo cells and its dependencies for
  !! interpolation/reconstruction
  type tem_levelDesc_type

    ! ---------new types for element description---------------------
    type(tem_element_type) :: elem
    !> This list includes treeIDs for which additionally
    !! neighbors have to be identified
    !! constructed in tem_init_elemLevels
    !! used in routine: identify_additionalNeigh
    type(dyn_longArray_type) :: require
    type(tem_haloList_type) :: haloList
    ! ---------------------------------------------------------------

    !> list of treeIDs for this level. Consists of ordered treeIDs of
    !! first fluid, then ghost, then halo elements.
    !! total:
    !! \[
    !! \newcommand\T{\Rule{0pt}{1em}{.3em}}
    !! \begin{array}{|c|c|c|c|}
    !! \hline fluid \T & ghostFromCoarser \T & ghostFromFiner \T & halo \\\hline
    !! \end{array}
    !! \]
    !! Array size: nElems ( =   nElems_fluid+nElems_ghostFromCoarser
    !!                        + nElems_ghostFromFiner+nElems_halo  )
    !! @todo: to be replaced by growing array
    integer(kind=long_k), allocatable :: total(:)

    !> Barycenter for all treeID in total list
    !! size: nElems in total, 3
    real(kind=rk), allocatable :: baryOfTotal(:,:)

    !> pointer to elem%tID list
    !! set in routine: identify_lists
    !! used in tem_build_listHorizontalDep, assemble_lists
    !! Array size: nElems
    !! @todo: to be replaced by growing array
    integer, allocatable :: totalPnt(:)

    !> list of property bits for this level. the same order as total list
    !! array size: nElems
    !! @todo: to be replaced by growing array
    integer(kind=long_k), allocatable :: property(:)

    !> pointer from the levelDescriptor to the original treeID list
    !! ( only for fluids )
    !! array size: nElems_fluid
    !! @todo: to be replaced by growing array
    integer, allocatable :: pntTID(:)

    !> neighbor relations for all fluid elements.
    !! Dimension: number of stencils
    !! We store the positions of the neighbor elements inside the total-list.
    !! If a fluid element does not have a neighbor in a direction
    !! (e.g. because of a boundary in that direction), we store the boundary ID
    !! as negative to indicate, that it is not a regular neighbor.
    type( tem_levelNeighbor_type ), allocatable :: neigh(:)

    !> Dependencies for ghost elements
    !! To reconstruct all the data you should
    !! iterate over this list and reconstruct the ghost elements
    !! with source element information from these data types
    !! data. Up = to coarser, down = to finer
    !! array size: nElems_ghostFromFiner
    type( depSource_type ), allocatable :: depFromFiner(:)
    !> In treelm, only the parent is stored.
    !! If more sources are needed, it has to be extend in the solver.
    !! array size: nElems_ghostFromCoarser
    type( depSource_type ), allocatable :: depFromCoarser(:)

    !> Store all the source elements for the ghostFromFiner
    !! Their positions in total list on source level
    type(dyn_intArray_type) :: sourceFromFiner
    !> Store all the source elements that needed for all ghostFromCoarser
    type(dyn_intArray_type) :: sourceFromCoarser

    !> Buffer storing intermediate values of the source elements for
    !! the interpolation
    !! @todo: move into solver?
    real(kind=rk), allocatable :: intpBufFromFiner(:,:)
    real(kind=rk), allocatable :: intpBufFromCoarser(:,:)

    !> List to store interpolation from coarser ghost elements
    !! How to use:
    !!  do indElem = 1, intpFromCoarser%nVals
    !!    posInDepFromCoarser = intpFromCoarser%val( indElem )
    !!    posInTotal = depFromCoarser%elem%val( posInDepFromCoarser )
    !!  end do
    !! Size of intpFromCoarser depends on interpolation order which intern
    !! depends on available number of source elements
    type(grw_intArray_type), allocatable :: intpFromCoarser(:)
    !> List to store interpolation from finer ghost elements
    type(grw_intArray_type) :: intpFromFiner

    !> pointing to the position of boundary elements
    !! in the levelDescriptor total list
    type(dyn_intArray_type) :: bc_elemBuffer

    !> Offsets in the assembled lists for
    !! fluid (1), ghostFromCoarser(2), ghostFromFiner(3) and halo(4) elements
    !! for the assembled lists, i.e the totalList, invSorted, ...
    !! gets the values (0, nElems_fluid,
    !!   nElems_fluid+nElems_ghostCoarser,
    !!   nElems_fluid+nElems_ghostCoarser+nELems_ghostFiner)
    integer :: offset( 2, eT_minRelevant:eT_maxRelevant ) = 0

    !> Local Fluids required by remote processes
    type( tem_communication_type ) :: sendBuffer
    !> Local ghostFromCoarser required by remote processes
    type( tem_communication_type ) :: sendBufferFromCoarser
    !> Local ghostFromFiner required by remote processes
    type( tem_communication_type ) :: sendBufferFromFiner
    !> My halos which are fluids on remote processes
    type( tem_communication_type ) :: recvBuffer
    !> My halos which are ghostFromCoarser on remote processes
    type( tem_communication_type ) :: recvBufferFromCoarser
    !> My halos which are ghostFromFiner on remote processes
    type( tem_communication_type ) :: recvBufferFromFiner

    !> total number of elements
    integer :: nElems

  end type tem_levelDesc_type

  contains

  ! ------------------------------------------------------------------------ !
  !> call this routine from your geometry initialization routine in the solver
  !! create all the necessary level-wise objects, such as element lists,
  !! dependencies
  !!
  !! 1.) build all dependencies for halos and ghost which
  !! are needed for interpolation/reconstruction (including MPI communications)
  !! 2.) build the pointers for each element to its neighbors/stencil elements.
  !! All this information is stored in tem_levelDesc_type
  !!
  subroutine tem_find_allElements( tree, levelDesc, levelPointer,              &
    &                              computeStencil, commPattern, cleanup,       &
    &                              reqNesting, proc )
    ! -------------------------------------------------------------------- !
    !> the global tree
    type(treelmesh_type), intent(inout)       :: tree
    !> the level descriptor to be filled
    type(tem_levelDesc_type), intent(inout)   :: levelDesc(tree%global%minLevel:)
    !> Pointer from treeIDlist entry to level-wise fluid part of total list
    integer, allocatable, intent(out)         :: levelPointer(:)
    !> array of all stencils used in the simulation
    type(tem_stencilHeader_type)              :: computeStencil(:)
    !> the communication pattern used
    type(tem_commPattern_type), intent(in)    :: commPattern
    !> nesting level
    integer, intent(in), optional             :: reqNesting
    !> cleanup arrays afterwards?
    logical, intent(in), optional             :: cleanup
    !> Process description to use.
    type(tem_comm_env_type), intent(in)       :: proc
    ! -------------------------------------------------------------------- !
    ! Geometry and Tree related variables
    integer :: iLevel
    type( tem_path_type ), allocatable :: pathFirst(:), pathLast(:)
    logical :: doAdditional
    logical :: clean_constructionArrays
    ! -------------------------------------------------------------------- !

    call tem_horizontalSpacer( fUnit = dbgUnit(1), before = 1 )
    write(dbgUnit(1),*) 'Inside routine: tem_find_allElements'

    if( present( cleanup )) then
      clean_constructionArrays = cleanup
    else
      clean_constructionArrays = .false.
    end if
    if( present( reqNesting )) then
      nestingLimit = reqNesting
    else
      nestingLimit = 1
    end if
    write(logUnit(3),*) 'Building up the element and total list.'

    ! Build the pathlist for improved performance when finding local element
    ! positions since the tree is the same for every scheme the allocation only
    ! has to be performed once
    if ( .not. allocated(tree%pathList))then
      allocate( tree%pathList(tree%nElems) )
    end if
    ! Fill the pathList for each element in the treeID list
    tree%pathList = tem_PathOf( tree%treeID )

    ! first and last path in every process
    allocate(pathFirst( tree%global%nParts ))
    allocate(pathLast(  tree%global%nParts ))
    pathFirst = tem_PathOf( tree%Part_First )
    pathLast  = tem_PathOf( tree%Part_Last  )

    ! Step 2: build levelDesc element list including identification of neighbor
    ! elements
    call build_levelElements( levelDesc      = levelDesc,         &
      &                       tree           = tree,              &
      &                       proc           = proc,              &
      &                       stencil        = computeStencil(1), &
      &                       pathFirst      = pathFirst,         &
      &                       pathLast       = pathLast           )

    ! Step 3: assign totalPnt to elem%tID in sorted fashion and prepare
    ! haloPrc list
    do iLevel = tree%global%minLevel, tree%global%maxLevel
      call identify_lists( levelDesc(iLevel) )
    end do

    ! Step 4: Communicate halo elements
    ! Here, first local halos are send to their corresponding process to
    ! check for existence of halo elements.
    ! Each process checks whether requested element exist and returns
    ! list of available elements.
    ! With this new information halo list is redefined.
    call communicate_elements( tree, proc, levelDesc, commPattern, &
      &                        pathFirst, pathLast, computeStencil )

    ! Step 5: do additional communication if there are elements in require list
    ! which are neighbors of higher order boundaries.
    call check_additionalComm( levelDesc, proc, doAdditional, &
      &                        tree%global%minlevel           )

    ! If doAdditional then identify neighbors of higher order boundary
    ! neighbor elements.
    ! After this identification, new halo elements might have to be added so
    ! we need to communicate all halo elements again.
    if( doAdditional ) then
      ! passing only first stencil as this is the required compute stencil
      call identify_additionalNeigh( tree, proc, levelDesc,                 &
        &                            pathFirst, pathLast, computeStencil(1) )
      call communicate_elements( tree, proc, levelDesc, commPattern, &
        &                        pathFirst, pathLast, computeStencil )
    end if

    ! Step 6: assemble levelDesc total(treeID) list, property list and
    ! pntTID list in sorted fashion (fluids+ghosts+halos)
    ! which are pre-assembled in element type
    call assemble_lists( levelDesc,                                  &
      &                  tree%global%minLevel, tree%global%maxLevel, &
      &                  tree                                        )

    ! Step 7:
    call tem_build_levelPointer( levelPointer, tree, levelDesc )
    call update_elemPosToTotalPos( levelDesc, levelPointer, tree, &
      &                            computeStencil                 )
    ! Warning: Truncation introduces a memory peak because of copy
    ! operations! Better not use...
    !call truncate_lists( levelDesc, tree%global%minLevel )

    if( clean_constructionArrays ) then
      call tem_cleanupDependencyArrays( levelDesc )
    endif

    write(logUnit(3),*) 'Finished building up element and total list. '
    deallocate(pathFirst)
    deallocate(pathLast)
    deallocate(hash)
    deallocate(hash_elempos)

    write(dbgUnit(1),*) 'Leave routine: tem_find_allElements'
    call tem_horizontalSpacer( fUnit = dbgUnit(1), after = 1 )

  end subroutine tem_find_allElements
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> subroutine to find neighbours of cells
  !!
  !! Typically every element requires information from its neighbors to perform
  !! an update of its state.
  !! All such required neighbors constitute a so called [[tem_stencil_module]].
  !! [[tem_build_horizontaldependencies]], or the
  !! connectivity of elements on the same refinement
  !! level in the tree, are therefore essential for the stencils in typical
  !! mesh-based numerical schemes.
  !! We now sketch the method to find the connectivity of elements in the mesh
  !! with the help of the described mesh layout.
  !! In a distributed mesh, the first distinction has to be made with respect to
  !! processor ownership.
  !! @todo: tem_find_depNeigh does not exist any more
  !! Is the element in question tem_find_depNeigh a local or remote
  !! element?
  !!
  !! For all [[identify_local_element]], their actual
  !! position in the sparse mesh has to be identified for a given
  !! [[treelmesh_type:treeID]].
  !!
  !! A [[tem_topology_module:tem_path_type]] comparison of two nodes to
  !! decide their relative position in the ordered list
  !! of elements has to take into account the hierarchy of the mesh.
  !!
  !! As in the actual simulation only leaves will be present, and no overlap of
  !! refinement levels is allowed, this is a sufficient determination.
  !! However, such overlapping regions might be created to enable the
  !! interpolation between different levels with virtual elements.
  !! These special elements are distinguished anyway, so this does not pose any
  !! problem in the search for neighborhood relations.
  !!
  !! An element can be either identified by its
  !! [[treelmesh_type:treeID]]
  !! or by a tuple with four integer entries for the spatial coordinates and the
  !! refinement level: \((x, y, z, L)\).
  !! This coordinate fully describes the spatial shape and position of an
  !! element and is important to determine spatial relations between elements.
  !! Conversion between treeIDs and coordinates can be achieved by the routines
  !! - [[tem_topology_module:tem_coordofid]] Coordinate from TreeID
  !! - [[tem_topology_module:tem_idofcoord]] TreeID from Coordinate
  !!
  !! The spatial indices are limited by the refinement level:
  !! \( x, y, z \in \mathbf{Z}_{\ge 0}: x, y, z < 2^L \).
  !! To avoid undefined coordinates, movements through the mesh by additions to
  !! indices are done in a modulo(\(2^L\)) safeguard resulting in an periodic
  !! universe cube.
  !! This enables the movement in the mesh in horizontal direction by index
  !! alteration and translation into a
  !! [[treelmesh_type:treeID]] again.
  !! The described procedure is completely reversible, enabling the construction
  !! of the [[treelmesh_type:treeID]] for any
  !! given coordinate.
  !! Thus, the conversion between this coordinate and the serialized
  !! [[treelmesh_type:treeID]]
  !! encoding is fully described by the topology of the octree and the chosen
  !! space-filling curve ordering of elements, as explained above.
  !!
  !! With this method we can describe the neighborhood of any given element
  !! independently of the actual solver by a simple list of relative offsets.
  !! In contrast to the generic horizontal relation, the vertical relation
  !! between child and parent nodes in the tree requires an interpolation
  !! operator between different refinement levels.
  !! This interpolation usually has to take into account solver specific
  !! requirements, but is otherwise quite isolated from the numerical operation
  !! on each refinement level.
  !! The <em>TreElM library</em> offers the solver a level-wise view, as
  !! suggested by the properties described above.
  !! To find all required neighbors in the distributed octree, the solver merely
  !! has to provide its horizontal dependencies.
  !! These are described with the help of an element specific
  !! [[tem_stencil_module]].
  !! A stencil is basically a set of element-offsets \((s_x, s_y, s_z)\),
  !! describing the relative positions of all required elements for a given
  !! element.
  !!
  !!
  !! In this routine, level descriptor is allocated,
  !! all elements in tree are added into `me%elem` as fluid type,
  !! including their `property, pntTID, stencil, neighID, sourceProc`
  !!
  subroutine tem_init_elemLevels( me, boundary, tree, stencils )
    ! -------------------------------------------------------------------- !
    !> neighbor list containing all the neighbours for the
    !! cells given in treeidsubset. Result of this routine
    type(tem_levelDesc_type), allocatable, intent(out)  :: me(:)
    !> boundaries for the elements with bnd property set
    type(tem_bc_prop_type), intent(in)                  :: boundary
    !> subset of tree ids for which the neighbours will be specified
    type(treelmesh_type), intent(in)                    :: tree
    !> the given stencil
    type(tem_stencilHeader_type), intent(in)            :: stencils(:)
    ! -------------------------------------------------------------------- !
    type(tem_stencilElement_type) :: tStencil
    integer :: posInTree, nElemsBnd, iQQN, iLevel, nProcs, hashpos
    integer :: x(4), nStencils, iStencil, elemPos, nStencilElems
    integer :: indElem, minLevel, maxLevel, QQN
    integer :: initlen
    integer(kind=long_k) :: treeID
    integer(kind=long_k), allocatable :: neighIDs(:)
    integer :: addedPos
    logical :: wasAdded
    integer :: posInBCID
    ! -------------------------------------------------------------------- !

    call tem_horizontalSpacer( fUnit = dbgUnit(1), before = 1 )
    write(dbgUnit(1),*) 'Inside routine: tem_init_elemLevels '

    ! precondition for this routine: the first stencil must be compute stencil
    if ( .not. stencils(1)%useAll ) then
      write(logUnit(1),*)'tem_init_elemLevels:'
      write(logUnit(1),*)'Error: The first stencil does not have the'//    &
        &                'useAll flag set'
      write(logUnit(1),*)'       The first one is considered the compute ' &
        &            //'stencil and'
      write(logUnit(1),*)'       needs to be set for all elements.'
      call tem_abort()
    end if

    nProcs   = tree%global%nParts
    minLevel = tree%global%minLevel
    maxLevel = tree%global%maxLevel
    nElemsBnd = 0
    nStencils = size( stencils )

    ! -------------   Prepare a hash table   -----------------------------------
    ! Find a suitable hash size.
    ! The hash lookup is O(1), while the identification grows with
    ! log(nElems) and log(nProcs), yet allocation, deallocation and
    ! filling will cause some overhead of O(nHashes), this is an attempt to
    ! find a compromise of these parameters.
    nHashes = int(max(tree%nElems, nProcs, 64), kind=long_k)

    ! @todo Find the closest smaller prime for the hash size, to minimize the
    !       number of collisions?
    ! Limit the size of the hash:
    nHashes = min(nHashes, int(io_buffer_size, kind=long_k))

    allocate(hash(0:int(nHashes-1)))
    allocate(hash_elemPos(0:int(nHashes-1)))
    ! Reset the hash
    hash = -1_long_k
    hash_elemPos = -1
    ! --------------------------------------------------------------------------


    initlen = ceiling(1.1 * real(tree%nElems)/real(maxlevel-minlevel+1))
    call tem_alloc_levelDesc( me, minLevel, maxLevel, initlen, nStencils )

    write(logUnit(5),*) 'Searching the neighbors for each element ...'

    ! Loop over all given stencils
    stencilLoop: do iStencil = 1, nStencils
      write(logUnit(6), "('   on stencil: ', I0)") iStencil
      QQN = stencils(iStencil)%QQN
      allocate( neighIDs( QQN ) )

      ! map user stencil to treelm definition
      call tem_stencil_map_toTreelmDef( me = stencils( iStencil ))

      if ( stencils(iStencil)%useAll ) then
        ! This stencil belongs to all elements in the tree (fluid stencil)
        nStencilElems = tree%nElems
      else
        ! This stencil belongs to only a subset of the tree
        ! (boundary and IBM stencils)
        nStencilElems = stencils( iStencil )%nElems
      end if

      write(logUnit(6),"(A,I0)") '   nElems: ', nStencilElems
      write(logUnit(6),"(A,I0)") '   QQN: ', QQN

      ! Loop over all the elements connected to the current stencil
      ! stencils( iStencil )
      elemLoop: do indElem = 1, nStencilElems
        ! Reset the empty stencil
        ! Create the temporary stencil for assigning the positions of the
        ! neighbor treeIDs which are added in elem%neighID
        ! also store the original stencil it depends on (depends = iStencil )
        call init( me        = tStencil, &
          &        QQN       = QQN,      &
          &        headerPos = iStencil )

        if( stencils( iStencil )%useAll ) then
          ! if the stencil is defined to be connected to all elements,
          ! then simply loop over all the fluid elements
          posInTree = indElem
        else
          ! if the stencil is defined for a subset, loop over the
          ! defined subset in %elem
          posInTree = stencils( iStencil )%elem%val( indElem )
        end if
        ! assign fluid element with info from the tree list and properties
        ! from disk
        treeID = tree%treeID( posInTree )
        ! Get topological info about coordinates+level
        x = tem_CoordOfId( treeID )

        ! Find position of this element in the elemList as it might be
        ! there already from previous stencil iterations
        hashpos = int(mod( treeID, nHashes))
        if (hash(hashpos) /= treeID) then ! cache miss
          call append( me         = me( x(4) )%elem,                  &
            &          tID        = treeID,                           &
            &          pntTID     = posInTree,                        &
            &          eType      = eT_fluid,                         &
            &          nNeighIDs  = QQN,                              &
            &          sourceProc = tree%global%myPart+1,             &
            &          property   = tree%ElemPropertyBits(posInTree), &
            &          pos        = elemPos                           )
          ! for differentiating between fluid and ghost nodes
          me( x(4) )%elem%property%val(elemPos) & 
            & = ibset( me( x(4) )%elem%property%val(elemPos), prp_fluid)
          !write(dbgunit(1),*) "prp elems = ", me( x(4) )%elem%property%val(elemPos)
          !flush(dbgUnit(1))
          !stop
          hash(hashpos) = treeID
          hash_elemPos( hashpos ) = elemPos
        else ! cache hit
          elemPos = hash_elemPos( hashpos )
        end if

        posInBCID = -1
        ! Only perform the boundary check based on the treeID boundary
        ! information for stencils which work on all elements (so we can
        ! actually read out the boundary information from the tree )
        if ( iStencil == 1 ) then
          ! First check, if I have prp_hasBnd
          ! If I don't take all the elements, then I can't analyze boundaries
          if( btest( me(x(4))%elem%property%val( elemPos ), prp_hasBnd )) then
            ! If yes, there might be no neighbors
            nElemsBnd = nElemsBnd +1
            posInBCID = nElemsBnd
          end if
        end if ! useAll elements? (boundary information)

        ! calculate possible neighbors or get BCID
        call tem_calc_neighbors( posInBCID   = posInBCID,            &
          &                      boundary_ID = boundary%boundary_ID, &
          &                      stencil     = stencils(iStencil),   &
          &                      x           = x,                    &
          &                      neighIDs    = neighIDs              )

        ! append neighbors into elem%neighID or require list
        do iQQN = 1, QQN

          call append( me       = me(x(4))%elem%neighID%val( elemPos ), &
            &          val      = neighIDs(iQQN),                       &
            &          pos      = addedPos,                             &
            &          wasAdded = wasAdded )
          tStencil%tIDpos(iQQN) = addedPos

          if ( stencils(iStencil)%requireNeighNeigh .and. neighIDs(iQQN)>0 ) then
            call append( me       = me(x(4))%require, &
              &          val      = neighIDs(iQQN),   &
              &          pos      = addedPos,         &
              &          wasAdded = wasAdded )
          end if

        end do

        ! Append the temporary stencil to the element
        call append( me  = me(x(4))%elem%stencil%val( elemPos ),             &
          &          val = tStencil )
      end do elemLoop ! indElem = 1, nStencilElems

      deallocate( neighIDs )
    end do stencilLoop ! iStencil = 1, nStencils
    write(logUnit(5),*) 'Finished searching the neighbors ID for each fluid.'

    if (tem_logging_isActive(main_debug%logger, 7)) then
      do iLevel = minLevel, maxLevel
        write(dbgUnit(1),"(A,I0)") 'Dump levelDesc%elem on level: ', iLevel
        call tem_elemList_dump( me      = me( iLevel )%elem,                  &
          &                     compact = .true.,                             &
          &                     nUnit   = dbgUnit(5),                         &
          &                     stencil = .true.,                             &
          &                     string  = 'after initializing level elements' &
          &                               //' i.e. only fluids')
        if( me( iLevel )%require%nVals > 0 ) then
          write(dbgUnit(1),"(A,I0)") 'Dump levelDesc%require on level: ', iLevel
          call tem_require_dump( me      = me( iLevel )%require,               &
            &                    nUnit   = dbgUnit(5),                         &
            &                    string  = 'after initializing level elements' )
        end if
      end do
    end if

    write(dbgUnit(1),*) 'Leave  routine: tem_init_elemLevels'
    call tem_horizontalSpacer( fUnit = dbgUnit(1), after = 1 )

  end subroutine tem_init_elemLevels
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Assemble the fluid list and identify neighbor relations
  !! start with building up the ghost and halo element collection as well
  !!
  subroutine build_levelElements( levelDesc, tree, proc, stencil,       &
    &                             pathFirst, pathLast )
    ! -------------------------------------------------------------------- !
    !> the global tree
    type(treelmesh_type),        intent(in) :: tree
    !> the level descriptor to be filled
    type(tem_levelDesc_type), intent(inout) :: levelDesc( tree%global%minLevel:)
    !> Process description to use.
    type(tem_comm_env_type), intent(in)     :: proc
    !> array of all stencils used in the simulation
    type(tem_stencilHeader_type)            :: stencil
    !> first and last treeID path in every process
    type(tem_path_type), intent(in)         :: pathFirst(:), pathLast(:)
    ! -------------------------------------------------------------------- !
    integer(kind=long_k) :: neighID   ! neighboring neighID
    integer :: elemPos, minLevel, maxLevel
    integer :: iElem, iNeighElem, iLevel, iStencil
    ! position of where to read the stencil neighbor neighID in the element
    integer :: neighPos
    ! -------------------------------------------------------------------- !

    minLevel = tree%global%minLevel
    maxLevel = tree%global%maxLevel

    write(logUnit(3),*) 'Building level-wise fluid list ...'

    call tem_horizontalSpacer( fUnit = dbgUnit(1), before = 1 )
    write(dbgUnit(3),*) 'Inside routine: build_levelElements'

    ! step 1. Identify the neighbors of each elements including ghost from finer
    ! and ghost from coarser
    ! Now iterate over all the neighbor lists and add the required elements if
    ! necessary
    do iLevel = minLevel, maxLevel
      write(dbgUnit(7),*) 'Level: ', iLevel
      do iElem = 1, levelDesc( iLevel )%elem%nElems( eT_fluid )
        do iStencil = 1, levelDesc( iLevel )%elem%stencil%val( iElem )%nVals
          do iNeighElem = 1, levelDesc( iLevel )%elem                          &
            &                         %stencil%val( iElem )%val( iStencil )%QQN
            ! Same as in identify_additionalNeigh
            neighPos = levelDesc( iLevel )%elem%stencil%val( iElem )           &
              &                           %val( iStencil )%tIDpos( iNeighElem )
            neighID = levelDesc(iLevel)%elem%neighID%val(iElem)%val( neighPos )

! write(dbgUnit(7),"(6(A,I0))") 'iElem: ', iElem, &
!   &                 ', treeID: ', levelDesc( iLevel )%elem%tID%val( iElem ), &
!   &                 ', iStencil: ', iStencil, &
!   &                 ', iNeighElem: ', iNeighElem, &
!   &                 ', neighPos: ', neighPos, &
!   &                 ', neighID: ', neighID

            if( neighID > 0_long_k ) then
              ! identify the process in which requested neighID is located
              !   - if remote process add to elem list as halo element.
              !   - if different level add to elem list as ghostFromFiner or
              !     ghostFromCoarser.
              ! Also return position (elemPos) of neighID in elem%tID list
              call identify_elements( treeID    = neighID,      &
                &                     tree      = tree,         &
                &                     pathFirst = pathFirst,    &
                &                     pathLast  = pathLast,     &
                &                     levelDesc = levelDesc,    &
                &                     elemPos   = elemPos,      &
                &                     proc      = proc,         &
                &                     nesting   = 0,            &
                &                     stencil   = stencil )
            else ! neighID < 0, i.e. it is a bcID
              elemPos = int( neighID )
            end if ! neighPos > 0
            levelDesc( iLevel )%elem%stencil%val( iElem )%val( iStencil )      &
              &                               %totalPos( iNeighElem ) = elemPos
          end do ! iNeighElem
        end do ! iStencil
      end do ! iElem
    end do ! iLevel

    ! JQ: halo and ghost (with empty stencil) are added to elem list.
    !     shall we calculate possible neighbors for them?

    ! Find neighbors of neighbors for elements in the require list.
    ! Same procedure as above but identify neighbors of elements in require
    ! list along 1st stencil directions.
    ! What exactly is the require list for?
    ! - Used ONLY for boundary stencil with higher order neighbors i.e
    !   only when require nVals > 0
    call identify_additionalNeigh( tree       = tree,             &
      &                            proc       = proc,             &
      &                            levelDesc  = levelDesc,        &
      &                            pathFirst  = pathFirst,        &
      &                            pathLast   = pathLast,         &
      &                            stencil    = stencil )

    ! dump elemList into debug unit
    if (tem_logging_isActive(main_debug%logger, 5)) then
      do iLevel = minLevel, maxLevel
          call tem_elemList_dump( me      = levelDesc( iLevel )%elem, &
            &                     nUnit   = dbgUnit(5),               &
            &                     stencil = .true.,                   &
            &                     string  = 'after build level elm'   )
      end do ! iLevel
    end if

    write(dbgUnit(1),*) 'Leave routine: build_levelElements'
    call tem_horizontalSpacer( fUnit = dbgUnit(1), after = 1 )

  end subroutine build_levelElements
  ! ------------------------------------------------------------------------ !
  

  ! ------------------------------------------------------------------------ !
  !> Check, on which partition a given element is located add required elements
  !! to corresponding lists:
  !! if remote, add to halo
  !! if ghost, add to resp. ghost list
  !!
  recursive subroutine identify_elements( treeID, tree, pathFirst, pathLast,  &
    &                                     levelDesc, elemPos, proc,           &
    &                                     Stencil, nesting,                   &
    &                                     skip_add_additionalGhost )
    ! -------------------------------------------------------------------- !
    !> treeID to identify
    integer(kind=long_k), intent(in) :: treeID
    !> tree information
    type(treelmesh_type), intent(in) :: tree
    !> first treeID path in every process
    type(tem_path_type), intent(in) :: pathFirst(:)
    !> last treeID path in every process
    type(tem_path_type), intent(in) :: pathLast(:)
    !> the level descriptor to be filled
    type(tem_levelDesc_type), intent(inout) :: levelDesc(tree%global%minLevel:)
    !> Process description to use.
    type(tem_comm_env_type), intent(in) :: proc
    !> nTreeID element position in the levelDesc % elem list
    integer, intent(out) :: elemPos
    !> current stencil definition
    type(tem_stencilHeader_type), intent(in) :: Stencil
    !> nesting level
    integer, intent(in) :: nesting
    !> logical, optional, if true no ghosts are added
    logical, intent(in), optional :: skip_add_additionalGhost
    ! -------------------------------------------------------------------- !
    integer(kind=long_k) :: children(8) ! child elements
    integer :: nDepProcs ! processes that this neigh depend on
    integer :: depProc
    integer :: iChild, neighLevel
    type(tem_path_type) :: elemPath
    type(tem_stencilElement_type) :: emptyStencil(1)
    integer :: childPos   ! position of child element
    integer :: nNesting ! new nesting level
    ! Position of the current element in the hash
    integer :: hashpos, elemNesting
    logical :: cacheHit
    logical :: updated ! was the element updated during identify_local_element
    logical :: l_skip_add_additionalGhost
    ! -------------------------------------------------------------------- !
    
    if (present(skip_add_additionalGhost)) then 
      l_skip_add_additionalGhost = skip_add_additionalGhost
    else 
      l_skip_add_additionalGhost = .false.
    end if 
    
    ! Increase the nesting
    nNesting = nesting + 1
    cacheHit = .false.
    elemNesting = -999

    if (treeID > 0) then ! it is a element, otherwise it is a bcID
      neighLevel = tem_LevelOf(TreeID)
      elemPath  = tem_PathOf(TreeID)

      hashpos = int(mod(TreeID, nHashes))

      hashmatch: if (hash(hashpos) == TreeID) then

        if ( nesting == levelDesc(neighlevel)         &
          &               %elem                       &
          &               %haloNesting                &
          &               %val(hash_elemPos(hashPos)) ) then
          cacheHit = .true.
        end if
        if ( levelDesc(neighlevel)         &
          &    %elem                       &
          &    %needsUpdate                &
          &    %val(hash_elemPos(hashPos)) ) then
          cacheHit = .false.
          ! Set the needs update to false as it was now updated
          levelDesc(neighLevel)            &
            &  %elem                       &
            &  %needsUpdate                &
            &  %val(hash_elemPos(hashPos)) = .false.
        end if

      end if hashmatch

      cachemiss: if (.not. cacheHit) then
        ! Cache miss
        ! Probably did not hit this element yet, put it in the hash now,
        ! and identify it.
        ! (If this treeID is already in the hash, we definitely have
        !  already processed this treeID, and need nothing to do anymore).
        ! This leaves elements to be identified multiple times only in
        ! the case of hash collisions, which should be pretty few, for
        ! a sufficiently large hash.
        ! -----------------------------------------------
        call tem_find_depProc( depProc   = depProc,   &
          &                    nDepProcs = nDepProcs, &
          &                    tree      = tree,      &
          &                    elemPath  = elemPath,  &
          &                    PathFirst = PathFirst, &
          &                    PathLast  = PathLast   )

        ! How many processes possess a part of the requested treeID
        if (nDepProcs == 1) then
          updated = .false.
          ! Might be a local or halo of same level or a ghostFromCoarser
          ! or a ghostFromFiner. If halo and ghost it is distributed
          ! in single process.
          ! elemPos here is position of TreeID in levelDesc elem list
          call single_process_element( targetID       = TreeID,               &
            &                          levelDesc      = levelDesc,            &
            &                          tree           = tree,                 &
            &                          proc           = proc,                 &
            &                          iProc          = depProc,              &
            &                          minLevel       = tree%global%minlevel, &
            &                          elemPos        = elemPos,              &
            &                          updated        = updated,              &
            &                          stencil        = Stencil,              &
            &                          nesting        = nesting,               &
            &              skip_add_additionalGhost = l_skip_add_additionalGhost)

          if (elemPos > 0) then
          
            elemNesting = min( nesting,              &
              &                levelDesc(neighLevel) &
              &                  %elem               &
              &                  %haloNesting        &
              &                  %val(elemPos)       )

            if ( nesting < levelDesc(neighLevel)%elem%haloNesting  &
              &                                      %val(elemPos) ) then
              ! element needs updating as the current nesting is smaller than the
              ! one in the element property
              levelDesc(neighLevel)%elem%needsUpdate%val(elemPos) = .true.
              levelDesc(neighLevel)%elem%haloNesting%val(elemPos) &
                &  = elemNesting
              updated = .true.
            end if
            
          endif

          ! it is ghost from coarser element in current process proc%rank
          if (updated) then
            if ( levelDesc( neighLevel )%elem%eType%val( elemPos )   &
              &                               == eT_ghostFromCoarser ) then
              if (nesting < nestingLimit) then
                ! Create the direct neighbors of the ghostFromCoarser

                ! identify all the compute neighbors of the current element
                call identify_stencilNeigh( iElem          = elemPos,        &
                  &                         iLevel         = neighLevel,     &
                  &                         tree           = tree,           &
                  &                         iStencil       = 1,              &
                  &                         pathFirst      = pathFirst,      &
                  &                         pathLast       = pathLast,       &
                  &                         levelDesc      = levelDesc,      &
                  &                         proc           = proc,           &
                  &                         stencil        = Stencil,        &
                  &                         nesting        = elemNesting + 1 )

              end if ! nesting < 1?

              ! Create all elements required up to the actual existing fluid
              ! element these include the neighbors of the parents. In a level
              ! jump >1, these intermediate levels have to provide valid
              ! quantities over two of their computation updates to account for
              ! the recursive algorithm.
              call create_allParentNeighbors( &
                &    targetID  = levelDesc(neighLevel)%elem%tID%val(elemPos), &
                &    level     = neighLevel,     &
                &    tree      = tree,           &
                &    stencil   = Stencil, &
                &    levelDesc = levelDesc,      &
                &    pathFirst = pathFirst,      &
                &    pathLast  = pathLast,       &
                &    proc      = proc )

            end if ! ghostFromCoarser?
          end if ! updated?

        else if ( nDepProcs > 1 ) then
          ! more than one depending processes
          ! create ghost and find children (can only have finer children)
          call init( me        = emptyStencil(1),   &
            &        QQN       = Stencil%QQN,       &
            &        headerPos = 1 )
          ! Add to the level-wise ghost list
          ! append and store position of element for return
          call append( me              = levelDesc( neighLevel )%elem, &
            &          tID             = TreeID,                       &
            &          property        = 0_long_k,                     &
            &          eType           = eT_distributedGhostFromFiner, &
            &          stencilElements = emptyStencil,                 &
            &          pos             = elemPos                       )

          !... and store the position in the ghost list !
          ! Now find all children of this ghost
          ! childPos is return but not used anymore
          children = tem_directChildren( TreeID )
          do iChild = 1,8
            call identify_elements( TreeID    = children( iChild ), &
              &                     tree      = tree,               &
              &                     pathFirst = pathFirst,          &
              &                     pathLast  = pathLast,           &
              &                     levelDesc = levelDesc,          &
              &                     proc      = proc,               &
              &                     elemPos   = childPos,           &
              &                     stencil   = Stencil,            &
              &                     nesting   = nNesting,            &
              &  skip_add_additionalGhost = l_skip_add_additionalGhost)
          end do
        else ! nDepProcs < 1
          ! cell not existing. stop.
          write(logUnit(6),*) '---------- WARNING!!! -------------------------'
          write(logUnit(6),*) 'cell', TreeID ,' not existing on proc ',       &
            &        proc%rank, '. Number of dependencies', nDepProcs
          call tem_tIDinfo( me = TreeID, tree = tree, nUnit = logUnit(6) )
          write(logUnit(6),*) 'WARNING!!! This should never occur and '        &
            &                 // 'points to a bug or a buggy mesh.'
        end if

        if (elemPos > 0) then
          ! Set the encountered element position to the hash
          hash(hashpos) = TreeID
          hash_elemPos( hashpos ) = elemPos
        end if
        ! end cache miss
        ! -----------------------------------------------
      else cachemiss
        ! -----------------------------------------------
        ! Cache hit, i.e. element already in levelDesc%elem, just update nesting
        elemPos = hash_elemPos( hashpos )
        levelDesc( neighLevel )%elem%haloNesting%val( elemPos )            &
          &  = min( levelDesc( neighLevel )%elem%haloNesting%val(elemPos), &
          &         nesting                                                )
        ! -----------------------------------------------
      end if cachemiss
    else ! treeID <= 0, i.e. it is a bcID
      elemPos = int( TreeID )
    end if ! TreeID > 0

  end subroutine identify_elements
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Find the partitions holding data on a given path
  !!
  !! Using a binary search over the processes first and last elements.
  !!
  !! Is the element in question a local or remote element?
  !! To look up a certain element by its
  !! [[treelmesh_type:treeID]] in the distributed
  !! list of elements, it is sufficient to know the splitting positions of all
  !! chunks.
  !! That is, the first and last
  !! [[treelmesh_type:treeID]] of each partition.
  !! With a binary search over the splitting positions any requested element can
  !! then be identified to be either outside the computational domain at all, or
  !! inside of one or several known partitions.
  subroutine tem_find_depProc( depProc, nDepProcs, tree, elemPath, PathFirst, &
    &                          PathLast )
    ! -------------------------------------------------------------------- !
    !> List of partitions
    integer, intent(out) :: depProc
    !> Number of partitions
    integer, intent(out) :: nDepProcs
    !> tree information
    type(treelmesh_type), intent(in) :: tree
    !> Element to look up
    type(tem_path_type), intent(in) :: elemPath
    !> Left partition bounds
    type(tem_path_type), intent(in) :: PathFirst(:)
    !> Right partition bounds
    type(tem_path_type), intent(in) :: PathLast(:)
    ! -------------------------------------------------------------------- !
    integer :: p_lb, p_ub ! process lower and upper bound
    integer :: relFirst, relLast
    integer :: myRank
    ! -------------------------------------------------------------------- !
    myRank = tree%global%myPart

    nDepProcs = 0
    ! First check if this neighbor is in my local tree range
    p_lb = 1
    ! rank starts from indice 0 where as fortran array pathFirst and
    ! pathLast starts from 1
    relFirst = tem_PathComparison(elemPath, PathFirst(myRank+1))
    relLast = tem_PathComparison(elemPath, PathLast(myRank+1))

    if (relFirst < 0) then
      ! The searched element is definitely left of myself
      ! Do only search up to my own range
      p_ub = myRank
    else
      p_ub = int( min( int(tree%global%nParts, kind=long_k), &
        &             tree%global%nElems ) )
      if (relLast > 0) then
        ! The searched element is definitely right of myself
        ! Do only search beyond my own range
        p_lb = myRank+2
      else

        ! The element might be (partly) on my own partition
        if (relLast < 0) then
          ! The element is at least left of my right border
          ! Set the upper bound to include myself
          p_ub = myRank + 1
        end if
        if (relFirst > 0) then
          ! The element is at least right of my left border
          ! Set the lower bound to include myself
          p_lb = myRank + 1
        end if
      end if ! relLast > 0
    end if ! relFirst < 0

    if ((p_lb == p_ub) .and. (p_lb == myRank+1)) then
      ! The element is local, do not go on to other processes
      nDepProcs = 1
      depProc = myRank+1
    else
      ! Possibly NON-LOCAL Element...
      ! Every process COULD hold a part of the current neighbor element
      call tem_find_depProc_globSearch( depProc   = depProc,   &
        &                               nDepProcs = nDepProcs, &
        &                               elemPath  = elemPath,  &
        &                               p_lb      = p_lb,      &
        &                               p_ub      = p_ub,      &
        &                               PathFirst = PathFirst, &
        &                               PathLast  = PathLast   )
    end if ! local or non-local

  end subroutine tem_find_depProc
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Find the remote partitions holding data on a given path
  !!
  !! Using a binary search over the processes first and last elements.
  !!
  !! Is the element in question a local or remote element?
  !! To look up a certain element by its
  !! [[treelmesh_type:treeID]] in the distributed
  !! list of elements, it is sufficient to know the splitting positions of all
  !! chunks.
  !! That is, the first and last
  !! [[treelmesh_type:treeID]] of each partition.
  !! With a binary search over the splitting positions any requested element can
  !! then be identified to be either outside the computational domain at all, or
  !! inside of one or several known partitions.
  subroutine tem_find_depProc_globSearch( depProc, nDepProcs, elemPath, p_lb, &
    &                                     p_ub, PathFirst, PathLast)
    ! -------------------------------------------------------------------- !
    !> List of partitions
    integer, intent(out) :: depProc
    !> Number of partitions
    integer, intent(out) :: nDepProcs
    !> Element to look up
    type(tem_path_type), intent(in) :: elemPath
    !> Left partition bounds
    type(tem_path_type), intent(in) :: PathFirst(:)
    !> Right partition bounds
    type(tem_path_type), intent(in) :: PathLast(:)
    !> Left interval bound to search in
    integer, intent(in) :: p_lb
    !> Right interval bound to search in
    integer, intent(in) :: p_ub
    ! -------------------------------------------------------------------- !
    integer :: lb, ub
    integer :: foundProc, lastProc
    integer :: curProc
    integer :: pComp
    ! -------------------------------------------------------------------- !

    lb = p_lb
    ub = p_ub

    ! Find the partition, to which elemPath is definitely right of
    ! the first element.
    do
      curProc = (lb + ub) / 2
      if( curProc <= 0 ) then
        ! Element not found on any proc
        foundProc = 0
        exit
      end if
      pComp = tem_PathComparison(elemPath, PathFirst(curProc))
      if (pComp > 0) then
        ! The first element of curProc is smaller than elemPath, therefore
        ! we do not have to search left of curProc. Can we even go to curProc+1?
        if (tem_PathComparison(elemPath, PathLast(curProc)) > 0) then
          ! The search element is not in the partition curProc at all, continue
          ! search with new lower bound of curProc+1.
          lb = min(curProc + 1, p_ub)
        else
          ! The left element is smaller than the search element and the right
          ! one is larger or equal to it, thus we found the first partition,
          ! which contains information on the element.
          foundProc = curProc
          exit
        end if
      else
        ! The first element of curProc is larger or equal to elemPath, thus
        ! we can restrict our search interval up to this partition.
        ub = curProc
        ! It is OK to not decrease the interval further here, as we always
        ! round down to find the middle, thus, if by this step we set
        ! ub = lb+1, the next iteration will test a new element, as we get
        ! curProc = (lb + lb + 1) / 2 = lb.
      end if
      if (ub <= lb) then
        ! Interval of size 0, found the first process to hold information on
        ! the element in neighpath.
        ! The lower bound was tested strictly, thus it is safe to set the first
        ! process to the lower bound here and exit.
        foundProc = lb
        pComp = tem_PathComparison(elemPath, PathFirst(foundProc))
        ! element is below lower bound of process 1.
        ! Therefor no process is holding the element
        if ( lb == 1 .and. pComp < 0 ) foundProc = 0
        ! element is below lower bound of process and the element is bigger
        ! than upper bound of previous process. -> element doesn't exist
        if ( ub == lb .and. pComp < 0 ) foundProc = 0
        exit
      end if
    end do

    ! Now linearly search for the last process spanned by the searched element.
    ! Usually this is just the single partition identified above or more
    lastProc = foundProc
    if (lastProc > 0) then
      if (lastProc < p_ub) then
        ! Found process that is not the last one, check whether we span
        ! more than this one (it might just be the first partition with that
        ! element).
        nDepProcs = 1
        depProc = lastProc
        do
          if (lastProc == p_ub) exit
          if (tem_PathComparison(elemPath, PathLast(lastProc)) < 0) exit
          if (tem_PathComparison(elemPath, PathFirst(lastProc+1)) < 0) exit
          lastProc = lastProc + 1
          nDepProcs = nDepProcs + 1
        end do
      else
        ! Last process, there can not be any more processes that would contain
        ! this element. Check whether this element actually is found in this
        ! last process.
        if (tem_PathComparison(elemPath, PathLast(lastProc)) <= 0) then
          ! Element is found on this process, return it accordingly.
          nDepProcs = 1
          depProc = lastProc
        else
          ! Element not found on last process, nowhere in all processes to be
          ! found, return accordingly 0 processes.
          nDepProcs = 0
          depProc = 0
        end if
      end if
    else
      ! Element is left of first process in given list, so no process contains
      ! it ...
      nDepProcs = 0
      depProc = 0
    end if

  end subroutine tem_find_depProc_globSearch
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> create all the neighbors of an element's parent
  !!
  !! Create all elements required up to the actual existing fluid
  !! element these include the neighbors of the parents. In a level
  !! jump >1, these intermediate levels have to provide valid
  !! quantities over two of their computation updates to account for
  !! the recursive algorithm.
  !!
  !! Here the fromCoarser interpolation should be handed in.
  !!
  recursive subroutine create_allParentNeighbors(                          &
    &                           targetID, level, stencil, tree, levelDesc, &
    &                           pathFirst, pathLast, proc                  )
    ! -------------------------------------------------------------------- !
    !> requested element position (child element) in LevelDesc elem list
    integer(kind=long_k), intent(in) :: targetID
    !> requested element level
    integer, intent(in) :: level
    !> first treeID path in every process
    type(tem_path_type), intent(in) :: pathFirst(:)
    !> last treeID path in every process
    type(tem_path_type), intent(in) :: pathLast(:)
    !> tree information
    type(treelmesh_type), intent(in) :: tree
    !> the level descriptor to be filled
    type(tem_levelDesc_type), intent(inout) :: levelDesc(tree%global%minLevel:)
    !> process
    type(tem_comm_env_type), intent(in) :: proc
    !> current stencil definition
    type(tem_stencilHeader_type), intent(in) :: stencil
    ! -------------------------------------------------------------------- !
    integer(kind=long_k) :: parentID, neighID  ! current tree ID
    integer :: coarserLevel, cPos, parentNesting, addedPos, iStencilElem
    integer :: neighIDpos
    ! -------------------------------------------------------------------- !

    ! exit if we have reached the minimal level
    if ( level == tree%global%minlevel ) return

    ! Get the parent of the current treeID
    parentID = tem_parentOf( targetID )

    ! ... and identify the parent
    parentNesting = -1
    call identify_elements( TreeID     = parentID,  &
      &                     tree       = tree,      &
      &                     pathFirst  = pathFirst, &
      &                     pathLast   = pathLast,  &
      &                     levelDesc  = levelDesc, &
      &                     elemPos    = cPos,      &
      &                     proc       = proc,      &
      &                     stencil    = stencil,   &
      &                     nesting    = -1         )

    if ( cPos <= 0 ) then
      write(dbgUnit(3),*) ' Element not found: ', parentID
      write(dbgUnit(3),*) '              cPos: ', cPos
    end if

    ! identify the stencil neighbors of the parent.
    ! Here we should identify the fromCoarser interpolation stencil neighbors
    ! instead of the compute stencil neighbors
    coarserLevel = level - 1
    call identify_stencilNeigh( iElem     = cPos,         &
      &                         iLevel    = coarserLevel, &
      &                         tree      = tree,         &
      &                         iStencil  = 1,            &
      &                         pathFirst = pathFirst,    &
      &                         pathLast  = pathLast,     &
      &                         levelDesc = levelDesc,    &
      &                         proc      = proc,         &
      &                         stencil   = stencil,      &
      &                         nesting   = parentNesting )
      
    ! adding neighs of neighs for interpolation at ML
    ! this is for the stencil of the fine ghost, composed by coarse fluid 
    do iStencilElem = 1, stencil%QQN
    
      neighIDpos = levelDesc(coarserLevel)%elem%stencil%val(cPos)    &
        &                           %val(1)%tIDpos(iStencilElem)
      if( neighIDpos > 0 ) then
        neighID = &
          & levelDesc( coarserLevel )%elem%neighID%val(cPos)%val(neighIDpos)
        ! This call might add new halo elements
        if ( neighID > 0_long_k ) then
          call append( me       = levelDesc( coarserLevel )%require, &
            &          val      = neighID,                           &
            &          pos      = addedPos )
        end if
      end if ! neighIDpos > 0
    enddo

  end subroutine create_allParentNeighbors
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Invoke the identify_elements for each neighbor of the stencil
  !! and store the positions of the encountered elements
  !!
  recursive subroutine identify_stencilNeigh( iElem, iLevel, iStencil, tree,  &
    &                                         pathFirst, pathLast, levelDesc, &
    &                                         proc, stencil, nesting          )
    ! -------------------------------------------------------------------- !
    !> element position in levelDesc to identify
    integer, intent(in)                      :: iElem
    !> element level
    integer, intent(in)                      :: iLevel
    !> stencil within the element to act on
    integer, intent(in)                      :: iStencil
    !> tree information
    type(treelmesh_type), intent(in)         :: tree
    !> first treeID path in every process
    type(tem_path_type), intent(in)          :: pathFirst(:)
    !> last treeID path in every process
    type(tem_path_type), intent(in)          :: pathLast(:)
    !> the level descriptor to be filled
    type(tem_levelDesc_type), intent(inout)  :: levelDesc(tree%global%minLevel:)
    !> process
    type(tem_comm_env_type), intent(in)      :: proc
    !> current stencil definition
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> nesting level
    integer, intent(in)                      :: nesting
    ! -------------------------------------------------------------------- !
    integer :: iStencilElem, elemPos
    integer :: neighIDpos
    integer(kind=long_k) :: neighID
    ! -------------------------------------------------------------------- !
    ! identify all the compute neighbors of the current element
    do iStencilElem = 1, stencil%QQN
      neighIDpos = levelDesc(iLevel)%elem%stencil%val(iElem)          &
        &                           %val(iStencil)%tIDpos(iStencilElem)
      if( neighIDpos > 0 ) then
        neighID = &
          & levelDesc( iLevel )%elem%neighID%val(iElem)%val(neighIDpos)
        ! This call might add new halo elements
        if ( neighID > 0_long_k ) then
          call identify_elements( TreeID     = neighID,   &
            &                     tree       = tree,      &
            &                     pathFirst  = pathFirst, &
            &                     pathLast   = pathLast,  &
            &                     levelDesc  = levelDesc, &
            &                     elemPos    = elemPos,   &
            &                     proc       = proc,      &
            &                     stencil    = stencil,   &
            &                     nesting    = nesting    )
        else ! neighID < 0
          elemPos = 0
        end if
      else ! neighIDpos < 0
        elemPos = 0
      end if ! neighIDpos > 0

      ! And add the encountered neighbor elements to the current element's
      ! stencil neighbors
      levelDesc(iLevel)%elem%stencil%val(iElem)              &
        &              %val(iStencil)%totalPos(iStencilElem) = elemPos
    end do

  end subroutine identify_stencilNeigh
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Determine the location (which process) of a requested element,
  !! which was identified to be located on one single process
  !! (can be local or remote)
  !! If it is located on a remote process: add to halo list
  !!                     local process: identify if ghost or fluid
  !!
  subroutine single_process_element( targetID, levelDesc, tree, proc, iProc, &
    &                                minLevel, elemPos, stencil, nesting,    &
    &                                updated, skip_add_additionalGhost )
    ! -------------------------------------------------------------------- !
    !> neighboring treeID
    integer(kind=long_k), intent(in)         :: targetID
    !> minimum level fluid element in the tree
    integer, intent(in)                      :: minLevel
    !> the level descriptor to be filled
    type(tem_levelDesc_type), intent(inout)  :: levelDesc(minLevel:)
    !> Process description to use.
    type(tem_comm_env_type), intent(in)      :: proc
    !> Process on which targetID is located
    integer, intent(in)                      :: iProc
    !> tree information
    type(treelmesh_type), intent(in)         :: tree
    !> current stencil definition
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> targetID element position in the levelDesc % elem list
    integer, intent(out)                     :: elemPos
    !> nesting level
    integer, intent(in)                      :: nesting
    !> was the element updated in this call?
    logical, intent(out)                     :: updated
    !> logical, optional, if true no ghosts are added
    logical, intent(in), optional :: skip_add_additionalGhost
    ! -------------------------------------------------------------------- !
    type(tem_stencilElement_type) :: emptyStencil(1)
    integer :: targetLevel
    logical :: wasAdded
    logical :: l_skip_add_additionalGhost
    ! -------------------------------------------------------------------- !
    
    if (present(skip_add_additionalGhost)) then 
      l_skip_add_additionalGhost = skip_add_additionalGhost
    else 
      l_skip_add_additionalGhost = .false.
    end if 

    targetLevel = tem_LevelOf(targetID) ! Has to be same as tLevel!?
    if ( (targetLevel < minLevel)                 &
      &  .or. (targetLevel > uBound(levelDesc,1)) ) then
      write(logUnit(1),*) ' ERROR: level which is not included in the fluid'// &
        &                 ' tree was demanded in singleProcNeigh'
      write(logUnit(1),"(2(A,I0))") 'treeID: ', targetID, ', level: ', &
        &                           targetLevel
      call tem_abort()
    end if

    ! Set the element updated flag as a default to false
    updated = .false.

    ! If it is a remote cell on only one process -> regular halo
    if (iProc /= proc%rank + 1) then
      ! REMOTE
      call init( me = emptyStencil(1), QQN=stencil%QQN )

      ! append this targetID as halo element to levelDesc elem list
      call append( me              = levelDesc(targetLevel)%elem, &
        &          tID             = targetID,                    &
        &          eType           = eT_halo,                     &
        &          property        = 0_long_k,                    &
        &          sourceProc      = iProc,                       &
        &          haloNesting     = nesting,                     &
        &          stencilElements = emptyStencil,                &
        &          pos             = elemPos,                     &
        &          wasAdded        = wasAdded                     )

      if (.not. wasAdded) then
        ! If this element was already there, make sure we use the minimal
        ! nesting level requested for this element.
        updated = ( nesting < levelDesc(targetLevel) &
          &                     %elem                &
          &                     %haloNesting         &
          &                     %val(elemPos)        )
        ! If the nesting has been updated (decreased, we need to revisit this
        ! element in search for its neighbors).
        levelDesc(targetLevel)%elem%needsUpdate%val(elemPos) = updated

        levelDesc(targetLevel)%elem%haloNesting%val(elemPos)              &
          &  = min( levelDesc(targetLevel)%elem%haloNesting%val(elemPos), &
          &         nesting                                               )
      else
        ! New halo element added
        updated = .true.
        write(dbgUnit(7),"(A,I0)") 'Added as a Halo: ', targetID, &
          &                        'to level: ', targetLevel
      end if ! wasAdded

    else ! iProc == proc%rank + 1
      ! LOCAL

      ! Either a local ghost or fluid cell
      call identify_local_element(                                 &
        &    targetID                 = targetID,                  &
        &    levelDesc                = levelDesc,                 &
        &    tree                     = tree,                      &
        &    elemPos                  = elemPos,                   &
        &    minLevel                 = minLevel,                  &
        &    nesting                  = nesting,                   &
        &    updated                  = updated,                   &
        &    stencil                  = stencil,                   &
        &    skip_add_additionalGhost = l_skip_add_additionalGhost )
    end if ! iProc /= proc%rank + 1

  end subroutine single_process_element
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Determine if the target element (local) targetID is fluid or ghost in the
  !! local process
  !! If fluid: do nothing, as it will be added later on anyway (or already is)
  !!    ghostFromFiner (coarser than requested):
  !!        add all virtual children, i.e. all levels between requested treeID
  !!        and found one.
  !!    ghostFromCoarser (finer than requested):
  !!    not existing( localPos=0): add to halo
  !!
  subroutine identify_local_element( targetID, levelDesc, tree, minLevel, &
    &                                elemPos, nesting, updated, stencil,  &
    &                                skip_add_additionalGhost             )
    ! -------------------------------------------------------------------- !
    !> neighboring treeID
    integer(kind=long_k), intent(in) :: targetID
    !> minimum level fluid element in the tree
    integer, intent(in) :: minLevel
    !> the level descriptor to be filled
    type(tem_levelDesc_type), intent(inout) :: levelDesc(minLevel:)
    !> tree information
    type(treelmesh_type), intent(in) :: tree
    !> nesting level
    integer, intent(in) :: nesting
    !> current stencil definition
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> targetID element position in the levelDesc % elem list
    integer, intent(out) :: elemPos
    !> was the element updated in this call?
    logical, intent(out) :: updated
    !> logical, optional, if true no ghosts are added
    logical, intent(in), optional :: skip_add_additionalGhost
    ! -------------------------------------------------------------------- !
    integer :: localPos, targetLevel, dPos, fluidLevel
    integer(kind=long_k) :: fluidID
    type(tem_path_type) :: targetPath
    logical :: l_skip_add_additionalGhost
    ! -------------------------------------------------------------------- !

    if (present(skip_add_additionalGhost)) then
      l_skip_add_additionalGhost = skip_add_additionalGhost
    else
      l_skip_add_additionalGhost = .false.
    end if

    ! Set the element updated flag as a default to false
    updated = .false.

    ! Position of neighbor treeID
    targetLevel = tem_LevelOf( targetID )
    targetPath  = tem_PathOf( targetID )

    ! Return position of targetID in the treeIDlist
    localPos = tem_PosOfPath( targetPath, tree%pathList )

    ! By localPos we can determine how the element (might) exist locally:
    ! - fluid
    ! - ghostFromCoarser
    ! - ghostFromFiner
    if (localPos > 0) then
      ! Path exist. It may be GhostFromCoarser or FLUID
      fluidID = tree%treeID( localPos )
      fluidLevel = tem_LevelOf( fluidID )

      if (fluidLevel == targetLevel) then
        ! It is a FLUID. Already exists in element list
        updated = .false.
        elemPos = PositionOfVal( me  = levelDesc( targetLevel )%elem%tID, &
          &                      val = targetID                           )

      else if (fluidLevel < targetLevel) then
        if (.not. l_skip_add_additionalGhost) then
          ! Target element is a GhostFromCoarser.
          ! Target element is a descendant of Fluid element.
          ! ---------------
          ! |             |
          ! |             |
          ! |             |
          ! |      F      |
          ! |-----        |
          ! | T  |        |
          ! |    |        |
          ! ---------------
          ! Add all the descendants of Fluid down to target( including
          ! intermediate levels).
          call add_all_virtual_children(                           &
            &    sourceID       = fluidID,                         &
            &    foundPos       = localPos,                        &
            &    elemPath       = targetPath,                      &
            &    sourceProperty = tree%ElemPropertyBits(localPos), &
            &    targetLevel    = targetLevel,                     &
            &    levelDesc      = levelDesc,                       &
            &    minlevel       = minLevel,                        &
            &    nesting        = nesting,                         &
            &    updated        = updated,                         &
            &    tree           = tree,                            &
            &    Stencil        = stencil                          )
        end if
      end if ! on same level?

    else if (localPos < 0) then
      if (.not. l_skip_add_additionalGhost) then
        ! ghostFromFiner
        ! Find all existing fluid cells within requested targetID position
        ! Add all the parents between requested targetID and available child ID in
        ! treeID list
        call add_ghostFromFiner( elemID     = targetID,  &
          &                      levelDesc  = levelDesc, &
          &                      minLevel   = minlevel,  &
          &                      tree       = tree,      &
          &                      foundPos   = dPos,      &
          &                      updated    = updated,   &
          &                      stencil    = stencil    )
      end if

    else ! localPos == 0

      write(dbgUnit(6),*) 'Warning: element not existing ', targetID,    &
        &                 'adding to nonexisting ...'
      call tem_tIDinfo( me = targetID, tree = tree, nUnit = dbgUnit(6) )
      ! This case occurs, when a remote halo was added, which was not existing
      ! before.
      ! Halos are added in the transfer_treeIDs routine
      call append( me         = levelDesc( targetLevel )%elem, &
        &          tID        = targetID,                      &
        &          eType      = eT_nonExisting,                &
        &          sourceProc = tree%global%myPart+1,          &
        &          pos        = dPos                           )

    endif ! localPos > 0? coarser ghost or fluid?

    ! position of added targetID in the levelDesc elem list
    elemPos = PositionOfVal( me  = levelDesc( targetLevel )%elem%tID, &
      &                      val = targetID                           )

  end subroutine identify_local_element
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Find all the virtual children of the sourceID down to the targetLevel
  !! and add to the level-wise ghostFromCoarser list in the level descriptor
  !!
  recursive subroutine add_all_virtual_children(                       &
    &                    sourceID, sourceProperty, foundPos, elemPath, &
    &                    targetLevel, levelDesc, minLevel, tree,       &
    &                    stencil, nesting, updated                     )
    ! -------------------------------------------------------------------- !
    !> source treeID (existing founded ID in tree%treeID list or children ID
    !! from recursion)
    integer(kind=long_k), intent(in)            :: sourceID
    !> property of source element
    integer(kind=long_k), intent(in)            :: sourceProperty
    !> element path
    type(tem_path_type), intent(in)             :: elemPath
    !> nesting level
    integer, intent(in)                         :: nesting
    !> position of this sourceID in elem%tID list
    integer, intent(in)                         :: foundPos
    !> level upto which virtual children must be created
    integer,intent(in)                          :: targetLevel
    !> minimum level in the tree
    integer, intent(in)                         :: minLevel
    !> the level descriptor to be filled
    type(tem_levelDesc_type), intent(inout)     :: levelDesc(minLevel:)
    !> tree information
    type(treelmesh_type), intent(in)            :: tree
    !> current stencil definition
    type(tem_stencilHeader_type), intent(in)    :: stencil
    !> was the element updated in this call?
    logical, intent(out)                        :: updated
    ! -------------------------------------------------------------------- !
    integer :: targetPos, iChild, iDir, nVals, sourceLevel
    ! position of the existing (source) tID in the elem list
    integer :: sourcePos
    logical :: wasAdded
    logical :: childUpdated
    integer(kind=long_k) :: cTreeID      ! current treeID to identify
    integer(kind=long_k) :: ichildID     ! child treeID
    integer(kind=long_k) :: curNeighborID
    integer(kind=long_k), allocatable :: tNeighID(:)
    type(tem_stencilElement_type) :: tStencil(1)
    integer :: iChildCoord(4), curLevel, offset(4), xc(4), childCoord(4)
    integer :: addedPos
    ! -------------------------------------------------------------------- !

    !Position of the coarser source element
    sourceLevel = tem_LevelOf(sourceID)
    sourcePos = PositionOfVal( me  = levelDesc( sourceLevel )%elem%tID, &
      &                        val = sourceID                           )
    allocate( tNeighID(stencil%QQN) )
    offset(4) = 0

    ! By default, set that no element was updated
    updated = .false.
    childUpdated = .false.

    ! create virual children until target level is reached
    if ( sourceLevel < targetLevel ) then
      curLevel = sourceLevel + 1
      call init( me = tStencil(1), QQN = stencil%QQN, headerPos = 1 )
      ! Add to the level-wise ghost list
      cTreeID = elemPath%node( targetLevel - sourceLevel )
      call append( me              = levelDesc( curLevel )%elem, &
        &          tID             = cTreeID,                    &
        &          property        = sourceProperty,             &
        &          eType           = eT_ghostFromCoarser,        &
        &          sourceProc      = tree%global%myPart+1,       &
        &          stencilElements = tStencil,                   &
        &          pos             = targetPos,                  &
        &          haloNesting     = nesting,                    &
        &          wasAdded        = wasAdded                    )

      if (wasAdded) then

        write(dbgUnit(7),"(2(A,I0))") 'Added as a GhostFromCoarser: ', &
          &                           cTreeID, ', to level: ', curLevel

        updated = .true.
        iChild = tem_childNumber(cTreeID)
        iChildCoord(:) = tem_coordOfId( int(iChild, long_k) )
        ! inherit the boundary infos from the parent
        tNeighID = 0_long_k
        ! Get neighIds of (source) coarse element
        do iDir = 1, 3
          call tem_find_BCs_fromCoarser( dir            = iDir,        &
            &                            childCoord     = iChildCoord, &
            &                            sourceLevel    = sourceLevel, &
            &                            sourcePos      = sourcePos,   &
            &                            neighID        = tNeighID,    &
            &                            minLevel       = minLevel,    &
            &                            levelDesc      = levelDesc,   &
            &                            computeStencil = Stencil      )
        end do

        ! loop over all directions to determine neighIDs for ghost child
        ! (targetPos)
        do iDir = 1, stencil%QQN
          ! compute virtual neighbor child from coarser neighbor in iDir
          ! direction
          if( tNeighID(iDir) > 0_long_k ) then
            ! find the corresponding children of the neighbor defined for my
            ! parent. Find the child in the corresponding direction
            ! get the corresponding neighbor in the offset direction, which
            ! is a child of the coarser neighbor element. Use +2 offset to get
            ! only positive numbers
            ! (child coordinates always range between {0..1})
            offset(1:3) = stencil%cxDir(1:3, iDir )
            xc(1:3) = mod( iChildCoord(1:3) + 2 + offset(1:3), 2 )
            xc(4) = 1
            ! get the child number of the corresponding neighbor within the
            ! coarser element. (is {1...8})
            ichildID  = tem_IdOfCoord( xc )
            ! calculate the child treeID from the coarser neighbor and the
            ! childID
            curNeighborID = tNeighID(iDir)*8_long_k + ichildID
          else if (tNeighID(iDir) == 0_long_k) then
            ! virtual neighbor child exist in current parent.
            ! If the neighbor is still 0, it must be a direct sibling,
            ! so we can directly compute its tID.
            ! find the neighbor according to stencil offset from me
            childCoord = tem_coordOfId( cTreeID )
            offset(1:3) = stencil%cxDir(1:3, iDir )
            curNeighborID = tem_IdOfCoord( childCoord + offset )
          else ! tNeighID(iDir) < 0_long_k
            ! inherit the boundary ID
            curNeighborID = tNeighID(iDir)
          end if

          ! append the neighbor ID ...
          call append( me  = levelDesc( curLevel )%elem%neighID     &
            &                                    %val( targetPos ), &
            &          val = curNeighborID,                         &
            &          pos = addedPos                               )
          ! ... and store this position in the stencil
          levelDesc( curLevel )%elem%stencil%val( targetPos )       &
            &                               %val(1)%tIDpos( iDir )  &
            &  = addedPos

          ! to append the neighbors of ghosts to required list
          !call append( me       = levelDesc( curLevel )%require, &
          !  &          val      = curNeighborID,                 &
          !  &          pos      = addedPos,                      &
          !  &          wasAdded = wasAdded )

        end do !iDir QQN

        ! Set prp_hasBnd if any of neighbors is boundary
        nVals = levelDesc( curLevel )%elem%neighID%val( targetPos )%nVals
        if ( minval( levelDesc(curLevel)%elem%neighID%val(targetPos)      &
          &          %val(1:nVals)                                  ) < 0 ) then
          ! Found boundary
          levelDesc( curLevel )%elem%property%val( targetPos )              &
            &  = ibset( levelDesc( curLevel )%elem%property%val(targetPos), &
            &           prp_hasBnd                                          )
        else
          ! Unset the property bit
          levelDesc( curLevel )%elem%property%val( targetPos )              &
            &  = ibclr( levelDesc( curLevel )%elem%property%val(targetPos), &
            &           prp_hasBnd                                          )
        end if

      else
        ! Existing element encountered.
        updated = .false.
      end if
      ! Overwrite the eventually existing nesting with the smallest value.
      ! The smallest nesting determines if further neighbors have to be
      ! retrieved in communicate_elements
      if ( nesting < levelDesc(curLevel)%elem%haloNesting    &
        &                                    %val(targetPos) ) then
        ! needs update
        updated = .true.
        levelDesc( curLevel )%elem%needsUpdate%val( targetPos ) = .true.
        levelDesc( curLevel )%elem%haloNesting%val( targetPos )            &
          &   = min( nesting,                                              &
          &          levelDesc( curLevel )%elem%haloNesting%val(targetPos) )
      end if
      ! In any case we have to recurse down to the target level
      ! lower levels might not yet exist.
      call add_all_virtual_children( sourceID       = cTreeID,        &
        &                            foundPos       = foundPos,       &
        &                            elemPath       = elemPath,       &
        &                            nesting        = nesting,        &
        &                            targetLevel    = targetLevel,    &
        &                            levelDesc      = levelDesc,      &
        &                            minLevel       = minLevel,       &
        &                            updated        = childUpdated,   &
        &                            tree           = tree,           &
        &                            sourceProperty = sourceProperty, &
        &                            Stencil        = Stencil         )
    end if ! sourceLevel < targetLevel
    updated = ( updated .or. childUpdated )

  end subroutine add_all_virtual_children
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Inherit the neighborhood from the sourceELem to the targetElem
  !!
  subroutine tem_find_BCs_fromCoarser( dir, childCoord, sourceLevel,       &
    &                                  sourcePos, neighID, computeStencil, &
    &                                  levelDesc, minLevel                 )
    ! -------------------------------------------------------------------- !
    !> coarse element level
    integer, intent(in) :: sourceLevel
    !> position of coarser element in original treeID list
    integer, intent(in) :: sourcePos
    !>
    integer, intent(in) :: dir
    !> minimum level in the tree
    integer, intent(in) :: minLevel
    !> neighbor treeIDs of child element
    integer(kind=long_k), intent(inout) :: neighID(:)
    !> coordinate of virtual(ghost) child
    integer, intent(in) :: childCoord(4)
    !> current stencil definition
    type(tem_stencilHeader_type), intent(in) :: computeStencil
    !> the level descriptor to be filled
    type(tem_levelDesc_type), intent(in) :: levelDesc(minLevel:)
    ! -------------------------------------------------------------------- !
    ! Tangential direction iterators
    integer :: iDirX, iDirY
    integer :: dirX ! first tangential direction
    integer :: dirY ! second tangential direction
    integer :: curDir ! child face normal direction, always {-1,1}
    integer :: myLink(4) ! child link direction
    integer :: parentLink(3) ! each entry must be either {-1,0,1}
    integer :: iChildStencil ! child stencil element iterator
    integer :: iStencilElem  ! stencil element iterator
    integer :: iStencil
    integer :: posInNeighID
    ! -------------------------------------------------------------------- !

    ! curDir is -1 if childCoord(dir) is 0
    ! curDir is 1 if childCoord(dir) is 1
    ! it is needed to find in which direction of child's neighIDs are to be
    ! found
    curDir = childToStencil( childCoord( dir ))
    myLink( dir ) = curDir
    myLink(4) = 0
    parentLink(dir) = curDir

    ! identify the tangential directions
    dirX = mod( dir, 3) + 1
    dirY = mod( dir+1, 3) + 1
    iStencil = 1

    ! Iterate over all nine child links of the surface with the normal pointing
    ! into dir
    do iDirY = -1, 1
      myLink( dirY ) = iDirY
      ! select the valid directions to look for the neighbors in the second
      ! tangential direction
      if( childCoord( dirY ) == 0) then
        parentLink( dirY ) = min( iDirY, 0 )
      else
        parentLink( dirY ) = max( iDirY, 0 )
      end if
      do iDirX = -1, 1
        myLink( dirX ) = iDirX
        if( childCoord( dirX ) == 0) then
          parentLink( dirX ) = min( iDirX, 0 )
        else
          parentLink( dirX ) = max( iDirX, 0 )
        end if
        ! matching direction of child to compute stencil
        iChildStencil = 0
        do iStencilElem = 1, computeStencil%QQN
          if (      computeStencil%cxDir(1, iStencilElem) == myLink(1)  &
            & .and. computeStencil%cxDir(2, iStencilElem) == myLink(2)  &
            & .and. computeStencil%cxDir(3, iStencilElem) == myLink(3)  &
            &                                                           ) then
            ! Found matching stencil entry for current child direction
            iChildStencil = iStencilElem
          end if
        end do

        if ( iChildStencil > 0 ) then
          ! matching direction of parent to compute stencil
          do iStencilElem = 1, computeStencil%QQN
            if (      computeStencil%cxDir(1, iStencilElem) == parentLink(1) &
              & .and. computeStencil%cxDir(2, iStencilElem) == parentLink(2) &
              & .and. computeStencil%cxDir(3, iStencilElem) == parentLink(3) &
              &                                                              ) then
              ! Found matching stencil entry for current child direction
              ! Set the parent's neighbor here. Later we replace the
              ! parent's neighbor by the current level's neighbor
              posInNeighID = levelDesc(sourceLevel) &
                &              %elem                &
                &              %stencil             &
                &              %val(sourcePos)      &
                &              %val(iStencil)       &
                &              %tIDpos(iStencilElem)
              neighID( iChildStencil )        &
                &  = levelDesc( sourceLevel ) &
                &      %elem                  &
                &      %neighID               &
                &      %val(sourcePos)        &
                &      %val(posInNeighID)
            end if
          end do ! iStencilElem
        end if ! iChildStencil > 0
      end do ! iDirX
    end do ! iDirY

  end subroutine tem_find_BCs_fromCoarser
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Add parentID as GhostFromFiner.
  !! Then set its BC from its children.
  !! If any children do NOT exist, recursively call this routine to add them as
  !! GhostFromFiner.
  !!
  recursive subroutine add_ghostFromFiner( elemID, levelDesc, minLevel,     &
    &                                      tree, updated, foundPos, stencil )
    ! -------------------------------------------------------------------- !
    !> requested treeID
    integer(kind=long_k), intent(in) :: elemID
    !> minimum level fluid element in the tree
    integer, intent(in) :: minLevel
    !> the level descriptor to be filled
    type(tem_levelDesc_type), intent(inout) :: levelDesc(minLevel:)
    !> tree information
    type(treelmesh_type), intent(in) :: tree
    ! position of elemID in elem%tID list
    integer, intent(out) :: foundPos
    !> was the current element updated in this call?
    logical, intent(out) :: updated
    !> current stencil definition
    type( tem_stencilHeader_type ), intent(in) :: stencil
    ! -------------------------------------------------------------------- !
    integer :: iChild, level
    integer(kind=long_k) :: children(8), property
    logical :: wasAdded, childUpdated
    integer :: childPos(8)
    type(tem_path_type) :: childPath
    ! -------------------------------------------------------------------- !
    ! Set as not updated by default
    updated = .false.

    ! Create the ghostFromFiner
    level = tem_LevelOf( elemID )
    call append( me         = levelDesc( level )%elem, &
      &          tID        = elemID,                  &
      &          eType      = eT_ghostFromFiner,       &
      &          property   = 0_long_k,                &
      &          sourceProc = tree%global%myPart+1,    &
      &          pos        = foundPos,                &
      &          wasAdded   = wasAdded                 )

    if( wasAdded ) then

      updated = .true.
      children = tem_directChildren( elemID )
      childPos = 0 ! reset child positions. non-existing children are 0
      ! reset property
      property = 0_long_k
      ! if added elemID is more than level coarser than available child treeID
      ! in original treeID list then add all children between level and
      ! neighLevel
      do iChild = 1, 8

        ! Return position in the treeIDlist
        childPath  = tem_PathOf( children( iChild ))
        childPos( iChild ) = tem_PosOfPath( childPath, tree%pathList )

        if( childPos( iChild ) < 0 ) then
           ! This child does NOT exists, recusively add it as a ghostFromFiner.
           call add_ghostFromFiner( elemID    = children( iChild ),    &
             &                      levelDesc = levelDesc,             &
             &                      minLevel  = minlevel,              &
             &                      tree      = tree,                  &
             &                      foundPos  = childPos( iChild ),    &
             &                      updated   = childUpdated,          &
             &                      stencil   = stencil                )
           ! Unify all properties of the children
           property = ieor( property,                  &
             &              levelDesc(level+1)         &
             &                %elem                    &
             &                %property                &
             &                %val( childPos(iChild) ) )
           updated = ( updated .or. childUpdated )
        else
          ! This child is a Fluid, i.e. already exists in element list
          childPos( iChild ) = PositionOfVal( &
            &   levelDesc( level+1 )%elem%tID, children(iChild) )
        end if
      end do ! iChild = 1, 8

      ! Now reconstruct current element's neighborhood
      ! based on the children's information
      call tem_find_BCs_fromFiner( childPos     = childPos,   &
        &                          sourceLevel  = level + 1,  &
        &                          targetLevel  = level,      &
        &                          targetPos    = foundPos,   &
        &                          levelDesc    = levelDesc,  &
        &                          minLevel     = minLevel,   &
        &                          stencil      = stencil     )

    else
      ! ghostFromFiner element was not added
      ! Hence, we found available information which now has to be returned
      updated = .false.
    end if

  end subroutine add_ghostFromFiner
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Inherit the neighborhood from the sourceELem to the targetElem
  !! Note that targetElem is inout, as it might have already values assigned.
  !!
  subroutine tem_find_BCs_fromFiner( childPos, sourceLevel, targetLevel, &
    &                                targetPos, levelDesc, minLevel,     &
    &                                stencil                             )
    ! -------------------------------------------------------------------- !
    !> position of all childs in the levelDesc elem tID list
    integer, intent(in) :: childPos(8)
    !> level of child
    integer, intent(in) :: sourceLevel
    !> level of parent
    integer, intent(in) :: targetLevel
    !> added position of parent in the levelDesc elem tID list
    integer, intent(in) :: targetPos
    !> minimum level in the tree
    integer, intent(in) :: minLevel
    !> the level descriptor to be filled
    type(tem_levelDesc_type ) :: levelDesc(minLevel:)
    !> current stencil definition
    type(tem_stencilHeader_type), intent(in) :: stencil
    ! -------------------------------------------------------------------- !
    integer :: dir
    ! Tangential direction iterators
    integer :: iDirX, iDirY, iDir
    integer :: dirX !< first tangential direction
    integer :: dirY !< second tangential direction
    integer :: iStencilElem  !< stencil element iterator
    integer :: childCoord(4)
    ! the type of the current stencil link (side, edge, corner)
    integer :: linktype
    integer :: iDirNormal
    integer :: iStencil, addedPos
    type(tem_stencilElement_type) :: tStencil
    integer(kind=long_k) :: tNeighID
    ! -------------------------------------------------------------------- !

    if ( .not. allocated( levelDesc( targetLevel )%elem%stencil%              &
      &                                              val(targetPos)%val )) then
      call init( me     = levelDesc( targetLevel )%elem%stencil%val(targetPos),&
        &        length = 1 )
    end if

    call init( me = tStencil, QQN = stencil%QQN, headerPos = 1 )
    call append( me  = levelDesc( targetLevel )%elem%stencil%val(targetPos),   &
      &          val = tStencil )
    if( levelDesc( targetLevel )%elem%neighID%val(targetPos)%nVals > 0 ) then
      write(dbgUnit(2),*) 'Warning: tem_find_BCs_fromFiner has already' &
        &                       //' encountered'
      write(dbgUnit(2),*) '         existing neighbor IDs for the stencil'
    end if
    iStencil = levelDesc( targetLevel )%elem%stencil%val( targetPos )%nVals

    childCoord(4) = 1
    do iStencilElem = 1, stencil%QQN
      ! Reset current direction
      ! only 0 can not be used for inheritance of BCs later because boundaries
      ! are stored as negative treeIDs. Transformation of 0 into a HUGE
      ! negative value enables easy inheritance.
      tNeighID = -huge( 0_long_k )
      ! Count the non-zero entries to decide the kind of the current
      ! link direction in the stencil. (only direct neighbors)
      ! Distinguishing three different kinds of neighbors:
      ! SIDE, EDGE and CORNER
      linktype = sum( abs(stencil%cxDir( :, iStencilElem )))
      select case(linktype)
      case(1)
        ! 1 non-zero
        ! SIDE, needs four children
        ! identify the normal direction of the side
        dir = maxloc( abs( stencil%cxDir(:, iStencilElem )), 1)
        ! Convert the stencil direction to child direction
        iDirNormal = stencilToChild( stencil%cxDir(dir, iStencilElem) )
        childCoord( dir ) = iDirNormal
        ! identify the tangential directions
        dirX = mod( dir,   3) + 1
        dirY = mod( dir+1, 3) + 1
        do iDirY = 0, 1
          childCoord( dirY ) = iDirY
          do iDirX = 0, 1
            childCoord( dirX ) = iDirX
            call update_childNeighborID( neighID      = tNeighID,            &
              &                          childCoord   = childCoord,          &
              &                          childPos     = childPos,            &
              &                          iStencil     = iStencil,            &
              &                          iStencilElem = iStencilElem,        &
              &                          elem  = levelDesc(sourceLevel)%elem )
          end do
        end do

      case(2)
        ! 2 non-zeroes
        ! EDGE, needs two children
        ! Determine the 0-direction (dimension with the zero entry)
        dir = minloc( abs( stencil%cxDir(:,iStencilElem )), 1)
        ! Get the edge directions
        dirX = mod( dir,   3) + 1
        dirY = mod( dir+1, 3) + 1
        childCoord( dirX ) = stencilToChild(stencil%cxDir(dirX,iStencilElem))
        childCoord( dirY ) = stencilToChild(stencil%cxDir(dirY,iStencilElem))
        do iDir = 0, 1
          childCoord( dir ) = iDir
          call update_childNeighborID( neighID      = tNeighID,            &
            &                          childCoord   = childCoord,          &
            &                          childPos     = childPos,            &
            &                          iStencil     = iStencil,            &
            &                          iStencilElem = iStencilElem,        &
            &                          elem  = levelDesc(sourceLevel)%elem )
        enddo

      case(3)
        ! No zero at all, all three directions have a offset
        ! CORNER, just a single child is connected to this link
        childCoord( 1:3 ) = stencilToChild( stencil%cxDir(:,iStencilElem) )
        call update_childNeighborID( neighID      = tNeighID,            &
          &                          childCoord   = childCoord,          &
          &                          childPos     = childPos,            &
          &                          iStencil     = iStencil,            &
          &                          iStencilElem = iStencilElem,        &
          &                          elem  = levelDesc(sourceLevel)%elem )
      end select

      ! Append the neighID of virtual parent
      call append( me  = levelDesc(targetLevel)%elem%neighID%val(targetPos), &
        &          val = tNeighID,                                           &
        &          pos = addedPos                                            )
      levelDesc(targetLevel)     &
        &  %elem                 &
        &  %stencil              &
        &  %val(targetPos)       &
        &  %val(iStencil)        &
        &  %tIDpos(iStencilElem) = addedPos

      ! to append the neighbors of ghosts to required list
      !call append( me       = levelDesc( targetLevel )%require, &
      !  &          val      = tNeighID,                         &
      !  &          pos      = addedPos )

    end do

  end subroutine tem_find_BCs_fromFiner
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Returns the absolute position in the total list of a given treeID
  !! opposed to PosOfId, where the relative position in one of the separate
  !! lists is returned. Herefore, total list has to be created beforehand.
  !!
  function tem_treeIDinTotal( tID, levelDesc, eType ) result( elemPos )
    ! -------------------------------------------------------------------- !
    !> the element you are looking for
    integer(kind=long_k), intent(in) :: tID
    !> the descriptor you use for searching
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> element type
    integer, intent(in), optional :: eType
    !> return position of tID in levelDesc%total list
    integer :: elemPos
    ! -------------------------------------------------------------------- !
    integer :: eType_loc
    ! -------------------------------------------------------------------- !

    if( present( eType )) then
      eType_loc = eType
    else
      eType_loc = tem_eTypeOfID( tID, levelDesc%elem )
    end if

    elemPos = 0
    if( eType_loc > 0 ) then
      elemPos = tem_PositionInSorted( me    = levelDesc%total,                 &
        &                             lower = levelDesc%offset(1, eType_loc)+1,&
        &                             upper = levelDesc%offset(2, eType_loc),  &
        &                             val   = tID )
    end if
    elemPos = max( elemPos, 0 )

  end function tem_treeIDinTotal
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> create the intermediate, static list totalPnt, which holds pointers to the
  !! elem%TID list, but in an ordered fashion. The order is the same as it will
  !! be in the total list later on, i.e.: fluid, ghostFC, ghostFF, halo.
  !! this four sub-lists are within sorted by their treeID.
  !! Additionally, the process-wise collections of halo elements are collected
  !! into haloList by grouping the treeIDs according to their belonging process
  !!
  subroutine identify_lists( me )
    ! -------------------------------------------------------------------- !
    !> the level descriptor to be filled
    type(tem_levelDesc_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !
    integer :: iElem, indElem
    integer :: iPnt( eT_minNumber:eT_maxNumber ), eType, iVal
    ! -------------------------------------------------------------------- !
    ! Destroy lists
    call tem_halo_destroy(me%haloList)
    ! init lists
    call tem_halo_init(me%haloList)

    ! --------------------------------------------------------------------------
    ! 1. count nElems
    me%nElems = sum( me%elem%nElems(eT_minRelevant:eT_maxRelevant) ) &
      &       + me%elem%nElems(et_distributedGhostFromFiner)

    call set_offsets( me       = me%offset(:,:),                               &
      &               nFluids  = me%elem%nElems( eT_fluid ),                   &
      &               nGhostFC = me%elem%nElems( eT_ghostFromCoarser ),        &
      &               nGhostFF = me%elem%nElems( eT_ghostFromFiner )           &
      &                        + me%elem%nElems(eT_distributedGhostfromFiner), &
      &               nHalos   = me%elem%nElems( eT_halo )                     )

    ! --------------------------------------------------------------------------

    if ( allocated( me%totalPnt )) deallocate( me%totalPnt )
    allocate( me%totalPnt(me%nElems) )
    me%totalPnt = -1

    ! Reset pointers to current eType element in the total list
    ! @todo: first add fluid. do not have to follow sorted order here.
    iPnt = 0
    do indElem = 1, me%elem%tID%nVals
      ! Access sorted list to maintain locality of space-filling curve order
      iElem = me%elem%tID%sorted( indElem )
      eType = me%elem%eType%val( iElem )
      if ( eType == eT_distributedGhostFromFiner) then
        eType = eT_ghostFromFiner
        call changeType( me%elem, iElem, eT_ghostFromFiner )
      end if

      ! increase counter for current element
      iPnt( eType ) = iPnt( eType ) + 1
      ! is the eType a required element (fluid, ghost or halo)?
      ! Get rid of nonExistent elements here.
      if ( eType >= eT_minRelevant .and. eType <= eT_maxRelevant ) then
        ! Add sorted position of tID to the pointer
        me%totalPnt( me%offset(1,eType) + iPnt(eType) ) = iElem
        ! And add the haloList for halo elements
        if ( eType == eT_halo ) then
          ! get the process from where to get the element
          ! Create an entry in the process list for the current source process
          call tem_halo_append( me     = me%haloList,                   &
            &                   proc   = me%elem%sourceProc%val(iElem), &
            &                   elemPos= iElem                          )
        end if ! eT_halo
      end if !  ( eType >= eT_minRelevant )
    end do ! indElem

    ! Security check
    ! Check if there are no entries < 1
    do iElem = 1, me%nElems
      if( me%totalPnt( iElem ) < 1 ) then
        write(dbgUnit(1),*) "Error: Found index < 1 in the totalPnt array."
        write(dbgUnit(1),*) 'Abort!'
        write(dbgUnit(1),*) 'offset: '&
          &     //trim(tem_toStr(me%offset(1,:),'; '))

        do iVal = 1, me%nElems
          write(dbgUnit(1),*)'totalPnt: ',me%totalPnt(iVal)
        end do
        call tem_abort
      end if
    enddo

  end subroutine identify_lists
! ****************************************************************************** !


! ****************************************************************************** !
  !> Set the offsets for accessing totallist, invsorted etc. arrays for
  !! fluids, ghosts and halos
  !!
  subroutine set_offsets( me, nFluids, nGhostFC, nGhostFF, nHalos )
    ! ---------------------------------------------------------------------------
    !> element type offsets
    integer, intent(out) :: me( 2, eT_minRelevant:eT_maxRelevant )
    !>
    integer, intent(in) :: nFluids, nGhostFC, nGhostFF, nHalos
    ! ---------------------------------------------------------------------------

    ! store the offset for the beginning ...
    me( 1, eT_fluid )             = 0
    me( 1, eT_ghostFromCoarser )  = nFluids
    me( 1, eT_ghostFromFiner )    = nFluids + nGhostFC
    me( 1, eT_halo )              = nFluids + nGhostFC + nGhostFF

    ! ... and the end of the list
    me( 2, eT_fluid )            = me( 1, eT_ghostFromCoarser )
    me( 2, eT_ghostFromCoarser ) = me( 1, eT_ghostFromFiner )
    me( 2, eT_ghostFromFiner )   = me( 1, eT_halo )
    me( 2, eT_halo )             = nFluids + nGhostFC + nGhostFF + nHalos

    ! me( :, eT_distributedGhostFromFiner ) = me( :, eT_GhostFromFiner )

  end subroutine set_offsets
! ****************************************************************************** !


! ****************************************************************************** !
  !> Create the level-wise total lists.
  !!
  !! They consist of
  !! fluid elements + ghost elements + halo elements
  !! and are sorted by the treeIDs
  !! Also create the property lists
  !!
  subroutine assemble_lists( me, minLevel, maxLevel, tree )
    ! ---------------------------------------------------------------------------
    !> Minimal level in the mesh
    integer, intent(in) :: minlevel, maxLevel
    !> Level descriptor to fill
    type( tem_levelDesc_type ), intent(inout) :: me(minlevel:maxLevel)
    !> tree information
    type(treelmesh_type), intent(in) :: tree
    ! ---------------------------------------------------------------------------
    integer :: iLevel, iIndex
    integer(kind=long_k) :: tID
    ! ---------------------------------------------------------------------------
    write(logUnit(5),*) 'Assembling total list ...'

    do iLevel = minlevel, ubound(me,1)
      allocate( me( iLevel )%total(   me( iLevel )%nElems ))
      allocate( me( iLevel )%baryOfTotal( me(iLevel)%nElems, 3 ) )
      allocate( me( iLevel )%property(me( iLevel )%nElems ))
      allocate( me( iLevel )%pntTID( me( iLevel )%elem%nElems( eT_fluid )))

      write(logUnit(5),"(A,I0)") '     Level: ', iLevel
      write(logUnit(5),"(A,I0)") '    nElems: ', me(iLevel)%nElems
      call print_nElems( me(iLevel)%elem%nElems(eT_minRelevant:eT_maxRelevant),&
        &                logUnit(5) )

      ! Counter across all lists
      do iIndex = 1, me( iLevel )%nElems
        tID = me( iLevel )%elem%tID%val( me( iLevel )%totalPnt(iIndex) )
        ! sorted treeID list
        me( iLevel )%total( iIndex ) = tID

        ! barycenter of treeID in total list
        me( iLevel )%baryOfTotal( iIndex, : ) = tem_BaryOfID( tree, tID )

        ! property
        me( iLevel )%property( iIndex ) = &
          &  me( iLevel )%elem%property%val( me( iLevel )%totalPnt(iIndex))
      end do

      do iIndex = 1, me( iLevel )%elem%nElems( eT_fluid )
        ! pointer to original tree%treeID list
        me( iLevel)%pntTID( iIndex ) = &
          & me( iLevel )%elem%pntTID%val( me( iLevel )%totalPnt( iIndex ))
      end do
      call destroy( me( iLevel )%elem%pntTID )

    end do ! iLevel

  end subroutine assemble_lists
! ****************************************************************************** !


! ****************************************************************************** !
  !> identify additionally required neighbor elements
  !! run over the 'require' list of elements, which was accumulated before in
  !! init_elemLevels. The list includes neighbor elements of stencil neighbors,
  !! for stencils with the requireNeighNeigh attribute set.
  !! This is needed for example for LBM boundary stencil elements, which in turn
  !! require their compute stencil neighborhood to allow PULL operations from
  !! there
  !!
  !! What exactly is the require list for?
  !! - Used ONLY for boundary stencil with higher order neighbors i.e
  !!   only when require nVals > 0
  !!
  subroutine identify_additionalNeigh( tree, proc, levelDesc, pathFirst,        &
    &                                  pathLast, stencil )
    ! ---------------------------------------------------------------------------
    !> the global tree
    type(treelmesh_type), intent(in) :: tree
    !> Process description to use.
    type(tem_comm_env_type), intent(in) :: proc
    !> the level descriptor to be filled
    type(tem_levelDesc_type), intent(inout) :: levelDesc(tree%global%minlevel:)
    !> first treeID path in every process
    type(tem_path_type), intent(in) :: pathFirst(:)
    !> last treeID path in every process
    type(tem_path_type), intent(in) :: pathLast(:)
    !> the compute stencil, for which the additional neighbors are reconstructed
    type(tem_stencilHeader_type), intent(in) :: stencil
    ! ---------------------------------------------------------------------------
    integer :: iLevel, posInElem, neighPos, elemPos, iNeighElem, iElem
    integer(kind=long_k) :: treeID
    ! ---------------------------------------------------------------------------

    call tem_horizontalSpacer( fUnit = dbgUnit(1), before = 1 )
    write(dbgUnit(3),*) 'Inside routine: identify_additionalNeigh'

    ! The position of the compute stencil
    do iLevel = tree%global%minlevel, tree%global%maxLevel
      ! Run over the additionally required element list
      do iElem = 1, levelDesc( iLevel )%require%nVals
        ! get the position of the treeID in the element list
        posInElem = PositionOfVal( me  = levelDesc( iLevel )%elem%tID,             &
          &                        val = levelDesc( iLevel )%require%val( iElem ))
        ! get the element position
        if( posInElem > 0 ) then
          if ( levelDesc( iLevel )%elem%eType%val( posInElem ) > 0) then
            ! Run over all the neighbors of the compute stencil
            do iNeighElem = 1, stencil%QQN
              ! position of neighbor treeID in dynamic array of neighID
              neighPos = levelDesc( iLevel )%elem%stencil%val( posInElem )         &
                &                         %val(1)%tIDpos( iNeighElem )
              if( neighPos > 0 ) then
                treeID = &
                  & levelDesc( iLevel )%elem%neighID%val( posInElem )%val( neighPos )
                call identify_elements( TreeID     = treeID,        &
                  &                     tree       = tree,           &
                  &                     pathFirst  = pathFirst,      &
                  &                     pathLast   = pathLast,       &
                  &                     levelDesc  = levelDesc,      &
                  &                     elemPos    = elemPos,        &
                  &                     proc       = proc,           &
                  &                     nesting    = 0,              &
                  &                     stencil    = stencil,        &
                  &       skip_add_additionalGhost = .true.             )
              else ! neighPos =< 0, i.e. no additional neighbor
                elemPos = 0
              end if
              levelDesc( iLevel )%elem%stencil%val( posInElem )%val(1)    &
                &                %totalPos( iNeighElem ) = elemPos
            end do ! neighPos > 0
          else ! eType <= 0
            write(logUnit(8),*) 'Can not find additional neighbor: ', treeID
          end if  ! valid elemType
        else ! posInElem <= 0
          write(logUnit(2),*) 'element which requires additional neighbor is not part of total list'
        end if ! posInElem > 0
      end do ! iElem
    end do ! iLevel

    write(dbgUnit(1),*) 'Leave routine: identify_additionalNeigh'
    call tem_horizontalSpacer( fUnit = dbgUnit(1), after = 1 )

  end subroutine identify_additionalNeigh
! ****************************************************************************** !


! ****************************************************************************** !
  !> exchange the requested treeIDs between all MPI processs
  !!
  subroutine communicate_elements( tree, proc, me, commPattern,          &
    &                              pathFirst, pathLast, computeStencil )
    ! ---------------------------------------------------------------------------
    !> the global tree
    type(treelmesh_type), intent(in) :: tree
    !> Process description to use.
    type(tem_comm_env_type), intent(in) :: proc
    !> the level descriptor to be filled
    type(tem_levelDesc_type), intent(inout) :: me(tree%global%minlevel:)
    !> the communication pattern used
    type(tem_commPattern_type), intent(in)    :: commPattern
    !> first and last treeID path in every process
    type(tem_path_type), intent(in) :: pathFirst(:), pathLast(:)
    !> stencil definition
    type(tem_stencilHeader_type), intent(in) :: computeStencil(:)
    ! ---------------------------------------------------------------------------
    integer :: iLevel, iErr, nProcs, iProc
    integer,allocatable :: nHalos(:)
    integer :: nIterations
    logical :: redo ! locally indicate if another iteration has to be performed
    logical :: redo_global    ! global indicator for another iteration
    ! ---------------------------------------------------------------------------

    call tem_horizontalSpacer( fUnit = logUnit(3) )
    write(logUnit(3),*) 'Communicating elements ...'

    call tem_horizontalSpacer( fUnit = dbgUnit(1), before = 1 )
    write(dbgUnit(1),*) 'Communicating elements ...'

    ! ----------- exchange number of requesting / requested treeIDs -------------
    nIterations = 0
    do  !iIter = 1, 1 ! exchange halos until no new elements are created
      nIterations = nIterations + 1
      write(logUnit(4),"(A,I0)") 'Halo count exchange iteration: ', nIterations
      write(dbgUnit(1),"(A,I0)") 'Halo count exchange iteration: ', nIterations
      redo = .false.
      ! communicate nHalos and allocate me%buffer
      call communicate_nElemsToTransfer( me, proc, tree%global%minLevel, &
        &                                tree%global%maxLevel )
      !
      ! ---done.--- exchange number of requesting / requested treeIDs ----------!
      write(dbgUnit(1),"(A)") 'Halo count exchange done!'
      write(logUnit(5),*) '  Done communicating nElems to Transfer'

      ! -----------                                                   ----------!
      !
      ! now we request the halo cells from the mpi processes,
      ! so they can send us their dependencies for these halo cells.
      ! We allow only dependencies with a difference of one level
      ! from lower to higher refinement levels and an arbitrary difference
      ! of refinement level from higher to lower refinement level.
      ! Since we get only the leaves as a dependency from the other mpi process
      ! we have to figure out later which cells we have to add locally
      !
      ! 1)  get the number of cells we will receive from the active processes,
      !     since we have
      !     to communicate only to the mpi processs which have/need informations
      ! 2)  exchange the treeIDs of the halos
      !
      ! First collapse the me%halos.
      ! We remove the processes where we don't have any elements to exchange
      ! with.
      ! Inverse Communciation (send to sourceProc, recv from targetProc)
      do iLevel = tree%global%minLevel, tree%global%maxLevel
        ! Send treeIDs
        write(logUnit(5),"(A,I0)") '  Requesting remote halos on level ', iLevel
        call request_remoteHalos( levelDesc = me,                  &
          &                       tree      = tree,                &
          &                       iLevel    = iLevel,              &
          &                       pathFirst = pathFirst,           &
          &                       pathLast  = pathLast,            &
          &                       stencil   = computeStencil(1),   &
          &                       proc      = proc )
        if (tem_logging_isActive(main_debug%logger, 7)) then
          call tem_elemList_dump( me = me( iLevel )%elem,      &
            &                     nUnit = dbgUnit(5),                 &
            &                     stencil = .true.,                   &
            &                     string = 'after request remoteHalos' )
        end if
        write(logUnit(5),*) '  Done requesting remote halos.'

        nProcs = me( iLevel )%haloList%partnerProc%nVals
        if( allocated( nHalos )) deallocate( nHalos )
        allocate( nHalos( nProcs ))
        nHalos(:nProcs) = me( iLevel )%haloList%halos%val(:nProcs )%nVals
        write(logUnit(5),*) '  Identifying lists'
        call identify_lists( me(iLevel) )
        ! If nHalos or nProcs changes, then do request again.
        ! As long as any level has change, do request again.
        redo = redo .or. ( any( me( iLevel )%haloList%halos%val(1:nProcs)%nVals   &
          &                     /= nHalos(1:nProcs) )  &
          &      .or. (nProcs /= me(iLevel)%haloList%partnerProc%nVals) )
      end do !iLevel

      write(logUnit(6),*) '  Allreduce to check if changes occurred on any' &
        &                 //' process'
      ! ------------------------------------------------------------------------
      ! JUROPA work-around for crash in the mpi_allreduce
      ! call mpi_barrier( proc%comm, iErr )
      ! ------------------------------------------------------------------------
      ! Determine among all neighbor processes, if further iterations required
      call mpi_allreduce( redo, redo_global, 1, mpi_logical, mpi_lor,          &
        &                 proc%comm, iErr )

      if ( .not. redo_global ) exit
    end do !exchange halos
    write(logUnit(3),*) 'Done exchanging number of elements to communicate.'

    !!   Now each process knows, which halos are requested.
    !!   Continue with identifying the actual leaf elements, which
    !!   are then communicated

    write(logUnit(5),"(A)") 'Return halo counts and Redefine halos ...'
    do iLevel = tree%global%minLevel, tree%global%maxLevel
      write(logUnit(5),"(A,I0)") '  Returning halo counts on level ', iLevel
      ! Receive the number of really existing halo elements
      call return_haloCounts( sendbuffer = me( iLevel )%sendbuffer,     &
        &                     recvbuffer = me( iLevel )%recvbuffer,     &
        &                     comm       = proc%comm )
      call return_haloCounts( &
        &    sendbuffer = me( iLevel )%sendbufferFromCoarser, &
        &    recvbuffer = me( iLevel )%recvbufferFromCoarser, &
        &    comm       = proc%comm )
      call return_haloCounts( &
        &    sendbuffer = me( iLevel )%sendbufferFromFiner,   &
        &    recvbuffer = me( iLevel )%recvbufferFromFiner,   &
        &    comm       = proc%comm )

      ! reset the halos in the elem list
      do iProc = 1, me( iLevel )%haloList%PartnerProc%nVals
        ! First declare all local halos as non-existent, and set only those
        ! actually provided by the remote process.
        call changeType( me(iLevel)%elem, &
          &              me(iLevel)%haloList%halos%val(iProc)%nVals, &
          &              me(iLevel)%haloList%halos%val(iProc)%val(:),&
          &              eT_nonExisting )
      end do

      write(logUnit(5),*) '  Redefining halos ... '
      ! Receive the number of really existing halo elements
      call redefine_halos( levelDesc      = me( iLevel ),               &
        &                  sendbuffer     = me( iLevel )%sendbuffer,    &
        &                  recvbuffer     = me( iLevel )%recvbuffer,    &
        &                  commPattern    = commPattern,                       &
        &                  computeStencil = computeStencil,                    &
        &                  proc           = proc )

      call redefine_halos( levelDesc      = me( iLevel ),               &
        &                  sendbuffer     = me( iLevel )%sendbufferFromCoarser,  &
        &                  recvbuffer     = me( iLevel )%recvbufferFromCoarser,  &
        &                  commPattern    = commPattern,                       &
        &                  computeStencil = computeStencil,                    &
        &                  proc           = proc )

      call redefine_halos( levelDesc      = me( iLevel ),               &
        &                  sendbuffer     = me( iLevel )%sendbufferFromFiner,   &
        &                  recvbuffer     = me( iLevel )%recvbufferFromFiner,   &
        &                  commPattern    = commPattern,                       &
        &                  computeStencil = computeStencil,                    &
        &                  proc           = proc )
      call identify_lists( me(iLevel) )
    end do

    ! dump debug output
    if (tem_logging_isActive(main_debug%logger, 7)) then
      do iLevel = tree%global%minLevel, tree%global%maxLevel
        call tem_elemList_dump( me      = me( iLevel )%elem,   &
          &                     nUnit   = dbgUnit(5),                 &
          &                     stencil = .true.,                     &
          &                     string  = 'after redefine remoteHalos' )
      end do
    end if

    write(logUnit(3),*) 'Done with communication of elements. '
    call tem_horizontalSpacer( fUnit = logUnit(3), after = 1 )

    write(dbgUnit(1),*) 'Done Communicating elements ...'
    call tem_horizontalSpacer( fUnit = dbgUnit(1), after = 1 )

  end subroutine communicate_elements
! ****************************************************************************** !


! ****************************************************************************** !
  !> Map requested halo to a position in my local fluid list or
  !! add recursively ghosts until I reach valid fluid elements
  !! return type of added element in levelPos(2)
  !! Also, non-existing elements are reported as such (levelPos(2))
  !!
  subroutine identify_halo( haloTreeID, elemPos, haloLevel, levelDesc, tree,   &
    &                       updated,   nesting, minLevel, stencil )
    ! ---------------------------------------------------------------------------
    !>
    type(treelmesh_type), intent(in) :: tree
    !>
    integer, intent(in) :: minLevel
    !>
    integer, intent(in) :: nesting
    !> neighboring treeID
    integer(kind=long_k), intent(in) :: haloTreeID
    !> type and position in list of found treeID
    integer, intent(out) :: elemPos
    !>
    integer, intent(out) :: haloLevel
    !>
    logical, intent(out) :: updated
    !>
    type(tem_levelDesc_type), intent(inout) :: levelDesc(minlevel:)
    !>
    type(tem_stencilHeader_type), intent(in) :: stencil
    ! ---------------------------------------------------------------------------
    haloLevel = tem_levelOf( haloTreeID )
    elemPos = PositionOfVal( me  = levelDesc( halolevel )%elem%tID,            &
      &                      val = haloTreeID )

    if( elemPos == 0 ) then
      ! requested halo does not exist and has to be created somehow
      call identify_local_element( targetID       = haloTreeID,                &
        &                          levelDesc      = levelDesc,                 &
        &                          tree           = tree,                      &
        &                          elemPos        = elemPos,                   &
        &                          nesting        = nesting,                   &
        &                          updated        = updated,                   &
        &                          minLevel       = minLevel,                  &
        &                          stencil = stencil )
    end if

  end subroutine identify_halo
! ****************************************************************************** !


! ****************************************************************************** !
  !> deallocate all dispensable dynamic data structures and arrays
  !!
  subroutine tem_cleanupDependencyArrays( levelDesc )
    ! ---------------------------------------------------------------------------
    !> the level descriptor
    type(tem_levelDesc_type) , intent(inout) :: levelDesc(:)
    ! ---------------------------------------------------------------------------
    integer :: iLevel
    ! ---------------------------------------------------------------------------

    do iLevel = 1,size(levelDesc)
      ! Kill empty fluid, ghost and halo lists
      call tem_halo_destroy( levelDesc(iLevel)%haloList )
      call destroy( levelDesc(iLevel)%require )
    end do

  end subroutine tem_cleanupDependencyArrays
! ****************************************************************************** !


! ****************************************************************************** !
  !> deallocate the stencil treeID neighbor arrays for each element
  !! This routine can only be called after build_horizontalDependencies
  !!
  subroutine tem_cleanup_arrays( levelDesc )
    ! ---------------------------------------------------------------------------
    !> the level descriptor
    type(tem_levelDesc_type) , intent(inout) :: levelDesc(:)
    ! ---------------------------------------------------------------------------
    integer :: iLevel
    ! ---------------------------------------------------------------------------

    ! loop over all levels and destroy the corresponding elem types
    do iLevel = 1,size(levelDesc)
      call destroy( me = levelDesc( iLevel )%elem )
    end do

  end subroutine tem_cleanup_arrays
! ****************************************************************************** !


! ****************************************************************************** !
  !> Build the vertical dependencies of ghost elements
  !!
  subroutine tem_build_verticalDependencies( levelDesc, minlevel, maxLevel )
    ! ---------------------------------------------------------------------------
    !> Level range
    integer, intent(in) :: minlevel, maxLevel
    !> the level descriptor
    type(tem_levelDesc_type), intent(inout) :: levelDesc(minlevel:maxLevel)
    ! ---------------------------------------------------------------------------
    integer(kind=long_k) :: ghostID, parentID, childTreeID(8)
    integer :: iLevel, iElem, iChild
    integer :: offset, posInTotal, childNum
    ! ---------------------------------------------------------------------------

    write(logUnit(3),*) 'Build vertical dep ...'
    write(dbgUnit(1),*) ''
    write(dbgUnit(1),*) ' into tem_build_verticalDependencies ... '

    ! Start with GhostFromCoarser
    do iLevel = minLevel+1, maxLevel

      allocate( levelDesc( iLevel )%depFromCoarser(         &
        &       levelDesc( iLevel )%elem%nElems( eT_ghostFromCoarser )))

      offset = levelDesc( iLevel )%offset(1,eT_ghostFromCoarser)

      do iElem = 1, levelDesc( iLevel )%elem%nElems( eT_ghostFromCoarser )

        ! ghost tree ID
        ghostID  = levelDesc( iLevel )%total( offset+iElem )
        parentID = tem_parentOf( ghostID )
        ! the order as a child for this target
        childNum = tem_ChildNumber( ghostID )
        ! Store the child number and coord
        ! calculate relative position of target to its parent:
        ! Bary of parent is considered as 0,0,0
        ! dx on coarser level is considered as 1
        !   childNum                 childPosition(childNum, 1:3)
        !          1                                   -1, -1, -1
        !          2                                    1, -1, -1
        !          3                                   -1,  1, -1
        !          4                                    1,  1, -1
        !          5                                   -1, -1,  1
        !          6                                    1, -1,  1
        !          7                                   -1,  1,  1
        !          8                                    1,  1,  1
        levelDesc( iLevel )%depFromCoarser( iElem )%childNum = childNum
        levelDesc( iLevel )%depFromCoarser( iElem )%coord(:) = &
          &                 0.25_rk * dble(childPosition(childNum, 1:3))

        ! position of ghost parent ID in total list in level-1
        posInTotal = tem_treeIDinTotal( parentID, levelDesc( iLevel-1 ) )

        ! Add the requested TreeID to the halos and add dependencies
        call appendGhostDependency(                                            &
          &           sourcePos   = posInTotal,                                &
          &           sourceLevel = iLevel-1,                                  &
          &           tgtDep      = levelDesc( iLevel )%depFromCoarser( iElem ))
      end do ! iElem
    end do ! iLevel

    ! Continue with GhostFromFiner
    do iLevel = minLevel, maxLevel
      allocate( levelDesc( iLevel )%depFromFiner(                              &
        &       levelDesc( iLevel )%elem%nElems( eT_ghostFromFiner ) ))
      offset = levelDesc( iLevel )%offset(1,eT_ghostFromFiner)
      do iElem = 1, levelDesc( iLevel )%elem%nElems( eT_ghostFromFiner )
        ! coarser ghost ID
        ghostID  = levelDesc( iLevel )%total( offset+iElem )
        ! children of coarser ghost
        childTreeID(:) = tem_directChildren( ghostID )
        ! append position of each children ID to depFromFiner list
        do iChild = 1, 8
          posInTotal = tem_treeIDinTotal( childTreeID( iChild ),               &
            &                             levelDesc( iLevel+1 ))
          if( posInTotal > 0 ) then
            ! if the child treeID exists, add dependency
            call appendGhostDependency(                                        &
              &         sourcePos   = posInTotal,                              &
              &         sourceLevel = iLevel+1,                                &
              &         tgtDep      = levelDesc( iLevel )%depFromFiner( iElem ))
          end if
        end do ! iChild = 1, 8
      end do ! iElem = 1, elem%nElems( eT_ghostFromFiner )
    end do ! iLevel = minLevel, maxLevel

    write(dbgUnit(1),*) ' leave tem_build_verticalDependencies ... '
    write(dbgUnit(1),*) ''

  end subroutine tem_build_verticalDependencies
! ****************************************************************************** !


! ****************************************************************************** !
  !> Building neighor array
  !! @todo: is neighID and stencil in element_type still used after this?
  !!
  subroutine tem_build_horizontalDependencies( iStencil, levelDesc, tree,      &
    &                                          computeStencil )
    ! ---------------------------------------------------------------------------
    !> Tree representation of your mesh.
    type(treelmesh_type), intent(in) :: tree
    !> Level descriptor for each level of your mesh (starting from minimum
    !! level).
    type(tem_levelDesc_type), intent(inout) :: levelDesc(tree%global%minLevel:)
    !> Index of your neighbor list??
    integer, intent(in) :: iStencil
    !> The stencil you build the horizontal dependencies for.
    type(tem_stencilHeader_type), intent(in) :: computeStencil
    ! ---------------------------------------------------------------------------
    integer :: iLevel, iIndex
    integer :: nElemsWithNeigh
    ! ---------------------------------------------------------------------------

    write(logUnit(4),"(A,I0)") 'Build neighbor array for stencil: ', iStencil

    do iLevel = tree%global%minLevel,tree%global%maxLevel

      if( computeStencil%useAll .and. iStencil == 1 ) then
        ! JQ: is it possible the 1st stencil does not use all?
        ! all fluid and ghost and halo elements
        nElemsWithNeigh = levelDesc(iLevel)%nElems
      else if( computeStencil%useAll ) then
        ! @todo: the logic of this IF condition is not clear
        ! all fluid and ghost elements no halos
        nElemsWithNeigh =   levelDesc(iLevel)%nElems     &
          &               - levelDesc(iLevel)%elem%nElems( eT_halo )
      else
        ! only the fluids and ghosts??? that have the stencil iStencil
        nElemsWithNeigh = computeStencil%elemLvl( iLevel )%nVals
      end if

      write(logUnit(4),"(A,I2,A,I0)") ' level: ', iLevel, &
        &                             ', nElems to treat: ', nElemsWithNeigh

      ! allocate the nghElems array with the number of elements this stencil is
      ! used by
      allocate( levelDesc( iLevel )%neigh( iStencil )%            &
        &          nghElems( computeStencil%QQN, nElemsWithNeigh ) )

      iIndex = 0
      if( computeStencil%QQN > 0 ) then
        ! we init everything by default with 0 to find a bug as fast as possible
        levelDesc( iLevel )%neigh( iStencil )%nghElems(:,:) = 0
        if (computeStencil%useAll) then
!          write(logUnit(5),*) '     fluid elements  ', iIndex
          call tem_build_listHorizontalDep(                                    &
            &                  levelDesc      = levelDesc( iLevel ),           &
            &                  iStencil       = iStencil,                      &
            &                  posInSortElem  = levelDesc( iLevel )%totalPnt,  &
            &                  nElems         = nElemsWithNeigh,               &
            &                  iIndex         = iIndex )

        else ! stencil not using all elements
!          write(logUnit(5),*) '     not all elems   ', iIndex
          call tem_build_treeHorizontalDep(                                    &
            &         levelDesc      = levelDesc( iLevel ),                    &
            &         computeStencil = computeStencil,                         &
            &         iStencil       = iStencil,                               &
            &         list           = computeStencil%elemLvl( iLevel )%val,   &
            &         nElems         = nElemsWithNeigh,                        &
            &         tree           = tree                                    )
        end if
      end if

    end do ! iLevel

  end subroutine tem_build_horizontalDependencies
! ****************************************************************************** !


! ****************************************************************************** !
  !> Building neighor array
  !!
  subroutine tem_debug_HorizontalDependencies( iStencil, levelDesc, tree,      &
    &                                          computeStencil )
    ! ---------------------------------------------------------------------------
    !> Tree representation of your mesh.
    type(treelmesh_type), intent(in) :: tree
    !> Level descriptor for each level of your mesh (starting from min level).
    type(tem_levelDesc_type), intent(inout) :: levelDesc(tree%global%minLevel:)
    !> Index of your neighbor list.
    integer, intent(in) :: iStencil
    !> The stencil you build the horizontal dependencies for.
    type(tem_stencilHeader_type) :: computeStencil
    ! ---------------------------------------------------------------------------
    integer :: iElem, iLevel, nElems, nEntries, iVal, elemPos, tIDpos
    integer :: neighPos, stencilSize, stencilPos
    integer(kind=long_k) :: tID, bufID, refID
    ! ---------------------------------------------------------------------------

    call tem_horizontalSpacer( fUnit = dbgUnit(1), before = 1 )
    write(dbgUnit(1),"(A   )") 'Get into tem_debug_horizontalDependencies routine'
    write(dbgUnit(1),"(A,I3)") 'Check horizontal dependencies stencil:',iStencil

    write(logUnit(1),*)' Checking, if the horizontal dependencies were set'// &
      &                ' correctly'

    do iLevel = tree%global%minLevel,tree%global%maxLevel
      nEntries = size( levelDesc( iLevel )%neigh(iStencil)%nghElems, 1)
      nElems   = size( levelDesc( iLevel )%neigh(iStencil)%nghElems, 2)

      write(dbgUnit(1),"(3(A,I0))") ' iLevel: ', iLevel, &
        &                    ' nEntries: ', nEntries, ', nElems: ', nElems
      if( nEntries > 0 ) then
        write(dbgUnit(1),"(A11, A10)") 'fluid, ', 'neighbors:'
        do iElem = 1, nElems
          ! First get the treeID
          if( computeStencil%useAll ) then
            tIDPos = iElem
            tID    = levelDesc( iLevel )%total( tIDpos )
          else
            if(  computeStencil%elemLvl( iLevel )%nVals <= 0 ) exit
            tIDPos = computeStencil%elemLvl( iLevel )%val( iElem )
            tID    = tree%treeID( tIDpos )
          end if
          ! Get the position of the treeID in the element list
          elemPos = PositionOfVal(levelDesc( iLevel )%elem%tID, tID )
          if( elemPos <= levelDesc( iLevel )%elem%nElems(eT_fluid) )then
            ! find the right stencil
            stencilPos = tem_stencil_getHeaderPos(                               &
              &           me  = levelDesc( iLevel )%elem%stencil%val( elemPos ), &
              &           val = iStencil )
          else
            ! set the stencil position to be the fluid stencil
            stencilPos = 1
          end if
!write(*,*)'in debug_HorDeps, level ', iLevel,' elemPos: ', elemPos
!write(*,*)'stencilPos: ', stencilPos, ' iStencil: ', iStencil
!write(*,*)'nVals: ', levelDesc( iLevel )%elem%stencil%val( elemPos )%nVals
!write(*,*)'val%pos', levelDesc( iLevel )%elem%stencil%val( elemPos )%val(:)%headerPos
          stencilSize = size( levelDesc( iLevel )%elem%stencil%val( elemPos )  &
            &                                       %val( stencilPos )%tIDpos )
!write(*,*)'Stencil entries: ', nEntries, ' stencilSize: ', stencilSize
          if (nEntries /= stencilSize) then
            write(dbgUnit(1),*)'stencil entries do not match for tID ',tID,    &
              &                ': ', nEntries, ' vs. ', stencilSize
          end if

          write( dbgUnit(3), '(I10,A)', advance='no' ) tID, ':'

          do iVal = 1, nEntries

            ! required neighbor ID
            tIDpos = levelDesc( iLevel )%elem%stencil%val( elemPos )%val(      &
              &                                     stencilPos )%tIDpos( iVal )
            if( tIDpos > 0 ) then
              refID = levelDesc( iLevel )%elem%neighID%val( elemPos )%val(     &
                &                                                      tIDpos )
            else
              refID = 0
            end if

            ! actually obtained neighbor ID
            neighPos = levelDesc( iLevel )%neigh(iStencil)%                    &
              &                                         nghElems( iVal, iElem )
            if( neighPos > 0 ) then
              bufID = levelDesc( iLevel )%total( neighPos )
            else
              bufID = neighPos
            endif
            write(dbgUnit(3),'(I10)',advance='no') bufID

            if( refID /= bufID ) then
              write(dbgUnit(3),'(a,i0,a)', advance='no') ' (',refID,')'
              if( bufID > 0_long_k ) then
                write(dbgUnit(1),*)' Non matching stencil entries!!'
                write(dbgUnit(1),*)'    original treeID:         ', tID
                write(dbgUnit(1),*)'    required stencil treeID: ', refID
                write(dbgUnit(1),*) '   encountered treeID:      ', bufID
                write(logUnit(1),*)' Found non-matching treeID in the neighborhood'
                write(logUnit(1),*)'    original treeID:         ', tID
                write(logUnit(1),*)'    required stencil treeID: ', refID
                write(logUnit(1),*)'    encountered treeID:      ', bufID
                call tem_abort()
              end if
            end if
          end do ! iVal = 1, nEntries
          write(dbgUnit(3),*)''

        end do ! iElem
        write(dbgUnit(3),*)''
      end if
    end do ! iLevel

    write(dbgUnit(1),*) 'Leave  tem_debug_HorizontalDependencies'
    call tem_horizontalSpacer( fUnit = dbgUnit(1), after = 1 )

  end subroutine tem_debug_HorizontalDependencies
! ****************************************************************************** !


! ****************************************************************************** !
  !> Update the neighor arrays depending on what is given in the element stencil
  !!
  !! The array levelDesc( iLevel )%neigh( iStenci )%nghElems( 1:QQN, 1:nElems )
  !! is being filled up here
  !!
  subroutine tem_build_treeHorizontalDep( iStencil, levelDesc, computeStencil, &
    &                                     list, nElems, tree )
    ! ---------------------------------------------------------------------------
    !> Level descriptor for each level of your mesh (starting from min level).
    type(tem_levelDesc_type), intent(inout) :: levelDesc
    !> tree information
    type(treelmesh_type), intent(in) :: tree
    !> Index of your neighbor list.
    integer, intent(in) :: iStencil
    !> The stencil you build the horizontal dependencies for.
    type(tem_stencilHeader_type), intent(in) :: computeStencil
    !> number of elements
    integer, intent(in) :: nElems
    !> stencil elemLvl points to sorted original treeID list
    integer, intent(in) :: list(:)
    ! ---------------------------------------------------------------------------
    integer :: iElem, iStencilElem
    integer :: elemPos, stencilPos, levelPos, totalPos
    integer(kind=long_k) :: tID
    logical :: missingNeigh
    real(kind=rk) :: bary(3) ! bary center of element
    ! ---------------------------------------------------------------------------

    ! run over the given element list to update
    do iElem = 1, nElems
      missingNeigh = .false.
      tID = tree%treeID( list( iElem ))
      ! Get the position of the treeID in the element list
      elemPos = PositionOfVal( me = levelDesc%elem%tID, val = tID )
      if( elemPos > 0 ) then
        totalPos = tem_treeIDinTotal(                                          &
          &                    tID       = tID,                                &
          &                    levelDesc = levelDesc,                          &
          &                    eType     = levelDesc%elem%eType%val( elemPos ))
        stencilPos = tem_stencil_getHeaderPos(                                 &
          &                       me  = levelDesc%elem%stencil%val( elemPos ), &
          &                       val = iStencil )
        if( stencilPos > 0 ) then
          do iStencilElem = 1, levelDesc%elem%stencil%val(elemPos)             &
            &                                        %val(stencilPos)%QQN
            ! Totalpos has before already been updated, so just read from
            ! totalPos
            levelPos = levelDesc%elem%stencil%val(elemPos)                     &
              &                              %val(stencilPos)                  &
              &                              %totalPos(iStencilElem)
            if( levelPos <= 0 .and. computeStencil%requireNeighNeigh ) then
              missingNeigh = .true.
            end if
            levelDesc%neigh( iStencil )%nghElems( iStencilElem, iElem )        &
              &                                                     = levelPos
          end do
        end if
        if( missingNeigh ) then
          bary = tem_BaryOfId(tree, tID)
          write(logUnit(10),"(A,I0,A,I0,A,3ES12.5)") &
            &                   'missing neighbor iElem: ', iElem,&
            &                   ', tID: ', tID, ', bary: ', bary
          ! Set the property of the current element to solid
          levelDesc%property( totalPos )                                       &
            &  = ibset( levelDesc%property( totalPos ), prp_solid )
        end if
      else ! elemPos <= 0, i.e. tID can not be found in levelDesc%elem%tID
        write(dbgUnit(1),*) 'Warning!!!'
        write(dbgUnit(1),*) tID, ' not found during build_treehorizontalDep'
        write(dbgUnit(1),*) 'That should not occur.'
      end if
    end do

  end subroutine tem_build_treeHorizontalDep
! ****************************************************************************** !


! ****************************************************************************** !
  !> Update the neighor arrays depending on what is given in the element stencil
  !!
  !! The array levelDesc( iLevel )%neigh( iStenci )%nghElems( 1:QQN, 1:nElems )
  !! is being filled up here
  subroutine tem_build_listHorizontalDep( iStencil, levelDesc, &
    &                                     posInSortElem, nElems, iIndex )
    ! ---------------------------------------------------------------------------
    !> Level descriptor for each level of your mesh (starting from min level).
    type(tem_levelDesc_type), intent(inout) :: levelDesc
    !> Index of your neighbor list.
    integer, intent(in) :: iStencil
    !>
    integer, intent(inout) :: iIndex
    !> number of elements
    integer, intent(in) :: nElems
    !> Positions in sorted elem%tID list
    integer, intent(in) :: posInSortElem(:)
    ! ---------------------------------------------------------------------------
    integer :: iElem, iStencilElem
    integer :: posInElem, levelPos, neighPos, neighVal
    integer(kind=long_k) :: nTreeID ! neighbor treeID value to treat
    logical :: invalid
    ! ---------------------------------------------------------------------------
    ! run over the given element list to update
    do iElem = 1, nElems
      posInElem = posInSortElem( iElem )
!write(*,*)'in build_listHorDep nVals: ', levelDesc%elem%stencil%val( posInElem )%nVals
!write(*,*)'in build_listHorDep iStencil: ', iStencil
      ! Check if there are enough stencils given in the element
      if( levelDesc%elem%stencil%val( posInElem )%nVals >= iStencil ) then
        ! Increase the index element counter, if a stencil was encountered
        iIndex = iIndex + 1
        ! Run over each entry of the stencil
        do iStencilElem = 1, levelDesc%elem%stencil%val( posInElem )%            &
          &                                                val( iStencil )%QQN
          neighPos = levelDesc%elem%stencil%val( posInElem )%val( iStencil )     &
            &                                           %tIDpos( iStencilElem )
          invalid = .false.
          if( neighPos > 0 ) then
            nTreeID = levelDesc%elem%neighID%val( posInElem )%val( neighPos )
            if( nTreeID < 0_long_k ) then
              ! If the entry in the stencil is < 0, it must be a boundary ID
              ! -> simply set
              neighVal = int(nTreeID)
            else
              ! If the entry is > 0, test if there is also an entry for the
              ! element position directly (%elem array)
              if( allocated( levelDesc%elem%stencil%val( posInElem )             &
                &                            %val( iStencil )%totalPos) ) then
                ! If there is an entry, test, if it is valid
                levelPos = levelDesc%elem%stencil%val( posInElem )               &
                  &                    %val( iStencil )%totalPos( iStencilElem )
                if( levelPos > 0 ) then
                  ! If there is a valid entry, get the position of this element
                  ! in the total list
                  neighVal = levelPos
                else
                  ! If the entry in %elem was not valid, use the tID to find
                  invalid = .true.
                end if ! levelPos > 0
              else
                ! if there was no %elem entry given at all, directly use the
                ! tID to find the element position in the total list
                invalid = .true.
              end if ! allocated totalPos
            end if ! nTreeID > 0

            ! If above anything went wrong, ...
            if( invalid ) then
              ! ... try to find it with the slowest method: just look up the
              ! tID in the total list
              neighVal = tem_treeIDinTotal( nTreeID, levelDesc )
            end if
          else ! neighPos < 0
            neighVal = neighPos
          end if ! neighPos > 0
          levelDesc%neigh(iStencil)%nghElems( iStencilElem, iIndex ) = neighVal
        end do ! stencil element
      end if ! nStencils > iStencil
    end do ! iElem

  end subroutine tem_build_listHorizontalDep
! ****************************************************************************** !


! ****************************************************************************** !
  !> Communicate with all existing process the number of requested halo elements
  !!
  !! After this routine, each process knows how many processes there are to
  !! communicate with and how many elements have to be transferred
  !!
  subroutine communicate_nElemsToTransfer( me, proc, minLevel, maxLevel )
    ! ---------------------------------------------------------------------------
    !> level range
    integer, intent(in) :: minLevel, maxLevel
    !> Process description to use.
    type(tem_comm_env_type), intent(in) :: proc
    !> the level descriptor to be filled
    type(tem_levelDesc_type), intent(inout) :: me(minlevel:maxLevel)
    ! ---------------------------------------------------------------------------
    ! the number of halos which will be sent to this mpi process for each level
    ! and for each process:
    ! first dimension is the level, second one is the mpi process
    integer, allocatable :: nHalosHis(:,:)
    integer, allocatable :: nHalosMine(:,:)
    integer, allocatable :: sources(:)
    integer, allocatable :: sourceHalos(:)
    integer :: iProc, iLevel, ierr
    integer :: nLevels, sourceProc, nProcs
    real(kind=rk) :: tStart, tEnd
    ! ---------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------
    ! Count number of halos on each level for each other process
    ! these values will be sent in the following
    nLevels = maxLevel - minLevel + 1

    tStart = mpi_wtime()

    if (use_sparse_alltoall) then

      write(logunit(4),*) 'Starting sparse(all-to-all) exchange of halo counts'

      ! Using levelwise sparse all-to-all here for now out of simplicity.
      ! If there is a need, we could first count the number of levels to
      ! exchange with each process, and then have only this in a single,
      ! sparse all-to-all.
      ! However, this is a little involved and probably requires us to
      ! employ a dynamic array for the processes...
      ! In most cases the levelwise sparse alltoall as done now is probably
      ! good enough.
      do iLevel = minLevel, maxLevel
        write(logunit(10),*) '    Level: ', iLevel
        nProcs = me(iLevel)%haloList%PartnerProc%nVals
        allocate( nHalosMine(nProcs,1) )
        nHalosMine(:,1) = me(iLevel)%haloList%halos%val(:nProcs)%nVals

        call tem_comm_init( me     = me(iLevel)%recvbuffer, &
          &                 nProcs = nProcs                 )
        me(iLevel)%recvBuffer%nElemsProc = nHalosMine(:,1)
        me(iLevel)%recvBuffer%proc &
          &  = me(iLevel)%halolist%partnerProc%val(:nProcs) - 1

        call tem_sparse_alltoall_int( targets     = me(iLevel)%recvBuffer &
          &                                                   %Proc,      &
          &                           send_buffer = nHalosMine(:,1),      &
          &                           sources     = sources,              &
          &                           recv_buffer = sourceHalos,          &
          &                           tag         = iLevel,               &
          &                           comm        = proc%comm             )

        call tem_comm_init( me     = me( iLevel )%sendbuffer, &
          &                 nProcs = size(sources)            )
        me(iLevel)%sendBuffer%nElemsProc = sourceHalos
        me(iLevel)%sendBuffer%proc = sources

        deallocate(sourceHalos)
        deallocate(sources)
        deallocate(nHalosMine)
      end do

    else

      write(logunit(4),*) 'Starting all-to-all exchange of halo counts'

      ! allocate the send buffer.
      ! Contains the number of elements for each level, which are requested
      ! from the specific process
      allocate( nHalosMine( minLevel:maxLevel,proc%comm_size ))
      ! allocate the recv buffer.
      ! Will contain the number of elements for each level, which other
      ! processes request from me. This one is filled in the mpi_alltoall call
      allocate( nHalosHis( minLevel:maxLevel, proc%comm_size ))

      ! Reset the number of requested and requesting halos
      nHalosMine(:,:) = 0
      nHalosHis(:,:)  = 0
      do iLevel = minLevel, maxLevel

        nProcs = me(iLevel)%haloList%PartnerProc%nVals
        do iProc = 1, nProcs
          sourceProc = me(iLevel)%haloList%partnerProc%val(iProc)
          nHalosMine(iLevel, sourceProc) =                  &
            &      me(iLevel)%haloList%halos%val(iProc)%nVals
        end do ! iProc

        call tem_comm_init( me(iLevel)%recvbuffer, nProcs )
        ! count valid procs and nElemsProc
        call tem_comm_count( me(iLevel)%recvBuffer, proc%comm_size, &
          &                  nHalosMine(iLevel,:)                   )

      end do   ! iLevel

      ! number of halos which are sent to a process
      ! might be different from number of halos
      ! received from the same process for that we use mpi_alltoall


      ! Exchange the number of requested treeIDs
      call mpi_alltoall( nHalosMine, nLevels, mpi_integer, &
        &                nHalosHis, nLevels, mpi_integer,  &
        &                proc%comm, ierr                   )
      do iLevel = minLevel, maxLevel
        ! now check which processs really ask for information and allocate
        ! buffer
        nProcs = count( nHalosHis( iLevel, : ) > 0 )
        call tem_comm_init( me( iLevel )%sendbuffer, nProcs )

        ! count valid procs and nElemsProc
        call tem_comm_count( me(iLevel)%sendBuffer, proc%comm_size, &
          &                  nHalosHis(iLevel,:))
      end do

      deallocate( nHalosMine )
      deallocate( nHalosHis  )

    end if

    tEnd = mpi_wtime()
    write(logunit(4),"(A)") 'Finished all-to-all exchanging number of elements'
    write(logunit(4),"(A,E12.6)") 'All to all cost: ', tEnd-tStart

    do iLevel = minLevel, maxLevel
      ! Copy all headers (information about target and source procs)
      ! to bufferFromCoarser and bufferFromFiner
      me( iLevel )%sendbufferFromCoarser = me( iLevel )%sendbuffer
      me( iLevel )%sendbufferFromFiner   = me( iLevel )%sendbuffer
      me( iLevel )%recvbufferFromCoarser = me( iLevel )%recvbuffer
      me( iLevel )%recvbufferFromFiner   = me( iLevel )%recvbuffer
    end do ! iLevel

    write(logUnit(5),*) ' Finished counting valid procs and nElemsProcs.'

  end subroutine communicate_nElemsToTransfer
! ****************************************************************************** !


! ****************************************************************************** !
  !> Inverse Communication: Communicate, which elements each process needs from
  !! me.
  !!
  !! In this routine, we send the treeIDs of the halo elements to the
  !! processes, where they are located. Later on, we fill these halos locally
  !! with information from these processes (sourceProcs). In this routine
  !! however, we now SEND information to these sourceProcs, so do not get
  !! confused here.
  !!
  subroutine request_remoteHalos( levelDesc, proc, tree, iLevel, stencil,&
    &                             pathFirst, pathLast )
    ! ---------------------------------------------------------------------------
    !> the global tree
    type(treelmesh_type), intent(in) :: tree
    !> the level descriptor to be filled
    type(tem_levelDesc_type), intent(inout) :: levelDesc(tree%global%minlevel: )
    !> Process description to use.
    type(tem_comm_env_type), intent(in) :: proc
    !> current level
    integer, intent(in) :: iLevel
    !> stencil definition
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> first and last treeID path in every process
    type(tem_path_type), intent(in) :: pathFirst(:), pathLast(:)
    ! ---------------------------------------------------------------------------
    integer :: iProc, iErr, iElem, elemPos, procPos
    integer :: haloLevel
    integer, allocatable :: rq_handle(:)
    integer, allocatable :: status(:,: )
    integer :: nCommunications, nesting
    type( grw_longArray_type ), allocatable ::  treeIDs_fromTarget(:)
    type( grw_intArray_type ),  allocatable :: nestings_fromTarget(:)
    type( grw_longArray_type ), allocatable ::  treeIDs_toSource(:)
    type( grw_intArray_type ),  allocatable :: nestings_toSource(:)
    logical :: updated
    integer, parameter :: message_flag_long = 24
    integer, parameter :: message_flag_int  = 25
    ! ---------------------------------------------------------------------------

    call tem_horizontalSpacer( fUnit = dbgUnit(1), before = 1 )
    write(dbgUnit(1),*) "Get into routine: request_remoteHalos"
    write(dbgUnit(1),*) 'Requesting remote halos on level: ', iLevel

    ! two communications: treeID and nesting
    nCommunications = 2 * (  levelDesc( iLevel )%sendbuffer%nProcs &
      &                    + levelDesc( iLevel )%recvbuffer%nProcs )
    allocate( rq_handle( nCommunications ) )
    allocate( status( mpi_status_size, nCommunications ) )
    rq_handle(:) = MPI_REQUEST_NULL

    ! Warning: Inverse Communication ! (send to source, recv from target)
    ! ---------------------------------------------------------------------
    ! I receive from target what elements are needed by them
    ! SendBuffer contains my elements required by remote targets
    allocate(  treeIDs_fromTarget( levelDesc( iLevel )%sendbuffer%nProcs ))
    allocate( nestings_fromTarget( levelDesc( iLevel )%sendbuffer%nProcs ))

    do iProc = 1, levelDesc( iLevel )%sendbuffer%nProcs

      ! Allocate the buffers
      call init( me     = treeIDs_fromTarget( iProc ),                      &
        &        length = levelDesc( iLevel )%sendbuffer%nElemsProc( iProc ))
      call init( me     = nestings_fromTarget( iProc ),                     &
        &        length = levelDesc( iLevel )%sendbuffer%nElemsProc( iProc ))

      ! Receive the element tree IDs
      call mpi_irecv( treeIDs_fromTarget( iProc )%val,                         &
        &             treeIDs_fromTarget( iProc )%ContainerSize,               &
        &             mpi_integer8,                                            &
        &             levelDesc( iLevel )%sendbuffer%proc(iProc),              &
        &             message_flag_long,                                       &
        &             proc%comm,                                               &
        &             rq_handle( iProc),                                       &
        &             iErr )
      ! Receive the element nestings
      call mpi_irecv( nestings_fromTarget( iProc )%val,                        &
        &             nestings_fromTarget( iProc )%ContainerSize,              &
        &             mpi_integer,                                             &
        &             levelDesc( iLevel )%sendbuffer%proc(iProc),              &
        &             message_flag_int,                                        &
        &             proc%comm,                                               &
        &             rq_handle( iProc+nCommunications/2 ),                    &
        &             iErr  )
      ! Update the number of elements inside the growing array of the recv
      ! buffer
       treeIDs_fromTarget( iProc )%nVals = treeIDs_fromTarget( iProc )%ContainerSize
      nestings_fromTarget( iProc )%nVals = treeIDs_fromTarget( iProc )%ContainerSize
    end do ! iProc
    ! ---------------------------------------------------------------------

    ! ---------------------------------------------------------------------
    ! I send to source what elements I need from them
    allocate(  treeIDs_toSource( levelDesc( iLevel )%recvbuffer%nProcs ) )
    allocate( nestings_toSource( levelDesc( iLevel )%recvbuffer%nProcs ) )

    do iProc = 1, levelDesc( iLevel )%recvbuffer%nProcs

      ! Get the position of the process in the dynamic halos list (might be
      ! unordered)
      procPos = PositionOfVal( me  = levelDesc(iLevel)%haloList%PartnerProc, &
        &                      val = levelDesc(iLevel)%recvbuffer            &
        &                                             %proc(iProc) + 1       )

      call init(  me     = treeIDs_toSource( iProc ),                        &
        &         length = levelDesc( iLevel )%recvbuffer%nElemsProc( iProc ))
      call init(  me     = nestings_toSource( iProc ),                       &
        &         length = levelDesc( iLevel )%recvbuffer%nElemsProc( iProc ))

      ! Collect the halo treeIDs into send buffers for all the processes to
      ! request from
      do iElem = 1, levelDesc( iLevel )%recvbuffer%nElemsProc( iProc )
        elemPos = levelDesc( iLevel )%haloList%halos%val(procPos)%val(iElem)
        call append( me  = treeIDs_toSource( iProc ),                  &
          &          val = levelDesc( iLevel )%elem%tID%val( elemPos ) )
        call append( me  = nestings_toSource( iProc ),                         &
          &          val = levelDesc( iLevel )%elem%haloNesting%val( elemPos ) )
      enddo

      ! Send treeIDs
      call mpi_isend(                                                           &
        &        treeIDs_toSource(iProc)%val,                                   &
        &        treeIDs_toSource(iProc)%nVals,                                 &
        &        mpi_integer8,                                                  &
        &        levelDesc( iLevel )%recvbuffer%proc( iProc ),                  &
        &        message_flag_long,                                             &
        &        proc%comm,                                                     &
        &        rq_handle( iProc + levelDesc( iLevel )%sendbuffer%nProcs),     &
        &        ierr  )
      ! Send nesting
      call mpi_isend(                                                           &
        &        nestings_toSource(iProc)%val,                                  &
        &        nestings_toSource(iProc)%nVals,                                &
        &        mpi_integer,                                                   &
        &        levelDesc( iLevel )%recvbuffer%proc( iProc ),                  &
        &        message_flag_int,                                              &
        &        proc%comm,                                                     &
        &        rq_handle( iProc + levelDesc( iLevel )%sendbuffer%nProcs       &
        &                         + nCommunications/2),                         &
        &        ierr  )

    end do ! iProc

    call mpi_waitall( nCommunications, rq_handle, status, ierr)
    write(logUnit(5),*) '  Received requested halo elements successfully.'
    deallocate(  treeIDs_toSource )
    deallocate( nestings_toSource )
    deallocate(         rq_handle )
    deallocate(         status    )
    ! Requested halo elements were received.
    ! ---------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------
    ! Now identify the requested halos.
    if( allocated( levelDesc( iLevel )%sendbufferFromCoarser%elemPos))         &
      &         deallocate( levelDesc( iLevel )%sendbufferFromCoarser%elemPos )
    allocate( levelDesc( iLevel )%sendbufferFromCoarser%elemPos(               &
      &                                levelDesc( iLevel )%sendbuffer%nProcs ))
    if( allocated( levelDesc( iLevel )%sendbufferFromFiner%elemPos ))          &
      &           deallocate( levelDesc( iLevel )%sendbufferFromFiner%elemPos )
    allocate( levelDesc( iLevel )%sendbufferFromFiner%elemPos(                 &
      &                                levelDesc( iLevel )%sendbuffer%nProcs ))

    ! Add elements of received buffers to elem
    do iProc = 1, levelDesc( iLevel )%sendbuffer%nProcs
      ! Allocate the buffer for the element position indices
      call init( me = levelDesc( iLevel )%sendbuffer%elemPos(iProc) )
      call init( me = levelDesc( iLevel )%sendbufferFromCoarser%elemPos(iProc) )
      call init( me = levelDesc( iLevel )%sendbufferFromFiner%elemPos(iProc) )

      do iElem = 1, treeIDs_fromTarget( iProc )%nVals
        nesting = nestings_fromTarget( iProc )%val( iElem )
        ! identify requested halo treeID in local process
        call identify_halo( haloTreeID = treeIDs_fromTarget( iProc )%val(iElem), &
          &                 elemPos    = elemPos,                          &
          &                 halolevel  = haloLevel,                        &
          &                 levelDesc  = levelDesc,                        &
          &                 nesting    = nesting,                          &
          &                 updated    = updated,                          &
          &                 tree       = tree,                             &
          &                 minLevel   = tree%global%minLevel,             &
          &                 stencil    = stencil )
        if( elemPos > 0 ) then
          ! if requested halo is ghostFromCoarser then find stencil neighbors of
          ! this halo element
          if ( (nesting < nestingLimit)                           &
            &  .and. (levelDesc( iLevel )%elem%eType%val(elemPos) &
            &         == eT_ghostFromCoarser)                     ) then
            ! identify all the compute neighbors of the current element
            call identify_stencilNeigh( iElem          = elemPos,    &
              &                         iLevel         = iLevel,     &
              &                         tree           = tree,       &
              &                         iStencil       = 1,          &
              &                         pathFirst      = pathFirst,  &
              &                         pathLast       = pathLast,   &
              &                         levelDesc      = levelDesc,  &
              &                         proc           = proc,       &
              &                         stencil        = stencil,    &
              &                         nesting        = nesting + 1 )
          end if

          ! if requested halo element haloNesting < found halo element (elemPos)
          ! haloNesting
          if ( nestings_fromTarget(iProc)%val(iElem)                  &
            &  < levelDesc( haloLevel )%elem%haloNesting%val(elemPos) ) then
            levelDesc(haloLevel)%elem%needsUpdate%val(elemPos) = .true.
            levelDesc(haloLevel)                                             &
              &  %elem                                                       &
              &  %haloNesting                                                &
              &  %val(elemPos) = min( nestings_fromTarget(iProc)%val(iElem), &
              &                       levelDesc(haloLevel)                   &
              &                         %elem                                &
              &                         %haloNesting                         &
              &                         %val(elemPos)                        )
          end if

          ! only add, if the element was added locally
          select case( levelDesc(iLevel)%elem%eType%val(elemPos) )
          ! Depending on the type of the element, add to the
          ! regular buffer, bufferFromCoarser, bufferFromFiner
          case( eT_fluid )
            call append( me  = levelDesc( iLevel )%sendbuffer%elemPos( iProc ),&
              &          val = elemPos )
          case( eT_ghostFromCoarser )
            call append( me  = levelDesc( iLevel )%sendbufferFromCoarser       &
              &                                             %elemPos( iProc ), &
              &          val = elemPos )

            ! for ghostFromCoarser determine neighbors of coarser element
            if( levelDesc( haloLevel )%elem%haloNesting%val( elemPos )         &
              &                                        < nestingLimit ) then
              call create_allParentNeighbors(                                  &
                &      targetID       = levelDesc(iLevel)%elem%tID%val( elemPos ),&
                &      level          = iLevel,                                &
                &      tree           = tree,                                  &
                &      stencil        = stencil,                        &
                &      levelDesc      = levelDesc,                             &
                &      pathFirst      = pathFirst,                             &
                &      pathLast       = pathLast,                              &
                &      proc           = proc )
            end if

          case( eT_ghostFromFiner )
            call append( me  = levelDesc( iLevel )%sendbufferFromFiner         &
              &                %elemPos( iProc ),                              &
              &          val = elemPos )
          case( eT_distributedGhostFromFiner)
             write(logUnit(1),*)' Found distributed ghost From Finer in '//    &
               &                'request remote Halos'
             write(logUnit(1),*)' This case should not occur!'
             call tem_abort()
          end select
        end if ! elemPos > 0
      end do ! ielem recvbuffer

      levelDesc( iLevel )%sendbuffer%nElemsProc( iProc )                       &
        &   = levelDesc( iLevel )%sendbuffer%elemPos( iProc )%nVals
      levelDesc( iLevel )%sendbufferFromCoarser%nElemsProc( iProc )            &
        &   = levelDesc( iLevel )%sendbufferFromCoarser%elemPos( iProc )%nVals
      levelDesc( iLevel )%sendbufferFromFiner%nElemsProc( iProc )              &
        &   = levelDesc( iLevel )%sendbufferFromFiner%elemPos( iProc )%nVals

      ! destroy temp variables
      call destroy( me =  treeIDs_fromTarget( iProc ) )
      call destroy( me = nestings_fromTarget( iProc ) )
    end do ! iProc
    deallocate(  treeIDs_fromTarget )
    deallocate( nestings_fromTarget )
    ! Now each Process knows, which elements to send to others
    ! ---------------------------------------------------------------------
    write(logUnit(5),*) 'Finished requesting remote halos'

    write(dbgUnit(1),*) "Leave  routine: request_remoteHalos"
    call tem_horizontalSpacer( fUnit = dbgUnit(1), after = 1 )

  end subroutine request_remoteHalos
! ****************************************************************************** !


! ****************************************************************************** !
  !> Report the actually existing elements, which were requested as halos
  !! from remote
  !!
  subroutine return_haloCounts( sendbuffer, recvbuffer, comm )
    ! ---------------------------------------------------------------------------
    !> send buffer
    type(tem_communication_type), intent(in) :: sendbuffer
    !> recv buffer
    type(tem_communication_type), intent(inout) :: recvbuffer
    !> MPI communicator
    integer, intent(in) :: comm
    ! ---------------------------------------------------------------------------
    integer :: iProc ! process iterator
    integer :: iErr  ! error flag
    integer, parameter :: message_tag = 23 ! arbitrary message flag
    integer, allocatable :: rq_handle(:)
    integer, allocatable :: status(:,:)
    integer :: nCommunications ! amount of communciations
    ! ---------------------------------------------------------------------------

    nCommunications = sendbuffer%nProcs + recvbuffer%nProcs
    allocate( rq_handle( nCommunications ) )
    allocate( status( mpi_status_size, nCommunications ) )
    ! ---------------------------------------------------------------------
    !! first report the number of actual halos
    rq_handle(:) = MPI_REQUEST_NULL

    do iProc = 1, recvbuffer%nProcs
      ! Receive the number of actual elements
      call mpi_irecv( recvbuffer%nElemsProc( iProc ),                          &
       &              1,                                                   &
       &              mpi_integer,                                             &
       &              recvbuffer%proc(iProc),                                  &
       &              message_tag,                                             &
       &              comm,                                               &
       &              rq_handle(iProc),                                        &
       &              iErr  )
    end do ! iProc

    ! I send the number of actual existing halo elements to the requester.
    do iProc = 1, sendbuffer%nProcs
      call mpi_isend( sendbuffer%nElemsProc( iProc ),                          &
       &              1,                                                   &
       &              mpi_integer,                                             &
       &              sendbuffer%proc( iProc ),                                &
       &              message_tag,                                             &
       &              comm,                                               &
       &              rq_handle(iProc + recvbuffer%nProcs),                    &
       &              iErr )
    end do ! iProc

    call mpi_waitall( nCommunications, rq_handle, status, iErr)
    ! Now we know the number of how actually existent cells, to receive
    ! allocate the recv buffers and communicate
    ! ---------------------------------------------------------------------

  end subroutine return_haloCounts
! ****************************************************************************** !


! ****************************************************************************** !
  !> Check if additional communications have to be performed
  !!
  subroutine check_additionalComm( levelDesc, proc, doAdditional, minlevel )
    ! ---------------------------------------------------------------------------
    !> minlevel in tree
    integer, intent(in) :: minlevel
    !> level descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc(minlevel:)
    !> Process description to use.
    type(tem_comm_env_type), intent(in) :: proc
    !> do addtional steps to identify neighbors of elems in require list
    ! and do overall communication again, if require%nVals > 0
    logical, intent(out) :: doAdditional
    ! ---------------------------------------------------------------------------
    integer :: iError, iLevel
    logical :: do_local
    ! ---------------------------------------------------------------------------
    doAdditional = .false.
    do_local = .false.
    do iLevel = minlevel, ubound( levelDesc, 1 )
      do_local = ( do_local .or. (levelDesc( iLevel )%require%nVals > 0))
    end do
    ! JUROPA work-around for crash in the mpi_allreduce
    call mpi_barrier( proc%comm, iError )
    call mpi_allreduce( do_local, doAdditional, 1, mpi_logical,                &
      &                 mpi_LOR, proc%comm, iError )
    if( doAdditional ) then
      write(dbgUnit(1),*)'Perform additional exchanges for required neighbors'
    end if

  end subroutine check_additionalComm
! ****************************************************************************** !


! ****************************************************************************** !
  !> Report the actually existing elements, which were requested as halos
  !! from remote
  !!
  subroutine redefine_halos( levelDesc, sendbuffer, recvbuffer, proc,           &
    &                        commPattern, computeStencil )
    ! ---------------------------------------------------------------------------
    !> the level descriptor of specific level
    type(tem_levelDesc_type), intent(inout) :: levelDesc
    !> send and receive communication buffer type
    type(tem_communication_type), intent(inout) :: sendbuffer, recvbuffer
    !> Process description to use.
    type(tem_comm_env_type), intent(in) :: proc
    !> communication pattern
    type(tem_commPattern_type)    :: commPattern
    !> array of all stencils used in the simulation
    type(tem_stencilHeader_type)  :: computeStencil(:)
    ! ---------------------------------------------------------------------------
    integer :: iProc, tPos, iElem, nElems
    integer :: elemCount
    integer(kind=long_k), allocatable :: recvBuf(:)
    integer :: message_flag
    logical :: wasAdded
    ! ---------------------------------------------------------------------------
    message_flag = 3

    allocate( recvBuf(sum(recvbuffer%nElemsProc)) )

    elemCount = 0
    do iProc = 1, recvbuffer%nProcs
      nElems = recvbuffer%nElemsProc(iProc)
      call init(  me     = recvBuffer%elemPos( iProc ),                        &
        &         length = nElems )
      do iElem=1,nElems
        elemCount = elemCount + 1
        call append( me  = recvbuffer%elemPos( iProc ),                        &
          &          val = elemCount )
      end do

      call commPattern%initBuf_long( me    = recvBuffer%buf_long( iProc ),     &
        &                            pos   = recvBuffer%elemPos( iProc )%val,  &
        &                            nVals = recvBuffer%nElemsProc( iProc ))
    end do

    if (tem_logging_isActive(main_debug%logger, 6)) then
      write(dbgUnit(0),*) '  SEND BUFFER '
      call tem_comm_dumpType( me = sendbuffer,    &
        &                     nUnit = dbgUnit(6) )
      write(dbgUnit(0),*) ' '
      write(dbgUnit(0),*) '  RECV BUFFER '
      call tem_comm_dumpType( me = recvbuffer,    &
        &                     nUnit = dbgUnit(6) )
    end if

    ! ---------------------------------------------------------------------
    ! initialize send buffer
    do iProc = 1, sendbuffer%nProcs
      call commPattern%initBuf_long( me    = sendBuffer%buf_long( iProc ),     &
        &                            pos   = sendBuffer%elemPos( iProc )%val,  &
        &                            nVals = sendBuffer%nElemsProc( iProc ))
    end do

    call commPattern%exchange_long( send         = sendbuffer,                 &
      &                             recv         = recvbuffer,                 &
      &                             state        = recvBuf,                    &
      &                             message_flag = message_flag,               &
      &                             comm         = proc%comm,                  &
      &                             send_state   = levelDesc%elem%tID%val )

    ! Now I received all the halos, which are really existing
    ! and required on the remote processes
    ! ---------------------------------------------------------------------

    ! ---------------------------------------------------------------------
    ! Add all received elements
    ! Set the eType accordingly or leave as nonExisting if it was not received
    elemCount = 0
    do iProc = 1, recvbuffer%nProcs
      nElems = recvBuffer%nElemsProc( iProc )
      call commPattern%finBuf_long(me = recvbuffer%buf_long(iProc))
      call init(  me     = recvBuffer%elemPos( iProc ),                        &
        &         length = nElems )
      do iElem = 1, nElems
        elemCount = elemCount + 1
        call append( me         = levelDesc%elem,             &
          &          tID        = recvBuf( elemCount ),       &
          &          sourceProc = recvbuffer%proc( iProc )+1, &
          &          eType      = eT_halo,                    &
          &          pos        = tPos,                       &
          &          wasAdded   = wasAdded )
        if ( .not. wasAdded ) then
          ! If it was not added, it was there before.
          ! Then simply set the element type to halo (was initialized with
          ! eT_nonExistingElem above )
          if ( levelDesc%elem%eType%val( tPos ) == eT_nonExisting ) then
            call changeType( levelDesc%elem, tPos, eT_halo )
          end if
        end if
        call append( me  = recvbuffer%elemPos( iProc ),                        &
          &          val = tPos )
      end do
      call commPattern%initBuf_long( me    = recvbuffer%buf_long( iProc ),     &
        &                            pos   = recvbuffer%elemPos( iProc )%val,  &
        &                            nVals = recvBuffer%nElemsProc( iProc ))
    end do

    deallocate(recvBuf)

    ! ---------------------------------------------------------------------
    ! communicate the properties of the halos
    call commPattern%exchange_long(                                            &
      &                         send         = sendbuffer,                     &
      &                         recv         = recvbuffer,                     &
      &                         state        = levelDesc%elem%property%val(:), &
      &                         message_flag = 5,                              &
      &                         comm         = proc%comm )

    ! communicate the stencil neighbors
    ! the stencil which includes the boundary information
    if ( recvbuffer%nProcs > 0 .or. sendbuffer%nProcs > 0 ) then
      ! communicate only fluid stencil which is defined for all
      ! elements in elem list including halo elements.
      ! In our case, fluid stencil is 1st Stencil so setting iStencil=1
      call tem_stencil_communicate( send           = sendbuffer,             &
        &                           recv           = recvbuffer,             &
        &                           elem           = levelDesc%elem,         &
        &                           proc           = proc,                   &
        &                           commPattern    = commPattern,            &
        &                           computeStencil = computeStencil(1),      &
        &                           iStencil       = 1 )
    end if

    ! Deallocate long type buffer
    do iProc = 1, sendBuffer%nProcs
      call commPattern%finbuf_long( sendBuffer%buf_long(iProc) )
    end do
    do iProc = 1, recvBuffer%nProcs
      call commPattern%finbuf_long( recvBuffer%buf_long(iProc) )
    end do

    write(logUnit(5),*) '  Done communicating stencil buffer'
    ! ---------------------------------------------------------------------

  end subroutine redefine_halos
! ****************************************************************************** !


! ****************************************************************************** !
  !> Communicate the complete stencil
  !!
  !! Currently, this assumes same stencils for all participating elements
  !!
  subroutine tem_stencil_communicate( send, recv, elem, computeStencil,        &
    &                                 proc, commPattern, iStencil )
    ! ---------------------------------------------------------------------------
    !> send and recv communication buffers
    type( tem_communication_type ), intent(inout) :: send, recv
    !> communication pattern
    type(tem_commPattern_type), intent(in)     :: commPattern
    !> levelDesc element list
    type(tem_element_type), intent(inout)  :: elem
    !> Process description to use.
    type(tem_comm_env_type), intent(in) :: proc
    !> array of all stencils used in the simulation
    type(tem_stencilHeader_type), intent(in) :: computeStencil
    !> amount of values to communicate
    integer, intent(in) :: iStencil
    ! ---------------------------------------------------------------------------
    integer :: iStencilElem, iElem, nElems
    integer(kind=long_k), allocatable :: buffer(:)
    type(tem_stencilElement_type ) :: tStencil
    integer :: addedPos, iProc, elemPos, neighPos, stencilPos
    logical :: wasAdded
    ! ---------------------------------------------------------------------------
    ! Number of local elements in element list
    nElems = elem%tID%nVals
    write(dbgUnit(10),*)'in stencil communicate, nElems: ', nElems

    ! buffer which sends stencil neighbors to remove process
    allocate( buffer( nElems ))

    ! Set the temporary (empty) stencil
    call init( me        = tStencil,                                           &
      &        QQN       = computeStencil%QQN ,                                &
      &        headerPos = iStencil )

    ! add stencil headerPos for all halo elements depends on requested stencil
    ! communication
    do iProc = 1, recv%nProcs
      do iElem = 1, recv%nElemsProc( iProc )
        elemPos = recv%elemPos( iProc )%val( iElem )
        if( elem%stencil%val( elemPos )%nVals > 0 ) then
          ! @todo: SZ: this does not make sense since all elements have a stencil
          !            elem%stencil%val( elemPos )%nVals .gt. 0 for all elements
          ! Find the stencil corresponding to the current one ( iStencil )
          stencilPos = tem_stencil_getHeaderPos(                               &
            &                               me  = elem%stencil%val( elemPos ), &
            &                               val = iStencil )
          ! if none was found, use the most appropriate one (headerPos = 0)
          ! this assigns the first stencil to the halo element
          ! why not assign it right away to be 1?
          if( stencilPos <= 0 )then
            stencilPos = tem_stencil_getHeaderPos(                             &
              &                             me  = elem%stencil%val( elemPos ), &
              &                             val = 0 )
            elem%stencil%val( elemPos )%val( stencilPos )%headerPos = iStencil
          end if
        else
          ! if stencil is not initialized for this element append tStencil
          call append( me = elem%stencil%val( elemPos ), val = tStencil )
        end if
      end do
    end do

    ! --------------------------------------------------------------------------
    ! Do the communication for each stencil tID entry
    do iStencilElem = 1, computeStencil%QQN
      ! Fill the send buffer
      do iProc = 1, send%nProcs
        do iElem = 1, send%nElemsProc( iProc )
          elemPos = send%elemPos( iProc )%val( iElem )
          neighPos = elem%stencil%val( elemPos )%val( iStencil )%              &
            &                                            tIDpos( iStencilElem )
          buffer( elemPos ) = elem%neighID%val( elemPos )%val( neighPos )
        end do ! iElem send
      end do ! iProc send

      ! Do the exchange between all participating processes
      call commPattern%exchange_long( send         = send,                     &
        &                             recv         = recv,                     &
        &                             state        = buffer,                   &
        &                             message_flag = 0,                        &
        &                             comm         = proc%comm )

      ! Read from the buffer
      do iProc = 1, recv%nProcs
        do iElem = 1, recv%nElemsProc( iProc )
          elemPos = recv%elemPos( iProc )%val( iElem )
          ! Find the stencil corresponding to the current one (iStencil )
          stencilPos = tem_stencil_getHeaderPos(                               &
            &                               me  = elem%stencil%val( elemPos ), &
            &                               val = iStencil )

          if( stencilPos > 0 ) then
            ! Only append elements, if there is only one stencil for the element
            ! Append the treeID to the element's neighIDs
            call append( me       = elem%neighID%val( elemPos ),               &
              &          val      = buffer( elemPos),                          &
              &          pos      = addedPos,                                  &
              &          wasAdded = wasAdded )
            ! ... and store the appended position to the stencil
            elem%stencil%val( elemPos )%val( stencilPos )%                     &
              &                               tIDpos( iStencilElem ) = addedPos
          end if
        end do ! iElem recv
      end do ! iProc recv

    end do ! iStencilElem

  end subroutine tem_stencil_communicate
! ****************************************************************************** !


! ****************************************************************************** !
  !> Update the found dependencies, which were built for non-ordered lists
  !! Out of fluid, ghost and halo lists, the totalList is constructed in an
  !! ordered fashion. The element order as in the TotalList is later passed on
  !! to the solver.
  !!
  subroutine update_elemPosToTotalPos( levelDesc, levelPointer, tree,          &
    &                                  computeStencil )
    ! ---------------------------------------------------------------------------
    !> the global tree
    type(treelmesh_type), intent(in)            :: tree
    !> Pointer from original treeID list to level wise fluid list
    integer, intent(in)                         :: levelPointer(:)
    !> the level descriptor to be filled
    type(tem_levelDesc_type), intent(inout)     :: levelDesc(tree%global%minlevel:)
    !> array of all stencils used in the simulation
    type(tem_stencilHeader_type), intent(inout) :: computeStencil(:)
    ! ---------------------------------------------------------------------------
    integer :: iLevel, iElem, iStencil, iStencilElem, pos
    integer :: nUpdates, iError
    integer(kind=long_k) :: nTreeID
    logical :: needUpdate
    ! ---------------------------------------------------------------------------
    write(logUnit(3),*) 'Updating dependencies ...'

    ! Update the element positions from the position in the original
    ! elem list to the total list
    do iLevel = tree%global%minlevel, tree%global%maxlevel
      call update_buffer_elemPos( buffer    = levelDesc( iLevel )%sendbuffer,  &
        &                         levelDesc = levelDesc( iLevel ),             &
        &                         iError    = iError )
      if (iError > 0) &
        write(dbgUnit(1),*)  'error in sendbuffer ',iError,' l',iLevel

      call update_buffer_elemPos( buffer    = levelDesc( iLevel )%recvbuffer,  &
        &                         levelDesc = levelDesc( iLevel ),             &
        &                         iError    = iError )
      if( iError > 0 ) &
        write(dbgUnit(1),*) 'error in recvbuffer ',iError,' l',iLevel

      call update_buffer_elemPos(                                              &
        &            buffer    = levelDesc( iLevel )%sendbufferFromCoarser,    &
        &            levelDesc = levelDesc( iLevel ),                          &
        &            iError    = iError )
      if( iError > 0 ) &
        write(dbgUnit(1),*) 'error in sendFCbuff ',iError,' l',iLevel

      call update_buffer_elemPos(                                              &
        &            buffer    = levelDesc( iLevel )%recvbufferFromCoarser,    &
        &            levelDesc = levelDesc( iLevel ),                          &
        &            iError    = iError )
      if( iError > 0 ) &
        write(dbgUnit(1),*) 'error in recvFCbuff ',iError,' l',iLevel

      call update_buffer_elemPos(                                              &
        &            buffer    = levelDesc( iLevel )%sendbufferFromFiner,      &
        &            levelDesc = levelDesc( iLevel ),                          &
        &            iError    = iError )
      if( iError > 0 ) &
        write(dbgUnit(1),*) 'error in sendFFbuff ',iError,' l',iLevel

      call update_buffer_elemPos(                                              &
        &            buffer    = levelDesc( iLevel )%recvbufferFromFiner,      &
        &            levelDesc = levelDesc( iLevel ),                          &
        &            iError    = iError )
      if( iError > 0 ) &
        write(dbgUnit(1),*) 'error in recvFFbuff ',iError,' l',iLevel
    end do ! iLevel

    write(logUnit(4),*) 'Updating the stencil entries'
    ! update stencil elem entry from original treeID list to levelwise fluid
    ! list
    do iStencil = 2, size( computeStencil )
      if( .not. computeStencil( iStencil )%useAll ) then
        do iElem = 1, computeStencil( iStencil )%nElems
          computeStencil( iStencil )%elem%val( iElem ) =                       &
            &   levelPointer( computeStencil( iStencil )%elem%val( iElem ))
        end do
      end if
    end do

    ! Update the totalPos in the element-stencils to the sorted total list
    write(logUnit(4),*) ' Updating entries in the element stencils...'
    nUpdates =  0
    do iLevel = tree%global%minlevel, tree%global%maxlevel
      do iElem = 1, levelDesc(iLevel)%elem%tID%nVals
        do iStencil = 1, levelDesc( iLevel )%elem%stencil%val( iElem )%nVals
          do iStencilElem = 1, levelDesc( iLevel )%elem%stencil%               &
            &                                  val( iElem )%val( iStencil )%QQN
            needUpdate = .false.
            pos = levelDesc( iLevel )%elem%stencil%val( iElem )%               &
              &                        val( iStencil )%totalPos( iStencilElem )
            if( pos > 0 ) then
              nTreeID = levelDesc( iLevel )%elem%tID%val( pos )
              if( pos > levelDesc( iLevel )%nElems ) needUpdate = .true.
              if( .not. needUpdate ) then
                if( levelDesc( iLevel )%total( pos ) /= nTreeID )            &
                  &                                         needUpdate = .true.
              end if
              if( needUpdate ) then
                ! Totalpos needs update!
                nUpdates = nUpdates + 1
                levelDesc( iLevel )%elem%stencil%val( iElem )%                 &
                  &                   val( iStencil )%totalPos( iStencilElem ) &
                  &  = tem_treeIDinTotal( nTreeID, levelDesc( iLevel ))
              end if
            end if ! pos > 0
          end do ! iStencilElem
        end do ! iStencil
      end do ! iElem
    end do ! iLevel
    write(logUnit(8),"(A,I0)") '   Updated nEntries: ', nUpdates
    write(logUnit(8),"(A   )") ' Finished updating dependencies'

   end subroutine update_elemPosToTotalPos
! ****************************************************************************** !


! ****************************************************************************** !
  !> Update the position of the elements inside the buffers from the original
  !! tID list to the later totalList
  !!
  subroutine update_buffer_elemPos( buffer, levelDesc, iError )
    ! ---------------------------------------------------------------------------
    !> communication buffer
    type(tem_communication_type), intent(inout) :: buffer
    !> levelDesc to be used
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> return encountered error
    integer, intent(out) :: iError
    ! ---------------------------------------------------------------------------
    integer :: iElem, iProc, elemPos, iVal
    ! ---------------------------------------------------------------------------
    iError = 0
    do iProc = 1, buffer%nProcs
      do iElem = 1, buffer%elemPos( iProc )%nVals
        elemPos = buffer%elemPos( iProc )%val( iElem )
        ! replace buffer elemPos from levelDesc%elem%tID list to
        ! levelDesc%total list
        buffer%elemPos( iProc )%val( iElem ) =                                 &
          &   tem_treeIDinTotal( tID       = levelDesc%elem%tID%val( elemPos ),&
          &                      eType     = levelDesc%elem%eType%val(elemPos),&
          &                      levelDesc = levelDesc )

        if( buffer%elemPos( iProc )%val( iElem ) < 1 ) then
          write(dbgUnit(2),*) ' tID ', levelDesc%elem%tID%val(elemPos),        &
            &                ' in buffer not found'
          iError = iProc*10000000+iElem
          do iVal = 1, levelDesc%elem%tID%nVals
            if( levelDesc%elem%tID%val( iVal ) == levelDesc%elem%tID%          &
              &                                            val( elemPos )) then
              write(dbgUnit(2),*)'found at', iVal
            end if
          end do
          do iVal = 1, size(levelDesc%total)
            if( levelDesc%total( iVal ) == levelDesc%elem%tID%val( elemPos )) &
              & write(dbgUnit(2),*) 'found in total at', iVal
          end do
        end if ! if elemPos < 1

      end do !iElem
    end do !iProc

  end subroutine update_buffer_elemPos
! ****************************************************************************** !


! ****************************************************************************** !
  !> add here the dependency for interpolation between the levels
  !! For each target cell, there are one or more source cells.
  !! The source cell can be of type fluid, ghost or halo.
  !! We save the type to update the correct element position later on,
  !! when the lists have been assembled.
  !!
  subroutine appendGhostDependency( sourcePos, sourceLevel, tgtDep )
    ! ---------------------------------------------------------------------------
    !> position of the source cell in total list
    integer, intent(in) :: sourcePos
    !> level of the source ghost cell
    integer, intent(in) :: sourceLevel
    !> dependent source elements for this target
    type(depSource_type), intent(inout) :: tgtDep
    ! ---------------------------------------------------------------------------

    tgtDep%dependencyLevel = sourceLevel

    ! append the new entry to the dependency list at the end+1
    call append( me     = tgtDep%elem, &
      &          val    = sourcePos,   &
      &          length = 8            )

  end subroutine appendGhostDependency
! ****************************************************************************** !


! ****************************************************************************** !
  !> Convert a non-zero stencil direction {-1,1} to child coordinate {0,1}
  !!
  elemental function stencilToChild( stencilCoord ) result( childCoord)
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: stencilCoord
    integer :: childCoord
    ! ---------------------------------------------------------------------------

    childCoord = (stencilCoord + 1)/2

  end function stencilToChild
! ****************************************************************************** !


! ****************************************************************************** !
  !> Convert a child coordinate {0,1} to non-zero stencil direction {-1,1}
  !!
  function childToStencil( childCoord ) result( stencilCoord)
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: childCoord
    integer :: stencilCoord
    ! ---------------------------------------------------------------------------

    stencilCoord = (childCoord*2) -1

  end function childToStencil
! ****************************************************************************** !


! ****************************************************************************** !
  !> Update the link into a given direction, based on the childs neighbor
  !! relations.
  !! Define here the trumping rule to decide, which of the neighbors or
  !! boundarie is taken for the ghostFromFiner element
  !!
  subroutine update_childNeighborID( neighID, childCoord, childPos,            &
    &                                iStencilElem, elem, iStencil)
    ! ---------------------------------------------------------------------------
    !> neighID for coarser
    integer(kind=long_k),intent(inout) :: neighID
    !> child coordinates
    integer, intent(in) :: childCoord(4)
    !> position of childIds in levelDesc elem tID list
    integer, intent(in) :: childPos(8)
    !> current stencil direction
    integer, intent(in) :: iStencilElem
    !>
    integer, intent(in) :: iStencil
    !>
    type(tem_element_type), intent(in) :: elem
    ! ---------------------------------------------------------------------------
    integer :: childID, posInElem, posInNeighID
    integer(kind=long_k) :: tNeighID
    ! ---------------------------------------------------------------------------

    childID = int(tem_idOfCoord( childCoord ))
    ! childPos holds the positions of the child treeIDs in levelDesc%elem
    ! only if the current child element really exists, get information from its
    ! stencil
    posInElem = childPos(childID)
    if ( posInElem > 0 ) then

      ! if child at current position exists,
      ! take the max tID neighbor parent or the
      ! min boundaryID -> trumping rule (bIDs are stored as negative integers)
      posInNeighID = elem%stencil%val( posInElem )%val(iStencil)%tIDpos( iStencilElem )
      tNeighID = elem%neighID%val( posInElem )%val( posInNeighID )
      if ( tNeighID > 0_long_k ) then
        ! neighbor exist, calculate its parent
        neighID = max( neighID, tem_parentOF(tNeighID) )
      else ! tNeighID <= 0
        ! neighbor is a BC ID
        neighID = max( neighID, tNeighID )
      end if

    end if

  end subroutine update_childNeighborID
! ****************************************************************************** !


! ****************************************************************************** !
  !> write out the complete list of elements of a given level
  subroutine tem_elemList_dump( me, nUnit, string, stencil, compact )
    ! ---------------------------------------------------------------------------
    type( tem_element_type ), intent(in) :: me
    integer, intent(in) :: nUnit
    character(len=*), intent(in) :: string
    logical, intent(in), optional :: stencil
    logical, intent(in), optional :: compact
    ! ---------------------------------------------------------------------------
    integer :: iElem
    ! ---------------------------------------------------------------------------
    write(nUnit, *) '========================================================='
    write(nUnit, *) '= ', trim(string), ', size: ', me%tID%nVals
    write(nUnit, *) '========================================================='
    do iElem = 1, me%tID%nVals
      call tem_element_dump( me         = me,                 &
        &                    elemPos    = iElem,              &
        &                    nUnit      = nUnit,              &
        &                    header     = (mod(iElem,32)==1), &
        &                    compact    = compact,            &
        &                    stencil    = stencil             )
    end do
    write(nUnit, *) '==  DONE! ', trim(string)
    write(nUnit, *) '========================================================='

  end subroutine tem_elemList_dump
! ****************************************************************************** !


! ****************************************************************************** !
  !> write out the complete list of elements of a given level
  !!
  subroutine tem_require_dump( me, nUnit, string )
    ! ---------------------------------------------------------------------------
    type( dyn_longArray_type ), intent(in) :: me
    integer, intent(in) :: nUnit
    character(len=*), intent(in) :: string
    ! ---------------------------------------------------------------------------
    integer :: iElem
    integer :: nEntries, iLoop, iEntry
    character(len=pathLen) :: buffer
    ! ---------------------------------------------------------------------------
    nEntries = 15
    write(nUnit, *) '========================================================='
    write(nUnit, *) '==      REQUIRE  ', trim(string)
    write(nUnit, *) '========================================================='
    iElem = 0
    do iLoop = 1, ceiling( real(me%nVals)/real(nEntries) )
      buffer = ''
      do iEntry = 1, nEntries
        iElem = iElem + 1
        if( iElem .ge. me%nVals ) exit
        write(buffer, '(a, i10)') trim(buffer), me%val( iElem )
      end do
      write(nUnit,*) trim(buffer)
      if( iElem .ge. me%nVals ) exit
    end do
    write(nUnit, *) '==  DONE         ', trim(string)
    write(nUnit, *) '========================================================='

  end subroutine tem_require_dump
! ****************************************************************************** !


! ****************************************************************************** !
  !> output all level-wise treeIDs in a clean way
  !!
  subroutine tem_dumpTreeIDlists( minLevel, maxLevel, LD )
    ! ---------------------------------------------------------------------------
    !> minimum level of the global tree
    integer, intent(in) :: minlevel, maxLevel
    !> level descriptor
    type(tem_levelDesc_type), intent(in) :: LD(minlevel:maxLevel)
    ! ---------------------------------------------------------------------------
    integer :: iLevel, first, last
    ! ---------------------------------------------------------------------------

    call tem_horizontalSpacer( dbgUnit(3), 2 )
    write(dbgUnit(3),*) '           Level Descriptor'
    call tem_horizontalSpacer( dbgUnit(3) )

    do iLevel = minLevel, maxLevel

      write(dbgUnit(3),"(A,I0)") '   Level: ', iLevel
      write(dbgUnit(3),"(A,I0)") '  nElems: ', LD(iLevel)%nElems
      call print_nElems( LD(iLevel)%elem%nElems(eT_minRelevant:eT_maxRelevant),&
        &                dbgUnit(5) )

      first = LD(iLevel)%offset(1,eT_fluid) + 1
      last  = LD(iLevel)%offset(2,eT_fluid)
      call tem_printArray( title = 'fluid',                     &
        &                  me = LD(iLevel)%total( first:last ), &
        &                  itemLength = 10,                     &
        &                  lineLength = 120,                    &
        &                  nUnit = dbgUnit(3)                   )

      if ( LD(iLevel)%elem%nElems( eT_ghostFromCoarser ) > 0 ) then
        first = LD(iLevel)%offset(1,eT_ghostFromCoarser) + 1
        last  = LD(iLevel)%offset(2,eT_ghostFromCoarser)
        call tem_printArray( title = 'ghostFromCoarser',          &
          &                  me = LD(iLevel)%total( first:last ), &
          &                  itemLength = 10,                     &
          &                  lineLength = 120,                    &
          &                  nUnit = dbgUnit(3)                   )
      end if

      if ( LD(iLevel)%elem%nElems( eT_ghostFromFiner ) > 0 ) then
        first = LD(iLevel)%offset(1,eT_ghostFromFiner) + 1
        last  = LD(iLevel)%offset(2,eT_ghostFromFiner)
        call tem_printArray( title = 'ghostFromFiner',          &
          &                  me = LD(iLevel)%total( first:last ), &
          &                  itemLength = 10,                     &
          &                  lineLength = 120,                    &
          &                  nUnit = dbgUnit(3)                   )
      end if

      if ( LD(iLevel)%elem%nElems( eT_halo ) > 0 ) then
        first = LD(iLevel)%offset(1,eT_halo) + 1
        last  = LD(iLevel)%offset(2,eT_halo)
        call tem_printArray( title = 'halo',          &
          &                  me = LD(iLevel)%total( first:last ), &
          &                  itemLength = 10,                     &
          &                  lineLength = 120,                    &
          &                  nUnit = dbgUnit(3)                   )
      end if
    end do ! iLevel

  end subroutine tem_dumpTreeIDlists
! ****************************************************************************** !


! ****************************************************************************** !
  !> This routine updates the property bits in the tree with those of the
  !! level descriptor.
  !!
  subroutine tem_updateTree_properties( levelDesc, tree )
    ! ---------------------------------------------------------------------------
    !> level descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> global tree
    type(treelmesh_type), intent(inout) :: tree
    ! ---------------------------------------------------------------------------
    ! counter
    integer :: iElem
    ! ---------------------------------------------------------------------------

    do iElem = 1, levelDesc%elem%nElems( eT_fluid )
      if( tree%ElemPropertyBits( levelDesc%pntTID(iElem)) /= &
        & levelDesc%property( iElem ))then
        tree%ElemPropertyBits( levelDesc%pntTID( iElem )) = &
          &             levelDesc%property( iElem )

        ! set the global mesh changed tag to be true such that the mesh is
        ! dumped
        ! @todo: why set meshchange here?
        tree%global%meshChange = .true.
      end if
    end do

  end subroutine tem_updateTree_properties
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine depSource_append( me, sourceList, mySources, n )
    ! ---------------------------------------------------------------------------
    type( depSource_type )    :: me
    type( dyn_intArray_type ) :: sourceList
    integer, intent(in) :: n
    integer, intent(in) :: mySources(n)
    ! ---------------------------------------------------------------------------
    integer :: iElem, posInSourceList
    ! ---------------------------------------------------------------------------

    ! initialize array for storing the positions of the source elements in
    ! the levelwise buffer
    call destroy( me = me%elem )
    call init( me = me%elem,        length = n )
    call init( me = me%elemBuffer,  length = n )

    do iElem = 1, n

      ! Add the found element to the dependency element list
      call append( me  = me%elem, &
        &          val = mySources(iElem) )

      ! Store to the levelwise source element list
      call append( me  = sourceList,     &
        &          val = mySources(iElem),     &
        &          pos = posInSourceList )

      ! ... and store the position in the levelwise source element
      ! list for each ghost dep
      call append( me  = me%elemBuffer,  &
        &          val = posInSourceList )

    end do

  end subroutine depSource_append
! ****************************************************************************** !

! ****************************************************************************** !
  !! Calculate nearest neighbors.
  !! if its fluid then identify its treeID from tem_IdOfCoord
  !! If neighbor is boundary then identify the boundary ID from boundary_ID list
  !!
  subroutine tem_calc_neighbors( posInBCID, boundary_ID, stencil, x, neighIDs )
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: posInBCID
    integer(kind=long_k), intent(in) :: boundary_ID(:,:)
    type( tem_stencilHeader_type ), intent(in) :: stencil
    integer, intent(in) :: x(4)
    integer(kind=long_k), intent(out) :: neighIDs(stencil%QQN)
    ! ---------------------------------------------------------------------------
    integer(kind=long_k) :: bcID, tOffset
    integer :: iQQN, iQQQ
    logical :: noBC
    ! ---------------------------------------------------------------------------

    neighIDs(:) = 0_long_k
    tOffset = tem_FirstIdAtLevel( x(4) )

    directionLoop: do iQQN = 1, stencil%QQN

      iQQQ = stencil%map( iQQN )
      noBC = .false.

      if ( posInBCID > 0 .and. iQQQ > 0 ) then
        ! there is a boundary in the some direction
        bcID = boundary_ID( iQQQ, posInBCID )
        if ( bcID > 0 ) then
          ! Set the negative boundary ID as neighbor, so we don't need
          ! more arrays for boundary information
          neighIDs(iQQN) = -bcID
        else if ( bcID < 0 ) then
          ! Is a periodic neighbor! value is -treeID_periodicNeighbor
          ! Set the id given as treeID and make positive
          neighIDs(iQQN) = abs( bcID )
        else ! bcID == 0
          noBC = .true.
        end if
      else ! posInBCID == 0
        noBC = .true.
      end if ! posInBCID > 0?

      if ( noBC ) then
        ! NO boundary condition for this direction.
        ! Just store the neighboring treeID which will have to be
        ! found or reconstructed (might be an actual fluid element
        ! or ghost / halo)
        neighIDs(iQQN) = tem_IdOfCoord(                       &
          &               [ x(1) + stencil%cxDir( 1, iQQN ),  &
          &                 x(2) + stencil%cxDir( 2, iQQN ),  &
          &                 x(3) + stencil%cxDir( 3, iQQN ),  &
          &                 x(4) ], tOffset)
      end if

    end do directionLoop

  end subroutine tem_calc_neighbors
! ****************************************************************************** !

! ****************************************************************************** !
  !> Allocate level descriptor and initilize its variables
  !!
  subroutine tem_alloc_levelDesc( me, minLevel, maxLevel, initlen, nStencils )
    ! ---------------------------------------------------------------------------
    type(tem_levelDesc_type), allocatable, intent(out)  :: me(:)
    integer, intent(in) :: minLevel, maxLevel, initlen, nStencils
    ! ---------------------------------------------------------------------------
    integer :: iLevel
    ! ---------------------------------------------------------------------------

    write(logUnit(5),*) 'Allocating level Descriptor ...'
    allocate( me(minLevel:maxLevel) )
    me(minLevel:maxLevel)%nElems = 0

    ! Initial length for the element array
    ! We don't know how many elements there will be per level, just using
    ! the average here as a starting point. This should at least for single
    ! level meshes reduce the number of expands, we need to do.
    ! Giving a small space for additional elements (ghosts and halos, might
    ! avoid immidiate doubling, for elements beyond the fluids).

    ! initialize dynamic require array for requireNeighNeigh
    ! and element type
    write(logUnit(7),*) 'Initializing elemList, require and neigh array'
    do iLevel = minLevel, maxLevel
      call init( me = me( iLevel )%require )
      call init( me = me( iLevel )%elem,    length = initlen )
      ! For each of the neighbor lists create the horizontal (neighbor)
      ! relations of the size of the total number of stencils used in this
      ! scheme
      allocate( me( iLevel )%neigh( nStencils ) )
    end do

  end subroutine tem_alloc_levelDesc
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine tem_build_levelPointer( levelPointer, tree, levelDesc )
    ! ---------------------------------------------------------------------------
    !> Pointer from original treeID list to level wise fluid list
    integer, allocatable, intent(out)     :: levelPointer(:)
    !> the global tree
    type(treelmesh_type), intent(in)      :: tree
    !> the level descriptor to be filled
    type(tem_levelDesc_type), intent(in)  :: levelDesc(tree%global%minlevel:)
    ! ---------------------------------------------------------------------------
    integer :: iLevel, iElem
    ! ---------------------------------------------------------------------------

    ! KM: NOTE: scheme independent levelPointer is build for each scheme
    ! if ( allocated( levelPointer ) ) then deallocate( levelPointer )
    allocate( levelPointer( tree%nElems ))

    do iElem = 1, tree%nElems
      iLevel = tem_levelOf( tree%treeID( iElem ))
      levelPointer( iElem ) = tem_PositionInSorted(                            &
        &                      me    = levelDesc( iLevel )%total,              &
        &                      val   = tree%treeID( iElem ),                   &
        &                      upper = levelDesc( iLevel )%elem%nElems( eT_fluid ) )
    end do

  end subroutine tem_build_levelPointer
! ****************************************************************************** !

end module tem_construction_module
! ****************************************************************************** !
!> # Example for the stencil construction
!!
!! To build the stencil for a given element, all required
!! [[treelmesh_type:treeID]]s can be obtained with
!! the following procedure:
!! First the integer coordinate tuple \((x,y,z,L)\) of the
!! [[treelmesh_type:treeID]], for which the stencil
!! has to be obtained, is computed.
!! Then for each stencil element the
!! [[treelmesh_type:treeID]] is obtained by the
!! following steps:
!! - Add the stencil offset to the coordinate of the central element to obtain
!!   the coordinate \((x + s_x, y + s_y, z + s_z, L)\) for the stencil element
!! - Convert the obtained coordinate back to a
!!   [[treelmesh_type:treeID]].
!!
!! To illustrate this procedure with a specific example, lets consider a two
!! dimensional mesh.
!! We would like to find the right neighbor of the element with
!! [[treelmesh_type:treeID]] \( 40 \).
!! An illustration of the geometrical layout for this example is given in the
!! following Figure.
!! As we are looking at a two dimensional mesh, there are \( 4^L \) elements on
!! each refinement level \( L \).
!!
!! ![Illustration of the right neighbor identification](|media|/find_rightNeighbor.png)
!!
!! The first information, that is needed, is the level of the given
!! [[treelmesh_type:treeID]].
!! To find this, and the Morton index on that level, we subtract the elements
!! on coarse levels, until the reduced index is smaller than the next number of
!! elements to subtract:
!!
!! - Level 0: \(40 - 4^0 = 39\)
!! - Level 1: \(39 - 4^1 = 35\)
!! - Level 2: \(35 - 4^2 = 19\)
!! - Level 3: \(4^3 > 19\)
!! - [[treelmesh_type:treeID]] \(40\) is located
!!   on Level 3 and has a Morton index of \(19\).
!!
!! The binary representation of \(19\) is \)010011_b\).
!! Due to the bit-interleaving rule, odd bits correspond to the X coordinate and
!! even bits correspond to the Y coordinate.
!! Thus we obtain for X \(101_b = 5\) and for Y \)001_b = 1\).
!! Now we found the coordinate tuple \((x, y, L) = (5, 1, 3)\) of
!! [[treelmesh_type:treeID]] \(40\) and can
!! simply add the offset \((1, 0, 0)\), that represents the right neighbor.
!! By this we obtain the coordinate of this neighbor to be \((6, 1, 3)\),
!! resulting in a binary representation of \((110_b, 001_b)\).
!! After interleaving the bits, the Morton index of the right neighbor is found
!! to be \(010110_b = 22\)
!! Finally, by adding the offset, we obtain the
!! [[treelmesh_type:treeID]] of the right neighbor
!! as \(22 + 4^0 + 4^1 + 4^2 = 43\).
!! Note, that there is no restriction on the offsets, that are to be applied in
!! the neighborhood search, and due to the periodic universe all offset
!! locations are well defined.
