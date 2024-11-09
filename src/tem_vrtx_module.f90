! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2012-2013, 2016, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014, 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Raphael Haupt <Raphael.Haupt@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2019 Daniel Fleischer <daniel.fleischer@student.uni-siegen.de>
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
!> author: Simon Zimny
!! author: Kartik Jain
!!
!! A module to provide informations on vertex property for elements.
!!
module tem_vrtx_module

  ! include treelm modules
  use mpi
  use env_module,              only: rk, long_k, eps, globalMaxLevels
  ! use tem_prophead_module,     only: tem_prophead_type
  use tem_property_module,     only: tem_property_type, prp_hasQVal
  use tem_grow_array_module,   only: grw_longArray_type, grw_real2dArray_type, &
    &                                init, append, expand, destroy
  use treelmesh_module,        only: treelmesh_type
  use tem_param_module,        only: qOffset, qBSW, qBSE, qBNW,                &
    &                                qBNE, qTSW, qTSE, qTNW, qTNE, q_NE, q_SE, &
    &                                q_NW, q_SW, q_TE, q_TW, q_BE, q_BW, q_TN, &
    &                                q_BN, q_TS, q_BS
  use tem_geometry_module,     only: tem_originOfID, tem_ElemSize,             &
    &                                tem_posOfID, tem_baryOfID
  use tem_topology_module,     only: tem_CoordOfID, tem_IDofCoord
  use tem_tools_module,        only: append
  ! use tem_vrtx_prop_module,    only: tem_vrtx_prop_type
  use tem_bc_prop_module,      only: tem_BC_prop_type
  use tem_subTree_type_module, only: tem_subTree_type, tem_treeIDfrom_subTree
  ! use tem_global_module,       only: tem_global_type
  use tem_logging_module,      only: logUnit

  implicit none

  private

  public :: tem_calc_vrtx_coord, tem_vrtx_type
  public :: tem_vrtx_finalize

  ! -----------------------------------------------------------------------------
  !> Datatype for the vrtx dependend information.
  !! A dynamic array and a growing array are 'coupled'.
  !! The comparison between the real coordinates is
  !! shifted to the integer comparison of the dynamic array.
  !! The size of the two arrays are the same.
  !! Additionally a map of the 8 vertices for each element
  !! to the global index is stored.
  type tem_vrtx_type
    !> total number of vertices
    integer :: nVertices
    !> growing array to store the coordinates
    type(grw_real2dArray_type) :: coord
    !> map of vertices for each element to global index
    !! vrtx_index_map(nelems, 8 vertices)
    integer, allocatable :: map2global(:,:)
    !> simulation time that coordinate info belongs to
    ! real(kind=rk) :: sim_time
    !> max number of vertices
    integer :: maxVertices
    !> array of elements with qValues
    logical, allocatable :: refine(:)
  end type tem_vrtx_type
  ! -----------------------------------------------------------------------------

  ! offsets for the vertices of an element
  ! entries 1:3 -> offsets
  ! entry 4 -> level higher than the element level (always 1)
  ! entry 5 -> position in the q-Val list
  integer, dimension( 5, 20 ),parameter  :: vrtxMap =                          &
  ! first 8 entries correspond to the corners
  reshape((/ 0, 0, 0, 1, qBSW, &          !  1 Bottom South West(BSW)
             1, 0, 0, 1, qBSE, &          !  2 Bottom South East(BSE)
             0, 0, 1, 1, qTSW, &          !  3 Top South West(TSW)
             1, 0, 1, 1, qTSE, &          !  4 Top South East(TSE)
             0, 1, 0, 1, qBNW, &          !  5 Bottom North West(BNW)
             1, 1, 0, 1, qBNE, &          !  6 Bottom North East(BNE)
             0, 1, 1, 1, qTNW, &          !  7 Top North West(TNW)
             1, 1, 1, 1, qTNE, &          !  8 Top North East(TNE)
  ! the following entries are valid for the refined elements ONLY
  ! 4 entries for the points on the middle of the edges for the front plane
             1, 0, 0, 1, q_BS, &          !  9 Bottom South(BS)
             0, 0, 1, 1, q_SW, &          ! 10 South West(SW)
             2, 0, 1, 1, q_SE, &          ! 11 South East(SE)
             1, 0, 2, 1, q_TS, &          ! 12 Top South(TS)
  ! 4 entries for the points on the middle of the edges for the center plane
             0, 1, 0, 1, q_BW, &          ! 13 Bottom West(BW)
             0, 1, 2, 1, q_TW, &          ! 14 Top West(TW)
             2, 1, 0, 1, q_BE, &          ! 15 Bottom East(BE)
             2, 1, 2, 1, q_TE, &          ! 16 Top East(TE)
  ! 4 entries for the points on the middle of the edges for the back plane
             1, 2, 0, 1, q_BN, &          ! 17 Bottom North(BN)
             0, 2, 1, 1, q_NW, &          ! 18 North West(NW)
             2, 2, 1, 1, q_NE, &          ! 19 North East(NE)
             1, 2, 2, 1, q_TN/),(/5,20/)) ! 20 Top North(TN)

contains


! ****************************************************************************** !
  !> Run over all 8 vertices for each element in the treeID list, calculate
  !! its coordinates and add its position to the map.
  !!
  subroutine tem_calc_vrtx_coord( tree, vrtx, subTree, boundary, useQVal )
    ! ---------------------------------------------------------------------------
    !> fluid mesh
    type(treelmesh_type),     intent(in) :: tree
    !> Vertex data
    type(tem_vrtx_type), intent(inout) :: vrtx
    !> optional subTree information
    type(tem_subTree_type), optional, intent(in) :: subTree
    !> boundary information incl. q-Values
    type(tem_BC_prop_type), optional, intent(in) :: boundary
    !> use the qValue information?
    logical, optional, intent(in) :: useQVal
    ! ---------------------------------------------------------------------------
    ! counters
    integer :: iVrtx, iElem
    integer :: local_nElems
    ! store all treeIDs for the vertices of each element in vrtxTreeID
    ! in the order
    ! -----------------------------------------------------------------
    ! | vrtx1 vrtx2 vrtx3 vrtx4 vrtx5 vrtx6 vrtx7 vrtx8 |     ...     |
    ! -----------------------------------------------------------------
    !                 iElem = 1                            iElem =...
    integer(kind=long_k), allocatable :: vrtxTreeID(:)
    integer(kind=long_k), allocatable :: sortedVrtxTreeID(:)
    integer(kind=long_k) :: vrtxID
    integer :: elemCoord(4)
    integer :: locVrtx(4)
    integer :: vrtxAnchor(4)
    integer :: iLevel
    ! tree with bounding cube length twice as big as in tree (treeID array will
    ! not be filled!!!!)
    type(treelmesh_type) :: bigTree
    type(tem_property_type), pointer :: tree_property(:) => NULL()
    integer(kind=long_k), allocatable :: treeID(:)
    ! type(tem_global_type) :: global
    ! counters
    integer :: globCounter
    integer :: uniqueCounter
    integer :: nElemsQVal
    integer :: iBCElem
    real(kind=rk) :: coord(3)
    logical :: local_useQVal
    ! ---------------------------------------------------------------------------

    if( present( useQVal ))then
      local_useQVal = useQVal
    else
      local_useQVal = .false.
    end if

    if( present( subTree ))then
      local_nElems = subTree%nElems
    else
      local_nElems = tree%nElems
    end if

    uniqueCounter = 0

    if (local_useQval) then
      vrtx%maxVertices = 0
      if (associated(tree_property)) deallocate(tree_property)
      if (allocated(treeID)) deallocate(treeID)
      if (allocated(vrtx%refine)) deallocate(vrtx%refine)

      if( present( subTree ))then
        ! global = subTree%global
        allocate( tree_property( subTree%global%nProperties ))
        tree_property = subTree%Property
        allocate( treeID( local_nElems ))
        allocate( vrtx%refine( local_nElems ))
        call tem_treeIDfrom_subTree( subTree, tree, treeID, (/1,local_nElems/) )
        ! possible q-Values attached calculate max number of vertices (8 per
        ! 'normal' element, 20 per element with q-Values) and allocate the array
        do iElem = 1, local_nElems
          if( btest(subTree%elemPropertyBits( iElem ), prp_hasQVal )) then
            vrtx%maxVertices = vrtx%maxVertices + 20
            vrtx%refine( iElem ) = .true.
          else
            vrtx%maxVertices = vrtx%maxVertices + 8
            vrtx%refine( iElem ) = .false.
          end if
        end do
      else
        ! global = tree%global
        allocate( tree_property( tree%global%nProperties ))
        tree_property = tree%Property
        allocate( treeID( local_nElems ))
        allocate( vrtx%refine( local_nElems ))
        treeID = tree%treeID
        ! possible q-Values attached calculate max number of vertices (8 per
        ! 'normal' element, 20 per element with q-Values) and allocate the array
        do iElem = 1, local_nElems
          if (btest(tree%elemPropertyBits( iElem ), prp_hasQVal) ) then
            vrtx%maxVertices = vrtx%maxVertices + 20
            vrtx%refine( iElem ) = .true.
          else
            vrtx%maxVertices = vrtx%maxVertices + 8
            vrtx%refine( iElem ) = .false.
          end if
        end do
      end if

      ! initialize the vertex type
      call tem_init_vrtx_prop( vrtx = vrtx)

      ! allocate the list of all vrtxTreeIDs including dublicates
      allocate( vrtxTreeID( vrtx%maxVertices ))

      write(logUnit(6),*) 'DEBUG: Filling the global vrtxTreeID ...'

      ! initialize counters
      globCounter = 0

      ! map the treeIDs to those of a tree with a bounding cube length
      ! twice as big -> treeIDs correspond to those 1 refinement level
      ! higher in the bigger tree
      !                       -------------------------------
      !                       |              |              |
      !                       |              |              |
      !                       |              |              |
      !                       |              |              |
      !                       |              |              |
      ! ----------------      -------------------------------
      ! |      |       |      |      |       |              |
      ! |      |       |      |      |       |              |
      ! ----------------  --> ----------------              |
      ! |      |       |      |      |       |              |
      ! |      |       |      |      |       |              |
      ! ----------------      -------------------------------
      !

      nElemsQVal = 0
      do iElem = 1, local_nElems
        ! if element has q-Values it has to be refined once
        if( vrtx%refine( iElem ))then
           if( present( subTree ) .and. .not. subTree%useGlobalMesh )then
             do iBCElem = 1, size( tree%property(2)%elemID )
               if( tree%property(2)%elemID( iBCElem ) .eq.                       &
                 &                              subTree%map2global( iElem )) then
                 nElemsQVal = iBCElem
                 exit
               end if
             end do
           else
             nElemsQVal = nElemsQVal + 1
           end if
          ! calculate the vertices for the element incl. q-Values
          do iVrtx = 1, 20
            ! check if the q-Value for vertex iVrtx is greater than 0.5 (this
            ! means the point might be shared between elements) or the q-Value
            ! is -1.0 (this means that no q-Value is set in the corresponding
            ! direction)
            if( ( abs(boundary%qVal( vrtxMap(5, iVrtx), nElemsQVal) -0.5_rk) .le.&
              &                                            eps)              .or.&
!              & ( abs(boundary%qVal( vrtxMap(5, iVrtx), nElemsQVal) +0.5_rk) .gt.&
!              &                                       1.0_rk + eps)          .or.&
              & ( abs(boundary%qVal( vrtxMap(5, iVrtx), nElemsQVal) +1.0_rk) .le.&
              &                                            eps) )then

              ! for the 8 corners get them from the treeID
              if( iVrtx .le. 8)then
                elemCoord = tem_coordOfID(treeID(iElem))
              else ! for the 12 intermediate get them from the refined treeID
                ! refine the element by 1 level
                elemCoord = tem_coordOfID( treeID(iElem)*8_long_k + 1_long_k )
              end if
              ! since the coordinates of the individual vertices are on level 1
              ! the level for the vrtxAnchor is increased by 1 matching the
              ! requirements of the new tree (bounding cube size)
              vrtxAnchor = elemCoord + vrtxMap( 1:4, iVrtx )
              ! retransforming the coords to the treeID on the 'new tree'
              vrtxID = tem_IDofCoord(vrtxAnchor)
              ! get the treeID on the highest refinement level possible as
              ! a unique identifier
              do iLevel=vrtxAnchor(4)+1,globalMaxLevels
                vrtxID = vrtxID*8_long_k + 1_long_k
              end do
            else ! q-Value .ne. 0.5 (assume vertex is unique)
              uniqueCounter = uniqueCounter - 1
              vrtxID = uniqueCounter
              ! calculate the vertex based on the q-Value and append it to the
              ! growing array of vertices
              coord = tem_calc_vrtxOf_qVal(                                      &
                &         treeID = treeID(iElem),                                &
                &         tree   = tree,                                         &
                &         qVal   = boundary%qVal( vrtxMap(5, iVrtx), nElemsQVal),&
                &         iVrtx  = iVrtx )
              call append( me = vrtx%coord, val = coord )
            end if
            vrtxTreeID( globCounter+iVrtx ) = vrtxID
          end do
          globCounter = globCounter + 20
        else ! no q-Values
          elemCoord = tem_coordOfID(treeID(iElem))
          do iVrtx=1,8
            locVrtx = tem_coordOfID(int(iVrtx, kind=long_k))
            ! since the coordinates of the individual vertices are on level 1
            ! the level for the vrtxAnchor is increased by 1 matching the
            ! requirements of the new tree (bounding cube size)
            vrtxAnchor = elemCoord + locVrtx
            ! retransforming the coords to the treeID on the 'new tree'
            vrtxID = tem_IDofCoord(vrtxAnchor)
            ! get the treeID on the highest refinement level possible as
            ! a unique identifier
            do iLevel=vrtxAnchor(4)+1,globalMaxLevels
              vrtxID = vrtxID*8_long_k + 1_long_k
            end do
            vrtxTreeID( globCounter+iVrtx ) = vrtxID
          end do
          globCounter = globCounter + 8
        end if
      end do

      write(logUnit(6),*) 'DEBUG: Filled it.'
    else
      call tem_calc_vrtx_coord_noqval( tree, vrtx, subTree, vrtxTreeID )
    end if

    allocate( sortedVrtxTreeID( vrtx%maxVertices ))
    sortedVrtxTreeID = vrtxTreeID

    write(logUnit(6),*) 'DEBUG: Start sorting ...'

    ! sort the treeID array
    call qsort_vrtx( sortedVrtxTreeID )

    write(logUnit(6),*) 'DEBUG: Done sorting, start inverting coords ...'

    ! in case q-Values are present reorganize the growing array of coords
    ! such that it is in the correct order
    if( (-1)*uniqueCounter .gt. 1 )then
      call tem_invertRealRkArray( me = vrtx%coord%val, nElems = vrtx%coord%nVals )
    end if

    write(logUnit(6),*) 'DEBUG: Done inverting, start to unify ...'

    ! redefine the tree bounding cube size
    bigTree%global%origin = tree%global%origin
    bigTree%global%BoundingCubeLength = 2.0_rk * tree%global%BoundingCubeLength

    ! make sorted array vrtxTreeID unique and map the elements to the right
    ! vertex real coordinates
    call tem_unify_vrtx( inList   = sortedVrtxTreeID,                          &
      &                  origList = vrtxTreeID,                                &
      &                  coord    = vrtx%coord,                                &
      &                  map      = vrtx%map2global,                           &
      &                  tree     = bigTree,                                   &
      &                  nElems   = local_nElems,                              &
      &                  nUnique  = (-1)*uniqueCounter,                        &
      &                  refine = vrtx%refine )

    write(logUnit(6),*) 'DEBUG: Done unifying.'

    ! update the number of calculated vertices
    vrtx%nVertices = vrtx%coord%nVals

  end subroutine tem_calc_vrtx_coord
! ****************************************************************************** !

! ****************************************************************************** !
  !> Run over all 8 vertices for each element in the treeID list, calculate
  !! its coordinates and add its position to the map.
  !!
  subroutine tem_calc_vrtx_coord_noqval( tree, vrtx, subTree, vrtxTreeID )
    ! ---------------------------------------------------------------------------
    !> fluid mesh
    type(treelmesh_type),     intent(in) :: tree
    !> Vertex data
    type(tem_vrtx_type), intent(inout) :: vrtx
    !> optional subTree information
    type(tem_subTree_type), optional, intent(in) :: subTree
    !> boundary information incl. q-Values
    integer(kind=long_k), allocatable, intent(out) :: vrtxTreeID(:)
    ! ---------------------------------------------------------------------------
    ! counters
    integer :: iVrtx, iElem
    integer :: local_nElems
    integer(kind=long_k) :: vrtxID
    integer :: elemCoord(4)
    integer :: locVrtx(4)
    integer :: vrtxAnchor(4)
    integer :: iLevel
    integer(kind=long_k), allocatable :: treeID(:)
    ! ---------------------------------------------------------------------------

    vrtx%maxVertices = 0

    if (allocated(treeID)) deallocate(treeID)
    if (allocated(vrtx%refine)) deallocate(vrtx%refine)

    if( present( subTree ))then
      local_nElems = subTree%nElems
      allocate( treeID( local_nElems ))
      call tem_treeIDfrom_subTree( subTree, tree, treeID, (/1,local_nElems/) )
    else
      local_nElems = tree%nElems
      allocate( treeID( local_nElems ))
      treeID = tree%treeID
    end if
    vrtx%maxVertices = 8*local_nElems
    allocate( vrtx%refine( local_nElems ))
    vrtx%refine = .false.

    ! initialize the vertex type
    call tem_init_vrtx_prop( vrtx = vrtx)

    ! allocate the list of all vrtxTreeIDs including dublicates
    allocate( vrtxTreeID( vrtx%maxVertices ))

    write(logUnit(6),*) 'DEBUG: Filling the global vrtxTreeID ...'

    ! map the treeIDs to those of a tree with a bounding cube length
    ! twice as big -> treeIDs correspond to those 1 refinement level
    ! higher in the bigger tree
    !                       -------------------------------
    !                       |              |              |
    !                       |              |              |
    !                       |              |              |
    !                       |              |              |
    !                       |              |              |
    ! ----------------      -------------------------------
    ! |      |       |      |      |       |              |
    ! |      |       |      |      |       |              |
    ! ----------------  --> ----------------              |
    ! |      |       |      |      |       |              |
    ! |      |       |      |      |       |              |
    ! ----------------      -------------------------------
    !

    do iVrtx=1,8
      locVrtx = tem_coordOfID(int(iVrtx, kind=long_k))
      do iElem = 1, local_nElems
        elemCoord = tem_coordOfID(treeID(iElem))

        ! since the coordinates of the individual vertices are on level 1
        ! the level for the vrtxAnchor is increased by 1 matching the
        ! requirements of the new tree (bounding cube size)
        vrtxAnchor = elemCoord + locVrtx
        ! retransforming the coords to the treeID on the 'new tree'
        vrtxID = tem_IDofCoord(vrtxAnchor)
        ! get the treeID on the highest refinement level possible as
        ! a unique identifier
        !NEC$ NOVECTOR
        do iLevel=vrtxAnchor(4)+1,globalMaxLevels
          vrtxID = vrtxID*8_long_k + 1_long_k
        end do
        vrtxTreeID( (iElem-1)*8 + iVrtx ) = vrtxID
      end do
    end do

    write(logUnit(6),*) 'DEBUG: Filled it.'


  end subroutine tem_calc_vrtx_coord_noqval
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine takes the sorted list as an input and unifies its entries
  !! the result is used to create a unique array of vertex coordinates and a map
  !! for the 8 vertices of each element.
  !!
  subroutine tem_unify_vrtx( inList, origList, coord, map, tree, nElems,       &
    &                        nUnique, refine )
    ! ---------------------------------------------------------------------------
    integer(kind=long_k), allocatable, intent(inout) :: inList(:)
    integer(kind=long_k), allocatable, intent(inout) :: origList(:)
    type(grw_real2dArray_type) :: coord
    integer, allocatable, intent(inout) :: map(:,:)
    type(treelmesh_type), intent(in) :: tree
    integer, intent(in) :: nElems
    !> number of unique vertices (from q-Values)
    integer, intent(in) :: nUnique
    logical, intent(in) :: refine(:)
    ! ---------------------------------------------------------------------------
    ! counters
    integer :: count1, count2
    integer :: iElem, iVrtx, pos, counter
    type( grw_longArray_type ) :: unique
    real(kind=rk) :: tmp_vrtx(3)
    ! ---------------------------------------------------------------------------
    count1 = nUnique+1
    count2 = nUnique+1

    call init( me = unique, length = 10 )

    ! append the unique treeIDs to the unique array
    do iElem = 1, nUnique
      call append( me  = unique,                                               &
        &          val = inList( iElem ))
    end do

    ! at first make the list of vertex treeIDs and their coordinates unique
    do while( count1 .lt. size( inList))
      count1 = count2
      call append( me  = unique,                                               &
        &          val = inList( count1 ))

      ! get the real coordinates of the unique treeID ...
      tmp_vrtx = tem_originOfId( tree, inList( count1 ))
      ! ... and store them in by definition unique array
      call append( me  = coord,                                                &
        &          val = tmp_vrtx )

      do while(( inList( count1 ) .eq. inList( count2 ))                       &
        &        .and. count2 .lt. size( inList ))
        count2 = count2+1
      end do
    end do

    deallocate( inList )

    counter = 0
    ! map the original treeID list to the unique and by this to the unique coord
    ! array
    do iElem = 1, nElems
      if( refine( iElem ))then
        do iVrtx = 1, 20
          ! if the element has a valid treeID search for it using the regular
          ! posOfID
          if( origList( counter + iVrtx ) .gt. 0 )then
            pos = tem_posOfID( origList( counter + iVrtx ),                    &
              &                unique%val( nUnique+1:unique%nVals ))
            pos = pos + nUnique
          ! if not do a simple search in the first part of the unique array
          ! (the negative entries)
          else
            pos = tem_posOfLong( origList( counter + iVrtx ),                  &
              &                  unique%val( 1:nUnique ))
          end if
          ! map the vertex
          call append( array     = map,                                        &
            &          position1 = iElem,                                      &
            &          position2 = iVrtx,                                      &
            &          value     = pos )
        end do
        counter = counter + 20
      else
        do iVrtx = 1, 8
          pos = tem_posOfID( origList( counter+iVrtx ),                        &
            &                unique%val( nUnique+1:unique%nVals ))
          pos = pos + nUnique
          call append( array     = map,                                        &
            &          position1 = iElem,                                      &
            &          position2 = iVrtx,                                      &
            &          value     = pos )
        end do
        counter = counter + 8
      end if
    end do

    deallocate(origList)
    call destroy( me = unique )

  end subroutine tem_unify_vrtx
! *******************************************************************************


! *******************************************************************************
  !> Quicksort for long integer kinds.
  !!
  recursive subroutine qsort_vrtx( list )
    ! ---------------------------------------------------------------------------
    !> list to be sorted
    integer( kind=long_k ), intent(inout)  :: list(:)
    ! ---------------------------------------------------------------------------
    integer :: split
    ! ---------------------------------------------------------------------------

    ! recursive call of qsort
    if( size( list ) .gt. 1)then
      call partition( list, split )
      call qsort_vrtx( list( :split-1 ))
      call qsort_vrtx( list( split: ))
    end if

  end subroutine qsort_vrtx
! *******************************************************************************


! *******************************************************************************
  !> This subroutine partitions the given list for the quicksort algorithm
  !!
  subroutine partition( list, marker )
    ! ---------------------------------------------------------------------------
    !> list to be partitioned
    integer( kind=long_k ), intent(inout)  :: list(:)
    !> marker where the list is partitioned
    integer, intent(out) :: marker
    ! ---------------------------------------------------------------------------
    integer :: left, right
    integer(kind=long_k) :: pivot, temp
    ! ---------------------------------------------------------------------------

    ! set the average of first and last entry as pivot element (bug for entries
    ! < 0)
!    pivot = ( list(1) + list( size( list )))/2
    ! choose the element in the middle as pivot element
    pivot = list(size(list)/2)

    ! initialize the pointers on the entries
    left = 0
    right = size( list ) + 1

    do while( left .lt. right)
      right = right - 1
      do while( list( right ) .gt. pivot)
        right = right - 1
      end do
      left = left + 1
      do while( list( left ) .lt. pivot)
        left = left + 1
      end do
      if( left .lt. right )then
        temp = list( left )
        list(left) = list(right)
        list(right) = temp
      end if
    end do

    if( left .eq. right )then
      marker = left + 1
    else
      marker = left
    end if

  end subroutine partition
! *******************************************************************************


! *******************************************************************************
  !> Invert a given array
  !! 1 2 3 4 5 -> 5 4 3 2 1
  !!
  subroutine tem_invertRealRkArray( me, nElems )
    ! ---------------------------------------------------------------------------
    !> array to invert
    real(kind=rk), intent(inout) :: me(:,:)
    !> number of elements in the array
    integer, intent(in) :: nElems
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: tmp_entry(3)
    integer :: maxPos, minPos
    ! ---------------------------------------------------------------------------

    minPos = 1
    maxPos = nElems
    do while( minPos .lt. maxPos )
      tmp_entry = me( :, minPos )
      me( :, minPos ) = me( :, maxPos )
      me( :, maxPos ) = tmp_entry
      minPos = minPos + 1
      maxPos = maxPos - 1
    end do

  end subroutine tem_invertRealRkArray
! *******************************************************************************


! ****************************************************************************** !
  !> Initialize the vertex property headers.
  !!
  subroutine tem_init_vrtx_prop(vrtx)
    ! ---------------------------------------------------------------------------
    !> vertex type
    type(tem_vrtx_type)         :: vrtx
    ! ---------------------------------------------------------------------------
    integer :: init_nelems = 1024
    ! ---------------------------------------------------------------------------
    if ( .not.allocated(vrtx%map2global) ) then
      allocate( vrtx%map2global(init_nelems,8) )
      ! vrtx%sim_time = 0.d0
      vrtx%nvertices = 0
    end if

    ! initialise the growing array of actual vrtx real coordinates
    call init( me = vrtx%coord, width = 3 )

  end subroutine tem_init_vrtx_prop
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine calculates the vertex coordinate for a given element
  !! depending on the treeID, the global tree, the q-Value and iVrtx
  !!
  function tem_calc_vrtxOf_qVal( treeID, tree, qVal, iVrtx) result( coord )
    ! ---------------------------------------------------------------------------
    integer(kind=long_k), intent(in) :: treeID
    type(treelmesh_type), intent(in) :: tree
    real(kind=rk), intent(in) :: qVal
    integer, intent(in) :: iVrtx
    real(kind=rk) :: coord(3)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: bary(3)
    real(kind=rk) :: dx
    ! ---------------------------------------------------------------------------

    ! get the barycenter of the current treeID
    bary = tem_BaryOfId( tree, treeID )
    ! get the length of the element
    dx = tem_ElemSize( tree, treeID )

    coord = real(qOffset( vrtxMap(5, iVrtx), : ), kind=rk)*qVal*dx + bary

  end function tem_calc_vrtxOf_qVal
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function detects the first position of an integer value of kind
  !! long_k in an array. When there is no match the return value is 0.
  !!
  pure function tem_posOfLong( long, array ) result( pos )
    ! ---------------------------------------------------------------------------
    integer(kind=long_k), intent(in) :: long
    integer(kind=long_k), intent(in) :: array(:)
    integer :: pos
    ! ---------------------------------------------------------------------------
    integer :: iEntry
    ! ---------------------------------------------------------------------------

    do iEntry = 1, size(array)
      if( long == array(iEntry))then
        pos = iEntry
        return
      end if
    end do

    pos = 0

  end function tem_posOfLong
! ****************************************************************************** !


  ! ************************************************************************ !
  !> Clean up allocated memory in vrtx structure
  subroutine tem_vrtx_finalize(vrtx)
    type(tem_vrtx_type), intent(inout) :: vrtx

    if (allocated(vrtx%refine)) deallocate(vrtx%refine)
  end subroutine tem_vrtx_finalize
  ! ************************************************************************ !


end module tem_vrtx_module
! ****************************************************************************** !
