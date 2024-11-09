! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2013, 2020-2021 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011 Aravindh Krishnamoorthy <aravindh28.4@gmail.com>
! Copyright (c) 2011 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011 Laura Didinger <l.didinger@grs-sim.de>
! Copyright (c) 2011 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2011 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2011-2012, 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2011-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012, 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2015-2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
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
! **************************************************************************** !
!> Geometric methods for the TreElM module to act on a octree based mesh.
!! Make sure to read the introduction to the distributed octree in the
!! Documentation.
!!
!! This module contains methods that work with the actual sparse mesh, and
!! therefore depends on the [[Treelmesh_module]].
module tem_geometry_module

  ! include treelm modules
  use mpi
  use env_module,              only: rk, long_k, output_unit, rk_mpi,          &
    &                                globalMaxLevels, labelLen
  use tem_float_module,        only: operator(.feq.)
  use tem_topology_module,     only: tem_LevelOf, tem_PathOf, tem_path_type,   &
    &                                tem_PathComparison,   &
    &                                tem_coordOfID
  use tem_param_module,        only: q__W, q__E, q__S, q__N, q__B, q__T, q_NE, &
    &                                q_SW, q_NW, q_SE, q_TE, q_BW, q_TW, q_BE, &
    &                                q_TN, q_BS, q_TS, q_BN, qTNE, qBSW, qBSE, &
    &                                qTNW, qTSE, qBNW, qBNE, qTSW,             &
    &                                childPosition
  use tem_tools_module,        only: tem_longList, append
  use treelmesh_module,        only: treelmesh_type
  use tem_property_module,     only: prp_hasBnd, prp_hasQVal
  use tem_subTree_type_module, only: tem_subTree_type
  use tem_debug_module,        only: tem_debug_type, tem_debug, main_debug, &
    &                                dbgUnit
  use tem_logging_module,      only: tem_toStr, logUnit

  implicit none

  private

  public :: tem_CoordOfReal
  public :: tem_tIDinfo
  public :: tem_BaryOfCoord, tem_BaryOfId
  public :: tem_originOfId
  public :: tem_endOfId
  public :: tem_ElemSize, tem_ElemSizeLevel
  public :: tem_vrtxCoordOfId
  public :: tem_eligibleChildren
  public :: tem_findElement, tem_findPath
  public :: tem_PosOfId, tem_PosOfPath
  public :: tem_direction_type
  public :: tem_neighbor_type
  public :: tem_GetLocalBoundingCube
  public :: tem_setEffBoundingBox
  public :: tem_determine_discreteVector
  public :: tem_intp_bilinear
  public :: tem_intp_trilinear
  public :: tem_build_treeToProp_pointer

  !> derived type for each direction
  type tem_direction_type
    !> first entry of the element list
    type(tem_longList),pointer :: first => null()
    !> number of entries in element list
    integer :: nElems
    !> which list JZ:(whether it is fluid, ghost or halo, see
    !! tem_levelDesc_type)
    integer :: list
    !> position in list
    integer :: pos
    !> the current neighbors are of a different level
    logical :: otherLevel = .false.
  end type tem_direction_type

  !> type neighbor_type includes the direct neighbors of each tree ID
  type tem_neighbor_type
    !> number of directions in which the neighbors are located
    !! e.g. for cubes and surface neighbors this will be 6
    !! and for cubes with surface, edge and vertex neighbors this will
    !! be 22
    !! because stencils will be mapped by this type as well, we are not able
    !! to use it as a constant
    integer :: nNeighborDirections
    !! directions for this neighbor
    !! JZ: in case of reconstructed DG this has to be extended to all
    !! JZ: cells inside the stencil
    type(tem_direction_type),allocatable :: dir(:)
  end type tem_neighbor_type

  interface tem_intp_bilinear
    module procedure tem_intp_bilinear_scalar
    module procedure tem_intp_bilinear_vec
  end interface
  interface tem_intp_trilinear
    module procedure tem_intp_trilinear_scalar
    module procedure tem_intp_trilinear_vec
  end interface

  interface tem_setEffBoundingBox
    module procedure tem_setEffBoundingBox_fromTree
    module procedure tem_setEffBoundingBox_fromSubTree
  end interface
  interface tem_GetLocalBoundingCube
    module procedure tem_GetLocalBoundingCube_fromTree
    module procedure tem_GetLocalBoundingCube_fromSubTree
  end interface


  contains


  ! ************************************************************************ !
  !> This function returns a coordinate in the given treelmesh for a physical
  !! point location on the finest possible level.
  !!
  pure function tem_CoordOfReal(mesh, point, level) result(coord)
    ! -------------------------------------------------------------------- !
    type(treelmesh_type), intent(in) :: mesh !< Mesh to locate the point in
    real(kind=rk), intent(in) :: point(3) !< Point to look up
    integer, intent(in), optional :: level !< optional level to return the
                                           !! coordinate on
    integer :: coord(4) !< x,y,z,level
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: locInCube(3)
    real(kind=rk) :: meshDensity
    integer :: dimLen
    integer :: coordlevel
    ! -------------------------------------------------------------------- !

    if (present(level)) then
      coordlevel = min(level, globalMaxLevels)
    else
      coordlevel = globalMaxLevels
    end if

    locInCube = point - mesh%global%Origin
    dimLen = 2**coordlevel
    meshDensity = real(dimLen, kind=rk) / mesh%global%BoundingCubeLength

    ! Look up the real coordinate on the finest possible resolution.
    coord(4) = coordlevel

    ! Coordinates range from 0 to dimLen-1 in each direction
    ! Do not use periodic domain here, instead use the first
    ! and last element to capture all points outside the
    ! domain to deal with numerical inaccuracies.
    ! That is all elements include their lower boundaries, except
    ! the last one in a given direction, which includes the lower
    ! as well as the upper.
    coord(1:3) = max( min( int(locInCube*meshDensity), dimLen-1 ), 0 )

  end function tem_CoordOfReal
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Prints information about a treeID to a file unit.
  !!
  subroutine tem_tIDinfo( me, tree, nUnit)
    ! -------------------------------------------------------------------- !
    !> given level to get the size for
    integer(kind=long_k), intent(in) :: me
    !> Mesh to locate the point in
    type(treelmesh_type), intent(in) :: tree
    !> the file unit to use for printing the information.
    integer, intent(in), optional :: nUnit
    ! -------------------------------------------------------------------- !
    integer :: writeUnit, x(4)
    ! -------------------------------------------------------------------- !
    if ( present( nUnit )) then
      writeUnit = nUnit
    else
      writeUnit = output_unit
    end if
    write(writeUnit,*) '  treeID info ', me
    x(:) = tem_coordOfId( me )
    write(writeUnit,*) '       level   ', x(4)
    write(writeUnit,*) '       coord   ', x(1:3)
    write(writeUnit,*) '       elemSize', tem_elemSize( tree, me )
    write(writeUnit,*) '       bary    ', tem_baryOfId( tree, me )

  end subroutine tem_tIDinfo
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Return the size of elements on a given levle in the mesh by taking into
  !! account the size of the bounding cube given in the global info of the tree
  !!
  pure function tem_ElemSizeLevel( tree, level ) result( dx )
    ! -------------------------------------------------------------------- !
    !> Mesh to locate the point in
    type(treelmesh_type), intent(in) :: tree
    !> given level to get the size for
    integer, intent(in) :: level
    !> size of Element
    real(kind=rk) :: dx
    ! -------------------------------------------------------------------- !

    dx = tree%global%BoundingCubeLength / real(2**level, kind=rk)

  end function tem_ElemSizeLevel
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Return the size of a given treeID in the mesh by taking into account the
  !! size of the bounding cube given in the global info of the tree
  !!
  pure function tem_ElemSize( tree, treeID ) result(dx)
    ! -------------------------------------------------------------------- !
    !> Mesh to locate point in
    type(treelmesh_type), intent(in) :: tree
    !> input elements
    integer(kind=long_k), intent(in) :: treeID
    !> size of element
    real(kind=rk) :: dx
    ! -------------------------------------------------------------------- !
    integer       :: level
    ! -------------------------------------------------------------------- !

    level = tem_levelOf(TreeID)
    dx = tem_elemSizeLevel( tree, level )

  end function tem_ElemSize
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Run through all the elements, check the vertices and return the fluid
  !! bounding cube
  !!
  function tem_GetLocalBoundingCube_fromTree( tree ) result(BoundingCube)
    ! -------------------------------------------------------------------- !
    !> global mesh information
    type(treelmesh_type), intent(in) :: tree
    !> xyz coordinate for min and max of bounding cube
    real(kind=rk) :: BoundingCube(3,2)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: vrtxCoord(3,8) ! coordinates for all eight vertices
    real(kind=rk) :: minX(3) ! xyz coordinate for min of bounding cube
    real(kind=rk) :: maxX(3) ! xyz coordinate for max of bounding cube
    integer(kind=long_k) :: tTreeID
    integer :: iElem, iVert
    ! -------------------------------------------------------------------- !

    minX =  huge( minX )
    maxX = -huge( maxX )

    do iElem = 1, tree%nElems
      tTreeID = tree%treeID( iElem )
      ! Calculate coordinates of vertices
      vrtxCoord(:,:) = tem_vrtxCoordOfId( tree, tTreeID )
      do iVert = 1, 8
        ! Compare with min max
        minX(:) = min( minX(:), vrtxCoord(:, iVert ))
        maxX(:) = max( maxX(:), vrtxCoord(:, iVert ))
      enddo
    enddo
    BoundingCube(:,1) = minX(:)
    BoundingCube(:,2) = maxX(:)

  end function tem_GetLocalBoundingCube_fromTree
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Run through all the elements, check the vertices and return the fluid
  !! bounding cube
  !!
  function tem_GetLocalBoundingCube_fromSubTree( subTree, globalTree ) &
    &        result(BoundingCube)
    ! -------------------------------------------------------------------- !
    !> subTree to locate point in
    type(tem_subTree_type), intent(in) :: subTree
    !> corresponding global tree
    type(treelmesh_type), intent(in) :: globalTree
    !> xyz coordinate for min and max of bounding cube
    real(kind=rk) :: BoundingCube(3,2)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: vrtxCoord(3,8) ! coordinates for all eight vertices
    real(kind=rk) :: minX(3) ! xyz coordinate for min of bounding cube
    real(kind=rk) :: maxX(3) ! xyz coordinate for max of bounding cube
    integer(kind=long_k) :: tTreeID
    integer :: iElem, iVert
    ! -------------------------------------------------------------------- !

    ! if the subTree equals to the global shape
    if( subTree%useGlobalMesh ) then
      ! ... the effective bounding cube equals to the one of the tree
      BoundingCube = globalTree%global%effboundingCube(3,2)
    else ! if this is not the case
      minX =  huge( minX )
      maxX = -huge( maxX )

      ! loop over all elements
      do iElem = 1, subTree%nElems
        ! ... and distinguish between local shapes and other shapes
        if( subTree%useLocalMesh ) then
          ! ... copy the treeIDs directly
          tTreeID = subTree%treeID( iElem )
        else
          ! ... the treeIDs from the globalTree via the positions in map2global
          tTreeID = globalTree%treeID( subTree%map2global( iElem ))
        end if
        ! Calculate coordinates of vertices
        vrtxCoord(:,:) = tem_vrtxCoordOfId( globalTree, tTreeID )
        do iVert = 1, 8
          ! Compare with min max
          minX(:) = min( minX(:), vrtxCoord(:, iVert ))
          maxX(:) = max( maxX(:), vrtxCoord(:, iVert ))
        enddo
      enddo
      BoundingCube(:,1) = minX(:)
      BoundingCube(:,2) = maxX(:)
    end if ! useGlobalShape

  end function tem_GetLocalBoundingCube_fromSubTree
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate the real bounding box around the fluid domain
  !!
  subroutine tem_setEffBoundingBox_fromTree( tree )
    ! -------------------------------------------------------------------- !
    !> Mesh
    type(treelmesh_type), intent(inout) :: tree
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: boundingBox(3,2)
    ! -------------------------------------------------------------------- !

    boundingBox(:,:) = tem_GetRealBoundingCube( tree )

    ! Set the effective origin and length in the global tree part
    tree%global%effOrigin(:) = boundingBox(:,1)
    tree%global%effLength(:) = boundingBox(:,2) -  boundingBox(:,1)

  end subroutine tem_setEffBoundingBox_fromTree
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate the real bounding box around the fluid domain
  !!
  subroutine tem_setEffBoundingBox_fromSubTree( subTree, globalTree )
    ! -------------------------------------------------------------------- !
    !> subTree to get effective bounding cube from
    type(tem_subTree_type)   :: subTree
    !> corresponding global tree
    type(treelmesh_type), intent(in) :: globalTree
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: boundingBox(3,2)
    real(kind=rk) :: tBounding(3)
    integer :: iErr
    ! -------------------------------------------------------------------- !

    ! if the subTree equals to the global tree take over the settings
    if( subTree%useGlobalMesh )then
      subTree%global%effOrigin = globalTree%global%effOrigin
      subTree%global%effLength = globalTree%global%effLength
    else ! subTree is not equal to the global tree
      ! Calculate process-local bounding cube
      boundingBox = tem_GetLocalBoundingCube( subTree, globalTree )
      ! Exchange with neighbors
      call mpi_allreduce( boundingBox(:,1), tBounding, 3, rk_mpi, mpi_min, &
        &                 subTree%global%comm, iErr  )
      boundingBox(:,1) = tBounding
      call mpi_allreduce( boundingBox(:,2), tBounding, 3, rk_mpi, mpi_max, &
        &                 subTree%global%comm, iErr  )
      boundingBox(:,2) = tBounding

      ! Set the effective origin and length in the global subTree part
      subTree%global%effOrigin(:) = boundingBox(:,1)
      subTree%global%effLength(:) = boundingBox(:,2) -  boundingBox(:,1)
    end if

  end subroutine tem_setEffBoundingBox_fromSubTree
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> The following function provides the coordinates of the barycenter for a
  !! given treeID in the complete mesh.
  !!
  pure function tem_BaryOfId(tree,TreeID) result(bary)
    ! -------------------------------------------------------------------- !
    !> mesh information
    type(treelmesh_type), intent(in) :: tree
    !> input Element ID
    integer(kind=long_k), intent(in) :: TreeID
    !> barycenter return value
    real(kind=rk) :: bary(3)
    ! -------------------------------------------------------------------- !
    integer       :: coord(4) ! spatial index triple for a given ID
    real(kind=rk) :: dx ! size of Element
    ! -------------------------------------------------------------------- !
    coord = tem_CoordOfId(TreeID)
    dx = tem_elemSizeLevel( tree, coord(4) )
    bary = tree%global%origin + (real(coord(:3), kind=rk) + 0.5_rk)*dx

  end function tem_BaryOfId
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> The following function provides the coordinates of the origin for a
  !! given ID in the complete mesh.
  !!
  pure function tem_originOfId(tree,TreeID) result(origin)
    ! -------------------------------------------------------------------- !
    !> mesh information
    type(treelmesh_type), intent(in) :: tree
    !> input element ID
    integer(kind=long_k), intent(in) :: TreeID
    !> origin return value
    real(kind=rk) :: origin(3)
    ! -------------------------------------------------------------------- !
    integer       :: coord(4) !< spatial index triple for a given ID
    real(kind=rk) :: dx !< size of Element
    ! -------------------------------------------------------------------- !
    coord = tem_CoordOfId( TreeID )
    dx = tem_elemSizeLevel( tree, coord(4) )
    origin = tree%global%origin + real(coord(:3), kind=rk)*dx

  end function tem_originOfId
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> The following function provides the coordinates of the end for a
  !! given ID in the complete mesh. The described element lies between the
  !! origin and this end point.
  !!
  !! The end point of the element with a given ID coincides with the origin of
  !! treeID with a coordinate shift of +1 in each dimension.
  pure function tem_endOfId(tree,TreeID) result(origin)
    ! -------------------------------------------------------------------- !
    !> mesh information
    type(treelmesh_type), intent(in) :: tree
    !> input element ID
    integer(kind=long_k), intent(in) :: TreeID
    !> origin return value
    real(kind=rk) :: origin(3)
    ! -------------------------------------------------------------------- !
    integer       :: coord(4) !< spatial index triple for a given ID
    real(kind=rk) :: dx !< size of Element
    ! -------------------------------------------------------------------- !
    coord = tem_CoordOfId( TreeID )
    dx = tem_elemSizeLevel( tree, coord(4) )
    origin = tree%global%origin + real(coord(:3)+1, kind=rk)*dx

  end function tem_endOfId
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides the coordinates of the element barycenters
  !! for a set of given element coordinates on the same refinement level.
  pure function tem_BaryOfCoord(coord, nPoints, origin, dx) result(bary)
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: nPoints !< Number of points to evaluate
    integer, intent(in) :: coord(nPoints,3)
    !> spatial index triple for a given ID
    real(kind=rk), intent(in) :: origin(3) !< origin of the universe cube
    real(kind=rk), intent(in) :: dx !< size of the elements
    real(kind=rk) :: bary(nPoints,3) !< barycenter to return
    ! -------------------------------------------------------------------- !
    integer :: ii
    real(kind=rk) :: c(3)
    ! -------------------------------------------------------------------- !

    ! bary = origin + ( coord + 0.5 ) * dx
    !      = origin + coord * dx + 0.5dx
    !      = c      + coord * dx
    c(:) = origin(:) + 0.5_rk * dx
    do ii = 1, nPoints
      bary(ii,1) = c(1) + real(coord(ii,1), kind=rk) * dx
      bary(ii,2) = c(2) + real(coord(ii,2), kind=rk) * dx
      bary(ii,3) = c(3) + real(coord(ii,3), kind=rk) * dx
    end do

  end function tem_BaryOfCoord
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate all eight vertex coordinates of a given tree ID
  !!
  !! the numbering of the vertices is according to
  !! [tem_param_module:childPosition]
  pure function tem_vrtxCoordOfId( tree, treeID) result(coord)
    ! -------------------------------------------------------------------- !
    !> complete tree for info about dimensions
    type(treelmesh_type), intent(in) :: tree
    !> input element ID
    integer(kind=long_k), intent(in) :: treeID
    !> all vertices coordinates function return value
    real(kind=rk) :: coord(3,8)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: elemBary(3) !< barycenter of current element
    real(kind=rk) :: length  !< element size
    integer :: x(4), iCoord, iVrtx
    ! -------------------------------------------------------------------- !

    x = tem_coordOfId( treeID )
    ! Get element size
    length = tree%global%BoundingCubeLength / real( 2**x(4), kind=rk)

    ! Get Barycenter
    elemBary = tem_BaryOfId( tree, treeID )

    ! ... and calculate from there the vertex coordinates
    do iVrtx = 1, 8
      do iCoord = 1, 3
        coord( iCoord, iVrtx ) = elemBary( iCoord )                         &
          &                    + length*0.5_rk*childPosition( iVrtx, iCoord )
      enddo
    end do

  end function tem_vrtxCoordOfId
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Identify all possible children local ids for each of the 27 direct
  !! neighbors results are saved in the ElemList
  !!
  !! How to use this information:
  !! There are two main ways of doing so.
  !!
  !! 1.) Use eligibleChildren array as a subset array of directChildren
  !! Let's say, you have a neighbor on the right border (=East border),
  !! which is on a coarser level.
  !! But you need the treeIDs of the neighbor on the same level.
  !! You then call tem_directChildren( coarseNEighborTreeID )
  !! to get all the children.
  !! Then, you can access the treeIDs on the east border (that is the
  !! treeIDs of the children on the West!! border of the coarser neighbor)
  !! by using the eligibleChildren as a subset to access the tem_directChildren
  !! childrenIDs( eligible_child( iChild ))
  !!
  !! 2.) Use elgibleChildren as an offset from the lower, left, bottom
  !! child tree ID of the parent.
  !!
  !! This is routine is Morton curve specific
  !! HK: and solver specific?
  !!
  pure subroutine tem_eligibleChildren(eligible_child, Direction)
    ! -------------------------------------------------------------------- !
    !> Candidate children, which might be considered as neighbors
    integer, intent(out), allocatable :: eligible_child(:)
    !> In which direction to search for neighbors
    integer, intent(in) :: direction
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    select case ( Direction )

    ! direct neighbor directions: 4 children
    case( q__W)
      allocate(eligible_child(4))
      eligible_child(1) = 1
      eligible_child(2) = 3
      eligible_child(3) = 5
      eligible_child(4) = 7
    case( q__E)
      allocate(eligible_child(4))
      eligible_child(1) = 2
      eligible_child(2) = 4
      eligible_child(3) = 6
      eligible_child(4) = 8
    case( q__S)
      allocate(eligible_child(4))
      eligible_child(1) = 1
      eligible_child(2) = 2
      eligible_child(3) = 5
      eligible_child(4) = 6
    case( q__N)
      allocate(eligible_child(4))
      eligible_child(1) = 3
      eligible_child(2) = 4
      eligible_child(3) = 7
      eligible_child(4) = 8
    case( q__B)
      allocate(eligible_child(4))
      eligible_child(1) = 1
      eligible_child(2) = 2
      eligible_child(3) = 3
      eligible_child(4) = 4
    case( q__T)
      allocate(eligible_child(4))
      eligible_child(1) = 5
      eligible_child(2) = 6
      eligible_child(3) = 7
      eligible_child(4) = 8

    ! edge neighbor directions: 2 children
    case( q_SW)
      allocate(eligible_child(2))
      eligible_child(1) = 1
      eligible_child(2) = 5
    case( q_NW)
      allocate(eligible_child(2))
      eligible_child(1) = 3
      eligible_child(2) = 7
    case( q_SE)
      allocate(eligible_child(2))
      eligible_child(1) = 2
      eligible_child(2) = 6
    case( q_NE)
      allocate(eligible_child(2))
      eligible_child(1) = 4
      eligible_child(2) = 8

    case( q_BW)
      allocate(eligible_child(2))
      eligible_child(1) = 1
      eligible_child(2) = 3
    case( q_TW)
      allocate(eligible_child(2))
      eligible_child(1) = 5
      eligible_child(2) = 7
    case( q_BE)
      allocate(eligible_child(2))
      eligible_child(1) = 2
      eligible_child(2) = 4
    case( q_TE)
      allocate(eligible_child(2))
      eligible_child(1) = 6
      eligible_child(2) = 8

    case( q_BS)
      allocate(eligible_child(2))
      eligible_child(1) = 1
      eligible_child(2) = 2
    case( q_TS)
      allocate(eligible_child(2))
      eligible_child(1) = 5
      eligible_child(2) = 6
    case( q_BN)
      allocate(eligible_child(2))
      eligible_child(1) = 3
      eligible_child(2) = 4
    case( q_TN)
      allocate(eligible_child(2))
      eligible_child(1) = 7
      eligible_child(2) = 8

   ! corner neighbor directions: only one child
    case( qBSW)
      allocate(eligible_child(1))
      eligible_child(1) = 1
    case( qTNE)
      allocate(eligible_child(1))
      eligible_child(1) = 8
    case( qTNW)
      allocate(eligible_child(1))
      eligible_child(1) = 7
    case( qBSE)
      allocate(eligible_child(1))
      eligible_child(1) = 2
    case( qBNW)
      allocate(eligible_child(1))
      eligible_child(1) = 3
    case( qTSE)
      allocate(eligible_child(1))
      eligible_child(1) = 6
    case( qTSW)
      allocate(eligible_child(1))
      eligible_child(1) = 5
    case( qBNE)
      allocate(eligible_child(1))
      eligible_child(1) = 4

    end select ! select ( Direction )

  end subroutine tem_eligibleChildren
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Recursive routine to find all actual (eligible) leave nodes in the
  !! local partition for a given treeID.
  !! Alternatively use tem_findPath, which uses precomputed paths in the tree
  !! and should speed up the search (at the expense of storing the paths
  !! beforehand).
  !!
  recursive subroutine tem_findElement( TreeID, eligible_child, ElemList, &
    &                                   treeIDlist, nElems, Part_First,   &
    &                                   Part_Last, otherLevel             )
    ! -------------------------------------------------------------------- !
    !> TreeID to find in the array of Elements
    integer(kind=long_k) :: TreeID
    !> Candidate childs, which might be considered as neighbors
    integer :: eligible_child(:)
    !> linked list of resulting elements building the neighbor
    type(tem_longList), pointer :: ElemList
    !> number of elements in list
    integer, intent(in)  :: nElems
    !> array of treeIDs
    integer(kind=long_k), intent(in)  :: treeIDlist(nElems)
    !> parts first entry
    integer(kind=long_k), intent(in) :: Part_First(:)
    !> parts last entry
    integer(kind=long_k), intent(in) :: Part_Last(:)
    !> entry is on another level
    logical,optional,intent(inout) :: otherLevel
    ! -------------------------------------------------------------------- !
    integer(kind=long_k) :: pos
    integer :: i
    integer(kind=long_k) :: childID, off
    ! -------------------------------------------------------------------- !

    ! binary search of the ID in the array of actual present elements
    ! Return pos < 0 if
    pos = tem_PosOfId(TreeID, treeIDlist)

    ! If the neighbor is on a level higher than myself, it should be
    ! delivered by binary search

    if (pos > 0 ) then
      ! Element actually exists, append it to the list
       call append(ElemList, pos)
    else
      if (pos < 0) then
        ! Element is a virtual neighbor, look for childs
        if( present( otherLevel ) ) otherLevel = .true.
        off = TreeID*8
        do i=1,size(eligible_child)
          childID = off + eligible_child(i)
          call tem_findElement( childID, eligible_child, ElemList,        &
            &                   treeIDlist, nElems, Part_First, Part_Last )
        end do
      end if
    end if

  end subroutine tem_findElement
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Recursive routine to find all actual (eligible) leave nodes in the
  !! local partition for a given treeID.
  !!@todo HK: when doing this for the complete domain it would probably better
  !!          just compute the ID on the finest level along with the information
  !!          on the level of each leaf, in order to speed up things a little.
  !!          This way comparison would be just a simple integer comparison.
  !!
  recursive subroutine tem_findPath( Path, eligible_child, ElemList, &
    &                                Pathlist, nElems, otherLevel    )
    ! -------------------------------------------------------------------- !
    !> Path to the leaf to find in the array of Elements
    type(tem_path_type), intent(in) :: Path
    !> Candidate childs, which might be considered as neighbors
    integer, intent(in) :: eligible_child(:)
    !> linked list of resulting elements building the neighbor
    type(tem_longList), pointer :: ElemList
    !> number of elements in list
    integer, intent(in)  :: nElems
    !> array of paths
    type(tem_path_type), intent(in)  :: pathlist(nElems)
    !> entry is on another level
    logical,optional,intent(inout) :: otherLevel
    ! -------------------------------------------------------------------- !
    integer(kind=long_k) :: pos
    integer :: i
    integer(kind=long_k) :: off
    type(tem_path_type) :: childPath
    ! -------------------------------------------------------------------- !

    ! binary search of the Path in the array of actual present elements
    ! Return pos < 0 if
    pos = tem_PosOfPath(Path, Pathlist)

    ! If the neighbor is on a level higher than myself, it should be
    ! delivered by binary search

    if (pos > 0 ) then
      ! Element actually exists, append it to the list
       call append(ElemList, pos)
    else if (pos < 0) then
      ! Element is a GhostFromFiner, look for childs
      if( present( otherLevel ) ) otherLevel = .true.
      off = Path%Node(1)*8
      do i=1,size(eligible_child)
        childPath = tem_pathOf(off + eligible_child(i))
        call tem_findPath( childPath, eligible_child, ElemList, &
          &                Pathlist, nElems, otherLevel)
      end do
    else ! pos == 0
    end if

  end subroutine tem_findPath
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This subroutine does a binary search on a given (sparse) list of elements.
  !! The result is the position of the given tree ID in the list,
  !! 0 if no corresponding node is found,
  !! or the negative of the found ID, if it is a virtual node.
  !!
  pure function tem_PosOfId(sTreeID, treeIDlist, lower, upper) result(IdPos)
    ! -------------------------------------------------------------------- !
    !> tree ID to search for
    integer(kind=long_k), intent(in)  :: sTreeID
    !> List to search in
    integer(kind=long_k), intent(in)  :: treeIDlist(:)
    !> lowerbound of search interval
    integer, intent(in), optional     :: lower
    !> upperbound of search interval
    integer, intent(in), optional     :: upper
    !> position of sTreeID in the list of elements
    integer                           :: IdPos
    ! -------------------------------------------------------------------- !
    integer :: lb, ub
    integer :: middleSearch
    type(tem_path_type) :: searched
    type(tem_path_type) :: current
    integer :: pathRelation
    ! -------------------------------------------------------------------- !

    if (present(lower)) then
      lb = lower
    else
      lb = lbound(treeIDList,1)
    end if
    if (present(upper)) then
      ub = upper
    else
      ub = ubound(treeIDList,1)
    end if

    !> Build the path to the searched TreeID from the leaf to the root.
    searched = tem_PathOf(sTreeID)

    ! Start the Binary search for the neighbor elements
    binSearchLoop: do
      middleSearch = (lb + ub) / 2
      ! Build the path to the currently investigated element from leaf to root.
      current = tem_PathOf(treeIDList(middleSearch))

      pathRelation = tem_PathComparison(searched, current)

      if ((pathRelation == 0) .or. (lb >= ub)) then
        ! Leave the loop, if element has been found, or this
        ! was the last element to investigate.
        exit binSearchLoop
      else
        halves: if (pathRelation == 1) then
          ! Continue the search in the higher half, as the looked up element is
          ! to small.
          lb = min(middleSearch + 1, ub)
        else
          ! Continue search in the lower half, as the looked up element is to
          ! large.
          ub = max(middleSearch - 1, lb)
        end if halves
      end if
    end do binSearchLoop

    if (pathRelation == 0) then
      if (current%Level <= searched%Level) then
        ! The found ID is actually a leaf
        IdPos = middleSearch
      else
        ! The found ID is a parent of the searched
        ! virtual treeID
        IdPos = -middleSearch
      end if
    else
      IdPos = 0 ! no matching element found.
    end if

  end function tem_PosOfId
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Find the position of a specific path in the list of all paths.
  !!
  pure function tem_PosOfPath(sPath, Pathlist, lower, upper) result(IdPos)
    ! -------------------------------------------------------------------- !
    !> Specific path to search.
    type(tem_path_type), intent(in) :: sPath
    !> List of paths for all elements to search in.
    type(tem_path_type), intent(in) :: Pathlist(:)
    !> Possibly only search on a subinterval, starting at PathList(lower:)
    integer, intent(in), optional :: lower
    !> Possibly only search on a subinterval, ending at PathList(:upper)
    integer, intent(in), optional :: upper
    !> Position where sPath is found in Pathlist
    integer  :: IdPos
    ! -------------------------------------------------------------------- !
    integer :: lb, ub
    integer :: middleSearch
    type(tem_path_type) :: current
    integer(kind=long_k) :: pathRelation
    integer :: maxPathLen
    ! -------------------------------------------------------------------- !

    if (present(lower)) then
      lb = lower
    else
      lb = lbound(PathList,1)
    end if
    if (present(upper)) then
      ub = upper
    else
      ub = ubound(PathList,1)
    end if

    ! Start the Binary search for the neighbor elements
    binSearchLoop: do
      middleSearch = (lb + ub) / 2

      ! Build the path to the currently investigated element from leaf to root.
      current = PathList(middleSearch)

      ! pathRelation = tem_PathComparison(sPath, current)
      ! Inlined tem_PathComparison:
      maxPathLen = min(sPath%level-1, current%level-1)

      pathRelation = sPath%Node(sPath%level - maxPathLen)   &
        &          - current%Node(current%level - maxPathLen)

      if ((pathRelation == 0_long_k) .or. (lb >= ub)) then
        ! Leave the loop, if element has been found, or this
        ! was the last element to investigate.
        exit binSearchLoop
      else
        halves: if (pathRelation > 0_long_k) then
          ! Continue the search in the higher half, as the looked up element is
          ! to small.
          lb = min(middleSearch + 1, ub)
        else
          ! Continue search in the lower half, as the looked up element is to
          ! large.
          ub = max(middleSearch - 1, lb)
        end if halves
      end if

    end do binSearchLoop

    if (pathRelation == 0_long_k) then
      if (current%Level <= sPath%Level) then
        ! The found ID is actually a leaf
        IdPos = middleSearch
      else
        ! The found ID is a child of the searched
        ! virtual treeID
        IdPos = -middleSearch
      end if
    else
      IdPos = 0 ! no matching element found.
    end if

  end function tem_PosOfPath
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Compare the incoming discrete vector against a set of prevailing vectors
  !! and return the found closest prevailing integer vector
  !!
  !! A set of real numbered unit vectors (compareVector) must be given.
  !! The integer numbered vector is compared against those and the closest
  !! integer vector overwrites the original vector.
  subroutine tem_determine_discreteVector( vector, compareVector, angle )
    ! -------------------------------------------------------------------- !
    !> The given vector, will be set to the vector in compareVector
    !! that is best aligned to it.
    integer, intent(inout) :: vector(:)
    !> Set of unit vectors to select from against.
    !! size is (3, nVectors)
    real( kind=rk ), intent(in) :: compareVector(:,:)
    !> angle between vectors
    real( kind=rk ), intent(out), optional :: angle
    ! -------------------------------------------------------------------- !
    integer :: iDir, bestIndex
    real(kind=rk) :: length, rv(3)
    real(kind=rk) :: rVectIn(3), angleLoc, dotProduct, maxplen, min_nz_comp
    ! -------------------------------------------------------------------- !

    angleLoc = 0.0_rk

    ! Only need to do anything for non-zero vectors.
    if ( maxval(abs(vector)) > 0 ) then
      rv = real(vector, kind=rk)

      length = sqrt(rv(1)**2 + rv(2)**2 + rv(3)**2)

      ! Normalize vector
      rVectIn = rv / length

      ! Keep track of the best fitting index.
      bestIndex = 0

      ! Keep track of the maximal projected length.
      ! As we are looking at two unit vectors, the smallest
      ! projected length should be -1 (opposing directions).
      ! Thus, starting with -2 here as a save comparison that is
      ! smaller than all possible outcomes from compareVectors.
      maxplen = -2.0_rk

      ! Go over all vectors in the compareVector
      do iDir = 1, size(compareVector, 2)
        ! Go through all prevailing directions and check the minimum
        dotProduct = rVectIn(1)*compareVector( 1, iDir ) &
          &        + rVectIn(2)*compareVector( 2, iDir ) &
          &        + rVectIn(3)*compareVector( 3, iDir )

        dotProduct = min( max(dotProduct, -1.0_rk), 1.0_rk )
        if ( dotProduct > maxplen ) then
          maxplen = dotProduct
          bestIndex = iDir
          if ( maxplen .feq. 1._rk ) EXIT
        end if
      end do

      ! Align the vector to the closest matching prevailing direction from
      ! the set of compareVector:

      ! Use the smallest nonzero component to scale the unit vector to
      ! components close to integer numbers.
      min_nz_comp = abs( minval( pack(compareVector(:, bestIndex),     &
        &                             abs(compareVector(:, bestIndex)) &
        &                             > epsilon(0.0_rk))               &
        &                      )                                       )
      ! Use the nearest integers to represent the projected vector.
      vector = nint(compareVector(:, bestIndex)/min_nz_comp)

      ! Angle between the original vector and the projected direction.
      angleLoc = acos(maxplen)
    end if ! not zero-vector

    if (present(angle)) angle = angleLoc

  end subroutine tem_determine_discreteVector
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function returns the bi-linearly interpolated values from the four
  !! source points to the target position located at targetCoord.
  !! The source points are arranged in a square from (0,0)x(1,1)
  !! The order of the source points are according to the morton curve
  !!  1      2      3      4
  !! (0,0); (1,0); (0,1); (1,1)
  !!
  function tem_intp_bilinear_scalar( srcVal, targetCoord ) result( phi )
    ! -------------------------------------------------------------------- !
    !> source values of the square corners
    real(kind=rk), intent(in) :: srcVal(4)
    !> interpolation location within the square
    real(kind=rk), intent(in) :: targetCoord(2)
    !> interpolated value
    real(kind=rk) :: phi
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: phi_north, phi_south
    ! -------------------------------------------------------------------- !
    ! Linear interpolation on the north nodes
    phi_north = (1._rk - targetCoord(1)) * srcVal(3) + targetCoord(1)*srcVal(4)
    ! Linear interpolation on the south nodes
    phi_south = (1._rk - targetCoord(1)) * srcVal(1) + targetCoord(1)*srcVal(2)
    ! Linear interpolation on the obtained north and south values
    phi       = (1._rk - targetCoord(2)) * phi_south + targetCoord(2)*phi_north

  end function tem_intp_bilinear_scalar
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function returns the bi-linearly interpolated values from the four
  !! source points to the target position located at targetCoord.
  !! The source points are arranged in a square from (0,0)x(1,1)
  !! The order of the source points are according to the morton curve
  !!  1      2      3      4
  !! (0,0); (1,0); (0,1); (1,1)
  !!
  function tem_intp_bilinear_vec( srcVal, targetCoord, nVals ) result( phi )
    ! -------------------------------------------------------------------- !
    !> number of values
    integer, intent(in) :: nVals
    !> source values of the square corners
    real(kind=rk), intent(in) :: srcVal(nVals,4)
    !> interpolation location within the square
    real(kind=rk), intent(in) :: targetCoord(2)
    !> interpolated values
    real(kind=rk) :: phi(nVals)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: phi_north, phi_south
    integer :: iVal
    ! -------------------------------------------------------------------- !
    do iVal = 1, nVals
      ! Linear interpolation on the north nodes
      phi_north = (1._rk - targetCoord(1)) * srcVal(iVal,3) &
        &         + targetCoord(1) * srcVal(iVal,4)
      phi_south = (1._rk - targetCoord(1)) * srcVal(iVal,1) &
        &         + targetCoord(1)*srcVal(iVal,2)
      phi(iVal) = (1._rk - targetCoord(2)) * phi_south &
        &         + targetCoord(2)*phi_north
    end do

  end function tem_intp_bilinear_vec
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function returns the tri-linearly interpolated values from the eight
  !! source points to the target position located at targetCoord.
  !! The source points are arranged in a square from (0,0,0)x(1,1,1)
  !! The order of the source points are according to the morton curve
  !!  1        2        3        4
  !! (0,0,0); (1,0,0); (0,1,0); (1,1,0)
  !!  5        6        7        8
  !! (0,0,1); (1,0,1); (0,1,1); (1,1,1)
  !!
  function tem_intp_trilinear_scalar( srcVal, targetCoord ) result( phi )
    ! -------------------------------------------------------------------- !
    !> source values of the square corners
    real(kind=rk), intent(in) :: srcVal(8)
    !> interpolation location within the square
    real(kind=rk), intent(in) :: targetCoord(3)
    !> interpolated value
    real(kind=rk) :: phi
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: phi_northFront, phi_southFront
    real(kind=rk) :: phi_northBack, phi_southBack
    real(kind=rk) :: phi_front, phi_back
    ! -------------------------------------------------------------------- !
    ! Linear interpolation on the cube front side (z = 0 )
    phi_northFront = (1._rk - targetCoord(1)) * srcVal(3) &
      &              + targetCoord(1) * srcVal(4)
    phi_southFront = (1._rk - targetCoord(1)) * srcVal(1) &
      &              + targetCoord(1) * srcVal(2)
    ! Linear interpolation on the cube back side (z = 1 )
    phi_northBack = (1._rk - targetCoord(1)) * srcVal(7) &
      &             + targetCoord(1) * srcVal(8)
    phi_southBack = (1._rk - targetCoord(1)) * srcVal(5) &
      &             + targetCoord(1) * srcVal(6)
    ! Linear interpolation on the cube front side (z = 0 )
    phi_front = (1._rk - targetCoord(2)) * phi_southFront &
      &         + targetCoord(2) * phi_northFront
    ! Linear interpolation on the cube back side (z = 1 )
    phi_back = (1._rk - targetCoord(2)) * phi_southBack &
      &        + targetCoord(2) * phi_northBack
    phi = (1._rk - targetCoord(3)) * phi_front &
      &   + targetCoord(3)*phi_back

  end function tem_intp_trilinear_scalar
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function returns the tri-linearly interpolated values from the eight
  !! source points to the target position located at targetCoord.
  !! The source points are arranged in a square from (0,0,0)x(1,1,1)
  !! The order of the source points are according to the morton curve
  !!  1        2        3        4
  !! (0,0,0); (1,0,0); (0,1,0); (1,1,0)
  !!  5        6        7        8
  !! (0,0,1); (1,0,1); (0,1,1); (1,1,1)
  !!
  function tem_intp_trilinear_vec( srcVal, targetCoord, nVals ) result( phi )
    ! -------------------------------------------------------------------- !
    !> number of values
    integer, intent(in) :: nVals
    !> source values of the square corners
    real(kind=rk), intent(in) :: srcVal(nVals,8)
    !> interpolation location within the square
    real(kind=rk), intent(in) :: targetCoord(3)
    !> interpolated value
    real(kind=rk) :: phi(nVals)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: phi_northFront, phi_southFront
    real(kind=rk) :: phi_northBack,  phi_southBack
    real(kind=rk) :: phi_front, phi_back
    integer :: iVal
    ! -------------------------------------------------------------------- !
    do iVal = 1, nVals
      ! Linear interpolation on the north nodes
      phi_northFront = (1._rk - targetCoord(1)) * srcVal(iVal,3) &
        &              + targetCoord(1) * srcVal(iVal,4)
      phi_southFront = (1._rk - targetCoord(1)) * srcVal(iVal,1) &
        &              + targetCoord(1) * srcVal(iVal,2)
      ! Linear interpolation on the cube back side (z = 1 )
      phi_northBack = (1._rk - targetCoord(1)) * srcVal(iVal,7) &
        &             + targetCoord(1) * srcVal(iVal,8)
      phi_southBack = (1._rk - targetCoord(1)) * srcVal(iVal,5) &
        &             + targetCoord(1) * srcVal(iVal,6)
      ! Linear interpolation on the cube front side (z = 0 )
      phi_front = (1._rk - targetCoord(2)) * phi_southFront &
        &         + targetCoord(2) * phi_northFront
      ! Linear interpolation on the cube back side (z = 1 )
      phi_back  = (1._rk - targetCoord(2)) * phi_southBack &
        &         + targetCoord(2) * phi_northBack
      phi(iVal) = (1._rk - targetCoord(3)) * phi_front &
        &         + targetCoord(3)*phi_back
    enddo

  end function tem_intp_trilinear_vec
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine builds mapping from tree%treeID to given property data array
  !! like boundary_ID or qVal
  subroutine tem_build_treeToProp_pointer( treeToProp, nElems,       &
    &                                      ElemPropertyBits, prp_bit )
    ! -------------------------------------------------------------------- !
    !> mesh elements
    integer :: nElems
    !> Pointer from tree to property data to be filled by this routine
    integer, intent(out) :: treeToProp(nElems)
    !> Elements Property Bits
    integer(kind=long_k), intent(in) :: ElemPropertyBits(nElems)
    !> property bit
    integer, intent(in) :: prp_bit
    ! -------------------------------------------------------------------- !
    character(len=labelLen) :: buffer
    integer :: iElem, counter
    ! -------------------------------------------------------------------- !

    select case(prp_bit)
    case (prp_hasBnd)
      buffer = 'boundary_ID'
    case (prp_hasQVal)
      buffer = 'qVal'
    end select

    write(logUnit(7),*) 'Building map from tree to '//trim(buffer)

    counter = 0
    do iElem = 1, nElems
      if ( btest(ElemPropertyBits( iElem ), prp_bit ) ) then
        counter = counter + 1
        treeToProp(iElem) = counter
      else
        treeToProp(iElem) = -1
      end if
    end do

  end subroutine tem_build_treeToProp_pointer
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate the real bounding box around the fluid domain
  !! and return only to root (efficiency reasons)
  !!
  function tem_GetRealBoundingCube( tree ) result( boundingCube )
    ! -------------------------------------------------------------------- !
    !>
    type(treelmesh_type), intent(in) :: tree
    !>
    real(kind=rk) :: boundingCube(3,2)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: tBounding(3)
    integer :: iErr
    ! -------------------------------------------------------------------- !
    ! Calculate process-local bounding cube
    boundingCube = tem_GetLocalBoundingCube( tree )
    ! Exchange with neighbors
    call mpi_reduce( boundingCube(:,1), tBounding, 3, rk_mpi, mpi_min, &
      &              0, tree%global%comm, iErr  )
    boundingCube(:,1) = tBounding
    call mpi_reduce( boundingCube(:,2), tBounding, 3, rk_mpi, mpi_max, &
      &              0, tree%global%comm, iErr  )
    boundingCube(:,2) = tBounding

  end function tem_GetRealBoundingCube
  ! ************************************************************************ !

end module tem_geometry_module
! **************************************************************************** !
