! Copyright (c) 2013-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> author: Simon Zimny
!! This module provides the functionality to read surface information from file
!! (stl, ...) and stores them in a surfaceData datatype.
!!
!!
module tem_surfaceData_module

!$ use omp_lib

  ! include treelm modules
  use env_module,               only: rk, long_k, PathLen
  use tem_aux_module,           only: tem_abort
  use tem_logging_module,       only: logUnit, tem_last_lu
  use tem_stlb_io_module,       only: tem_read_stlb, tem_size_stlb
  use treelmesh_module,         only: treelmesh_type
  use tem_dyn_array_module,     only: dyn_longArray_type,                      &
    &                                 init, append, destroy
  use tem_grow_array_module,    only: grw_real2darray_type, grw_intArray_type, &
    &                                 init, append, placeAt, destroy
  use tem_geometry_module,      only: tem_CoordOfReal
  use tem_topology_module,      only: tem_IdOfCoord, tem_PathComparison,       &
    &                                 tem_path_type, tem_PathOf
  use tem_time_module,          only: tem_time_type
  use tem_spacetime_fun_module, only: tem_spacetime_fun_type, tem_spacetime_for
  use tem_property_module,      only: prp_hasIBM
  use tem_math_module,          only: cross_product3D
  use tem_construction_module,  only: tem_levelDesc_type, tem_treeIDinTotal
  use tem_element_module,       only: eT_fluid, eT_halo, eT_ghostFromCoarser, &
    &                                 eT_ghostFromFiner
  use tem_timeControl_module,   only: tem_timeControl_type,              &
    &                                 tem_timeControl_load,              &
    &                                 tem_timeControl_dump,              &
    &                                 tem_timeControl_out
  use tem_comm_env_module,      only: tem_comm_env_type

  ! include aotus modules
  use aot_table_module, only: aot_table_open, aot_table_close, aot_table_length
  use aotus_module,     only: aot_get_val, flu_State

  implicit none

  private

  public :: tem_surfData_type
  public :: tem_load_surfData
  public :: tem_readAndUnify_surfData
  public :: tem_init_surfData
  public :: tem_update_surfPos
  public :: tem_calcTriaAreas
  public :: tem_freeSurfData

  type tem_surfaceData_stlHead_type
    !> filename of the data to be read from
    character(len=PathLen) :: filename
  end type tem_surfaceData_stlHead_type

  type tem_parentIDs_type
    !> levelwise pointers to the parent eulerian elements in the levelDesc
    !! size: pointCoords%nVals
    ! @todo: SZ: Check the size of this array!!
    integer, allocatable :: ptrs(:)
  end type tem_parentIDs_type

  !> Datatype to store the surface information in. The surface data consists of
  !! an array of unique points (XYZ coordinates) and their connectivity list
  !! (triangles).
  type tem_surfData_type
    !> data (filename) for the surface data header
    type( tem_surfaceData_stlHead_type ), allocatable :: stlHead(:)
    !> output prefix
    character(len=PathLen) :: outprefix
    !> dump min and max force to a seperate file (debug output)
    logical :: dumpForce
    !> time control type for controlling the dumping of the stl file
    type(tem_timeControl_type) :: timeControl
    !> number of unique point coordinates
    integer :: nUniquePoints_total
    !> linearized array of point coordinates (X,Y,Z)
    !! the coordinates are stored one after another
    !!        ---------------------------
    !!        | X1,Y1,Z1, ... , Xn,Yn,Zn|
    !!        ---------------------------
    !! size: 3*nUniquePoints_total
    real(kind=rk), allocatable :: pointCoords(:)
    !> array of surface areas attached to this point
    real(kind=rk), allocatable :: surfArea(:)
    !> array of levelwise pointers to the parent eulerian elements
    !! of the lagrangian points in the levelDesc (size: nLevels)
    type( tem_parentIDs_type ), allocatable :: parentIDs(:)
    !> connectivity array of the points
    !! size: 3, nTrias
    integer, allocatable :: trias(:,:)
    !> total number of triangles stored
    integer :: nTrias
    !> backup for linearized array of point coordinates (X,Y,Z)
    !! needed for defining offsets based on the initial position
    !! the coordinates are stored one after another
    !!        ---------------------------
    !!        | X1,Y1,Z1, ... , Xn,Yn,Zn|
    !!        ---------------------------
    !! size: 3*nUniquePoints_total
    real(kind=rk), allocatable :: backPointCoords(:)
  end type tem_surfData_type

contains


! ****************************************************************************** !
  !> @todo: SZ: Add parts for catching aotus errors
  !!
  !!
  subroutine tem_load_surfData( me, conf, sd_handle )
    ! ---------------------------------------------------------------------------
    !> datatype to store the surface information
    type( tem_surfData_type ), intent(inout) :: me
    !> handle of the lua config file
    type( flu_state ) :: conf
    !> handle for the surfaceData table
    integer :: sd_handle
    ! ---------------------------------------------------------------------------
    integer :: iError
    integer :: stlFile_handle
    integer :: stlDump_handle
    integer :: nSTLFiles
    integer :: iSTLFile
    ! ---------------------------------------------------------------------------

    ! try to open the filename identifier as a table
    call aot_table_open( L       = conf,                                       &
      &                  parent  = sd_handle,                                  &
      &                  thandle = stlFile_handle,                             &
      &                  key     = 'stlfiles')

    ! if there is no table named filename ...
    if( stlFile_handle == 0 ) then
      ! ... close the table again
      call aot_table_close( L       = conf,                                    &
        &                   thandle = stlFile_handle )
      ! allocate the stlHead with 1
      allocate( me%stlHead( 1 ))
      ! ... and read the filename directly
      call aot_get_val( L       = conf,                                        &
        &               thandle = sd_handle,                                   &
        &               val     = me%stlHead(1)%filename,                      &
        &               ErrCode = iError,                                      &
        &               key     = 'stlfiles' )

    else ! there is a table
      ! ... get the table length
      nSTLFiles = aot_table_length( L       = conf,                            &
        &                           thandle = stlFile_handle )
      ! ... allocate the list of stlHeads
      allocate( me%stlHead( nSTLFiles ))
      ! ... and read them one by one
      do iSTLFile = 1, nSTLFiles
        call aot_get_val( L       = conf,                                      &
          &               thandle = stlFile_handle,                            &
          &               val     = me%stlHead(iSTLFile)%filename,             &
          &               ErrCode = iError,                                    &
          &               pos     = iSTLFile )
      end do
      ! ... and close the table in the end
      call aot_table_close( L       = conf,                                    &
        &                   thandle = stlFile_handle)
    end if


    ! open the stl dump table
    call aot_table_open( L       = conf,                                       &
      &                  parent  = sd_handle,                                  &
      &                  thandle = stlDump_handle,                             &
      &                  key     = 'dump_stl')

    ! if the table exists ...
    if( stlFile_handle /= 0 ) then

      ! read the output path
      call aot_get_val( L       = conf,                                        &
        &               thandle = stlDump_handle,                              &
        &               val     = me%outprefix,                                &
        &               ErrCode = iError,                                      &
        &               key     = 'outprefix',                                 &
        &               default = '' )

      ! whether to dump the force values or not
      call aot_get_val( L       = conf,                                        &
        &               thandle = stlDump_handle,                              &
        &               val     = me%dumpForce,                                &
        &               ErrCode = iError,                                      &
        &               key     = 'dumpForce',                                 &
        &               default = .false. )

      ! read the time control information for dumping the stls
      ! load time control to output tracking
      call tem_timeControl_load( conf           = conf,                        &
        &                        parent         = stlDump_handle,              &
        &                        me             = me%timeControl )
      write(logUnit(2),*) 'Writing stls using the prefix: '// trim(me%outprefix)
      write(logUnit(2),*) 'at the following timings: '
      call tem_timeControl_dump(me%timeControl, logUnit(2))
    end if

    call aot_table_close( L       = conf,                                      &
      &                   thandle = stlDump_handle )

  end subroutine tem_load_surfData
! ****************************************************************************** !


! ****************************************************************************** !
  !> This routine reads the surface data from a set of stl files and stores
  !! it in the surfaceData_type.
  !!
  subroutine tem_readAndUnify_surfData( me, useInitPos )
    ! ---------------------------------------------------------------------------
    !> datatype to store the surface information
    type( tem_surfData_type ), intent(inout) :: me
    !> shall the initial points be stored and used for updating the points
    !! later on ???
    logical, optional, intent(in) :: useInitPos
    ! ---------------------------------------------------------------------------
    ! number of points to be read from the mesh
    integer :: nPoints( size( me%stlHead ))
    ! number of triangles to be read from the mesh
    integer :: nTrias( size( me%stlHead ))
    ! total number of points to be read from file
    integer :: nPoints_total
    ! tmp point coordinates
    real(kind=rk), allocatable :: tmp_pointCoords(:,:)
    ! Offset in the temporary pointCoords list for multiple STLs
    integer :: nodeOffset
    ! Offset in the triangle list for multiple STLs
    integer :: triOffset
    integer :: iFile
    integer :: ierr
    logical :: tmp_useInitPos
    ! min and max of all X,Y,Z coordinates
    real(kind=rk) :: minX, minY, minZ, maxX, maxY, maxZ
    ! ---------------------------------------------------------------------------
    if( present( useInitPos ))then
      tmp_useInitPos = useInitPos
    else
      tmp_useInitPos = .false.
    end if

    write(logUnit(10),*) " Reading in STL Headers ..."
    ! Read headerfiles to determine number of nodes and tris
    do iFile = 1, size(me%stlHead)
      call tem_size_stlb( filename = me%stlHead(iFile)%filename,               &
        &                 nNodes   = nPoints( iFile ),                         &
        &                 nTris    = nTrias( iFile ))
    end do

    nPoints_total = sum( nPoints( : ))
    me%nTrias     = sum(  nTrias( : ))

    write(logUnit(2),*) " Total number of triangles: ", me%nTrias
    write(logUnit(2),*) " Total number of points:    ", nPoints_total

    ! allocate the temporary array of coordinates
    allocate(tmp_pointCoords( 3, nPoints_total ))
    ! allocate the global number of triangles (this is constant regardless of
    ! unification)
    allocate(me%trias( 3, me%nTrias ))

    write(logUnit(10),*) " Reading in STL Files ..."

    ! Read in node values from STL files
    nodeOffset = 1
    triOffset = 1
    do iFile = 1, size(me%stlHead)
      ! read in the nodes and triangles to the right positions in the temporary
      ! pointCoords array and in the global triangle array
      call tem_read_stlb( filename   = me%stlHead( iFile )%filename,           &
        &                 nNodesRead = nPoints( iFile ),                       &
        &                 nTrisRead  = nTrias( iFile ),                        &
        &                 nodes      = tmp_pointCoords( 1:3,                   &
        &                           nodeOffset:nodeOffset+nPoints( iFile )-1), &
        &                 tri_node   = me%trias(1:3,                           &
        &                              triOffset:triOffset+nTrias( iFile )-1), &
        &                 ierror     = ierr)

      if( iErr /= 0) then
        write(logUnit(0),*) "An error appeared when reading the surface data " &
          &               //"file "//trim( me%stlHead( iFile )%filename )      &
          &               //". Stopping!!!"
        call tem_abort()
      end if

      ! update the positions in the triangle array according to the nodes read
      me%trias(1:3,triOffset:triOffset+nTrias( iFile )-1) =                    &
        &  me%trias(1:3,triOffset:triOffset+nTrias( iFile )-1) + (nodeOffset-1)

      ! update the triangle and node offsets
      nodeOffset = nodeOffset + nPoints( iFile )
      triOffset = triOffset + nTrias( iFile )

    end do

    write(logUnit(10),*) " Done Reading. Starting to unify ..."

    ! unify the coordinates
    call tem_unify_surfaceData( me, tmp_pointCoords )

    write(logUnit(10),*) " Done Unifying. Starting to calulate the initial "//  &
      &                  "surface areas ..."

    ! allocate the surface areas array with the number of unique surface points
    allocate( me%surfArea( me%nUniquePoints_total ))

    deallocate( tmp_pointCoords )

    write(logUnit(10),*) " Done calculating the surface areas."
    write(logUnit(10),*) " Number of unique surface points: ",                 &
      &                  me%nUniquePoints_total
    write(logUnit(10),*) " Number of triangles: ",                             &
      &                  me%nTrias

    if( tmp_useInitPos )then
      allocate( me%backPointCoords( 3*me%nUniquePoints_total))
      me%backPointCoords = me%pointCoords
    end if

    ! get the min and max coordinates
    minX = minval(me%pointCoords(1:((me%nUniquePoints_total-1))*3+1:3))
    minY = minval(me%pointCoords(2:((me%nUniquePoints_total-1)*3+2):3))
    minZ = minval(me%pointCoords(3:me%nUniquePoints_total*3:3))
    maxX = maxval(me%pointCoords(1:((me%nUniquePoints_total-1)*3+1):3))
    maxY = maxval(me%pointCoords(2:((me%nUniquePoints_total-1)*3+2):3))
    maxZ = maxval(me%pointCoords(3:me%nUniquePoints_total*3:3))

    write(logUnit(3),*)'minX: ', minX
    write(logUnit(3),*)'minY: ', minY
    write(logUnit(3),*)'minZ: ', minZ
    write(logUnit(3),*)'maxX: ', maxX
    write(logUnit(3),*)'maxY: ', maxY
    write(logUnit(3),*)'maxZ: ', maxZ

  end subroutine tem_readAndUnify_surfData
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine makes the temporary of pointCoordinates unique, updates the
  !! triangle connectivity and sets the actual pointCoordinates to be the
  !! barycenters of the elements on the highest refinement level possible.
  !!
  !! @todo: IBM: Think about making the triangles unique as well!!!
  !!             Currently only points are unique!
  subroutine tem_unify_surfaceData( me, all_pointCoords )
    ! ---------------------------------------------------------------------------
    !> datatype to store the surface information
    type( tem_surfData_type ), intent(inout) :: me
    !> tmp point coordinates to be unified and stored in me
    real(kind=rk), intent(in) :: all_pointCoords(:,:)
    ! ---------------------------------------------------------------------------
    ! temporary map from the all_pointCoords to the unique ones
    integer :: map2unique(size(all_pointCoords,2))
    ! array of unique point coordinates in physical coordinates
    ! size: 3, nUniquePoints_total
    type( grw_real2darray_type ) :: pointCoords
    ! temporary dynamic array of unique treeIDs (for unifying the points)
    type( dyn_longArray_type ) :: uniqueTreeIDs
    ! growing array of integers counting the number of identical points
    type( grw_intArray_type ) :: grw_counter
    ! temporary coordinates of an element (iX,iY,iZ, level)
    integer :: tmp_coord(4)
    ! temporary treeID
    integer(kind=long_k) :: treeID
    ! position of the element in the dynamic array of unique treeIDs, in the
    ! growing array of PointCoords and in the growing array of counters
    integer :: pos
    ! was this treeID appended to the unique list
    logical :: wasAdded
    ! local tree with no treeIDs needed for unification of the points
    type( treelmesh_type ) :: loc_tree
    ! minimum x, y, z coordinate
    real(kind=rk) :: minX, minY, minZ
    ! maximum x, y, z coordinate
    real(kind=rk) :: maxX, maxY, maxZ
    ! counters
    integer :: iCoord
    integer :: iTria
    ! ---------------------------------------------------------------------------

    write(logUnit(10),*) " Unifying the coordinates ..."

    ! initialize the growing and dynamic arrays
    call init( me = pointCoords, width = 3 )

    ! get the min and max of the coordinates
    minX = minval(all_pointCoords(1,:))
    minY = minval(all_pointCoords(2,:))
    minZ = minval(all_pointCoords(3,:))

    maxX = maxval(all_pointCoords(1,:))
    maxY = maxval(all_pointCoords(2,:))
    maxZ = maxval(all_pointCoords(3,:))

    ! initialize the local tree with an origin and a bc_length
    loc_tree%global%Origin = (/minX,minY,minZ/)
    loc_tree%global%BoundingCubeLength = max( maxX - minX,                     &
      &                                       maxY - minY,                     &
      &                                       maxZ - minZ )

    ! get the treeIDs on the highest level to unify the coordinates and map
    ! the coordinates on these unique treeIDs
    do iCoord = 1, size( all_pointCoords, 2)
      ! get the coordinates
      tmp_coord = tem_CoordOfReal( loc_tree, all_pointCoords(:,iCoord ) )
      treeID = tem_IdOfCoord( tmp_coord )
      call append( me       = uniqueTreeIDs,                                   &
        &          val      = treeID,                                          &
        &          pos      = pos,                                             &
        &          wasAdded = wasAdded )
      map2unique(iCoord) = pos
      ! if the value was added the first time ...
      if( wasAdded )then
        ! ... append the coordinate
        call append( me  = pointCoords,                                        &
          &          val = all_pointCoords(:,iCoord ),                         &
          &          pos = pos )
        ! ... set the counter to 1 for this element
        call placeAt( me  = grw_counter, &
          &           val = 1,           &
          &           pos = pos          )
      else ! if the value was existing before ...
        ! ...read the coordinate and multiply it with the counter (number of
        ! times this element was added before) add the new coordinate average
        ! it by counter + 1 and write it back
        pointCoords%val( :, pos ) = ( pointCoords%val( :, pos )                &
          &     * grw_counter%val( pos ) + all_pointCoords( :, iCoord ))       &
          &     / ( grw_counter%val( pos ) + 1 )
        ! ... and increase the counter by one
        grw_counter%val( pos ) = grw_counter%val( pos ) + 1
      end if
    end do

    ! initialize the number of unique points and allocate the pointCoords
    ! array
    me%nUniquePoints_total = pointCoords%nVals
    allocate( me%pointCoords( 3*me%nUniquePoints_total ))
    do iCoord = 1, me%nUniquePoints_total
      me%PointCoords((iCoord-1)*3+1:(iCoord-1)*3+3) = pointCoords%val(:,iCoord)
    end do

    ! loop over the triangles and triangle vertices and store the right
    ! coordinate (barycenter of the parent element on the highest level)
    ! and update the position in the trias connectivity array
    do iTria = 1, me%nTrias
      do iCoord = 1, 3
        ! update the position in the triangle connectivity array
        me%trias( iCoord, iTria ) = map2unique( me%trias( iCoord, iTria ))
      end do
    end do

    call destroy( uniqueTreeIDs )
    call destroy( grw_counter )
    call destroy( pointCoords )

    write(logUnit(1),*) " nTriangles: ", me%nTrias
    write(logUnit(1),*) " nPoints: ", me%nUniquePoints_total

  end subroutine tem_unify_surfaceData
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine identifies the parent treelm elements of the surface data
  !! points.
  !!
  subroutine tem_init_surfData( me, levelDesc, globTree, iLevel )
    ! ---------------------------------------------------------------------------
    !> datatype to store the surface information
    type( tem_surfData_type ), intent(inout) :: me
    !> the level descriptor incl. ghost and halo elements as well as the
    !! communicator information on the level iLevel
    type( tem_levelDesc_type ), intent(inout) :: levelDesc
    !> global Tree information
    type( treelmesh_type ), intent(in) :: globTree
    !> the current level
    integer, intent(in) :: iLevel
    ! ---------------------------------------------------------------------------
    ! counters
    integer :: iCoord, iVal
    integer :: pos
    integer( kind=long_k ) :: tmpTreeID
    type( dyn_longArray_type ) :: tmpTreeIDs
    integer, allocatable :: posInTotal(:)
    integer, allocatable :: tmpPos(:)
    real( kind=rk ) :: huge_real
    type(tem_path_type) :: tmpPath
    type(tem_path_type) :: minHaloPath
    type(tem_path_type) :: maxHaloPath
    integer :: nFluids
    integer :: firstHalo
    integer :: lastHalo
    logical :: wasAdded
    ! ---------------------------------------------------------------------------

    huge_real = huge(huge_real)
    if( .not. allocated(me%parentIDs(iLevel)%ptrs))then
      ! allocate the array of parent positions in the local tree id list
      allocate( me%parentIDs(iLevel)%ptrs( me%nUniquePoints_total ))
    end if

    nFluids = globTree%nElems
    ! correct the halo position (in case of serial runs no halos available)
    if( globTree%global%nParts > 1 )then
      firstHalo =   levelDesc%elem%nElems( eT_fluid ) &
        &         + levelDesc%elem%nElems( eT_ghostFromCoarser )   &
        &         + levelDesc%elem%nElems( eT_ghostFromFiner ) &
        &         + 1
      lastHalo  =   levelDesc%elem%nElems( eT_fluid ) &
        &         + levelDesc%elem%nElems( eT_ghostFromCoarser )   &
        &         + levelDesc%elem%nElems( eT_ghostFromFiner ) &
        &         + levelDesc%elem%nElems( eT_halo )

      ! get the path of the min and max halo
      minHaloPath = tem_PathOf( levelDesc%total( firstHalo ))
      maxHaloPath = tem_PathOf( levelDesc%total( lastHalo  ))
    end if

    allocate( tmpPos( me%nUniquePoints_total ))

    ! for all coordinates
    do iCoord = 1, me%nUniquePoints_total
      pos = 0
!      if( any( me%pointCoords( (iCoord-1)*3+1:(iCoord-1)*3+3) <                &
!        &                                                  huge_real ))then
        ! ... get the parent treeID (first integer from real coordinates then
        !                            treeID)
        tmpTreeID = tem_IdOfCoord(                                             &
          &     tem_CoordOfReal( globTree,                                     &
          &               me%pointCoords((iCoord-1)*3+1:(iCoord-1)*3+3 ),      &
          &                              iLevel ))

        call append( me  = tmpTreeIDs,                                         &
          &          val = tmpTreeID,                                          &
          &          pos = tmpPos( iCoord ),                                   &
          &          wasAdded = wasAdded )
    end do

    allocate( posInTotal( tmpTreeIDs%nVals ))

!$omp parallel

!$omp do private( tmpPath, pos )

    do iVal=1, tmpTreeIDs%nVals
      pos = 0
      tmpPath = tem_PathOf( tmpTreeIDs%val(iVal) )
      ! compare the paths to check wether the surface point is located
      ! on this proc
      if( tem_PathComparison( tmpPath, globTree%pathList( 1 )) > -1 .and.      &
        & tem_PathComparison( tmpPath, globTree%pathList( nFluids )) < 1) then
        ! search in the fluid elements
        pos = tem_treeIDinTotal( tmpTreeIDs%val(iVal), levelDesc, eT_fluid )
      else if( tem_PathComparison( tmpPath, minHaloPath ) > -1 .and.           &
        &      tem_PathComparison( tmpPath, maxHaloPath ) <  1 )then
        ! surf point is not located on this proc but might be a halo
        pos = tem_treeIDinTotal( tmpTreeIDs%val(iVal), levelDesc, eT_halo)
      end if
      posInTotal( iVal ) = pos
    end do

!$omp end do

!      end if

!$omp do
    do iCoord = 1, me%nUniquePoints_total
      ! ... store the new coordinate
      me%parentIDs(iLevel)%ptrs(iCoord) = posInTotal( tmpPos(iCoord) )
      if( me%parentIDs(iLevel)%ptrs(iCoord) > 0 )then
        ! ... and set the corresponding property bits
        levelDesc%property( me%ParentIDs(iLevel)%ptrs( iCoord )) =             &
          &   ibset( levelDesc%property( me%ParentIDs(iLevel)%ptrs( iCoord )), &
          &          prp_hasIBM )
      end if
    end do
!$omp end do

!$omp end parallel

  end subroutine tem_init_surfData
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine updates the surface points and the parentIDs array as well
  !! as sets the correct property bits.
  !!
  subroutine tem_update_surfPos( me, levelDesc, globTree, movement, time, &
    &                            iLevel, IBMUnit, useInitPos, movPredef )
    ! ---------------------------------------------------------------------------
    !> datatype to store the surface information
    type( tem_surfData_type ), intent(inout) :: me
    !> the level descriptor incl. ghost and halo elements as well as the
    !! communicator information on the level iLevel
    type( tem_levelDesc_type ), intent(inout) :: levelDesc
    !> global Tree information
    type( treelmesh_type ), intent(inout) :: globTree
    !> spacetime function to define the motion of the surface points
    type( tem_spacetime_fun_type ) :: movement
    !> timing information
    type(tem_time_type) :: time
    !> the current level
    integer, intent(inout) :: iLevel
    !> optional output log unit other than the global logUnit
    integer, optional, intent(in) :: IBMUnit(0:tem_last_lu)
    !> shall the initial points be stored and used for updating the points
    !! later on ???
    logical, optional, intent(in) :: useInitPos
    !> logical to define wether the motion is predefined or not
    !! if not: initialize the values differently
    logical, intent(in) :: movPredef
    ! ---------------------------------------------------------------------------
    ! counters
    integer :: iPoint
    real(kind=rk) :: pos(1,3)
    real(kind=rk) :: huge_real
    integer :: IBMUnit_loc(0:tem_last_lu)
    logical :: tmp_useInitPos
    ! ---------------------------------------------------------------------------

    if( present( useInitPos ))then
      tmp_useInitPos = useInitPos
    else
      tmp_useInitPos = .false.
    end if

    if( present( IBMUnit ))then
      IBMUnit_loc = IBMUnit
    else
      IBMUnit_loc = logUnit
    end if

    huge_real = huge( huge_real )

    ! loop over the surface points and ...
    do iPoint = 1, me%nUniquePoints_total
      if( movPredef )then
        ! only update points which have fluid elements as parents
        if( me%ParentIDs(iLevel)%ptrs( iPoint ) > 0 )then
          ! ... clean the IBM property bits from the element
          levelDesc%property( me%ParentIDs(iLevel)%ptrs( iPoint)) =              &
            &     ibclr( levelDesc%property( me%ParentIDs(iLevel)%ptrs( iPoint)),&
            &            prp_hasIBM )
        end if
        ! check wether the initial point coordinates shall be used ...
        if( tmp_useInitPos )then
          ! ... yes: use array backPointCoords and ...
          ! ... store the coordinates in a temporary variable
          pos(1,1:3) = me%backPointCoords( (iPoint-1)*3+1:(iPoint-1)*3+3 )
        else
          ! ... yes: use array pointCoords and ...
          ! ... store the coordinates in a temporary variable
          pos(1,1:3) = me%pointCoords( (iPoint-1)*3+1:(iPoint-1)*3+3 )
        end if
        ! ... apply the movement and store the new positions
        pos = tem_spacetime_for( me    = movement, &
          &                      coord = pos,      &
          &                      time  = time,     &
          &                      n     = 1,        &
          &                      nComp = 3         )
        ! ... copy back the new positions
        me%pointCoords( (iPoint-1)*3+1:(iPoint-1)*3+3 ) = pos(1,1:3)
      else ! .not. movPredef -> initialize Xk with non local fluid parents with
           ! huge
        ! only update points which have fluid elements as parents
        if( me%ParentIDs(iLevel)%ptrs( iPoint ) > 0 .and.                        &
          & me%ParentIDs(iLevel)%ptrs( iPoint ) <= levelDesc%elem%nElems( eT_fluid ) )then
          ! ... clean the IBM property bits from the element
          levelDesc%property( me%ParentIDs(iLevel)%ptrs( iPoint)) =              &
            &     ibclr( levelDesc%property( me%ParentIDs(iLevel)%ptrs( iPoint)),&
            &            prp_hasIBM )
          ! check wether the initial point coordinates shall be used ...
          if( tmp_useInitPos )then
            ! ... yes: use array backPointCoords and ...
            ! ... store the coordinates in a temporary variable
            pos(1,1:3) = me%backPointCoords( (iPoint-1)*3+1:(iPoint-1)*3+3 )
          else
            ! ... yes: use array pointCoords and ...
            ! ... store the coordinates in a temporary variable
            pos(1,1:3) = me%pointCoords( (iPoint-1)*3+1:(iPoint-1)*3+3 )
          end if
          ! ... apply the movement and store the new positions
          pos = tem_spacetime_for( me    = movement, &
            &                      coord = pos,      &
            &                      time  = time,     &
            &                      n     = 1,        &
            &                      nComp = 3         )
          ! ... copy back the new positions
          me%pointCoords( (iPoint-1)*3+1:(iPoint-1)*3+3 ) = pos(1,1:3)
        else
          ! set the coordinates to be infinity such that they are not recognized
          ! by this proc any more
          me%pointCoords( (iPoint-1)*3+1:(iPoint-1)*3+3 ) = (/huge_real,         &
            &                                                 huge_real,         &
            &                                                 huge_real/)
        end if
      end if
    end do ! iPoint

    ! update the parentIDs array and set the correct property bits
    call tem_init_surfData( me        = me,                                    &
      &                     levelDesc = levelDesc,                             &
      &                     globTree  = globTree,                              &
      &                     iLevel    = iLevel )

  end subroutine tem_update_surfPos
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine calculates the surface area attached to each point
  !!
  subroutine tem_calcTriaAreas( me )
    ! ---------------------------------------------------------------------------
    !> datatype to store the surface information
    type( tem_surfData_type ), intent(inout) :: me
    ! ---------------------------------------------------------------------------
    ! counters
    integer :: iTria
    ! temporary variable to store the triangle area
    real(kind=rk) :: area
    ! temporary vectors of two sides of the triangle
    real(kind=rk) :: vec1(3), vec2(3)
    ! temporary variable storing the start and end positions in the pointCoords
    ! array
    integer :: minPos1, maxPos1, minPos2, maxPos2, minPos3, maxPos3
    ! ---------------------------------------------------------------------------

    ! reset the array of surface areas
    me%surfArea = 0._rk

    ! loop over the triangles and ...
    do iTria = 1, me%nTrias
      ! ... set the minimum and maximum positions in the pointCoords array
      minPos1 = (me%trias(1, iTria) - 1 )*3 + 1
      maxPos1 = (me%trias(1, iTria) - 1 )*3 + 3
      minPos2 = (me%trias(2, iTria) - 1 )*3 + 1
      maxPos2 = (me%trias(2, iTria) - 1 )*3 + 3
      minPos3 = (me%trias(3, iTria) - 1 )*3 + 1
      maxPos3 = (me%trias(3, iTria) - 1 )*3 + 3
!      if( me%parentIDs(iLevel)%ptrs(me%trias(1, iTria)) > 0 .and.              &
!        & me%parentIDs(iLevel)%ptrs(me%trias(2, iTria)) > 0 .and.              &
!        & me%parentIDs(iLevel)%ptrs(me%trias(3, iTria)) > 0 )then
        ! ... calculate two vectors of the triangle
        vec1 = me%pointCoords(minPos1:maxPos1)                                 &
          &  - me%pointCoords(minPos2:maxPos2)
        vec2 = me%pointCoords(minPos1:maxPos1)                                 &
          &  - me%pointCoords(minPos3:maxPos3)
        ! ... calculate the surface area of the triangle
        area = sqrt( sum( cross_product3D( vec1, vec2 )**2))
        ! ... distribute it among the attached surface points
        me%surfArea( me%trias(1, iTria)) = me%surfArea( me%trias(1, iTria))    &
          &                              + area/3._rk
        me%surfArea( me%trias(2, iTria)) = me%surfArea( me%trias(2, iTria))    &
          &                              + area/3._rk
        me%surfArea( me%trias(3, iTria)) = me%surfArea( me%trias(3, iTria))    &
          &                              + area/3._rk
!      end if
    end do

  end subroutine tem_calcTriaAreas
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine deallocates all arrays in the tem_surfaceData_type.
  !! This is used when unloading and reloading the stl surface mesh during
  !! dynamic load balancing.
  !! General information like outputprefix, timeControl and backPointCoords
  !! are NOT touched!!!
  !!
  subroutine tem_freeSurfData( me, minLevel, maxLevel )
    ! ---------------------------------------------------------------------------
    !> datatype to store the surface information
    type( tem_surfData_type ), intent(inout) :: me
    !> Level range
    integer, intent(in) :: minLevel, maxLevel
    ! ---------------------------------------------------------------------------
    integer :: iLevel
    ! ---------------------------------------------------------------------------

    me%nUniquePoints_total = 0
    me%nTrias = 0

    if( allocated( me%pointCoords ))     deallocate( me%pointCoords )
    if( allocated( me%surfArea ))        deallocate( me%surfArea )
    if( allocated( me%parentIDs ))then
      do iLevel = minLevel, maxLevel
        if( allocated( me%parentIDs(iLevel)%ptrs ))then
          deallocate( me%parentIDs(iLevel)%ptrs )
        end if
      end do
    end if
    if( allocated( me%trias ))           deallocate( me%trias )
    if( allocated( me%stlHead ))         deallocate( me%stlHead )

  end subroutine tem_freeSurfData
! ****************************************************************************** !


end module tem_surfaceData_module
! ****************************************************************************** !
