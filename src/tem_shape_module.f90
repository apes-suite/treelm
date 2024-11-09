! Copyright (c) 2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2014, 2016, 2019, 2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2013, 2017, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012, 2014, 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013 Nikhil Anand <nikhil.anand@uni-siegen.de>
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
! ****************************************************************************** !
!> Geometrical shape definitions for pre-defined bodies
!!
!! This module provides data types and definitions for pre-defined geometrical
!! objects such as points, lines, planes, and cubes.
!!
!! Shapes
!! ---
!!
!! - For points, lines, planes and cubes see [[tem_canonicalnd_module]]
!! - A global shape to select all elements
!!
module tem_shape_module

  use mpi
  ! include treelm modules
  use env_module,             only: rk, long_k, labelLen, pathLen, &
    &                               globalMaxLevels
  use tem_aux_module,         only: tem_abort
  use tem_bc_prop_module,     only: tem_bc_prop_type
  use treelmesh_module,       only: treelmesh_type
  use tem_topology_module,    only: tem_levelOf, tem_firstIdAtLevel, &
    &                               tem_lastIdAtLevel
  use tem_pointData_module,   only: tem_grwPoints_type, init, append
  use tem_geometry_module,    only: tem_CoordOfReal, tem_PosOfId, &
    &                               tem_BaryOfId, tem_ElemSize
  use tem_dyn_array_module,   only: dyn_intArray_type, append
  use tem_property_module,    only: prp_hasBnd, prp_fluidify, prp_hasQVal
  use tem_param_module,       only: qQQQ, qDir, qOffset
  use tem_stencil_module,     only: tem_stencilHeader_type
  use tem_canonicalND_module, only: tem_canonicalND_type, tem_canonicalND_out, &
    &                               tem_load_canonicalND,                      &
    &                               tem_cano_initSubTree,                      &
    &                               tem_cano_checkNeigh,                       &
    &                               tem_cano_storePntsInSubTree
  use tem_triangle_module,    only: tem_triangle_type, tem_load_triangle,      &
    &                               tem_triangle_out, tem_triangleCubeOverlap
  use tem_stl_module,         only: tem_load_stl, tem_stlData_type, &
    &                               tem_stlCubeOverlap, tem_stlHead_out
  use tem_sphere_module,      only: tem_sphere_type, tem_load_sphere,          &
    &                               tem_sphere_out, tem_sphereCubeOverlap
  use tem_cylinder_module,    only: tem_cylinder_type, tem_load_cylinder,      &
    &                               tem_cylinder_out, tem_cylinderCubeOverlap
  use tem_ellipsoid_module,   only: tem_ellipsoid_type, tem_load_ellipsoid,    &
    &                               tem_ellipsoid_out, tem_ellipsoidCubeOverlap
  use tem_logging_module,     only: logUnit, tem_log, tem_toStr
  use tem_cube_module,        only: tem_cube_type, tem_convertTreeIDtoCube
  use tem_transformation_module, only: tem_transformation_type, &
    &                                  tem_load_transformation
  use tem_debug_module,       only: dbgUnit

  ! include aotus modules
  use aot_table_module, only: aot_table_open, aot_table_close,                 &
    &                         aot_table_length, aot_get_val
  use aotus_module,     only: flu_State, aoterr_NonExistent
  use aot_out_module,   only: aot_out_type, aot_out_val,                       &
    &                         aot_out_open_table, aot_out_close_table

  implicit none

  private

  public :: tem_shape_type
  public :: tem_load_shape
  public :: tem_global_shape
  public :: tem_geometrical_shape
  public :: tem_property_shape
  public :: tem_local_shape
  public :: tem_boundary_shape
  public :: tem_level_shape
  public :: tem_shape_out
  public :: tem_shape2subTree

  !> Parameters for different tracking shapes
  !> Global mesh
  integer, parameter :: tem_global_shape = 0
  !> treelm geometrical object like canoND, sphere, cylinder, eclipse, triangle
  !! STL (STLs are converted further into triangles)
  integer, parameter :: tem_geometrical_shape = 1
  !> elements that has a certain property
  integer, parameter :: tem_property_shape = 2
  !> elements on the local partition
  integer, parameter :: tem_local_shape = 3
  !> elements of one or more boundaries
  integer, parameter :: tem_boundary_shape = 4
  !> elements of certain levels
  integer, parameter :: tem_level_shape = 5

  !> Complete shape definitions
  type tem_shape_type

    !> a kind of the shape defined.
    character(len=labelLen) :: kind

    !> a identification for the shape
    integer :: shapeID = 0

    !> canonical definition
    type(tem_canonicalND_type), allocatable :: canoND(:)

    !> triangle definition
    type(tem_triangle_type), allocatable :: triangle(:)

    !> STL definition
    type(tem_stlData_type) :: stl_data

    !> spheres definition
    type(tem_sphere_type), allocatable :: sphere(:)

    !> ellipsoid definition
    type(tem_ellipsoid_type), allocatable :: ellipsoid(:)

    !> cylinder definition
    type(tem_cylinder_type), allocatable :: cylinder(:)

    !> property bits
    integer(kind=long_k) :: propBits = 0_long_k

    !> boundary labels, used to identify elements belong to these boundaries
    !! It is allocated and set in routine: tem_shape_load_bcLabels
    character(len=labelLen), allocatable :: bcLabels(:)

    !> boundary elements below this threshold are omitted
    real(kind=rk) :: cutOffQVal

    !> level range for level shape type
    integer :: minLevel = 1
    integer :: maxLevel = globalMaxLevels

    !> If true then subTree is created for inverted shape i.e nonintersected
    logical :: inverted = .false.
  end type tem_shape_type

  interface tem_shape_out
    module procedure tem_shape_out_scal
    module procedure tem_shape_out_vec
  end interface tem_shape_out

  interface tem_load_shape
    module procedure tem_load_shapes
    module procedure tem_load_shape_single
  end interface tem_load_shape

  contains

! ****************************************************************************** !
  !> Read in an arbitrary shapes from a lua file defined as multiple tables
  !!
  !! read a shape like for example inside a tracking table
  !!```lua
  !! tracking = {
  !!             {variable = { 'velocity' },
  !!              shape = { kind = 'canoND',
  !!                        object = { origin = { 1.0, 1.0, 1.0 },
  !!                                   vec = { 2.0, 2.0, 2.0 },
  !!                                   segments = { 10, 20, 30 } }
  !!                      }}
  !!            } -- tracking table
  !!```
  !! elements that has a certain property can also be tracked.
  !! This feature enables us to track boundary elements.
  !!```lua
  !! tracking = {
  !!              { variable = { 'velocity' },
  !!                shape = { kind = 'property',
  !!                          property = {'boundary'} },
  !!                      }
  !!              }
  !!            } -- tracking table
  !!```
  !!
  subroutine tem_load_shapes( me, conf, parent, key, iError, reqSegments )
    !---------------------------------------------------------------------------
    !> array of shape type defined in a lua file
    type(tem_shape_type), allocatable, intent(out) :: me(:)
    !> lua config file to load shape from
    type(flu_state) :: conf
    !> optional parent handle
    integer, optional, intent(in) :: parent
    !> optional key to load from
    character(len=*), optional, intent(in) :: key
    !> error flag
    integer, intent(out), optional :: iError
    !> Is true if use_get_point is true in output table
    logical, optional, intent(in) :: reqSegments
    !---------------------------------------------------------------------------
    character(len=32) :: localKey
    integer :: nShapes  ! number of shape table entries
    integer :: iShape, shape_table, sub_table
    !---------------------------------------------------------------------------
    if( present( key )) then
      localKey = key
    else
      localKey = 'shape'
    endif

    ! open the table
    ! shape = {}
    call aot_table_open( L       = conf,                                       &
      &                  thandle = shape_table,                                &
      &                  parent  = parent,                                     &
      &                  key     = trim( localKey ))

    call aot_table_open( L       = conf,                                       &
      &                  parent  = shape_table,                                &
      &                  thandle = sub_table,                                  &
      &                  pos     = 1 )

    ! no shape table is defined return 0-sized array
    if ( shape_table == 0 ) then
      write(logUnit(2),*) ' Shape table is not defined'
      write(logunit(2),*) ' ... using global mesh'
      allocate(me(1))
      me(1)%kind = 'all'
      me(1)%shapeID = tem_global_shape
      if(present(iError)) iError = ibset(0, aoterr_NonExistent)
    else if ( sub_table == 0 ) then ! shape is a single table
      ! load table from parent shape_table
      call aot_table_close(L=conf, thandle=sub_table)
      allocate( me( 1 ))
      call tem_load_shape_single( me          = me(1),       &
        &                         conf        = conf,        &
        &                         sub_table   = shape_table, &
        &                         iError      = iError,      &
        &                         reqSegments = reqSegments  )
    else ! multiple table
      call aot_table_close( L=conf, thandle=sub_table )
      ! open the first entry in the shape table
      ! shape = {{this entry}, {second entry}}
      ! get the number of entries in  the shape table
      nShapes = aot_table_length( L=conf, thandle=shape_table )
      write(logUnit(1),*) 'Number of shapes defined: ', nShapes
      allocate( me( nShapes ))
      do iShape = 1, nShapes
        call aot_table_open( L       = conf,                                   &
          &                  thandle = sub_table,                              &
          &                  parent  = shape_table,                            &
          &                  pos     = iShape )
        call tem_load_shape_single( me          = me(iShape), &
          &                         conf        = conf,       &
          &                         sub_table   = sub_table,  &
          &                         iError      = iError,     &
          &                         reqSegments = reqSegments )
        call aot_table_close( L = conf, thandle = sub_table )
      end do
    end if

    call aot_table_close( L=conf, thandle=shape_table )

  end subroutine tem_load_shapes
! ****************************************************************************** !


! ****************************************************************************** !
  !> Read in an arbitrary shape from a lua file defined in a single table
  !!
  !! read a shape like for example inside a tracking table
  !!```lua
  !! tracking = { variable = { 'velocity' },
  !!              shape = { kind = 'canoND',
  !!                        object = { origin = { 1.0, 1.0, 1.0 },
  !!                                   vec = { 2.0, 2.0, 2.0 },
  !!                                   segments = { 10, 20, 30 } }
  !!                      }
  !!            } -- tracking table
  !!```
  !!
  subroutine tem_load_shape_single( me, conf, key, parent, sub_table, iError, &
    &                               reqSegments                               )
    !---------------------------------------------------------------------------
    !> shape type defined in a lua file
    type(tem_shape_type), intent(out) :: me
    !> lua state
    type(flu_state) :: conf
    !> optional key to load from
    character(len=*), optional, intent(in) :: key
    !> optional parent handle
    integer, optional, intent(in) :: parent
    !> shape table handle
    integer, optional, intent(in) :: sub_table
    !> error flag
    integer, optional, intent(out) :: iError
    !> Is true if use_get_point is true in output table
    logical, optional, intent(in) :: reqSegments
    !---------------------------------------------------------------------------
    character(len=PathLen) :: buffer
    integer :: iError_loc
    integer :: shape_table
    character(len=32) :: localKey
    type(tem_transformation_type) :: transform
    !---------------------------------------------------------------------------

    if( present( iError )) iError = 0
    if( present( key )) then
      localKey = key
    else
      localKey = 'shape'
    endif

    ! open the table
    ! shape = {}
    if( present( sub_table )) then
      shape_table = sub_table
    else
      call aot_table_open( L       = conf,          &
        &                  thandle = shape_table,   &
        &                  parent  = parent,        &
        &                  key     = trim(localKey) )
    endif

    if (shape_table > 0 ) then

      ! load transformation table if defined and pass the transformation
      ! to geom table to apply the transformation to geometry object
      call tem_load_transformation(transform = transform, conf = conf, &
        &                          thandle = shape_table)

      ! Get the entry 'kind'
      ! shape = {{kind='test'}}
      call aot_get_val( L       = conf,        &
        &               thandle = shape_table, &
        &               val     = buffer,      &
        &               ErrCode = iError_loc,  &
        &               key     = 'kind',      &
        &               default = 'none'       )

      me%kind = trim(buffer)
      write(logUnit(1),*)' Loading shape kind: '//trim(me%kind)
      ! Read in the type and configuration of tracking shape
      ! (point, line, plane, ...)
      select case( trim(buffer) )
      case('global', 'all', 'globalmesh')
        me%shapeID = tem_global_shape

      case('canonicalND', 'canoND')
        me%shapeID = tem_geometrical_shape
        me%kind = 'canoND'
        call tem_load_canonicalND( me          = me%canoND,   &
          &                        conf        = conf,        &
          &                        thandle     = shape_table, &
          &                        transform   = transform,   &
          &                        reqSegments = reqSegments  )

      case('triangle')
        me%shapeID = tem_geometrical_shape
        call tem_load_triangle( me        = me%triangle, &
          &                     conf      = conf,        &
          &                     thandle   = shape_table, &
          &                     transform = transform    )

      case('stl')
        me%shapeID = tem_geometrical_shape
        call tem_load_stl( stl_data  = me%stl_data, &
          &                conf      = conf,        &
          &                thandle   = shape_table, &
          &                transform = transform    )

      case('sphere')
        me%shapeID = tem_geometrical_shape
        call tem_load_sphere( me        = me%sphere,   &
          &                   conf      = conf,        &
          &                   thandle   = shape_table, &
          &                   transform = transform    )

      case('ellipsoid')
        me%shapeID = tem_geometrical_shape
        call tem_load_ellipsoid( me        = me%ellipsoid, &
          &                      conf      = conf,         &
          &                      thandle   = shape_table,  &
          &                      transform = transform     )

      case('cylinder')
        me%shapeID = tem_geometrical_shape
        call tem_load_cylinder( me        = me%cylinder, &
          &                     conf      = conf,        &
          &                     thandle   = shape_table, &
          &                     transform = transform    )

      case('property')
        me%shapeID = tem_property_shape
        call tem_shape_load_propLabel( me%propBits, conf, shape_table )

      case( 'level' )
        me%shapeID = tem_level_shape
        call tem_shape_load_level( me%minLevel, me%maxLevel, conf, shape_table )

      case('boundary')
        me%shapeID = tem_boundary_shape
        call tem_shape_load_bcLabels( me%bcLabels, conf, shape_table )
        call aot_get_val( L       = conf,            &
          &               thandle = shape_table,     &
          &               val     = me%cutOffQVal,   &
          &               ErrCode = iError_loc,      &
          &               key     = 'cutoff_qvalue', &
          &               default = 1.0_rk           )
        write(logUnit(1),*)' Cutoff qValue: ', me%cutOffQVal


      ! If find unknown shape, code will stop!
      case default
        write(logUnit(1),*) '  Unknown shape specified: '//trim(buffer)
        write(logUnit(1),*) '  Choose: '
        write(logUnit(1),*) '    shape = { kind = "canoND", '
        write(logUnit(1),*) '              object = {origin = {-0.5,0.,0.}, '
        write(logUnit(1),*) '              vec = {1.0, 0.0, 0.0},'
        write(logUnit(1),*) '              segments = 1000 } } or'
        write(logUnit(1),*) '    shape = { kind = "all" } or'
        write(logUnit(1),*) '    shape = { kind = "property",'
        write(logUnit(1),*) '              property = {"boundary"} }'
        write(logUnit(1),*) '    shape = { kind = "boundary",'
        write(logUnit(1),*) '              boundary = {"bc1", "bc2"} }'
        write(logUnit(1),*) '  Stopping ...'
        call tem_abort()
      end select

      call aot_get_val( L       = conf,        &
        &               thandle = shape_table, &
        &               val     = me%inverted, &
        &               ErrCode = iError_loc,  &
        &               key     = 'inverted',  &
        &               default = .false.      )

      write(logUnit(1),*)' Inverted shape: ',me%inverted

      if( .not. present( parent ))                                             &
        &       call aot_table_close( L = conf, thandle = shape_table )

    else
      write(logUnit(2),*) 'Shape table is not defined'
      write(logunit(2),*) '... using global mesh'
      me%kind = 'all'
      me%shapeID = tem_global_shape
      if( present( iError )) iError = ibset( iError, aoterr_NonExistent )
    endif

  end subroutine tem_load_shape_single
! ****************************************************************************** !


  ! ************************************************************************** !
  !> This routine creates subTree from geometry intersection
  subroutine tem_shape_subTreeFromGeomInters( me, inTree, countElems, &
    &                                         countPoints, grwPnts,   &
    &                                         storePnts, map2global   )
    ! --------------------------------------------------------------------------
    !> shape objects on which to work
    type(tem_shape_type ),intent(in) :: me
    !> Global tree
    type(treelmesh_type), intent(in) :: inTree
    !> How many elements there will be for each level in the track
    integer, intent( inout ) :: countElems( globalMaxLevels )
    !> How many points there will be
    integer, intent( inout ) :: countPoints
    !> growing array for the map2global
    type(dyn_intArray_type), intent(inout) :: map2global
    !> growing array to store tracking points
    type(tem_grwPoints_type), intent(inout) :: grwPnts
    !> to Store points in grwPnts
    logical, intent(in) :: storePnts
    ! --------------------------------------------------------------------------
    real(kind=rk) :: tStart, tEnd
    integer :: iElem, dPos, tLevel, iObj
    logical :: wasAdded, intersects, addToSubTree
    logical :: foundAny, foundGlobally
    type(tem_cube_type) :: cube
    integer(kind=long_k) :: treeID
    integer :: iError
    ! --------------------------------------------------------------------------
    write(logUnit(5),"(A)") '   Extracting subTree from geometrical shape' &
      &                     // ' - tree intersection'
    tStart = mpi_wtime()

    select case (trim(me%kind))
    case ('canoND')
      ! treat canoND seperately to efficiently handle cano kind = points
      call tem_cano_initSubTree( me            = me%canoND(:), &
        &                        inTree        = inTree,       &
        &                        countElems    = countElems,   &
        &                        map2global    = map2global,   &
        &                        shapeInverted = me%inverted   )

      foundany = (sum(countElems) > 0)
      call MPI_Allreduce( foundAny, foundGlobally, 1,                      &
        &                 MPI_LOGICAL, MPI_LOR, inTree%global%comm, iError )
      if (storePnts .and. (.not. foundGlobally)) then
        write(logUnit(5),"(A)") '  WARNING: did not find any element for canonical shape'
        write(logUnit(5),"(A)") '           looking at potential neighboring elements'
        ! No elements were found for the canoND objects, check whether
        ! any neighboring element is close, to extrapolate nearby points.
        call tem_cano_checkNeigh( me            = me%canoND(:), &
          &                       inTree        = inTree,       &
          &                       countElems    = countElems,   &
          &                       map2global    = map2global    )
        write(logUnit(5),*) '  --> after considering neighbors found ', sum(countElems)
        write(logUnit(5),*) '      elements on this process (', inTree%global%myPart, ')'
      end if

    case ('triangle', 'stl', 'sphere', 'ellipsoid', 'cylinder')
      do iElem = 1, inTree%nElems
        treeID = inTree%treeID(iElem)
        call tem_convertTreeIDtoCube(cube, inTree, treeID)
        intersects = .false.
        ! Treat each object of same geometry kind
        select case (trim(me%kind))
        case ('triangle')
          do iObj = 1, size(me%triangle)
            intersects = intersects                                 &
              & .or. tem_triangleCubeOverlap(me%triangle(iObj), cube)
          end do
        case ('stl')
          intersects = tem_stlCubeOverlap(me%stl_data, cube)
        case ('sphere')
          do iObj = 1, size(me%sphere)
            intersects = intersects                             &
              & .or. tem_sphereCubeOverlap(me%sphere(iObj), cube)
          end do
        case ('ellipsoid')
          do iObj = 1, size(me%ellipsoid)
            intersects = intersects                                   &
              & .or. tem_ellipsoidCubeOverlap(me%ellipsoid(iObj), cube)
          end do
        case ('cylinder')
          do iObj = 1, size(me%cylinder)
            intersects = intersects                                 &
              & .or. tem_cylinderCubeOverlap(me%cylinder(iObj), cube)
          end do
        end select

        addToSubTree = .false.
        if (.not. me%inverted .and. intersects) then
          ! Shape intersects with current element and not inverted
          addToSubTree = .true.
        else if (me%inverted .and. .not. intersects) then
          ! shape not intersected and is inverted shape so add this to subTree
          addToSubTree = .true.
        end if

        if (addToSubTree) then
          ! append iElem in inTree to the map (note that already existing
          ! ones are omitted)
          call append( me       = map2global, &
            &          pos      = dpos,       &
            &          val      = iElem ,     &
            &          wasAdded = wasAdded    )

          ! Count up if it was added
          if( wasAdded ) then
            tLevel   = tem_levelOf( treeID )
            countElems( tLevel ) = countElems( tLevel ) + 1
          end if ! wasAdded
        end if !intersects
      end do !iElem
    case default
      call tem_abort('In tem_shape_subTreeFromGeomInters: Unknown shape kind')
    end select

    countPoints = 0
    call init(me=grwPnts)
    if ( storePnts .and. (map2global%nVals > 0) ) then
      select case (trim(me%kind))
      case ('canoND')
        call tem_cano_storePntsInSubTree( me          = me%canoND(:), &
          &                               inTree      = inTree,       &
          &                               countPoints = countPoints,  &
          &                               grwPnts     = grwPnts,      &
          &                               map2global  = map2global    )
      case default
        call tem_abort('Use get_points supported only for canoND')
      end select
    end if
    tEnd = mpi_wtime()
    write(logunit(4),"(A,E12.6)") 'Done. This process cost: ', tEnd-tStart

  end subroutine tem_shape_subTreeFromGeomInters
  ! ************************************************************************** !

! ****************************************************************************** !
  !> Write a array of shapes to lua file
  !!
  subroutine tem_shape_out_vec( me, conf )
    !---------------------------------------------------------------------------
    !> shape types to write out
    type( tem_shape_type ), intent(in) :: me(:)
    !> Aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    !---------------------------------------------------------------------------
    integer :: iVal
    !---------------------------------------------------------------------------

    ! create a table with name shape if not exist
    if( conf%level .eq. 0 ) then
      call aot_out_open_table( put_conf = conf, tname = 'shape' )
    else
      call aot_out_open_table( put_conf = conf )
    endif

    do iVal = 1, size(me)
      call tem_shape_out_scal( me(iVal), conf )
    end do

    call aot_out_close_table( put_conf = conf )

  end subroutine tem_shape_out_vec
! ****************************************************************************** !


! ****************************************************************************** !
  !> Write a shape to lua file
  !!
  subroutine tem_shape_out_scal( me, conf )
    !---------------------------------------------------------------------------
    !> shape types to write out
    type( tem_shape_type ), intent(in) :: me
    !> Aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    !---------------------------------------------------------------------------

    ! create a table with name shape if not exist
    if( conf%level .eq. 0 ) then
      call aot_out_open_table( put_conf = conf, tname = 'shape' )
    else
      call aot_out_open_table( put_conf = conf )
    end if

    call aot_out_val( put_conf = conf, vname = 'kind', &
      &               val = trim(me%kind) )
    ! choose output shape kind
    select case( trim(me%kind) )
    case('canoND')
      call tem_canonicalND_out(me%canoND, conf)
    case('triangle')
      call tem_triangle_out(me%triangle, conf)
    case('stl')
      call tem_stlHead_out(me%stl_data%head, conf)
    case('sphere')
      call tem_sphere_out(me%sphere, conf)
    case('ellipsoid')
      call tem_ellipsoid_out(me%ellipsoid, conf)
    case('cylinder')
      call tem_cylinder_out(me%cylinder, conf)
    case('property')
      call tem_shape_propLabel_out(me%propBits, conf)
    case('level')
      call tem_shape_level_out(me%minLevel, me%maxLevel, conf)
    case('boundary')
      call tem_shape_bcLabel_out(me%bcLabels, conf)
    end select

    call aot_out_close_table( put_conf = conf )

  end subroutine tem_shape_out_scal
! ****************************************************************************** !


! ****************************************************************************** !
  !> Loading property labels from the config file,
  !! set the property bits accordingly
  !!
  subroutine tem_shape_load_propLabel( propBits, conf, thandle )
    !---------------------------------------------------------------------------
    !> propBits
    integer( kind=long_k ) :: propBits
    !> lua config handle
    type(flu_state) :: conf
    !> table handle from which to read
    integer, intent(in) :: thandle
    !---------------------------------------------------------------------------
    ! lua handles
    integer :: propLabel_handle, nPropLabels
    character(len=labelLen) :: labelBuff
    integer :: iLabel, iErr
    !---------------------------------------------------------------------------

    call aot_table_open( L       = conf,                                       &
      &                  parent  = thandle,                                    &
      &                  thandle = propLabel_handle,                           &
      &                  key     = 'property' )
    ! get the number of property labels
    nPropLabels = aot_table_length( L = conf, thandle = propLabel_handle )

    do iLabel = 1, nPropLabels

      ! Now read in property labels
      call aot_get_val( L       = conf,                                        &
        &               thandle = propLabel_handle,                            &
        &               val     = labelBuff,                                   &
        &               ErrCode = iErr,                                        &
        &               pos     = iLabel )
      write(logUnit(1),*) '  name of property label: '//labelBuff
      ! set propBits according to property labels
      if( iErr == 0 ) then
        select case( trim(labelBuff) )
        case( 'boundary' )
          propBits = ibset( propBits, prp_hasBnd )
        case( 'solidified' )
          propBits = ibset( propBits, prp_fluidify )
        case default
          write(logUnit(1),*) '  Unknown property label in a shape specified! '&
            &             //trim(labelBuff)
          write(logUnit(1),*) '  Ignored it. '
        endselect
      endif
    enddo

    call aot_table_close( L = conf, thandle = propLabel_handle )

  end subroutine tem_shape_load_propLabel
! ****************************************************************************** !


! ****************************************************************************** !
  !> Loading bc labels from the config file,
  !! save those labels for further use.
  !!
  subroutine tem_shape_load_bcLabels( bcLabels, conf, thandle )
    !---------------------------------------------------------------------------
    !> bc labels
    character(len=labelLen), allocatable :: bcLabels(:)
    !> lua config handle
    type(flu_state) :: conf
    !> table handle from which to read
    integer, intent(in) :: thandle
    !---------------------------------------------------------------------------
    ! lua handles
    integer :: bcLabel_handle, nBCLabels
    character(len=labelLen) :: labelBuff
    integer :: iLabel, iErr
    !---------------------------------------------------------------------------

    call aot_table_open( L       = conf,            &
      &                  parent  = thandle,         &
      &                  thandle = bcLabel_handle,  &
      &                  key     = 'boundary' )

    ! get the number of property labels
    nBCLabels = aot_table_length( L = conf, thandle = bcLabel_handle )

    if ( nBCLabels == 0 ) then
      write(logUnit(3),*) ' Boundary label table is empty!'
      allocate( bcLabels(0) )
    else
      allocate( bcLabels( nBCLabels ) )

      do iLabel = 1, nBCLabels

        ! Now read in BC labels
        call aot_get_val( L       = conf,           &
          &               thandle = bcLabel_handle, &
          &               val     = labelBuff,      &
          &               ErrCode = iErr,           &
          &               pos     = iLabel          )


        ! set propBits according to property labels
        if( iErr == 0 ) then
          write(logUnit(3),*) '  name of boundary label: '//trim(labelBuff)
          bcLabels( iLabel ) = trim(labelBuff)
        else
          write(logUnit(3),*) '  failed to read boundary label!'
          bcLabels( iLabel ) = ''
        end if

      end do

    end if ! nBCLabels == 0

    call aot_table_close( L = conf, thandle = bcLabel_handle )

  end subroutine tem_shape_load_bcLabels
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine tem_shape_load_level( minLevel, maxLevel, conf, thandle )
    !---------------------------------------------------------------------------
    !> level range
    integer :: minLevel, maxLevel
    !> lua config handle
    type(flu_state) :: conf
    !> table handle from which to read
    integer, intent(in) :: thandle
    !---------------------------------------------------------------------------
    ! lua handles
    integer :: this_handle, iErr, nVals
    !---------------------------------------------------------------------------

    call aot_table_open( L       = conf,        &
      &                  parent  = thandle,     &
      &                  thandle = this_handle, &
      &                  key     = 'level'     )

    ! get the number of property labels
    nVals = aot_table_length( L = conf, thandle = this_handle )

    if ( nVals == 0 ) then
      minLevel = 1
      maxLevel = globalMaxLevels
    else if( nVals == 1 ) then
      call aot_get_val( L       = conf,       &
        &               thandle = this_handle,&
        &               val     = minLevel,   &
        &               ErrCode = iErr,       &
        &               pos     = 1 )
      maxLevel = globalMaxLevels
    else
      ! Now read in property labels
      call aot_get_val( L       = conf,       &
        &               thandle = this_handle,&
        &               val     = minLevel,   &
        &               ErrCode = iErr,       &
        &               pos     = 1 )

      call aot_get_val( L       = conf,       &
        &               thandle = this_handle,&
        &               val     = maxLevel ,  &
        &               ErrCode = iErr,       &
        &               pos     = 2 )
    end if

    call aot_table_close( L = conf, thandle = this_handle )
    call tem_log(5, ' Loaded shape level range: '//trim(tem_toStr(minLevel))&
      &             //' to '//trim(tem_toStr(maxLevel)))

  end subroutine tem_shape_load_level
! ****************************************************************************** !

  ! ************************************************************************** !
  !> Write out a shape level in lua format
  !!
  subroutine tem_shape_level_out( minLevel, maxLevel, conf )
    ! --------------------------------------------------------------------------
    !> Minlevel and maxlevel
    integer, intent(in) :: minLevel, maxLevel
    !> Aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    ! --------------------------------------------------------------------------

    call aot_out_open_table( put_conf = conf, tname = 'level' )

    call aot_out_val( put_conf = conf, val = minlevel )
    call aot_out_val( put_conf = conf, val = maxlevel )

    call aot_out_close_table( put_conf = conf )

  end subroutine tem_shape_level_out
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Write out a shape property label in lua format
  !!
  subroutine tem_shape_propLabel_out( propBits, conf )
    ! --------------------------------------------------------------------------
    !> property bits
    integer(kind=long_k), intent(in) :: propBits
    !> Aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    ! --------------------------------------------------------------------------

    call aot_out_open_table( put_conf = conf, tname = 'property' )

    if (btest(propBits, prp_hasBnd)) &
      call aot_out_val( put_conf = conf, val = 'boundary' )

    if (btest(propBits, prp_fluidify)) &
      call aot_out_val( put_conf = conf, val = 'solidified' )

    call aot_out_close_table( put_conf = conf )

  end subroutine tem_shape_propLabel_out
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Write out a shape boundary label in lua format
  !!
  subroutine tem_shape_bcLabel_out( bcLabels, conf )
    ! --------------------------------------------------------------------------
    !> Boundary labels
    character(len=labelLen), intent(in) :: bcLabels(:)
    !> Aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    ! --------------------------------------------------------------------------
    integer :: iLabel
    ! --------------------------------------------------------------------------

    call aot_out_open_table( put_conf = conf, tname = 'boundary' )

    do iLabel = 1, size(bcLabels)
      call aot_out_val( put_conf = conf, val = bcLabels(iLabel) )
    end do

    call aot_out_close_table( put_conf = conf )

  end subroutine tem_shape_bcLabel_out
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine identifies elements that belong to certain bounaries.
  !! Labels of required boundaries are given by bcLabels.
  !! bc_prop contains boudnary_ID of all local elements.
  !! Firstly, bcLabels are converted into bcIDs.
  !! Then all elements in bc_prop are looped over to check if it matches one of
  !! the required bcID. If match, its position is save in map2global.
  !! Number of elements found on each level is saved in countElems.
  subroutine tem_shape_findElemByBCLabels( bcLabels, cutOffQVal, bc_prop,     &
    &                                      foundAny, map2global, inTree,      &
    &                                      countPoints, grwPnts, storePnts,   &
    &                                      bcIDs, stencil )
    ! --------------------------------------------------------------------------
    !> bcLabels
    character(len=labelLen), intent(in) :: bcLabels(:)
    !> Cut off qValue
    real(kind=rk), intent(in) :: cutOffQVal
    !> bc property
    type(tem_bc_prop_type), intent(in) :: bc_prop
    !> if any element be identified
    logical,                 intent(out) :: foundAny
    !> dynamic array. Elements positions in bc_prop%property%elemID
    type(dyn_intArray_type), intent(inout) :: map2global
    !> Global tree
    type(treelmesh_type), intent(in) :: inTree
    !> How many points there will be
    integer, intent( inout ) :: countPoints
    !> growing array to store tracking points
    type(tem_grwPoints_type), intent(inout) :: grwPnts
    !> to Store points in grwPnts
    logical, intent(in) :: storePnts
    !> id of boundary condition to be tracked
    integer, allocatable, intent(out) :: bcIDs(:)
    !> stencil required to get useful links
    type( tem_stencilHeader_type ), optional, intent(in) :: stencil
    ! --------------------------------------------------------------------------
    integer :: dPos, iElem, iBC, nBCs, iBCtype, posInTree, QQN
    logical :: wasAdded, found
    integer, allocatable :: map(:)
    real(kind=rk), allocatable :: qValWeight(:), weightedQVal(:)
    real(kind=rk) :: dx, bary(3), point(3), minQVal
    integer :: iDir, length, iElem_qVal, minQValiDir
    ! --------------------------------------------------------------------------

    if( present( stencil ) ) then
      allocate( map( size(stencil%map) ) )
      map = stencil%map
      qqn = stencil%QQN
    else
      allocate( map(qQQQ) )
      map = qDir
      qqn = qQQQ
    end if

    ! calculate weights on qValue for each direction. 
    ! It depends on the link length.
    allocate(weightedQVal(qqn))
    allocate(qValWeight(qqn))
    do iDir = 1, qqn
      length   = qOffset( map(iDir), 1 )**2 &
        &      + qOffset( map(iDir), 2 )**2 &
        &      + qOffset( map(iDir), 3 )**2

      qValWeight(iDir) = sqrt(real(length, kind=rk))
    end do

    foundAny = .false.
    ! convert bcLabels to bcIDs
    nBCs = size(bcLabels)
    allocate( bcIDs( nBCs ) )

    ! loop over all required bcLabels
    do iBC = 1, nBCs

      found = .false.

      ! loop over all BCs in the mesh
      do iBCtype = 1, bc_prop%nBCtypes
        if ( trim(bcLabels(iBC)) == bc_prop%bc_label( iBCtype ) ) then
          bcIDs( iBC ) = iBCtype
          found = .true.
          exit
        end if
      end do

      if ( .not. found ) then
        write(logUnit(1),*) 'Required BC label: '//trim(bcLabels(iBC))
        write(logUnit(1),*) 'can not be found in given mesh!'
        stop
      end if
    end do ! iBC

    iElem_qVal = 0
    ! Loop over all element with boundary property
    do iElem = 1, bc_prop%property%nElems
      posInTree = bc_prop%property%elemID(iElem)
      ! Count the elem with qVal property to access the qVal from
      ! bc_prop%qVal
      if ( btest( inTree%elemPropertyBits( posInTree ), prp_hasQVal) ) then
        iElem_qVal = iElem_qVal + 1
      end if

      do iBC = 1, nBCs
        if ( any(int(bc_prop%boundary_ID( map(:qqn), iElem )) == bcIDs(iBC) ) ) then

          ! If element has qValues then consider element which satisfy
          ! cutOffQVal
          if (btest(inTree%ElemPropertyBits(posInTree), prp_hasQVal)) then
            do iDir = 1, qqn
              weightedQVal(iDir) = qValWeight(iDir)                  &
                &                * bc_prop%qVal(map(iDir), iElem_qVal)
            end do
            minQValiDir = minloc(weightedQVal, DIM=1, MASK=(weightedQVal>0))
            minQVal = bc_prop%qVal(map(minQValiDir), iElem_qVal)
          else
            minQVal = cutOffQVal
          end if

          if (minQVal <= cutOffQVal) then

            ! Append to treeID list (note that already existing ones are
            ! omitted)
            call append( me       = map2global,   &
              &          pos      = dPos,         &
              &          val      = posInTree,    &
              &          wasAdded = wasAdded )

            ! Count up if it was added
            if( wasAdded ) then
              ! tLevel = tem_levelOf( inTree%treeID(posInTree) )
              ! countElems( tLevel ) = countElems( tLevel ) + 1
              foundAny = .true.

              if (storePnts) then
                ! Create growing array of points for boundary elements with 
                ! qValues.
                ! Point is created along minimum qValue direction.
                bary = tem_BaryOfId( inTree, inTree%treeID(posInTree) )
                if (btest(inTree%ElemPropertyBits(posInTree), prp_hasQVal)) then
                  dx = tem_ElemSize( inTree, inTree%treeID(posInTree) ) 
                  point = bary + dx * qOffset(map(minQValiDir), :)
                  ! append the physical points to the growing array of points
                  call append( me  = grwPnts, &
                    &          val = point    )
                else
                  call append( me  = grwPnts, &
                    &          val = bary     )
                end if
                countPoints = countPoints + 1
              end if
            end if !storePnts

            exit ! continue with next element
          end if ! wasAdded

        end if ! boundary_ID == bcIDs
      end do

    end do ! iElem

  end subroutine tem_shape_findElemByBCLabels
  ! ************************************************************************** !


  ! ************************************************************************** !
  subroutine tem_shape_initByLevels( inTree, minLevel, maxLevel, countElems, &
    &                                map2global )
    ! ---------------------------------------------------------------------------
    !> Global mesh from which the elements are identified and then stored to
    type( treelmesh_type ), intent(in) :: inTree
    !> level range of target elements
    integer, intent(in)  :: minLevel, maxLevel
    !> How many elements there will be for each level in the track
    integer, intent(out) :: countElems( globalMaxLevels )
    !> growing array. Elements positions in inTree%treeID
    type(dyn_intArray_type), intent(inout) :: map2global
    ! ---------------------------------------------------------------------------
    integer(kind=long_k)  :: myID, minID, maxID
    integer :: tLevel, dPos, iElem, loc_min, loc_max
    logical :: wasAdded
    ! ---------------------------------------------------------------------------

    loc_min = minLevel
    loc_max = maxLevel
    if ( minLevel > maxLevel ) then
      ! take inverse
      loc_min = maxLevel
      loc_max = minLevel
    end if

    if ( minLevel < 1 )               loc_min = 1
    if ( maxLevel > globalMaxLevels ) loc_max = globalMaxLevels

    call tem_log(3, 'Initializing shapes by elements between level '&
      &        //trim(tem_toStr(loc_min))//' and '//trim(tem_toStr(loc_max)) )

    ! the treeID range is the first ID on min level and the last ID on max level
    minID = tem_firstIdAtLevel( loc_min )
    maxID =  tem_lastIdAtLevel( loc_max )

    ! Loop over all elements in inTree
    do iElem = 1, inTree%nElems

      myID = inTree%treeID(iElem)

      if( (myID >= minID) .and. (myID <= maxID) ) then
        ! Append to treeID list (note that already existing ones are
        ! omitted)
        call append( me       = map2global, &
          &          pos      = dPos,       &
          &          val      = iElem,      &
          &          wasAdded = wasAdded )

        ! Count up if it was added
        if( wasAdded ) then
          tLevel   = tem_levelOf( inTree%treeID(iElem) )
          countElems( tLevel ) = countElems( tLevel ) + 1
        end if ! wasAdded

      end if

    end do ! iElem

  end subroutine tem_shape_initByLevels
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine identify all the elements in inTree that has a certain property
  !! bit, save their positions in inTree into array: map2global,
  !! save the number of these elements into level wise array: countElems
  !! (e.g. for shape kind='property').
  !!
  subroutine tem_shape_initPropElements( propBits, inTree, countElems, map2global )
    ! ---------------------------------------------------------------------------
    !> shape objects on which to work
    integer( kind=long_k ), intent(in) :: propBits
    !> Global mesh from which the elements are identified and then stored to
    type( treelmesh_type ), intent(in) :: inTree
    !> How many elements there will be for each level in the track
    integer, intent( out ) :: countElems( globalMaxLevels )
    !> growing array. Elements positions in inTree%treeID
    type(dyn_intArray_type), intent(inout) :: map2global
    ! ---------------------------------------------------------------------------
    integer(kind=long_k)  :: elemProp, match
    integer :: tLevel, dPos, iElem
    logical :: wasAdded
    ! ---------------------------------------------------------------------------

    write(logUnit(1),*) 'Initializing shapes of elements that have a ' &
      &                 // 'certain property'

    ! Loop over all elements in inTree
    do iElem = 1, inTree%nElems

      ! get its property bits
      elemProp = inTree%elemPropertyBits( iElem )

      ! check whether its property
      match =  iand( propBits, elemProp )

      if( match > 0_long_k ) then

        ! Append to treeID list (note that already existing ones are
        ! omitted)
        call append( me       = map2global, &
          &          pos      = dPos,       &
          &          val      = iElem,      &
          &          wasAdded = wasAdded    )

        ! Count up if it was added
        if( wasAdded ) then
          tLevel   = tem_levelOf( inTree%treeID(iElem) )
          countElems( tLevel ) = countElems( tLevel ) + 1
        end if ! wasAdded

      end if ! match > 0

    end do ! iElem

  end subroutine tem_shape_initPropElements
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Identify elements matching a given shape to build a subTree.
  subroutine tem_shape2subTree( me, iShape, inTree, storePnts, map2global, &
    &                           grwPnts, countElems, countPnts, bcIDs,     &
    &                           bc_prop, stencil                           )
    ! ---------------------------------------------------------------------- !
    !> The shape to identify elements for
    type(tem_shape_type), intent(in) :: me

    !> Numbering of the shape (only for logging)
    integer, intent(in) :: iShape

    !> The tree to look in for elements that match the shape definition
    type(treelmesh_type), intent(in) :: inTree

    !> Whether to store point values
    logical, intent(in) :: storePnts

    !> Mapping to global elements in the tree
    type(dyn_intArray_type), intent(inout) :: map2global

    !> Growing list of Points to be observed
    type(tem_grwPoints_type), intent(inout) :: grwPnts

    !> Number of elements on each level matching the shape
    integer, intent(out) :: countElems(globalMaxLevels)

    !> Number of points to be observed
    integer, intent(out) :: countPnts

    !> Field of boundary ids that are to be tracked
    integer, allocatable, intent(out) :: bcIDs(:)

    !> Boundary condition property to identify boundary elements
    type(tem_bc_prop_type), optional, intent(in) :: bc_prop

    !> Stencil associated with the boundary to find respective neighbors
    type(tem_stencilHeader_type), optional, intent(in) :: stencil
    ! ---------------------------------------------------------------------- !
    logical :: foundAny
    integer :: nShapeElems(globalMaxLevels)
    ! ---------------------------------------------------------------------- !

    foundAny = .false.
    nShapeElems = 0

    select case( me%shapeID )
    case( tem_geometrical_shape )
      ! Use elements intersecting a geometrical object
      write(logUnit(5),*) 'iShape ', iShape, ' is a geometrical shape.'
      call tem_shape_subTreeFromGeomInters( me          = me,          &
        &                                   inTree      = inTree,      &
        &                                   countElems  = nShapeElems, &
        &                                   countPoints = countPnts,   &
        &                                   grwPnts     = grwPnts,     &
        &                                   storePnts   = storePnts,   &
        &                                   map2global  = map2global   )

    case( tem_property_shape )
      ! Only use elements with a certain property
      write(logUnit(5),*) 'iShape ', iShape, ' is a property shape.'
      call tem_shape_initPropElements( me%propBits, inTree,   &
        &                              nShapeElems, map2global )

    case( tem_boundary_shape )
      ! Only use elements belong to certain boundaries
      write(logUnit(5),*) 'iShape ', iShape, ' is a boundary shape.'
      if (present(bc_prop) .and. present(stencil)) then
        call tem_shape_findElemByBCLabels( bcLabels    = me%bcLabels,   &
          &                                cutOffQVal  = me%cutOffQVal, &
          &                                bc_prop     = bc_prop,       &
          &                                foundAny    = foundAny,      &
          &                                map2global  = map2Global,    &
          &                                inTree      = inTree,        &
          &                                countPoints = countPnts,     &
          &                                grwPnts     = grwPnts,       &
          &                                storePnts   = storePnts,     &
          &                                bcIDs       = bcIDs,         &
          &                                stencil     = stencil        )
        if (foundAny) nShapeElems = 1
      else
        call tem_abort('In tem_shape2subTree: Stencil or bc_prop not passed!')
      end if

    case( tem_level_shape )
      write(logUnit(5),*) 'iShape ', iShape, ' is a level shape.'
      call tem_shape_initByLevels( inTree, me%minLevel,   &
        &                          me%maxLevel,           &
        &                          nShapeElems, map2global )
    end select ! shapeID

    countElems = countElems + nShapeElems

  end subroutine tem_shape2subTree
  ! ************************************************************************** !

end module tem_shape_module
! ****************************************************************************** !
