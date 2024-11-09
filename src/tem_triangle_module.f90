! Copyright (c) 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2019 Harald Klimach <harald.klimach@uni-siegen.de>
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
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012, 2015-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013 Daniel Harlacher <d.harlacher@grs-sim.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2015-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2015-2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Daniel Fleischer <daniel.fleischer@student.uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
!
! Parts of this file were written by Harald Klimach, Simon Zimny and Manuel
! Hasert for German Research School for Simulation Sciences GmbH.
!
! Parts of this file were written by Harald Klimach, Kannan Masilamani,
! Daniel Harlacher, Kartik Jain, Verena Krupp, Jiaxing Qi, Peter Vitt,
! Daniel Fleischer, Tobias Girresser and Daniel Petró for University Siegen.
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

! This file contains the source code for growing and dynamic arrays.
! This is used for arrays of primitives (int, long_int, real, ...) as well as
! for arrays of derived datatypes (tem_variable_type,...).
!
! To use these macros include the following to your source file.
!
! Smart growing array (GA) for ?tstring?
! Growing Arrays:
!
! declaration
!
!
! implementation
!

! -----------------------------------------------------------------
! 2d Array, which can grow in second dimension only (GA2d)
! tname ... indicates type of dynamic array (long, int, real, ...)


!
!------------------------------------------------------------------------------
!
! dynamic Arrays:
!
! declaration
!
!
! implementation
!
!> \brief Module to describe triangles. Contains triangle type
!! with vertices of the triangle
module tem_triangle_module
  use env_module,                only: rk, minLength, zeroLength, eps
  use tem_aux_module,            only: tem_abort
  use tem_logging_module,        only: logunit
  use tem_math_module,           only: cross_product3D
  use tem_cube_module,           only: tem_cube_type
  use tem_transformation_module, only: tem_transformation_type

  ! include aotus modules
  use aotus_module,     only: aot_get_val, aoterr_Fatal, aoterr_WrongType,     &
    &                         flu_State
  use aot_table_module, only: aot_table_open, aot_table_close, aot_table_length
  use aot_out_module,   only: aot_out_type, aot_out_val,                       &
    &                         aot_out_open_table, aot_out_close_table

  implicit none

  private

  !> Indicating that the proximity evaluation for this triangle is invalid
  integer, public, parameter :: tem_triangle_close_invalid = 0

  !> Indicating that the closest point in the proximity test resides
  !! within the triangle
  integer, public, parameter :: tem_triangle_close_face = 1

  !> Indicating that the closest point in the proximity test resides
  !! on an edge of the triangle
  integer, public, parameter :: tem_triangle_close_edge = 2

  !> Indicating that the closest point in the proximity test resides
  !! on a vertex of the triangle
  integer, public, parameter :: tem_triangle_close_vertex = 3

  public :: grw_triangleArray_type
  public :: tem_triangle_type
  public :: init, append, truncate, destroy, empty, placeAt
  public :: tem_triangleCubeOverlap
  public :: tem_load_triangle
  public :: tem_triangle_out
  public :: tem_evaluate_normal_triangle
  public :: tem_triangle_normal_proximity

  !> This type contains three vertices of the triangle
  type tem_triangle_type
    !> 1st index contains x,y,z coordinates
    !! and 2nd index the vertex number.
    real(kind=rk) :: nodes(3,3)
    real(kind=rk) :: normal(3)
  end type tem_triangle_type

  !> growing array type for type(tem_triangle_type)
  type grw_trianglearray_type
    integer :: nvals = 0
    integer :: containersize = 0
    type(tem_triangle_type), allocatable :: val(:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga_triangle
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_ga_triangle
  end interface

  !> empty the entries  without changing arrays
  interface empty
    module procedure empty_ga_triangle
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga_triangle
  end interface

  !> insert an element at a given position
  interface placeat
    module procedure placeat_ga_triangle
    module procedure placeat_ga_triangle_vec
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_ga_triangle
    module procedure append_ga_triangle_vec
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga_triangle
  end interface


  !> interface to write out triangles in lua format to a file
  interface tem_triangle_out
    module procedure tem_triangle_out_scal
    module procedure tem_triangle_out_vec
  end interface tem_triangle_out

  !> interface to load triangles
  interface tem_load_triangle
    module procedure tem_load_triangle
    module procedure tem_load_triangle_single
  end interface tem_load_triangle

contains

  ! ****************************************************************************
  !> Load triangle information from config file.
  subroutine tem_load_triangle(me, transform, conf, thandle )
    !--------------------------------------------------------------------------!
    !inferface variables
    !> array of triangles
    type(tem_triangle_type), allocatable, intent(out) :: me(:)
    !> transformation for spatial object
    type(tem_transformation_type), intent(in) :: transform
    !> lua state
    type(flu_state) :: conf
    integer, intent(in) :: thandle !< handle for canonical objects
    !--------------------------------------------------------------------------!
    ! local varaibles
    integer :: tri_handle, tri_subHandle
    integer :: iObj, nObjects
    !--------------------------------------------------------------------------!

    write(logunit(1),*) 'Loading triangle: '

    call aot_table_open(L = conf, parent = thandle, thandle = tri_handle, &
      &                 key = 'object')
    call aot_table_open(L=conf, parent = tri_handle, thandle = tri_subHandle, &
      & pos = 1 )

    if ( tri_subHandle .eq. 0) then
      !object is a single table
      allocate(me(1))
      call aot_table_close(L=conf, thandle=tri_subHandle)
      call tem_load_triangle_single(me(1), transform, conf, tri_handle)
    else
      !object is a multiple table
      call aot_table_close(L=conf, thandle=tri_subHandle)
      nObjects = aot_table_length(L=conf, thandle=tri_handle)
      allocate(me(nObjects))
      do iObj=1,nObjects
        write(logunit(2),*) '  triangle ', iObj
        call aot_table_open(L=conf, parent=tri_handle, thandle=tri_suBHandle,&
          & pos=iObj)
        call tem_load_triangle_single(me(iObj), transform, conf, tri_Subhandle)
        call aot_table_close(L=conf, thandle=tri_subHandle)
        write(logunit(2),*) ''
      end do
    end if

    call aot_table_close(L=conf, thandle=tri_Handle)

  end subroutine tem_load_triangle
  ! ****************************************************************************


  ! ****************************************************************************
  !> Load single triangle from config file
  subroutine tem_load_triangle_single(me, transform, conf, thandle )
    !--------------------------------------------------------------------------!
    !inferface variables
    !> single triangle
    type(tem_triangle_type), intent(out) :: me
    !> transformation for spatial object
    type(tem_transformation_type), intent(in) :: transform
    !> lua state
    type(flu_state) :: conf
    integer, intent(in) :: thandle !< handle for canonical objects
    !--------------------------------------------------------------------------!
    ! local varaibles
    integer :: vError(3), errFatal(3)
    integer :: iNode
    integer :: sub_handle
    !--------------------------------------------------------------------------!
    errFatal = aoterr_fatal

    call aot_table_open(L=conf, parent=thandle, thandle=sub_handle, &
      &                 key='nodes')
    do iNode=1,3
      call aot_get_val(L=conf, thandle=sub_handle, &
        &              val=me%nodes(:,iNode), &
        &              pos=iNode, ErrCode=vError)

      if (any(btest(vError, errFatal))) then
        write(logunit(0),*) ' Error in configuration: triangle node nr ', iNode
        write(logunit(0),*) ' is not defined'
        call tem_abort()
      end if
    end do
    call aot_table_close(L=conf, thandle=sub_handle)

    write(logunit(2),*) '   node 1: ', me%nodes(:,1)
    write(logunit(2),*) '   node 2: ', me%nodes(:,2)
    write(logunit(2),*) '   node 3: ', me%nodes(:,3)

    !apply transformation
    if(transform%active) then
      if(transform%deform%active) then
        do iNode=1,3
          me%nodes(:,iNode) = matmul( transform%deform%matrix, &
            &                                  me%nodes(:,iNode) )
        enddo
      endif
      if(transform%translate%active) then
        do iNode=1,3
          me%nodes(:,iNode) = transform%translate%vec &
            &                         + me%nodes(:,iNode)
        enddo
      endif
    endif

  end subroutine tem_load_triangle_single
  ! ****************************************************************************


  ! ****************************************************************************
  !> Compute, if the triangle intersects the cube.
  function tem_triangleCubeOverlap(triangle, cube) result(overlaps)
    !--------------------------------------------------------------------------!
    type(tem_triangle_type), intent(in) :: triangle
    type(tem_cube_type), intent(in) :: cube
    logical :: overlaps
    !--------------------------------------------------------------------------!
    real(kind=rk) :: halfwidth
    !--------------------------------------------------------------------------!

    ! halfwidth is increased by eps to avoid precision problem in
    ! triangöe box intersection.
    halfwidth = cube%halfwidth + eps
    overlaps = triBoxOverlap_loc( cube%center,                       &
      &                           [halfwidth, halfwidth, halfwidth], &
      &                           triangle%nodes )

  end function tem_triangleCubeOverlap
  ! ****************************************************************************

  ! ****************************************************************************
  !> Compute, the outward normal of a triangle. To work nodes must be in
  !> counter-clockwise order. For STL this is the standard.
  !> It follows pg 136 of Jiri Blaze, CFD Principles and Applications
  subroutine tem_evaluate_normal_triangle(triangle)
    !--------------------------------------------------------------------------!
    type(tem_triangle_type), intent(inout) :: triangle
    !--------------------------------------------------------------------------!
    real(kind=rk) :: outward_vector(3), delta_x(3), delta_y(3), delta_z(3)
    real(kind=rk) :: magnitude_normal
    !--------------------------------------------------------------------------!

    ! Eq.s 5.8
    delta_x(1) = ( triangle%nodes(1,1) - triangle%nodes(1,2) ) * &
      &          ( triangle%nodes(2,1) + triangle%nodes(2,2) )
    delta_x(2) = ( triangle%nodes(1,2) - triangle%nodes(1,3) ) * &
      &          ( triangle%nodes(2,2) + triangle%nodes(2,3) )
    delta_x(3) = ( triangle%nodes(1,3) - triangle%nodes(1,1) ) * &
      &          ( triangle%nodes(2,3) + triangle%nodes(2,1) )
    delta_y(1) = ( triangle%nodes(2,1) - triangle%nodes(2,2) ) * &
      &          ( triangle%nodes(3,1) + triangle%nodes(3,2) )
    delta_y(2) = ( triangle%nodes(2,2) - triangle%nodes(2,3) ) * &
      &          ( triangle%nodes(3,2) + triangle%nodes(3,3) )
    delta_y(3) = ( triangle%nodes(2,3) - triangle%nodes(2,1) ) * &
      &          ( triangle%nodes(3,3) + triangle%nodes(3,1) )
    delta_z(1) = ( triangle%nodes(3,1) - triangle%nodes(3,2) ) * &
      &          ( triangle%nodes(1,1) + triangle%nodes(1,2) )
    delta_z(2) = ( triangle%nodes(3,2) - triangle%nodes(3,3) ) * &
      &          ( triangle%nodes(1,2) + triangle%nodes(1,3) )
    delta_z(3) = ( triangle%nodes(3,3) - triangle%nodes(3,1) ) * &
      &          ( triangle%nodes(1,3) + triangle%nodes(1,1) )

    ! Eq.s 5.9
    outward_vector(1) = 0.5_rk * ( delta_y(1) + delta_y(2) + delta_y(3) )
    outward_vector(2) = 0.5_rk * ( delta_z(1) + delta_z(2) + delta_z(3) )
    outward_vector(3) = 0.5_rk * ( delta_x(1) + delta_x(2) + delta_x(3) )

    ! Eq.s 5.12
    magnitude_normal = sqrt( outward_vector(1)**2 + outward_vector(2)**2 + outward_vector(3)**2 )

    ! normalization of the outward_vector
    triangle%normal(:) = outward_vector(:) / magnitude_normal

  end subroutine tem_evaluate_normal_triangle
  ! ****************************************************************************

  ! ***************************************************************************
  !> This routine checks for triangle box overlap
  !!
  !! this routine is conversion of c-code tribox3.c triBoxOverlap function.
  !! use separating axis theorem to test overlap between triangle and box
  !!  need to test for overlap in these directions:
  !!  1) the {x,y,z}-directions (actually, since we use the AABB of the triangle
  !!     we do not even need to test these)
  !!  2) normal of the triangle
  !!  3) separating axis test. crossproduct(edge from tri, {x,y,z}-directin)
  !!     this gives 3x3=9 more 7tests
  !! This code is available at:
  !! http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tribox3.txt
  function triBoxOverlap_loc( boxCenter, boxHalfwidth, triNodes ) &
    &                         result(overlaps)
    !--------------------------------------------------------------------------!
    !> box center
    real(kind=rk), intent(in) :: boxCenter(3)
    !> halfwidth of the box
    real(kind=rk), intent(in) :: boxHalfwidth(3)
    !> nodes of the triangle
    !! 1st index denote x,y,z coordinates and
    !! 2nd index denote nodes
    real(kind=rk), intent(in) :: triNodes(3,3)
    logical :: overlaps
    !--------------------------------------------------------------------------!
    real(kind=rk) :: nodes(3,3)
    real(kind=rk) :: dirVec1(3), dirVec2(3), dirVec3(3)
    real(kind=rk) :: mmin, mmax
    real(kind=rk) :: normal(3), edge1(3), edge2(3), edge3(3)
    !--------------------------------------------------------------------------!
    overlaps = .false.

    dirVec1 = (/ 1.0_rk, 0.0_rk, 0.0_rk /)
    dirVec2 = (/ 0.0_rk, 1.0_rk, 0.0_rk /)
    dirVec3 = (/ 0.0_rk, 0.0_rk, 1.0_rk /)

    ! move everything so that the box center is in (0,0,0)
    nodes(:,1) = triNodes(:,1) - boxCenter
    nodes(:,2) = triNodes(:,2) - boxCenter
    nodes(:,3) = triNodes(:,3) - boxCenter

    ! compute triangle edges
    edge1 = nodes(:,2) - nodes(:,1)
    edge2 = nodes(:,3) - nodes(:,2)
    edge3 = nodes(:,1) - nodes(:,3)


    ! separating axis test
    if(axistest(dirVec1, edge1, nodes, boxhalfwidth)) return
    if(axistest(dirVec1, edge2, nodes, boxhalfwidth)) return
    if(axistest(dirVec1, edge3, nodes, boxhalfwidth)) return

    if(axistest(dirVec2, edge1, nodes, boxhalfwidth)) return
    if(axistest(dirVec2, edge2, nodes, boxhalfwidth)) return
    if(axistest(dirVec2, edge3, nodes, boxhalfwidth)) return

    if(axistest(dirVec3, edge1, nodes, boxhalfwidth)) return
    if(axistest(dirVec3, edge2, nodes, boxhalfwidth)) return
    if(axistest(dirVec3, edge3, nodes, boxhalfwidth)) return


    ! Bullet 1:
    !  first test overlap in the {x,y,z}-directions
    !  find min, max of the triangle each direction, and test for overlap in
    !  that direction -- this is equivalent to testing a minimal AABB around
    !  the triangle against the AABB
    ! test in X-direction
    mmin = min( nodes(1,1), nodes(1,2), nodes(1,3) )
    mmax = max( nodes(1,1), nodes(1,2), nodes(1,3) )
    if ( mmin > boxhalfwidth(1) .or. mmax < -boxhalfwidth(1) ) return

    ! test in Y-direction
    mmin = min( nodes(2,1), nodes(2,2), nodes(2,3) )
    mmax = max( nodes(2,1), nodes(2,2), nodes(2,3) )
    if ( mmin > boxhalfwidth(2) .or. mmax < -boxhalfwidth(2) ) return

    ! test in Z-direction
    mmin = min( nodes(3,1), nodes(3,2), nodes(3,3) )
    mmax = max( nodes(3,1), nodes(3,2), nodes(3,3) )
    if ( mmin > boxhalfwidth(3) .or. mmax < -boxhalfwidth(3) ) return


    ! Bullet 2:
    !  test if the box intersects the plane of the triangle
    !  compute plane equation of triangle: normal*x+d=0
    normal = cross_product3D(edge1,edge2)

    if( .not. planeBoxOverlap(normal, nodes(:,1), boxhalfwidth)) return

    overlaps = .true.

  end function triBoxOverlap_loc
  ! ***************************************************************************

  !> This function check whether there is separating axis between triangle
  !! and box.
  !! This function can be optimized by explicitly providing cross_product result
  !! see example:
  !! http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tribox3.txt
  function Axistest( dirVec, edge, nodes, boxhalfwidth ) result (success)
    !--------------------------------------------------------------------------!
    real(kind=rk), intent(in) :: dirVec(3)
    real(kind=rk), intent(in) :: edge(3)
    real(kind=rk), intent(in) :: nodes(3,3)
    real(kind=rk), intent(in) :: boxhalfwidth(3)
    logical :: success
    !--------------------------------------------------------------------------!
    real(kind=rk) :: vecA(3), p1, p2, p3, pmin, pmax, rad
    !--------------------------------------------------------------------------!

    success = .true.

    ! a_(i,j) = e_i x f_j
    vecA = cross_product3D( dirVec, edge )

    p1 = dot_product( vecA, nodes(:,1) )
    p2 = dot_product( vecA, nodes(:,2) )
    p3 = dot_product( vecA, nodes(:,3) )

    pmin = min(p1, p2, p3)
    pmax = max(p1, p2, p3)

    rad = abs(vecA(1)) * boxhalfwidth(1) &
      & + abs(vecA(2)) * boxhalfwidth(2) &
      & + abs(vecA(3)) * boxhalfwidth(3)

    if (pmin > rad .or. pmax < -rad) return

    success = .false.

  end function Axistest
  ! ***************************************************************************

  ! ***************************************************************************
  !> This routine checks for plane box overlap
  !> this routine is conversion of c-code tribox3.c planeBoxOverlap function
  function planeBoxOverlap( normal, origin, boxhalfwidth ) result(overlaps)
    !--------------------------------------------------------------------------!
    !> normal direction of the plane
    real(kind=rk), intent(in) :: normal(3)
    !> origin of the plane
    real(kind=rk), intent(in) :: origin(3)
    !> halfwidth of the box
    real(kind=rk), intent(in) :: boxhalfwidth(3)
    logical :: overlaps
    !--------------------------------------------------------------------------!
    integer :: iDir
    real(kind=rk) :: vmin(3), vmax(3), tmp
    !--------------------------------------------------------------------------!
    overlaps = .false.
    ! find min and max of x,y,z of distance between boxhalfwidth and origin
    ! depends on the direction of the plane normal
    do iDir=1,3
      tmp = origin(iDir)
      if(normal(iDir) > 0.0_rk) then
        vmin(iDir) = - boxhalfwidth(iDir) - tmp
        vmax(iDir) =   boxhalfwidth(iDir) - tmp
      else
        vmin(iDir) =   boxhalfwidth(iDir) - tmp
        vmax(iDir) = - boxhalfwidth(iDir) - tmp
      endif
    end do

    if(dot_product(normal, vmin) > 0.0_rk) return

    if(dot_product(normal, vmax) >= 0.0_rk) overlaps = .true.

  end function planeBoxOverlap
  ! ***************************************************************************

  ! ************************************************************************** !
  !> Compute the point of the triangle closest to the given `point`.
  !!
  !! This function computes the point of the triangle with the smallest distance
  !! to the provided point.
  !! Input parameters:
  !!
  !! * `me` - the triangle to check
  !! * `point` - the point for which the closest neighbor point within the
  !!             triangle is to be found
  !!
  !! Output parameters:
  !!
  !! * `closest` - point of the triangle which is the closest to the given
  !!               `point`
  !! * `closekind` - there are three different kinds of closest point:
  !!                 - it may be inside the triangle
  !!                 - it may be on one of the edges
  !!                 - it may be one of the three vertices
  !!                 A fourth case arises with degenerate triangles for which
  !!                 no proper distance can be computed, if this is set, the
  !!                 other values have no proper values assigned to them.
  !!                 This parameter returns which kind of closest point we have.
  !! * `distance` - the smallest distance between `point` and the triangle (`me`)
  !! * `normal` - depending on the `closekind` this returns:
  !!              - the normal vector on the triangle pointing from `closest`
  !!                to point (has the length `distance`) if closest is withing
  !!                the triangle
  !!              - the unit normal vector to the side of `point` if `closest`
  !!                is on an edge or vertex of the triangle
  subroutine tem_triangle_normal_proximity(me, point, closest, closekind, &
    &                                      distance, normal               )
    ! ---------------------------------------------------------------------- !
    type(tem_triangle_type), intent(in) :: me
    real(kind=rk), intent(in) :: point(3)
    real(kind=rk), intent(out) :: closest(3)
    integer, intent(out) :: closekind
    real(kind=rk), intent(out) :: distance
    real(kind=rk), intent(out) :: normal(3)
    ! ---------------------------------------------------------------------- !
    integer :: min_edge, min_vertex

    real(kind=rk) :: edge(3,3) ! The three edges of the triangle
    real(kind=rk) :: edgelen_squared(3) ! Squared length the edges
    real(kind=rk) :: edgescale(3)
    real(kind=rk) :: edgeprod1_2 ! Dot Product of first two edges
    real(kind=rk) :: tri2p(3) ! Vector pointing from the triangles first vertex
                              ! to the given point.
    real(kind=rk) :: plane_point(3) ! Point projected onto plane of triangle
    real(kind=rk) :: edge_point(3,3) ! Projection of point onto the three edges
    real(kind=rk) :: plane_vector(3) ! Vector from first node to projected point
    real(kind=rk) :: tri_prod(3) ! Point projections onto edges
    real(kind=rk) :: tri_coord(2) ! Point in coordinates of first and second
                                  ! edge
    real(kind=rk) :: vertex_distance(3) ! Distances to each of the three
                                        ! vertices
    real(kind=rk) :: edge_distance(3) ! Distances to each of the three edges
    real(kind=rk) :: normlen_squared
    real(kind=rk) :: normlen
    real(kind=rk) :: det
    ! ---------------------------------------------------------------------- !

    closekind = tem_triangle_close_invalid
    closest = 0.0_rk
    distance = huge(distance)

    ! First: project the point onto the plane of the triangle
    !        if the projected point is inside the triangle the projected point
    !        is the `closest`.
    edge(:,1) = me%nodes(:,2) - me%nodes(:,1)
    edge(:,2) = me%nodes(:,3) - me%nodes(:,1)
    edge(:,3) = me%nodes(:,3) - me%nodes(:,2)

    normal = cross_product3D(edge(:,1), edge(:,2))
    normlen_squared = normal(1)**2 + normal(2)**2 + normal(3)**2

    validtri: if (normlen_squared > epsilon(normlen_squared)) then

      normlen = sqrt(normlen_squared)
      normal = normal / normlen

      tri2p = point - me%nodes(:,1)

      edgelen_squared(1) = dot_product(edge(:,1), edge(:,1))
      edgelen_squared(2) = dot_product(edge(:,2), edge(:,2))
      edgeprod1_2 = dot_product(edge(:,1), edge(:,2))

      plane_point = point - dot_product(tri2p, normal)
      plane_vector = plane_point - me%nodes(:,1)

      tri_prod(1) = dot_product(plane_vector, edge(:,1))
      tri_prod(2) = dot_product(plane_vector, edge(:,2))

      det = edgeprod1_2**2 - (edgelen_squared(1)*edgelen_squared(2))

      tri_coord(1) = (  (edgeprod1_2 * tri_prod(2))        &
        &             - (edgelen_squared(2) * tri_prod(1)) )
      tri_coord(2) = (  (edgeprod1_2 * tri_prod(1))        &
        &             - (edgelen_squared(1) * tri_prod(2)) )

      if (tri_coord(1) > tiny(0.0_rk) .and. tri_coord(2) > tiny(0.0_rk) &
        & .and. (tri_coord(1) + tri_coord(2)) < 1.0_rk                  ) then

        ! The projected point lies within the triangle
        closekind = tem_triangle_close_face
        closest = plane_point
        normal = point - plane_point
        distance = sqrt(dot_product(normal, normal))

      else

        edgelen_squared(3) = dot_product(edge(:,3), edge(:,3))
        tri_prod(3) = dot_product(plane_point - me%nodes(:,2), edge(:,3))

        ! The projected point lies outside the triangle (or on its boundary).
        ! Thus we check for the closest point on the boundary
        ! (vertices and edges).
        vertex_distance(1) = sqrt(sum((point - me%nodes(:,1))**2))
        vertex_distance(2) = sqrt(sum((point - me%nodes(:,2))**2))
        vertex_distance(3) = sqrt(sum((point - me%nodes(:,3))**2))
        min_vertex = minloc(vertex_distance, 1)

        edgescale = tri_prod / edgelen_squared

        ! Limit the points to the length of the edge
        edgescale = min(edgescale, 1.0_rk)
        edgescale = max(edgescale, 0.0_rk)
        edge_point(:,1) = me%nodes(:,1) + edgescale(1) * edge(:,1)
        edge_point(:,2) = me%nodes(:,1) + edgescale(2) * edge(:,2)
        edge_point(:,3) = me%nodes(:,2) + edgescale(3) * edge(:,3)

        edge_distance(1) = sqrt(sum((point - edge_point(:,1))**2))
        edge_distance(2) = sqrt(sum((point - edge_point(:,2))**2))
        edge_distance(3) = sqrt(sum((point - edge_point(:,3))**2))
        min_edge = minloc(vertex_distance, 1)

        closekind = tem_triangle_close_edge
        closest = edge_point(:, min_edge)
        distance = edge_distance(min_edge)

        if (vertex_distance(min_vertex) < distance + epsilon(distance)) then
          closekind = tem_triangle_close_vertex
          closest = me%nodes(:, min_vertex)
          distance = vertex_distance(min_vertex)
        end if

        if (dot_product(normal, point-closest) < 0.0_rk) then
          ! flip the normal to point to the side of the given point
          normal = -normal
        end if

      end if

    end if validtri

  end subroutine tem_triangle_normal_proximity
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Write out an array of triangles in lua format
  !! Only if nTriangles <= 10
  !!
  subroutine tem_triangle_out_vec( me, conf )
    ! --------------------------------------------------------------------------
    !> triangle types to write out
    type( tem_triangle_type ), intent(in) :: me(:)
    !> Aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    ! --------------------------------------------------------------------------
    ! counter
    integer :: i
    ! --------------------------------------------------------------------------

    ! create a table with name triangle
    call aot_out_open_table( put_conf = conf, tname = 'triangle' )

    if (size(me)<=10) then
      do i = 1, size(me)
        call tem_triangle_out_scal( me(i), conf )
      end do
    else
      call aot_out_val( put_conf = conf, vname = 'nTriangles', val = size(me) )
      call aot_out_val( put_conf = conf,                                  &
        &               val = 'Write triangle is limited to 10 triangles' )
    end if

    call aot_out_close_table( put_conf = conf )

  end subroutine tem_triangle_out_vec
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Write out a triangle shape in lua format
  !!
  subroutine tem_triangle_out_scal( me, conf )
    ! --------------------------------------------------------------------------
    !> triangle types to write out
    type( tem_triangle_type ), intent(in) :: me
    !> Aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    ! --------------------------------------------------------------------------
    integer :: iNode
    ! ---------------------------------------------------------------------------

    ! create a table with name triangle if not exist
    if( conf%level == 0 ) then
      call aot_out_open_table( put_conf = conf, tname = 'triangle' )
    else
      call aot_out_open_table( put_conf = conf )
    end if

    !write nodes
    call aot_out_open_table( put_conf = conf, tname = 'nodes' )
    do iNode=1,3
      call aot_out_val( put_conf = conf, val = me%nodes(:,iNode) )
    enddo
    call aot_out_close_table( put_conf = conf )

    call aot_out_close_table( put_conf = conf )

  end subroutine tem_triangle_out_scal
  ! ************************************************************************** !

  subroutine init_ga_triangle(me, length)
    type(grw_trianglearray_type), intent(out) :: me !< dynamic array to init
    integer, intent(in), optional :: length !< initial length of the container

    if (present(length)) then
      me%containersize = length
    else
      me%containersize = zerolength
    end if
    ! deallocate ...
    if( allocated( me%val ))     &
      deallocate(me%val)
    ! ... and reallocate
    allocate(me%val(me%containersize))
    me%nvals = 0

  end subroutine init_ga_triangle

  subroutine destroy_ga_triangle(me)
    type(grw_trianglearray_type), intent(inout) :: me !< dynamic array to destroy

    me%containersize = 0
    me%nvals = 0
    if( allocated( me%val ) ) deallocate(me%val)

  end subroutine destroy_ga_triangle


  subroutine truncate_ga_triangle(me)
    !------------------------------------------------------------------------
    type(grw_trianglearray_type) :: me !< array to truncate
    !------------------------------------------------------------------------
    type(tem_triangle_type), allocatable :: tarray(:)
    !------------------------------------------------------------------------
    integer :: ii
    !------------------------------------------------------------------------

    ! nothing to do if container size is not larger than the number of values
    ! in the array.
    if (me%containersize > me%nvals) then
      allocate(tarray(me%nvals))
      do ii = 1, me%nvals
        tarray(ii) = me%val(ii)
      end do
      call move_alloc(tarray, me%val)
      me%containersize = me%nvals
    end if

  end subroutine truncate_ga_triangle


  subroutine empty_ga_triangle(me)
    !------------------------------------------------------------------------
    type(grw_trianglearray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------

    me%nvals = 0

  end subroutine empty_ga_triangle

  !> adds the value to a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! element at the requested position will be replaced.
  subroutine placeat_ga_triangle(me, val, pos, length)
    type(grw_trianglearray_type) :: me !< array to place the value into
    type(tem_triangle_type), intent(in) :: val !< value to place at the given position
    integer, intent(in) :: pos !< predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length


    ! value to append is larger than all existing ones,
    ! just put it to the end of the list, this captures
    ! also the case of empty lists.
    ! in this case foundpos = me%nvals + 1 holds.
    if (pos > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand(me = me, pos = pos, length = length)
    end if

    me%nvals = max( pos, me%nvals )
    me%val(pos) = val

  end subroutine placeat_ga_triangle


  !> adds the values starting from a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! elements starting from the requested position will be replaced up to
  !! the element at position `pos + size(val) - 1`.
  subroutine placeat_ga_triangle_vec(me, val, pos, length)
    type(grw_trianglearray_type) :: me !< array to append the value to
    type(tem_triangle_type), intent(in) :: val(:) !< values to append
    integer, intent(in) :: pos !< predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    integer :: ub, ii

    if (me%nvals == huge(me%nvals)) then
      write(*,*) "reached end of integer range for growing array!"
      write(*,*) "aborting!!"
      stop
    end if

    ub = pos + size(val) - 1

    if (ub > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand(me = me, pos = ub, length = length)
    end if

    me%nvals = max( ub, me%nvals )
    do ii = pos, ub
      me%val(ii) = val(1+ii-pos)
    end do

  end subroutine placeat_ga_triangle_vec


  subroutine append_ga_triangle(me, val, length)
    type(grw_trianglearray_type) :: me !< array to append the value to
    type(tem_triangle_type), intent(in) :: val !< value to append
    !> optional length to expand the array
    integer, intent(in), optional :: length

    ! value to append is larger than all existing ones,
    ! just put it to the end of the list, this captures
    ! also the case of empty lists.
    ! in this case foundpos = me%nvals + 1 holds.
    if (me%nvals+1 > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand(me = me, length = length)
    end if

    me%nvals = me%nvals+1
    me%val(me%nvals) = val

  end subroutine append_ga_triangle

  subroutine append_ga_triangle_vec(me, val, length)
    type(grw_trianglearray_type) :: me !< array to append the value to
    type(tem_triangle_type), intent(in) :: val(:) !< values to append
    !> optional length to expand the array
    integer, intent(in), optional :: length

    integer :: lb, ub, ii

    if (me%nvals == huge(me%nvals)) then
      write(*,*) "reached end of integer range for growing array!"
      write(*,*) "aborting!!"
      stop
    end if

    lb = me%nvals + 1
    ub = lb + size(val) - 1

    if (ub > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand(me = me, pos = ub, length = length)
    end if

    me%nvals = max( ub, me%nvals )
    do ii = lb, ub
      me%val(ii) = val(1+ii-lb)
    end do

  end subroutine append_ga_triangle_vec


  subroutine expand_ga_triangle(me, pos, length)
    type(grw_trianglearray_type) :: me !< array to resize
    integer, intent(in), optional :: pos !< optional predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    type(tem_triangle_type), allocatable :: swpval(:)
    integer :: explen, ii

    explen = 0
    ! increase the container by the requested length of double it
    if( present(length) ) then
      explen = max( length, minlength )
    else
      ! set the global minimum length, if doubling would be smaller than that
      explen = max(me%containersize, minlength)
    end if

    ! if a position is given, increase the container to at least the size to
    ! fit the position.
    if( present(pos) ) explen = max(explen, pos-me%containersize)

    ! if the current size plus explen exceeds the max container size,
    ! reduce the size to the max container size.
    if( (huge(me%containersize) - explen) <= me%containersize) then
      ! set max container size
      me%containersize = huge(me%containersize)
    else
      ! set the new container size
      me%containersize = me%containersize + explen
    end if

    if ( me%nvals > 0 ) then
      allocate(swpval(me%containersize))
      do ii = 1, me%nvals
        swpval(ii) = me%val(ii)
      end do
      call move_alloc( swpval, me%val )
    else ! me%nvals == 0
      if ( allocated(me%val) ) deallocate( me%val )
      allocate( me%val(me%containersize) )
    end if

  end subroutine expand_ga_triangle


end module tem_triangle_module

!> \page triangle Triangle
!! Triangles are defined by three nodes. It is possible to define triangles
!! as a geometry kind. Also, the planes and hollow boxes are converted to
!! triangles internally.
!!
!! Valid definition:
!! \li Single triangle
!! \verbatim
!! geometry={
!!   kind='triangle',
!!     object={
!!       nodes = {             -- triangle 1
!!         {0.0, 0.0, 0.5},    -- node 1
!!         {-0.5, -0.5, -0.5}, -- node 2
!!         {0.5, -0.5, -0.5}      -- node 3
!!       }
!!     }
!! }
!! \endverbatim
!!
!! \li Multiple triangles
!! \verbatim
!! geometry={
!!   kind='triangle',
!!     object={
!!       {
!!         nodes = {             -- triangle 1
!!           {0.0, 0.0, 0.5},    -- node 1
!!           {-0.5, -0.5, -0.5}, -- node 2
!!           {0.5, -0.5, -0.5}      -- node 3
!!         }
!!       },
!!       {
!!         nodes = {             -- triangle 2
!!           {0.0, 0.0, 0.5},    -- node 1
!!           {-0.5, -0.5, -0.5}, -- node 2
!!           {-0.5, 0.5, -0.5}      -- node 3
!!         }
!!       }
!!     }
!! }
!! \endverbatim
!! \n\n
!! A pyramid has created with four triangles and a plane. The corresponding
!! seeder code is given below:
!! \include testsuite/triangle/seeder.lua
!! \n\n
!! Mesh with pyramid generated by the seeder file:
!! \image html triangle.png
!! \image html triangle_withedges.png
