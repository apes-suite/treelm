! Copyright (c) 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> summary: This module contains lines description and growing array
!! author: Kannan Masilamani
!! for lines

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

module tem_line_module
  use env_module,          only: rk, minLength, zeroLength, eps
  use tem_math_module,     only: cross_product3D
  use tem_cube_module,     only: tem_cube_type
  use tem_triangle_module, only: tem_triangle_type
  use tem_point_module,    only: tem_point_type
  use tem_logging_module,  only: logUnit, tem_toStr
  use tem_float_module,    only: operator(.flt.), &
    &                            operator(.fge.), &
    &                            operator(.feq.)

  implicit none

  private

  public :: grw_lineArray_type
  public :: tem_line_type
  public :: tem_lineCubeOverlap
  public :: init, append, truncate, destroy, empty, placeAt
  public :: intersect_RayTriangle
  public :: fraction_PointLine

  !> This type contains line definition i.e origin and
  !! vector which defines the direction of the line
  type tem_line_type
    real(kind=rk) :: origin(3) !< line origin
    real(kind=rk) :: vec(3) !< vector which defines direction of the line
  end type tem_line_type

  !> growing array type for type(tem_line_type)
  type grw_linearray_type
    integer :: nvals = 0
    integer :: containersize = 0
    type(tem_line_type), allocatable :: val(:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga_line
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_ga_line
  end interface

  !> empty the entries  without changing arrays
  interface empty
    module procedure empty_ga_line
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga_line
  end interface

  !> insert an element at a given position
  interface placeat
    module procedure placeat_ga_line
    module procedure placeat_ga_line_vec
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_ga_line
    module procedure append_ga_line_vec
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga_line
  end interface


contains

! *****************************************************************************!
  !> Function computes intersection of line with cube
  !!
  !! If optional argument pntIntersect contains the intersection point
  !! of the line with cube
  function tem_lineCubeOverlap( line, cube, pntIntersect ) result(overlap)
    ! ---------------------------------------------------------------------------!
    !> line segment to check for intersection
    type(tem_line_type), intent(in) :: line
    !> cube to intersect with
    type(tem_cube_type), intent(in) :: cube
    !> intersection point if there is intersection
    real(kind=rk), optional, intent(out) :: pntIntersect(3)
    logical :: overlap
    ! ---------------------------------------------------------------------------!
    real(kind=rk) :: proj
    real(kind=rk) :: loc_pntIntersect(3)
    ! ---------------------------------------------------------------------------!

    overlap = .false.
    !check whether line is intersect the cube by rayCubeOverlap test
    !then check whether intersected point is within the line segment
    if(rayCubeOverlap( line, cube, loc_pntIntersect )) then
      !project the intersected point on the line
      !and return true only if intersected point is
      !within line segment length
      !The point is inside the line segment if the
      !projected value is >= 0 and < 1.
      proj = dot_product(loc_pntIntersect-line%origin, line%vec) &
        &  / dot_product(line%vec, line%vec)
      overlap = (proj >= 0.0_rk) .and. (proj < 1.0_rk)
    endif

    if(present(pntIntersect))  pntIntersect = loc_pntIntersect

  end function tem_lineCubeOverlap
! ******************************************************************************!

! *****************************************************************************!
  !> Function computes intersection of ray with cube
  !!
  !! The algorithm for lineCubeOverlap used in this function is
  !! taken from
  !! http://www.siggraph.org/education/materials/HyperGraph/raytrace/
  !! rtinter3.htm
  !! http://gamedev.stackexchange.com/questions/18436/
  !! most-efficient-aabb-vs-ray-collision-algorithms
  function rayCubeOverlap( line, cube, pntIntersect ) result(overlap)
    ! ---------------------------------------------------------------------------!
    !> line segment to check for interection
    type(tem_line_type), intent(in) :: line
    !> cube to check intersection of line
    type(tem_cube_type), intent(in) :: cube
    !> intersection point if there is intersection
    real(kind=rk), optional, intent(out) :: pntIntersect(3)
    logical :: overlap
    ! ---------------------------------------------------------------------------!
    integer :: i
    real(kind=rk) :: t_near, t_far
    real(kind=rk) :: T_1, T_2, tmp
    ! ---------------------------------------------------------------------------!

    !initialize near point and var point
    t_near = 0.0_rk
    t_far = huge(t_far)

    dirLoop: do i=1,3 !x,y,z
      if (line%vec(i) .feq. 0._rk) then
        !line parallel to planes in this direction.
        !Line exactly on the cube origin is considered as overlap
        if ( (line%origin(i) < cube%origin(i)) &
          &  .or. (line%origin(i) >= cube%endPnt(i)) ) then
          !parallel and outside cube : no intersection possible
          overlap = .false.
          return
        end if
      else
        !line not parallel to cube

        !1st intersection point on one side of the cube plane
        T_1 = (cube%origin(i) - line%origin(i)) / line%vec(i)

        !2nd intersection point on one side of the cube plane
        T_2 = (cube%endPnt(i) - line%origin(i)) / line%vec(i)

        if (T_1 > T_2) then
          ! we want T_1 to hold values for intersection with near plane
          tmp = T_2
          T_2 = T_1
          T_1 = tmp
        end if

        if (T_1 > t_near) t_near = T_1

        if (T_2 < t_far) t_far = T_2

        if ( (t_near > t_far) .or. (t_far < 0) ) then
          overlap = .false.
          return
        end if
      end if
    end do dirLoop

    !point of intersection
    if(present(pntIntersect)) then
      pntIntersect = line%origin + t_near * line%vec
    endif
    !If we made it here, there is an intesection
    overlap = .true.

  end function rayCubeOverlap
! ******************************************************************************!

! *****************************************************************************!
  !> Function computes intersection of ray with triangle
  !!
  !! http://geomalgorithms.com/a06-_intersect-2.html
  !! intersect_RayTriangle(): intersect a ray with a 3D triangle
  !!    Input:  a ray R, and a triangle T
  !!    Output: *I = intersection point (when it exists)
  !!    Return: -1 = triangle is degenerate (a segment or point)
  !!             0 = disjoint (no intersect)
  !!             1 = intersect in unique point I1
  !!             2 = are in the same planeint
  !! todo: when line lies in triangle, need to treat properly
  function intersect_RayTriangle( line, triangle, intersect_p ) result(isIntersect)
    ! ---------------------------------------------------------------------------!
    !> line segment to check for interection
    type(tem_line_type), intent(in) :: line
    !> cube to check intersection of line
    type(tem_triangle_type), intent(in) :: triangle
    !> intersection point if there is intersection
    type( tem_point_type), optional, intent(out) :: intersect_p
    logical :: isIntersect
    ! ---------------------------------------------------------------------------!
    real(kind=rk) :: u(3), v(3), n(3)  ! triangle vectors and normal vector
    real(kind=rk) :: dir(3), w0(3), w(3)    ! ray vectors
    real(kind=rk) :: r, a, b       ! params to calc ray-plane intersect
    real(kind=rk) :: uu, uv, vv, wu, wv, D
    real(kind=rk) :: s, t
    real(kind=rk) :: temp_p(3)
    ! ---------------------------------------------------------------------------!
    isIntersect = .false. ! set not intersect as default

    ! get triangle edge vectors u & v, and plane normal n
    u(:) = triangle%nodes(:,2) - triangle%nodes(:,1)
    v(:) = triangle%nodes(:,3) - triangle%nodes(:,1)
    n = cross_product3D( u, v )

    ! triangle is degenerate
    ! do not deal with this case
    if (all(n .feq. 0._rk)) return

    dir = line%vec   ! ray direction vector
    w0 = line%origin - triangle%nodes(:,1)
    a = -dot_product( n, w0);
    b =  dot_product( n, dir);

    ! if ray parallel to triangle plane, treat it as no intersecting
    if (abs(b) < eps .and. abs(a) > tiny(a)) return
!    {   if (a == 0.0_rk)         ! ray lies in triangle plane
!          return 2;
!        else return 0;             ! ray disjoint from plane
!    }

    if (abs(b) < eps) then
      ! origin is very close to triangle plane, but parallel.
      ! Just use the origin itself as point to check.
      r = 0.0_rk
    else
      ! get intersect point of ray with triangle plane
      r = a / b
      if (r < 0._rk ) then     ! ray goes away from triangle
        return
      endif
    end if
    ! for a segment, also test if (r > 1.0) => no intersect

    ! intersect point of ray and plane
    temp_p = line%origin + r * dir

    ! is I inside T?
    uu = dot_product(u,u)
    uv = dot_product(u,v)
    vv = dot_product(v,v)
    w = temp_p - triangle%nodes(:,1)
    wu = dot_product(w,u)
    wv = dot_product(w,v)
    D = uv * uv - uu * vv

    ! get and test parametric coords
    s = (uv * wv - vv * wu) / D
    ! point is outside triangle
    if (s < 0._rk-eps .or. s > 1._rk+eps) then
      return
    endif

    t = (uv * wu - uu * wv) / D
    ! point is outside triangle
    if (t < 0._rk-eps .or. (s + t) > 1._rk+eps) then
      return
    endif

    isIntersect = .true. ! point is inside triangle
    if( present( intersect_p ) ) then
      intersect_p%coord = temp_p
    endif

  end function intersect_RayTriangle
! ******************************************************************************!

! ******************************************************************************!
  !> This evaluates relative distance of given point on line
  function fraction_PointLine( point, line ) result( frac )
    type( tem_point_type ), intent(in) :: point
    type( tem_line_type  ), intent(in) :: line
    real(kind=rk) :: numerator, denominator, frac

    numerator =   ( point%coord(1) - line%origin(1) ) &
      &         * ( point%coord(1) - line%origin(1) ) &
      &         + ( point%coord(2) - line%origin(2) ) &
      &         * ( point%coord(2) - line%origin(2) ) &
      &         + ( point%coord(3) - line%origin(3) ) &
      &         * ( point%coord(3) - line%origin(3) )

    denominator =   line%vec(1) * line%vec(1) &
      &           + line%vec(2) * line%vec(2) &
      &           + line%vec(3) * line%vec(3)

    frac = sqrt( numerator / denominator )
  end function fraction_PointLine
! ******************************************************************************!

  subroutine init_ga_line(me, length)
    type(grw_linearray_type), intent(out) :: me !< dynamic array to init
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

  end subroutine init_ga_line

  subroutine destroy_ga_line(me)
    type(grw_linearray_type), intent(inout) :: me !< dynamic array to destroy

    me%containersize = 0
    me%nvals = 0
    if( allocated( me%val ) ) deallocate(me%val)

  end subroutine destroy_ga_line


  subroutine truncate_ga_line(me)
    !------------------------------------------------------------------------
    type(grw_linearray_type) :: me !< array to truncate
    !------------------------------------------------------------------------
    type(tem_line_type), allocatable :: tarray(:)
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

  end subroutine truncate_ga_line


  subroutine empty_ga_line(me)
    !------------------------------------------------------------------------
    type(grw_linearray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------

    me%nvals = 0

  end subroutine empty_ga_line

  !> adds the value to a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! element at the requested position will be replaced.
  subroutine placeat_ga_line(me, val, pos, length)
    type(grw_linearray_type) :: me !< array to place the value into
    type(tem_line_type), intent(in) :: val !< value to place at the given position
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

  end subroutine placeat_ga_line


  !> adds the values starting from a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! elements starting from the requested position will be replaced up to
  !! the element at position `pos + size(val) - 1`.
  subroutine placeat_ga_line_vec(me, val, pos, length)
    type(grw_linearray_type) :: me !< array to append the value to
    type(tem_line_type), intent(in) :: val(:) !< values to append
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

  end subroutine placeat_ga_line_vec


  subroutine append_ga_line(me, val, length)
    type(grw_linearray_type) :: me !< array to append the value to
    type(tem_line_type), intent(in) :: val !< value to append
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

  end subroutine append_ga_line

  subroutine append_ga_line_vec(me, val, length)
    type(grw_linearray_type) :: me !< array to append the value to
    type(tem_line_type), intent(in) :: val(:) !< values to append
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

  end subroutine append_ga_line_vec


  subroutine expand_ga_line(me, pos, length)
    type(grw_linearray_type) :: me !< array to resize
    integer, intent(in), optional :: pos !< optional predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    type(tem_line_type), allocatable :: swpval(:)
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

  end subroutine expand_ga_line


end module tem_line_module

!> \page line Line
!! Lines are defined in the configuration file through canonical
!! geometry kind with an origin and vector defining the length
!! and direction of the line. \n
!! Valid definition:
!! \li Single line
!! \verbatim
!! geometry = {
!!   kind = 'canoND',
!!   object = {
!!     origin = { 0.0,0.0,0.0 },
!!     vec = { 2.0,0.0,0.0 }
!!   }
!! }
!! \endverbatim
!! \li Multiple line
!! \verbatim
!! geometry = {
!!   kind = 'canoND',
!!   object = {
!!     {
!!     origin = { 0.0,0.0,0.0 },
!!     vec = { 2.0,0.0,0.0 }
!!     },
!!     {
!!     origin = { 1.0,0.0,0.0 },
!!     vec = { 0.0,2.0,0.0 }
!!     },
!!   }
!! }
!! \endverbatim
!! \n\n
!! Seeder file to generate the mesh with line is generated using above canonical
!! geometry kind and the code is given below:
!! \include testsuite/line/seeder.lua
!! \n\n
!! The mesh generated with line inside mesh:
!! \image html line.png
!! \n
!! \image html line_withedges.png
!! \n\n
!! Example lua file is available at \link testsuite/line/seeder.lua
!! \example testsuite/line/seeder.lua
