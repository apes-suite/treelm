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
!> summary: This module provides the description of geometrical object 'box' and
!! routines related to them
module tem_box_module
  use env_module,       only: rk, minLength, zeroLength
  use tem_cube_module,  only: tem_cube_type
  use tem_plane_module, only: tem_plane_type, tem_planeCubeOverlap, &
    &                         tem_createPlane

  implicit none

  private

  public :: grw_boxArray_type
  public :: tem_box_type
  public :: init, append, truncate, destroy, empty, placeAt
  public :: tem_boxCubeOverlap
  public :: tem_createBox

  !> This type contains origin and vec of box in each direction
  type tem_box_type
    !> center of the box
    real(kind=rk) :: center(3)
    !> vector which defines length of the box in x,y,z direction
    !! 1st dimension contains x,y,z coord and
    !! 2nd dimension defines three direction of vector.
    !! This vector is need to do transformation of box
    real(kind=rk) :: halfvec(3,3)
    !> halfwidth of the box from center in each direction.
    !! It is computed from halfvec during initialization
    !! since it is needed for boxboxoverlap
    real(kind=rk) :: halfwidth(3)
    !> unit normal vectors which defines the box orientation.
    !! It is computed from halfvec during initialization
    real(kind=rk) :: normal(3,3)
    !> To choose what to do with intersection of box
    !! if only_surface = true than the only the surface of the object
    !! is intersected
    !! if only_surface = false then the whole object is intersected
    !! default is set to false
    logical :: only_surface
    !> For only_surface, box is converted into 6 planes and each
    !! plane is converted further into triangles
    type(tem_plane_type) :: plane(6)
  end type tem_box_type

  !> growing array type for type(tem_box_type)
  type grw_boxarray_type
    integer :: nvals = 0
    integer :: containersize = 0
    type(tem_box_type), allocatable :: val(:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga_box
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_ga_box
  end interface

  !> empty the entries  without changing arrays
  interface empty
    module procedure empty_ga_box
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga_box
  end interface

  !> insert an element at a given position
  interface placeat
    module procedure placeat_ga_box
    module procedure placeat_ga_box_vec
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_ga_box
    module procedure append_ga_box_vec
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga_box
  end interface


contains

  subroutine init_ga_box(me, length)
    type(grw_boxarray_type), intent(out) :: me !< dynamic array to init
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

  end subroutine init_ga_box

  subroutine destroy_ga_box(me)
    type(grw_boxarray_type), intent(inout) :: me !< dynamic array to destroy

    me%containersize = 0
    me%nvals = 0
    if( allocated( me%val ) ) deallocate(me%val)

  end subroutine destroy_ga_box


  subroutine truncate_ga_box(me)
    !------------------------------------------------------------------------
    type(grw_boxarray_type) :: me !< array to truncate
    !------------------------------------------------------------------------
    type(tem_box_type), allocatable :: tarray(:)
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

  end subroutine truncate_ga_box


  subroutine empty_ga_box(me)
    !------------------------------------------------------------------------
    type(grw_boxarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------

    me%nvals = 0

  end subroutine empty_ga_box

  !> adds the value to a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! element at the requested position will be replaced.
  subroutine placeat_ga_box(me, val, pos, length)
    type(grw_boxarray_type) :: me !< array to place the value into
    type(tem_box_type), intent(in) :: val !< value to place at the given position
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

  end subroutine placeat_ga_box


  !> adds the values starting from a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! elements starting from the requested position will be replaced up to
  !! the element at position `pos + size(val) - 1`.
  subroutine placeat_ga_box_vec(me, val, pos, length)
    type(grw_boxarray_type) :: me !< array to append the value to
    type(tem_box_type), intent(in) :: val(:) !< values to append
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

  end subroutine placeat_ga_box_vec


  subroutine append_ga_box(me, val, length)
    type(grw_boxarray_type) :: me !< array to append the value to
    type(tem_box_type), intent(in) :: val !< value to append
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

  end subroutine append_ga_box

  subroutine append_ga_box_vec(me, val, length)
    type(grw_boxarray_type) :: me !< array to append the value to
    type(tem_box_type), intent(in) :: val(:) !< values to append
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

  end subroutine append_ga_box_vec


  subroutine expand_ga_box(me, pos, length)
    type(grw_boxarray_type) :: me !< array to resize
    integer, intent(in), optional :: pos !< optional predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    type(tem_box_type), allocatable :: swpval(:)
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

  end subroutine expand_ga_box


  ! ***************************************************************************!
  !> This routine creates box from canoND definition i.en origin and
  !! three vectors. If only_surface is defined then box is converted further
  !! to plane and then to triangles.
  !! \verbatim
  !! vecB______vecC_
  !! /\       -    |
  !! |      -      |
  !! |    -        |
  !! |  -          |
  !! |-            |
  !! |------------>|
  !! origin      vecA
  !! \endverbatim
  subroutine tem_createBox(me, origin, vecA, vecB, vecC, only_surface)
    !--------------------------------------------------------------------------!
    type(tem_box_type), intent(out) :: me !< plane to return
    real(kind=rk), intent(in) :: origin(3) !< plane origin
    real(kind=rk), intent(in) :: vecA(3) !< 1st vector
    real(kind=rk), intent(in) :: vecB(3) !< 2nd vector
    real(kind=rk), intent(in) :: vecC(3) !< 3nd vector
    logical, intent(in) :: only_surface
    !--------------------------------------------------------------------------!
    real(kind=rk) :: extent
    real(kind=rk) :: secondOrigin(3)
    ! -------------------------------------------------------------------------!
    me%only_surface = only_surface
    if (only_surface) then
      ! convert box to planes
      ! planes 1,2,3
      call tem_createPlane(me = me%plane(1), origin = origin, &
        &                  vecA = vecA, vecB = vecB           )

      call tem_createPlane(me = me%plane(2), origin = origin, &
        &                  vecA = vecA, vecB = vecC           )

      call tem_createPlane(me = me%plane(3), origin = origin, &
        &                  vecA = vecB, vecB = vecC           )

      ! planes 4,5,6
      secondOrigin = origin + vecA + vecB + vecC
      call tem_createPlane(me = me%plane(4), origin = secondOrigin, &
        &                  vecA = -vecA, vecB = -vecB               )

      call tem_createPlane(me = me%plane(5), origin = secondOrigin, &
        &                  vecA = -vecA, vecB = -vecC               )

      call tem_createPlane(me = me%plane(6), origin = secondOrigin, &
        &                  vecA = -vecB, vecB = -vecC               )

    else ! convert me to box
      !vector is normalized first and then converted to half vector
      !length of the vector
      extent = sqrt(dot_product(vecA(:),vecA(:)))
      !normal of the vector(unitNormal)
      me%normal(:,1) = vecA(:)/extent
      !halfwidth of the vector
      me%halfwidth(1) = extent*0.5_rk
      !halfvector
      me%halfvec(:,1) = 0.5_rk*vecA(:)

      !length of the vector
      extent = sqrt(dot_product(vecB(:),vecB(:)))
      !normal of the vector(unitNormal)
      me%normal(:,2) = vecB(:)/extent
      !halfwidth of the vector
      me%halfwidth(2) = extent*0.5_rk
      !halfvector
      me%halfvec(:,2) = 0.5_rk*vecB(:)

      !length of the vector
      extent = sqrt(dot_product(vecC(:),vecC(:)))
      !normal of the vector(unitNormal)
      me%normal(:,3) = vecC(:)/extent
      !halfwidth of the vector
      me%halfwidth(3) = extent*0.5_rk
      !halfvector
      me%halfvec(:,3) = 0.5_rk*vecC(:)

      me%center = origin          &
        &       + me%halfvec(:,1) &
        &       + me%halfvec(:,2) &
        &       + me%halfvec(:,3)

    end if

  end subroutine tem_createBox
  ! ***************************************************************************!

  ! ****************************************************************************!
  !> This function checks for intersection of box and cube
  !!
  !! Currently support only axis aligned box
  !!KM: @todo implemented oriented box - cube overlap
  function tem_boxCubeOverlap( box, cube ) result(overlap)
    type(tem_box_type), intent(in) :: box
    type(tem_cube_type), intent(in) :: cube
    logical :: overlap
    integer :: iPlane

    if (box%only_surface) then
      overlap = .false.
      do iPlane = 1, 6
        overlap = overlap .or. tem_planeCubeOverlap( box%plane(iPlane), cube)
      end do
    else
      overlap = boxBoxOverlap( cube%center, &
        &                      [cube%halfwidth, &
        &                       cube%halfwidth, &
        &                       cube%halfwidth], &
        &                      box%center,  &
        &                      [box%halfwidth(1), &
        &                       box%halfwidth(2), &
        &                       box%halfwidth(3)], &
        &                      box%normal)
    end if

  end function tem_boxCubeOverlap
  ! ****************************************************************************!

! ******************************************************************************!
  !> This function checks for intersection of a axis aligned box and a
  !! parallelepiped.
  function boxBoxOverlap(center_a, dim_a, center_b, dim_b, norm_b) &
    &      result(overlap)

    real(kind=rk), intent(in) :: center_a(3) !< Center of the axis aligned box
    real(kind=rk), intent(in) :: dim_a(3) !< Length of the AAB in each direction
    real(kind=rk), intent(in) :: center_b(3) !< Center of the parallelepiped

    !> Halflength of the parallelepiped in each direction.
    !!@todo HK: Dim is a bad naming here, halflength would be better.
    real(kind=rk), intent(in) :: dim_b(3)

    real(kind=rk), intent(in) :: norm_b(3,3) !< normal vector of parallelepiped
    logical :: overlap


    ! Local Variables

    real(kind=rk)             :: vec_ab(3) ! Vector connecting center_a and
                                           ! center_b
!HK!    real(kind=rk)             :: coeffs(3,3)  ! Matrix for Coefficients needed
!HK!                                              ! for the algorithm.
!HK!
!HK!    real(kind=rk)             :: norm_a(3,3) ! edge directions of the axis
    !                                          aligned box

    ! Initialize variable
    overlap = .true.

    vec_ab = center_b - center_a

!HK!    norm_a(:,1) = (/1., 0., 0./)
!HK!    norm_a(:,2) = (/0., 1., 0./)
!HK!    norm_a(:,3) = (/0., 0., 1./)

    ! coeffs Matrix is structured as follows:
    ! --> second dimension
    ! xx xy xz
    ! yx yy yz
    ! zx zy zz


    ! Case 1:
    ! Calculation of so far needed Coefficients:
!HK!    coeffs(1,1) = dot_product(norm_a(:,1),norm_b(:,1))
!HK!    coeffs(1,2) = dot_product(norm_a(:,1),norm_b(:,2))
!HK!    coeffs(1,3) = dot_product(norm_a(:,1),norm_b(:,3))

    if ( abs(vec_ab(1)) &
      &  .gt. dim_a(1) + abs(dim_b(1)*norm_b(1,1)) &
      &                + abs(dim_b(2)*norm_b(1,2)) &
      &                + abs(dim_b(3)*norm_b(1,3))) then
      overlap = .false.
      return
    end if


    ! Case 2:
    ! Calculation of so far needed Coefficients:
!HK!    coeffs(2,1) = dot_product(norm_a(:,2),norm_b(:,1))
!HK!    coeffs(2,2) = dot_product(norm_a(:,2),norm_b(:,2))
!HK!    coeffs(2,3) = dot_product(norm_a(:,2),norm_b(:,3))

    if ( abs(vec_ab(2)) &
      &  .gt. dim_a(2) + abs(dim_b(1)*norm_b(2,1)) &
      &                + abs(dim_b(2)*norm_b(2,2)) &
      &                + abs(dim_b(3)*norm_b(2,3))) then
      overlap = .false.
      return
    end if

    ! Case 3:
    ! Calculation of so far needed Coefficients:
!HK!    coeffs(3,1) = dot_product(norm_a(:,3),norm_b(:,1))
!HK!    coeffs(3,2) = dot_product(norm_a(:,3),norm_b(:,2))
!HK!    coeffs(3,3) = dot_product(norm_a(:,3),norm_b(:,3))

    if ( abs(vec_ab(3)) &
      &  .gt. dim_a(3) + abs(dim_b(1)*norm_b(3,1)) &
      &                + abs(dim_b(2)*norm_b(3,2)) &
      &                + abs(dim_b(3)*norm_b(3,3))) then
      overlap = .false.
      return
    end if


    ! Case 4:
    if ( abs(dot_product(vec_ab, norm_b(:,1))) &
      &  .gt. dim_b(1) + abs(dim_a(1)*norm_b(1,1)) &
      &                + abs(dim_a(2)*norm_b(2,1)) &
      &                + abs(dim_a(3)*norm_b(3,1))) then
      overlap = .false.
      return
    end if

    ! Case 5:
    if ( abs(dot_product(vec_ab, norm_b(:,2))) &
      &  .gt. dim_b(2) + abs(dim_a(1)*norm_b(1,2)) &
      &                + abs(dim_a(2)*norm_b(2,2)) &
      &                + abs(dim_a(3)*norm_b(3,2))) then
      overlap = .false.
      return
    end if

    ! Case 6:
    if ( abs(dot_product(vec_ab, norm_b(:,3))) &
      &  .gt. dim_b(3) + abs(dim_a(1)*norm_b(1,3)) &
      &                + abs(dim_a(2)*norm_b(2,3)) &
      &                + abs(dim_a(3)*norm_b(3,3))) then
      overlap = .false.
      return
    end if


    ! Case 7:
    if ( abs(vec_ab(3))*norm_b(2,1) - vec_ab(2)*norm_b(3,1) &
      &  .gt. abs(dim_a(2)*norm_b(3,1)) + abs(dim_a(3)*norm_b(2,1)) &
      &       + abs(dim_b(2)*norm_b(1,3)) + abs(dim_b(3)*norm_b(1,2)) ) then
      overlap = .false.
      return
    end if

    ! Case 8:
    if ( abs(vec_ab(3))*norm_b(2,2) - vec_ab(2)*norm_b(3,2) &
      &  .gt. abs(dim_a(2)*norm_b(3,2)) + abs(dim_a(3)*norm_b(2,2)) &
      &       + abs(dim_b(1)*norm_b(1,3)) + abs(dim_b(3)*norm_b(1,1)) ) then
      overlap = .false.
      return
    end if

    ! Case 9:
    if ( abs(vec_ab(3))*norm_b(2,3) - vec_ab(2)*norm_b(3,3) &
      &  .gt. abs(dim_a(2)*norm_b(3,3)) + abs(dim_a(3)*norm_b(2,3)) &
      &       + abs(dim_b(1)*norm_b(1,2)) + abs(dim_b(2)*norm_b(1,1)) ) then
      overlap = .false.
      return
    end if

    ! Case 10:
    if ( abs(vec_ab(1))*norm_b(3,1) - vec_ab(3)*norm_b(1,1) &
      &  .gt. abs(dim_a(1)*norm_b(3,1)) + abs(dim_a(3)*norm_b(1,1)) &
      &       + abs(dim_b(2)*norm_b(2,3)) + abs(dim_b(3)*norm_b(2,2)) ) then
      overlap = .false.
      return
    end if

    ! Case 11:
    if ( abs(vec_ab(1))*norm_b(3,2) - vec_ab(3)*norm_b(1,2) &
      &  .gt. abs(dim_a(1)*norm_b(3,2)) + abs(dim_a(3)*norm_b(1,2)) &
      &       + abs(dim_b(1)*norm_b(2,3)) + abs(dim_b(3)*norm_b(2,1)) ) then
      overlap = .false.
      return
    end if

    ! Case 12:
    if ( abs(vec_ab(1))*norm_b(3,3) - vec_ab(3)*norm_b(1,3) &
      &  .gt. abs(dim_a(1)*norm_b(3,3)) + abs(dim_a(3)*norm_b(1,3)) &
      &       + abs(dim_b(1)*norm_b(2,2)) + abs(dim_b(2)*norm_b(2,1)) ) then
      overlap = .false.
      return
    end if

    ! Case 13:
    if ( abs(vec_ab(2))*norm_b(1,1) - vec_ab(1)*norm_b(2,1) &
      &  .gt. abs(dim_a(1)*norm_b(2,1)) + abs(dim_a(2)*norm_b(1,1)) &
      &       + abs(dim_b(2)*norm_b(3,3)) + abs(dim_b(3)*norm_b(3,2)) ) then
      overlap = .false.
      return
    end if

    ! Case 14:
    if ( abs(vec_ab(2))*norm_b(1,2) - vec_ab(1)*norm_b(2,2) &
      &  .gt. abs(dim_a(1)*norm_b(2,2)) + abs(dim_a(2)*norm_b(1,2)) &
      &       + abs(dim_b(1)*norm_b(3,3)) + abs(dim_b(3)*norm_b(3,1)) ) then
      overlap = .false.
      return
    end if

    ! Case 15:
    if ( abs(vec_ab(2))*norm_b(1,3) - vec_ab(1)*norm_b(2,3) &
      &  .gt. abs(dim_a(1)*norm_b(2,3)) + abs(dim_a(2)*norm_b(1,3)) &
      &       + abs(dim_b(1)*norm_b(3,2)) + abs(dim_b(2)*norm_b(3,1)) ) then
      overlap = .false.
      return
    end if

  end function boxBoxOverlap

end module tem_box_module

!> \page box Box
!! Boxes are defined by an origin and three vectors.
!! Box is considered to be solid as default i.e. all the cubes inside the
!! box are marked as intersected cubes.
!! It is possible to created hollow boxes by setting only_surface = true,
!! it will mark only the cubes intersect with sphere surface as intersected
!! cubes. Hollow box is created by converting box definion to triangles
!! internally and do triangle cube intersection.
!!
!! Seeder supports non-axis aligned oriented boxes.
!!
!! Valid definition:
!! \li Single box
!! \verbatim
!! geometry={
!!   kind='box',
!!     object={
!!       origin={0.0,0.0,0.0},
!!       vec = {
!!         { 0.0,0.5,0.0},
!!         { 0.2,0.0,0.0},
!!         { 0.0,0.0,0.8}
!!       }
!!       only_surface = true, -- If not defined default is set to false
!!     }
!! }
!! \endverbatim
!!
!! \li Rotated box
!! \verbatim
!! geometry={
!!   kind='sphere',
!!     object={
!!       {
!!       origin={0.0,0.0,0.0},
!!        vec = { { 0.0,
!!                  0.4*math.cos(45*math.pi/180),
!!                  -0.4*math.sin(45*math.pi/180) },
!!                { 0.8, 0.0, 0.0 },
!!                { 0.0,
!!                  0.4*math.sin(45*math.pi/180),
!!                  0.4*math.cos(45*math.pi/180) },
!!        }, -- box rotated along x-axis in anti-clockwise direction by
!!           -- 45 degrees
!!       }
!!     }
!! }
!! \endverbatim
!! \n\n
!! Seeder file to generate the mesh with box
!! include testsuite/box/seeder.lua
!! \n\n
!! Mesh with hollow box (Hollow => only_surface = true)
!! \image html box.png
!! \n\n
!! \image html box_withedges.png
!! \n\n
!! Cutview of Mesh with hollow box
!! \image html box_hollow.png
!! \n\n
!! Mesh generated with solid box (Solid => only_surface = false)
!! \nCutview of Mesh with solid box
!! \image html box_solid.png
!! \n\n
!! Example lua file is available at \link testsuite/box/seeder.lua
!! \exmaple testsuite/box/seeder.lua
