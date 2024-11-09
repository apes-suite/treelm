! Copyright (c) 2011-2014 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2014-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
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
!> Module to provide simple growing data structures.
!! The dynamic arrays provided by this module are
!! capable of handling lists of values, which might
!! need to grow over time.
!! Removal of entries is not possible directly.
!! The complete module might be put into a CoCo Text
!! template, to create new modules of this object
!! for different types. For now, two different
!! templates are used for the declaration part and
!! the implementation part.
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
!!
module tem_grow_array_module

  ! include treelm modules
  use env_module, only: long_k, rk, minLength, zeroLength, labelLen

  implicit none

  type intArray2d_type
    integer, allocatable :: val(:,:)
  end type

! -----------------------------------------------------------------
! Growing array (GA)
! tname ... indicates type of dynamic array (long, int, real, ...)

  !> growing array type for logical
  type grw_logicalarray_type
    integer :: nvals = 0
    integer :: containersize = 0
    logical, allocatable :: val(:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga_logical
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_ga_logical
  end interface

  !> empty the entries  without changing arrays
  interface empty
    module procedure empty_ga_logical
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga_logical
  end interface

  !> insert an element at a given position
  interface placeat
    module procedure placeat_ga_logical
    module procedure placeat_ga_logical_vec
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_ga_logical
    module procedure append_ga_logical_vec
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga_logical
  end interface

  !> growing array type for integer(kind=long_k)
  type grw_longarray_type
    integer :: nvals = 0
    integer :: containersize = 0
    integer(kind=long_k), allocatable :: val(:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga_long
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_ga_long
  end interface

  !> empty the entries  without changing arrays
  interface empty
    module procedure empty_ga_long
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga_long
  end interface

  !> insert an element at a given position
  interface placeat
    module procedure placeat_ga_long
    module procedure placeat_ga_long_vec
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_ga_long
    module procedure append_ga_long_vec
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga_long
  end interface

  !> growing array type for integer
  type grw_intarray_type
    integer :: nvals = 0
    integer :: containersize = 0
    integer, allocatable :: val(:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga_int
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_ga_int
  end interface

  !> empty the entries  without changing arrays
  interface empty
    module procedure empty_ga_int
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga_int
  end interface

  !> insert an element at a given position
  interface placeat
    module procedure placeat_ga_int
    module procedure placeat_ga_int_vec
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_ga_int
    module procedure append_ga_int_vec
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga_int
  end interface

  !> growing array type for real(kind=rk)
  type grw_realarray_type
    integer :: nvals = 0
    integer :: containersize = 0
    real(kind=rk), allocatable :: val(:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga_real
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_ga_real
  end interface

  !> empty the entries  without changing arrays
  interface empty
    module procedure empty_ga_real
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga_real
  end interface

  !> insert an element at a given position
  interface placeat
    module procedure placeat_ga_real
    module procedure placeat_ga_real_vec
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_ga_real
    module procedure append_ga_real_vec
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga_real
  end interface

  !> growing array type for type(intarray2d_type)
  type grw_dtint2darray_type
    integer :: nvals = 0
    integer :: containersize = 0
    type(intarray2d_type), allocatable :: val(:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga_dtint2d
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_ga_dtint2d
  end interface

  !> empty the entries  without changing arrays
  interface empty
    module procedure empty_ga_dtint2d
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga_dtint2d
  end interface

  !> insert an element at a given position
  interface placeat
    module procedure placeat_ga_dtint2d
    module procedure placeat_ga_dtint2d_vec
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_ga_dtint2d
    module procedure append_ga_dtint2d_vec
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga_dtint2d
  end interface

  !> growing array type for character(len=labellen)
  type grw_labelarray_type
    integer :: nvals = 0
    integer :: containersize = 0
    character(len=labellen), allocatable :: val(:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga_label
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_ga_label
  end interface

  !> empty the entries  without changing arrays
  interface empty
    module procedure empty_ga_label
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga_label
  end interface

  !> insert an element at a given position
  interface placeat
    module procedure placeat_ga_label
    module procedure placeat_ga_label_vec
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_ga_label
    module procedure append_ga_label_vec
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga_label
  end interface

  !> growing array type for character
  type grw_chararray_type
    integer :: nvals = 0
    integer :: containersize = 0
    character, allocatable :: val(:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga_char
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_ga_char
  end interface

  !> empty the entries  without changing arrays
  interface empty
    module procedure empty_ga_char
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga_char
  end interface

  !> insert an element at a given position
  interface placeat
    module procedure placeat_ga_char
    module procedure placeat_ga_char_vec
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_ga_char
    module procedure append_ga_char_vec
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga_char
  end interface


! Growing array
! -----------------------------------------------------------------


! -----------------------------------------------------------------
! Growing array (GA2d)
! tname ... indicates type of dynamic array (long, int, real, ...)
! 2d array, but only the second dimension can grow.
!

  type grw_char2darray_type
    integer :: nvals = 0
    integer :: containersize = 0
    integer :: containerwidth = 0
    character, allocatable :: val(:,:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga2d_char
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga2d_char
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_singlega2d_char
    module procedure append_arrayga2d_char
  end interface

  interface placeat
    module procedure placeat_singlega2d_char
    module procedure placeat_arrayga2d_char
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga2d_char
  end interface

  type grw_logical2darray_type
    integer :: nvals = 0
    integer :: containersize = 0
    integer :: containerwidth = 0
    logical, allocatable :: val(:,:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga2d_logical
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga2d_logical
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_singlega2d_logical
    module procedure append_arrayga2d_logical
  end interface

  interface placeat
    module procedure placeat_singlega2d_logical
    module procedure placeat_arrayga2d_logical
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga2d_logical
  end interface

  type grw_long2darray_type
    integer :: nvals = 0
    integer :: containersize = 0
    integer :: containerwidth = 0
    integer(kind=long_k), allocatable :: val(:,:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga2d_long
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga2d_long
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_singlega2d_long
    module procedure append_arrayga2d_long
  end interface

  interface placeat
    module procedure placeat_singlega2d_long
    module procedure placeat_arrayga2d_long
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga2d_long
  end interface

  type grw_int2darray_type
    integer :: nvals = 0
    integer :: containersize = 0
    integer :: containerwidth = 0
    integer, allocatable :: val(:,:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga2d_int
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga2d_int
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_singlega2d_int
    module procedure append_arrayga2d_int
  end interface

  interface placeat
    module procedure placeat_singlega2d_int
    module procedure placeat_arrayga2d_int
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga2d_int
  end interface

  type grw_real2darray_type
    integer :: nvals = 0
    integer :: containersize = 0
    integer :: containerwidth = 0
    real(kind=rk), allocatable :: val(:,:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga2d_real
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga2d_real
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_singlega2d_real
    module procedure append_arrayga2d_real
  end interface

  interface placeat
    module procedure placeat_singlega2d_real
    module procedure placeat_arrayga2d_real
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga2d_real
  end interface


! Growing array 2D
! -----------------------------------------------------------------




contains

! -----------------------------------------------------------------
! Growing array (GA)
! tname ... indicates type of dynamic array (long, int, real, ...)

  subroutine init_ga_logical(me, length)
    type(grw_logicalarray_type), intent(out) :: me !< dynamic array to init
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

  end subroutine init_ga_logical

  subroutine destroy_ga_logical(me)
    type(grw_logicalarray_type), intent(inout) :: me !< dynamic array to destroy

    me%containersize = 0
    me%nvals = 0
    if( allocated( me%val ) ) deallocate(me%val)

  end subroutine destroy_ga_logical


  subroutine truncate_ga_logical(me)
    !------------------------------------------------------------------------
    type(grw_logicalarray_type) :: me !< array to truncate
    !------------------------------------------------------------------------
    logical, allocatable :: tarray(:)
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

  end subroutine truncate_ga_logical


  subroutine empty_ga_logical(me)
    !------------------------------------------------------------------------
    type(grw_logicalarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------

    me%nvals = 0

  end subroutine empty_ga_logical

  !> adds the value to a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! element at the requested position will be replaced.
  subroutine placeat_ga_logical(me, val, pos, length)
    type(grw_logicalarray_type) :: me !< array to place the value into
    logical, intent(in) :: val !< value to place at the given position
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

  end subroutine placeat_ga_logical


  !> adds the values starting from a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! elements starting from the requested position will be replaced up to
  !! the element at position `pos + size(val) - 1`.
  subroutine placeat_ga_logical_vec(me, val, pos, length)
    type(grw_logicalarray_type) :: me !< array to append the value to
    logical, intent(in) :: val(:) !< values to append
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

  end subroutine placeat_ga_logical_vec


  subroutine append_ga_logical(me, val, length)
    type(grw_logicalarray_type) :: me !< array to append the value to
    logical, intent(in) :: val !< value to append
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

  end subroutine append_ga_logical

  subroutine append_ga_logical_vec(me, val, length)
    type(grw_logicalarray_type) :: me !< array to append the value to
    logical, intent(in) :: val(:) !< values to append
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

  end subroutine append_ga_logical_vec


  subroutine expand_ga_logical(me, pos, length)
    type(grw_logicalarray_type) :: me !< array to resize
    integer, intent(in), optional :: pos !< optional predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    logical, allocatable :: swpval(:)
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

  end subroutine expand_ga_logical

  subroutine init_ga_long(me, length)
    type(grw_longarray_type), intent(out) :: me !< dynamic array to init
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

  end subroutine init_ga_long

  subroutine destroy_ga_long(me)
    type(grw_longarray_type), intent(inout) :: me !< dynamic array to destroy

    me%containersize = 0
    me%nvals = 0
    if( allocated( me%val ) ) deallocate(me%val)

  end subroutine destroy_ga_long


  subroutine truncate_ga_long(me)
    !------------------------------------------------------------------------
    type(grw_longarray_type) :: me !< array to truncate
    !------------------------------------------------------------------------
    integer(kind=long_k), allocatable :: tarray(:)
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

  end subroutine truncate_ga_long


  subroutine empty_ga_long(me)
    !------------------------------------------------------------------------
    type(grw_longarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------

    me%nvals = 0

  end subroutine empty_ga_long

  !> adds the value to a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! element at the requested position will be replaced.
  subroutine placeat_ga_long(me, val, pos, length)
    type(grw_longarray_type) :: me !< array to place the value into
    integer(kind=long_k), intent(in) :: val !< value to place at the given position
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

  end subroutine placeat_ga_long


  !> adds the values starting from a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! elements starting from the requested position will be replaced up to
  !! the element at position `pos + size(val) - 1`.
  subroutine placeat_ga_long_vec(me, val, pos, length)
    type(grw_longarray_type) :: me !< array to append the value to
    integer(kind=long_k), intent(in) :: val(:) !< values to append
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

  end subroutine placeat_ga_long_vec


  subroutine append_ga_long(me, val, length)
    type(grw_longarray_type) :: me !< array to append the value to
    integer(kind=long_k), intent(in) :: val !< value to append
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

  end subroutine append_ga_long

  subroutine append_ga_long_vec(me, val, length)
    type(grw_longarray_type) :: me !< array to append the value to
    integer(kind=long_k), intent(in) :: val(:) !< values to append
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

  end subroutine append_ga_long_vec


  subroutine expand_ga_long(me, pos, length)
    type(grw_longarray_type) :: me !< array to resize
    integer, intent(in), optional :: pos !< optional predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    integer(kind=long_k), allocatable :: swpval(:)
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

  end subroutine expand_ga_long

  subroutine init_ga_int(me, length)
    type(grw_intarray_type), intent(out) :: me !< dynamic array to init
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

  end subroutine init_ga_int

  subroutine destroy_ga_int(me)
    type(grw_intarray_type), intent(inout) :: me !< dynamic array to destroy

    me%containersize = 0
    me%nvals = 0
    if( allocated( me%val ) ) deallocate(me%val)

  end subroutine destroy_ga_int


  subroutine truncate_ga_int(me)
    !------------------------------------------------------------------------
    type(grw_intarray_type) :: me !< array to truncate
    !------------------------------------------------------------------------
    integer, allocatable :: tarray(:)
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

  end subroutine truncate_ga_int


  subroutine empty_ga_int(me)
    !------------------------------------------------------------------------
    type(grw_intarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------

    me%nvals = 0

  end subroutine empty_ga_int

  !> adds the value to a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! element at the requested position will be replaced.
  subroutine placeat_ga_int(me, val, pos, length)
    type(grw_intarray_type) :: me !< array to place the value into
    integer, intent(in) :: val !< value to place at the given position
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

  end subroutine placeat_ga_int


  !> adds the values starting from a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! elements starting from the requested position will be replaced up to
  !! the element at position `pos + size(val) - 1`.
  subroutine placeat_ga_int_vec(me, val, pos, length)
    type(grw_intarray_type) :: me !< array to append the value to
    integer, intent(in) :: val(:) !< values to append
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

  end subroutine placeat_ga_int_vec


  subroutine append_ga_int(me, val, length)
    type(grw_intarray_type) :: me !< array to append the value to
    integer, intent(in) :: val !< value to append
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

  end subroutine append_ga_int

  subroutine append_ga_int_vec(me, val, length)
    type(grw_intarray_type) :: me !< array to append the value to
    integer, intent(in) :: val(:) !< values to append
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

  end subroutine append_ga_int_vec


  subroutine expand_ga_int(me, pos, length)
    type(grw_intarray_type) :: me !< array to resize
    integer, intent(in), optional :: pos !< optional predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    integer, allocatable :: swpval(:)
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

  end subroutine expand_ga_int

  subroutine init_ga_real(me, length)
    type(grw_realarray_type), intent(out) :: me !< dynamic array to init
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

  end subroutine init_ga_real

  subroutine destroy_ga_real(me)
    type(grw_realarray_type), intent(inout) :: me !< dynamic array to destroy

    me%containersize = 0
    me%nvals = 0
    if( allocated( me%val ) ) deallocate(me%val)

  end subroutine destroy_ga_real


  subroutine truncate_ga_real(me)
    !------------------------------------------------------------------------
    type(grw_realarray_type) :: me !< array to truncate
    !------------------------------------------------------------------------
    real(kind=rk), allocatable :: tarray(:)
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

  end subroutine truncate_ga_real


  subroutine empty_ga_real(me)
    !------------------------------------------------------------------------
    type(grw_realarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------

    me%nvals = 0

  end subroutine empty_ga_real

  !> adds the value to a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! element at the requested position will be replaced.
  subroutine placeat_ga_real(me, val, pos, length)
    type(grw_realarray_type) :: me !< array to place the value into
    real(kind=rk), intent(in) :: val !< value to place at the given position
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

  end subroutine placeat_ga_real


  !> adds the values starting from a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! elements starting from the requested position will be replaced up to
  !! the element at position `pos + size(val) - 1`.
  subroutine placeat_ga_real_vec(me, val, pos, length)
    type(grw_realarray_type) :: me !< array to append the value to
    real(kind=rk), intent(in) :: val(:) !< values to append
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

  end subroutine placeat_ga_real_vec


  subroutine append_ga_real(me, val, length)
    type(grw_realarray_type) :: me !< array to append the value to
    real(kind=rk), intent(in) :: val !< value to append
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

  end subroutine append_ga_real

  subroutine append_ga_real_vec(me, val, length)
    type(grw_realarray_type) :: me !< array to append the value to
    real(kind=rk), intent(in) :: val(:) !< values to append
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

  end subroutine append_ga_real_vec


  subroutine expand_ga_real(me, pos, length)
    type(grw_realarray_type) :: me !< array to resize
    integer, intent(in), optional :: pos !< optional predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    real(kind=rk), allocatable :: swpval(:)
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

  end subroutine expand_ga_real

  subroutine init_ga_dtint2d(me, length)
    type(grw_dtint2darray_type), intent(out) :: me !< dynamic array to init
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

  end subroutine init_ga_dtint2d

  subroutine destroy_ga_dtint2d(me)
    type(grw_dtint2darray_type), intent(inout) :: me !< dynamic array to destroy

    me%containersize = 0
    me%nvals = 0
    if( allocated( me%val ) ) deallocate(me%val)

  end subroutine destroy_ga_dtint2d


  subroutine truncate_ga_dtint2d(me)
    !------------------------------------------------------------------------
    type(grw_dtint2darray_type) :: me !< array to truncate
    !------------------------------------------------------------------------
    type(intarray2d_type), allocatable :: tarray(:)
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

  end subroutine truncate_ga_dtint2d


  subroutine empty_ga_dtint2d(me)
    !------------------------------------------------------------------------
    type(grw_dtint2darray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------

    me%nvals = 0

  end subroutine empty_ga_dtint2d

  !> adds the value to a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! element at the requested position will be replaced.
  subroutine placeat_ga_dtint2d(me, val, pos, length)
    type(grw_dtint2darray_type) :: me !< array to place the value into
    type(intarray2d_type), intent(in) :: val !< value to place at the given position
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

  end subroutine placeat_ga_dtint2d


  !> adds the values starting from a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! elements starting from the requested position will be replaced up to
  !! the element at position `pos + size(val) - 1`.
  subroutine placeat_ga_dtint2d_vec(me, val, pos, length)
    type(grw_dtint2darray_type) :: me !< array to append the value to
    type(intarray2d_type), intent(in) :: val(:) !< values to append
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

  end subroutine placeat_ga_dtint2d_vec


  subroutine append_ga_dtint2d(me, val, length)
    type(grw_dtint2darray_type) :: me !< array to append the value to
    type(intarray2d_type), intent(in) :: val !< value to append
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

  end subroutine append_ga_dtint2d

  subroutine append_ga_dtint2d_vec(me, val, length)
    type(grw_dtint2darray_type) :: me !< array to append the value to
    type(intarray2d_type), intent(in) :: val(:) !< values to append
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

  end subroutine append_ga_dtint2d_vec


  subroutine expand_ga_dtint2d(me, pos, length)
    type(grw_dtint2darray_type) :: me !< array to resize
    integer, intent(in), optional :: pos !< optional predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    type(intarray2d_type), allocatable :: swpval(:)
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

  end subroutine expand_ga_dtint2d

  subroutine init_ga_label(me, length)
    type(grw_labelarray_type), intent(out) :: me !< dynamic array to init
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

  end subroutine init_ga_label

  subroutine destroy_ga_label(me)
    type(grw_labelarray_type), intent(inout) :: me !< dynamic array to destroy

    me%containersize = 0
    me%nvals = 0
    if( allocated( me%val ) ) deallocate(me%val)

  end subroutine destroy_ga_label


  subroutine truncate_ga_label(me)
    !------------------------------------------------------------------------
    type(grw_labelarray_type) :: me !< array to truncate
    !------------------------------------------------------------------------
    character(len=labellen), allocatable :: tarray(:)
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

  end subroutine truncate_ga_label


  subroutine empty_ga_label(me)
    !------------------------------------------------------------------------
    type(grw_labelarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------

    me%nvals = 0

  end subroutine empty_ga_label

  !> adds the value to a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! element at the requested position will be replaced.
  subroutine placeat_ga_label(me, val, pos, length)
    type(grw_labelarray_type) :: me !< array to place the value into
    character(len=labellen), intent(in) :: val !< value to place at the given position
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

  end subroutine placeat_ga_label


  !> adds the values starting from a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! elements starting from the requested position will be replaced up to
  !! the element at position `pos + size(val) - 1`.
  subroutine placeat_ga_label_vec(me, val, pos, length)
    type(grw_labelarray_type) :: me !< array to append the value to
    character(len=labellen), intent(in) :: val(:) !< values to append
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

  end subroutine placeat_ga_label_vec


  subroutine append_ga_label(me, val, length)
    type(grw_labelarray_type) :: me !< array to append the value to
    character(len=labellen), intent(in) :: val !< value to append
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

  end subroutine append_ga_label

  subroutine append_ga_label_vec(me, val, length)
    type(grw_labelarray_type) :: me !< array to append the value to
    character(len=labellen), intent(in) :: val(:) !< values to append
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

  end subroutine append_ga_label_vec


  subroutine expand_ga_label(me, pos, length)
    type(grw_labelarray_type) :: me !< array to resize
    integer, intent(in), optional :: pos !< optional predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    character(len=labellen), allocatable :: swpval(:)
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

  end subroutine expand_ga_label

  subroutine init_ga_char(me, length)
    type(grw_chararray_type), intent(out) :: me !< dynamic array to init
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

  end subroutine init_ga_char

  subroutine destroy_ga_char(me)
    type(grw_chararray_type), intent(inout) :: me !< dynamic array to destroy

    me%containersize = 0
    me%nvals = 0
    if( allocated( me%val ) ) deallocate(me%val)

  end subroutine destroy_ga_char


  subroutine truncate_ga_char(me)
    !------------------------------------------------------------------------
    type(grw_chararray_type) :: me !< array to truncate
    !------------------------------------------------------------------------
    character, allocatable :: tarray(:)
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

  end subroutine truncate_ga_char


  subroutine empty_ga_char(me)
    !------------------------------------------------------------------------
    type(grw_chararray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------

    me%nvals = 0

  end subroutine empty_ga_char

  !> adds the value to a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! element at the requested position will be replaced.
  subroutine placeat_ga_char(me, val, pos, length)
    type(grw_chararray_type) :: me !< array to place the value into
    character, intent(in) :: val !< value to place at the given position
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

  end subroutine placeat_ga_char


  !> adds the values starting from a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! elements starting from the requested position will be replaced up to
  !! the element at position `pos + size(val) - 1`.
  subroutine placeat_ga_char_vec(me, val, pos, length)
    type(grw_chararray_type) :: me !< array to append the value to
    character, intent(in) :: val(:) !< values to append
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

  end subroutine placeat_ga_char_vec


  subroutine append_ga_char(me, val, length)
    type(grw_chararray_type) :: me !< array to append the value to
    character, intent(in) :: val !< value to append
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

  end subroutine append_ga_char

  subroutine append_ga_char_vec(me, val, length)
    type(grw_chararray_type) :: me !< array to append the value to
    character, intent(in) :: val(:) !< values to append
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

  end subroutine append_ga_char_vec


  subroutine expand_ga_char(me, pos, length)
    type(grw_chararray_type) :: me !< array to resize
    integer, intent(in), optional :: pos !< optional predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    character, allocatable :: swpval(:)
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

  end subroutine expand_ga_char

! Growing array
! -----------------------------------------------------------------



! ****************************************************************************** !
! -----------------------------------------------------------------
! 2d Array, which can grow in second dimension only (GA2d)
! tname ... indicates type of dynamic array (long, int, real, ...)

  ! ****************************************************************************
  subroutine init_ga2d_char(me, width, length)
    ! --------------------------------------------------------------------------
    !> dynamic array to init
    type(grw_char2darray_type), intent(out) :: me
    !> width of the container
    integer, intent(in)           :: width
    !> initial length of the container
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------

    me%containerwidth = width
    if (present(length)) me%containersize = length

    allocate(me%val(me%containerwidth,me%containersize))
    ! reset all values (only positive values are valid)
    me%val(:,:) = ''
    me%nvals = 0

  end subroutine init_ga2d_char
  ! ****************************************************************************


  ! ****************************************************************************
  !> destroy the 2d growing array
  !!
  subroutine destroy_ga2d_char(me)
    ! --------------------------------------------------------------------------
    !> dynamic array to init
    type(grw_char2darray_type), intent(inout) :: me
    ! --------------------------------------------------------------------------

    deallocate(me%val)
    me%nvals  = 0
    me%containersize = 0

  end subroutine destroy_ga2d_char
  ! ****************************************************************************


  ! ****************************************************************************
  !> append a single value to the growing 2d array.
  !!
  subroutine append_singlega2d_char(me, val, pos1, pos2, length)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_char2darray_type) :: me
    !> value to append
    character, intent(in) :: val
    !> position in first dimension (cannot grow)
    integer, intent(in) :: pos1
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> the position in second dimension the element were added to
    integer, intent(out), optional :: pos2
    ! --------------------------------------------------------------------------
    integer :: newpos
    ! --------------------------------------------------------------------------

    newpos = me%nvals + 1
    if (newpos > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand( me = me, length = length )
    end if

    me%val(pos1, newpos) = val

    if( present( pos2 ) ) pos2 = newpos

  end subroutine append_singlega2d_char
  ! ****************************************************************************


  ! ****************************************************************************
  !> append an array of values to the growing 2d array.
  !!
  subroutine append_arrayga2d_char(me, val, length, pos)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_char2darray_type) :: me
    !> array of values to append
    character, intent(in) :: val(:)
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> the position in second dimension the elements were added to
    integer, intent(out), optional :: pos
    ! --------------------------------------------------------------------------
    integer :: newpos
    ! --------------------------------------------------------------------------

    if (me%nvals == huge(me%nvals)) then
       write(*,*) "reached end of integer range for growing array!"
       write(*,*) "aborting!!"
       stop
    end if

    newpos = me%nvals + 1
    if (newpos > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand( me = me, length = length )
    end if

    me%nvals = max(newpos, me%nvals)
    me%val( : , newpos ) = val(:)

    if( present( pos ) ) pos = newpos

  end subroutine append_arrayga2d_char
  ! ****************************************************************************


  ! ****************************************************************************
  !> append a single value to the growing 2d array.
  !!
  subroutine placeat_singlega2d_char(me, val, pos1, pos2, length)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_char2darray_type) :: me
    !> value to append
    character, intent(in) :: val
    !> position in first dimension (cannot grow)
    integer, intent(in) :: pos1
    !> position in second dimension (can grow)
    integer, intent(in) :: pos2
    !> optional length to expand the array
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------
    integer :: expandlength
    ! --------------------------------------------------------------------------

    if (pos2 > me%containersize) then
      expandlength = pos2 - me%containersize
      if( present( length ) ) expandlength = max(expandlength, length)
      ! expand the array, if its boundary is reached
      call expand( me = me, length = expandlength )
    end if
    me%nvals = max(pos2, me%nvals)
    me%val(pos1,pos2) = val

  end subroutine placeat_singlega2d_char
  ! ****************************************************************************


  ! ****************************************************************************
  !> append an array of values to the growing 2d array.
  !!
  subroutine placeat_arrayga2d_char(me, val, pos2, length)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_char2darray_type) :: me
    !> array of values to append
    character, intent(in) :: val(:)
    !> position in second dimension (can grow)
    integer, intent(in) :: pos2
    !> optional length to expand the array
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------
    integer :: expandlength
    ! --------------------------------------------------------------------------

    if (me%nvals == huge(me%nvals)) then
       write(*,*) "reached end of integer range for growing array!"
       write(*,*) "aborting!!"
       stop
    end if
    if (pos2 > me%containersize) then
      expandlength = pos2 - me%containersize
      if( present( length ) ) expandlength = max(expandlength, length)
      ! expand the array, if its boundary is reached
      call expand( me = me, length = expandlength )
    end if

    me%nvals = max(pos2, me%nvals)
    me%val( : , pos2 ) = val(:)

  end subroutine placeat_arrayga2d_char
  ! ****************************************************************************


  ! ****************************************************************************
  !> expand the growing 2d array
  subroutine expand_ga2d_char( me, length )
    ! --------------------------------------------------------------------------
    !> array to resize
    type(grw_char2darray_type) :: me
    !> optional length to expand the array
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------
    character, allocatable :: swpval(:,:)
    integer :: explen
    ! --------------------------------------------------------------------------

    explen = 0
    ! increase the container by the requested length of double it
    if( present(length) ) then
      explen = length
    else
      ! set the global minimum length, if doubling would be smaller than that
      explen = max(me%containersize, minlength)
    end if


    ! if the current size plus explen exceeds the max container size,
    ! reduce the size to the max container size.
    if( (huge(me%containersize) - explen) <= me%containersize) then
      ! set max container size
      me%containersize = huge(me%containersize)
    else
      ! set the calculated size
      me%containersize = me%containersize + explen
    end if


    if ( me%nvals > 0 ) then
      allocate(swpval(me%containerwidth,me%containersize))
      !> first reset all entries
      swpval(:,:) = ''
      swpval(:,:me%nvals) = me%val(:,:me%nvals)
      call move_alloc( swpval, me%val )
    else ! me%nvals == 0
      if ( allocated(me%val) ) deallocate( me%val )
      allocate( me%val(me%containerwidth,me%containersize) )
    end if

  end subroutine expand_ga2d_char
  ! ****************************************************************************

! growing array 2d
! -----------------------------------------------------------------


  ! ****************************************************************************
  subroutine init_ga2d_logical(me, width, length)
    ! --------------------------------------------------------------------------
    !> dynamic array to init
    type(grw_logical2darray_type), intent(out) :: me
    !> width of the container
    integer, intent(in)           :: width
    !> initial length of the container
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------

    me%containerwidth = width
    if (present(length)) me%containersize = length

    allocate(me%val(me%containerwidth,me%containersize))
    ! reset all values (only positive values are valid)
    me%val(:,:) = .false.
    me%nvals = 0

  end subroutine init_ga2d_logical
  ! ****************************************************************************


  ! ****************************************************************************
  !> destroy the 2d growing array
  !!
  subroutine destroy_ga2d_logical(me)
    ! --------------------------------------------------------------------------
    !> dynamic array to init
    type(grw_logical2darray_type), intent(inout) :: me
    ! --------------------------------------------------------------------------

    deallocate(me%val)
    me%nvals  = 0
    me%containersize = 0

  end subroutine destroy_ga2d_logical
  ! ****************************************************************************


  ! ****************************************************************************
  !> append a single value to the growing 2d array.
  !!
  subroutine append_singlega2d_logical(me, val, pos1, pos2, length)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_logical2darray_type) :: me
    !> value to append
    logical, intent(in) :: val
    !> position in first dimension (cannot grow)
    integer, intent(in) :: pos1
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> the position in second dimension the element were added to
    integer, intent(out), optional :: pos2
    ! --------------------------------------------------------------------------
    integer :: newpos
    ! --------------------------------------------------------------------------

    newpos = me%nvals + 1
    if (newpos > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand( me = me, length = length )
    end if

    me%val(pos1, newpos) = val

    if( present( pos2 ) ) pos2 = newpos

  end subroutine append_singlega2d_logical
  ! ****************************************************************************


  ! ****************************************************************************
  !> append an array of values to the growing 2d array.
  !!
  subroutine append_arrayga2d_logical(me, val, length, pos)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_logical2darray_type) :: me
    !> array of values to append
    logical, intent(in) :: val(:)
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> the position in second dimension the elements were added to
    integer, intent(out), optional :: pos
    ! --------------------------------------------------------------------------
    integer :: newpos
    ! --------------------------------------------------------------------------

    if (me%nvals == huge(me%nvals)) then
       write(*,*) "reached end of integer range for growing array!"
       write(*,*) "aborting!!"
       stop
    end if

    newpos = me%nvals + 1
    if (newpos > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand( me = me, length = length )
    end if

    me%nvals = max(newpos, me%nvals)
    me%val( : , newpos ) = val(:)

    if( present( pos ) ) pos = newpos

  end subroutine append_arrayga2d_logical
  ! ****************************************************************************


  ! ****************************************************************************
  !> append a single value to the growing 2d array.
  !!
  subroutine placeat_singlega2d_logical(me, val, pos1, pos2, length)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_logical2darray_type) :: me
    !> value to append
    logical, intent(in) :: val
    !> position in first dimension (cannot grow)
    integer, intent(in) :: pos1
    !> position in second dimension (can grow)
    integer, intent(in) :: pos2
    !> optional length to expand the array
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------
    integer :: expandlength
    ! --------------------------------------------------------------------------

    if (pos2 > me%containersize) then
      expandlength = pos2 - me%containersize
      if( present( length ) ) expandlength = max(expandlength, length)
      ! expand the array, if its boundary is reached
      call expand( me = me, length = expandlength )
    end if
    me%nvals = max(pos2, me%nvals)
    me%val(pos1,pos2) = val

  end subroutine placeat_singlega2d_logical
  ! ****************************************************************************


  ! ****************************************************************************
  !> append an array of values to the growing 2d array.
  !!
  subroutine placeat_arrayga2d_logical(me, val, pos2, length)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_logical2darray_type) :: me
    !> array of values to append
    logical, intent(in) :: val(:)
    !> position in second dimension (can grow)
    integer, intent(in) :: pos2
    !> optional length to expand the array
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------
    integer :: expandlength
    ! --------------------------------------------------------------------------

    if (me%nvals == huge(me%nvals)) then
       write(*,*) "reached end of integer range for growing array!"
       write(*,*) "aborting!!"
       stop
    end if
    if (pos2 > me%containersize) then
      expandlength = pos2 - me%containersize
      if( present( length ) ) expandlength = max(expandlength, length)
      ! expand the array, if its boundary is reached
      call expand( me = me, length = expandlength )
    end if

    me%nvals = max(pos2, me%nvals)
    me%val( : , pos2 ) = val(:)

  end subroutine placeat_arrayga2d_logical
  ! ****************************************************************************


  ! ****************************************************************************
  !> expand the growing 2d array
  subroutine expand_ga2d_logical( me, length )
    ! --------------------------------------------------------------------------
    !> array to resize
    type(grw_logical2darray_type) :: me
    !> optional length to expand the array
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------
    logical, allocatable :: swpval(:,:)
    integer :: explen
    ! --------------------------------------------------------------------------

    explen = 0
    ! increase the container by the requested length of double it
    if( present(length) ) then
      explen = length
    else
      ! set the global minimum length, if doubling would be smaller than that
      explen = max(me%containersize, minlength)
    end if


    ! if the current size plus explen exceeds the max container size,
    ! reduce the size to the max container size.
    if( (huge(me%containersize) - explen) <= me%containersize) then
      ! set max container size
      me%containersize = huge(me%containersize)
    else
      ! set the calculated size
      me%containersize = me%containersize + explen
    end if


    if ( me%nvals > 0 ) then
      allocate(swpval(me%containerwidth,me%containersize))
      !> first reset all entries
      swpval(:,:) = .false.
      swpval(:,:me%nvals) = me%val(:,:me%nvals)
      call move_alloc( swpval, me%val )
    else ! me%nvals == 0
      if ( allocated(me%val) ) deallocate( me%val )
      allocate( me%val(me%containerwidth,me%containersize) )
    end if

  end subroutine expand_ga2d_logical
  ! ****************************************************************************

! growing array 2d
! -----------------------------------------------------------------


  ! ****************************************************************************
  subroutine init_ga2d_long(me, width, length)
    ! --------------------------------------------------------------------------
    !> dynamic array to init
    type(grw_long2darray_type), intent(out) :: me
    !> width of the container
    integer, intent(in)           :: width
    !> initial length of the container
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------

    me%containerwidth = width
    if (present(length)) me%containersize = length

    allocate(me%val(me%containerwidth,me%containersize))
    ! reset all values (only positive values are valid)
    me%val(:,:) = -1_long_k
    me%nvals = 0

  end subroutine init_ga2d_long
  ! ****************************************************************************


  ! ****************************************************************************
  !> destroy the 2d growing array
  !!
  subroutine destroy_ga2d_long(me)
    ! --------------------------------------------------------------------------
    !> dynamic array to init
    type(grw_long2darray_type), intent(inout) :: me
    ! --------------------------------------------------------------------------

    deallocate(me%val)
    me%nvals  = 0
    me%containersize = 0

  end subroutine destroy_ga2d_long
  ! ****************************************************************************


  ! ****************************************************************************
  !> append a single value to the growing 2d array.
  !!
  subroutine append_singlega2d_long(me, val, pos1, pos2, length)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_long2darray_type) :: me
    !> value to append
    integer(kind=long_k), intent(in) :: val
    !> position in first dimension (cannot grow)
    integer, intent(in) :: pos1
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> the position in second dimension the element were added to
    integer, intent(out), optional :: pos2
    ! --------------------------------------------------------------------------
    integer :: newpos
    ! --------------------------------------------------------------------------

    newpos = me%nvals + 1
    if (newpos > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand( me = me, length = length )
    end if

    me%val(pos1, newpos) = val

    if( present( pos2 ) ) pos2 = newpos

  end subroutine append_singlega2d_long
  ! ****************************************************************************


  ! ****************************************************************************
  !> append an array of values to the growing 2d array.
  !!
  subroutine append_arrayga2d_long(me, val, length, pos)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_long2darray_type) :: me
    !> array of values to append
    integer(kind=long_k), intent(in) :: val(:)
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> the position in second dimension the elements were added to
    integer, intent(out), optional :: pos
    ! --------------------------------------------------------------------------
    integer :: newpos
    ! --------------------------------------------------------------------------

    if (me%nvals == huge(me%nvals)) then
       write(*,*) "reached end of integer range for growing array!"
       write(*,*) "aborting!!"
       stop
    end if

    newpos = me%nvals + 1
    if (newpos > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand( me = me, length = length )
    end if

    me%nvals = max(newpos, me%nvals)
    me%val( : , newpos ) = val(:)

    if( present( pos ) ) pos = newpos

  end subroutine append_arrayga2d_long
  ! ****************************************************************************


  ! ****************************************************************************
  !> append a single value to the growing 2d array.
  !!
  subroutine placeat_singlega2d_long(me, val, pos1, pos2, length)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_long2darray_type) :: me
    !> value to append
    integer(kind=long_k), intent(in) :: val
    !> position in first dimension (cannot grow)
    integer, intent(in) :: pos1
    !> position in second dimension (can grow)
    integer, intent(in) :: pos2
    !> optional length to expand the array
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------
    integer :: expandlength
    ! --------------------------------------------------------------------------

    if (pos2 > me%containersize) then
      expandlength = pos2 - me%containersize
      if( present( length ) ) expandlength = max(expandlength, length)
      ! expand the array, if its boundary is reached
      call expand( me = me, length = expandlength )
    end if
    me%nvals = max(pos2, me%nvals)
    me%val(pos1,pos2) = val

  end subroutine placeat_singlega2d_long
  ! ****************************************************************************


  ! ****************************************************************************
  !> append an array of values to the growing 2d array.
  !!
  subroutine placeat_arrayga2d_long(me, val, pos2, length)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_long2darray_type) :: me
    !> array of values to append
    integer(kind=long_k), intent(in) :: val(:)
    !> position in second dimension (can grow)
    integer, intent(in) :: pos2
    !> optional length to expand the array
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------
    integer :: expandlength
    ! --------------------------------------------------------------------------

    if (me%nvals == huge(me%nvals)) then
       write(*,*) "reached end of integer range for growing array!"
       write(*,*) "aborting!!"
       stop
    end if
    if (pos2 > me%containersize) then
      expandlength = pos2 - me%containersize
      if( present( length ) ) expandlength = max(expandlength, length)
      ! expand the array, if its boundary is reached
      call expand( me = me, length = expandlength )
    end if

    me%nvals = max(pos2, me%nvals)
    me%val( : , pos2 ) = val(:)

  end subroutine placeat_arrayga2d_long
  ! ****************************************************************************


  ! ****************************************************************************
  !> expand the growing 2d array
  subroutine expand_ga2d_long( me, length )
    ! --------------------------------------------------------------------------
    !> array to resize
    type(grw_long2darray_type) :: me
    !> optional length to expand the array
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------
    integer(kind=long_k), allocatable :: swpval(:,:)
    integer :: explen
    ! --------------------------------------------------------------------------

    explen = 0
    ! increase the container by the requested length of double it
    if( present(length) ) then
      explen = length
    else
      ! set the global minimum length, if doubling would be smaller than that
      explen = max(me%containersize, minlength)
    end if


    ! if the current size plus explen exceeds the max container size,
    ! reduce the size to the max container size.
    if( (huge(me%containersize) - explen) <= me%containersize) then
      ! set max container size
      me%containersize = huge(me%containersize)
    else
      ! set the calculated size
      me%containersize = me%containersize + explen
    end if


    if ( me%nvals > 0 ) then
      allocate(swpval(me%containerwidth,me%containersize))
      !> first reset all entries
      swpval(:,:) = -1_long_k
      swpval(:,:me%nvals) = me%val(:,:me%nvals)
      call move_alloc( swpval, me%val )
    else ! me%nvals == 0
      if ( allocated(me%val) ) deallocate( me%val )
      allocate( me%val(me%containerwidth,me%containersize) )
    end if

  end subroutine expand_ga2d_long
  ! ****************************************************************************

! growing array 2d
! -----------------------------------------------------------------


  ! ****************************************************************************
  subroutine init_ga2d_int(me, width, length)
    ! --------------------------------------------------------------------------
    !> dynamic array to init
    type(grw_int2darray_type), intent(out) :: me
    !> width of the container
    integer, intent(in)           :: width
    !> initial length of the container
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------

    me%containerwidth = width
    if (present(length)) me%containersize = length

    allocate(me%val(me%containerwidth,me%containersize))
    ! reset all values (only positive values are valid)
    me%val(:,:) = -1
    me%nvals = 0

  end subroutine init_ga2d_int
  ! ****************************************************************************


  ! ****************************************************************************
  !> destroy the 2d growing array
  !!
  subroutine destroy_ga2d_int(me)
    ! --------------------------------------------------------------------------
    !> dynamic array to init
    type(grw_int2darray_type), intent(inout) :: me
    ! --------------------------------------------------------------------------

    deallocate(me%val)
    me%nvals  = 0
    me%containersize = 0

  end subroutine destroy_ga2d_int
  ! ****************************************************************************


  ! ****************************************************************************
  !> append a single value to the growing 2d array.
  !!
  subroutine append_singlega2d_int(me, val, pos1, pos2, length)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_int2darray_type) :: me
    !> value to append
    integer, intent(in) :: val
    !> position in first dimension (cannot grow)
    integer, intent(in) :: pos1
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> the position in second dimension the element were added to
    integer, intent(out), optional :: pos2
    ! --------------------------------------------------------------------------
    integer :: newpos
    ! --------------------------------------------------------------------------

    newpos = me%nvals + 1
    if (newpos > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand( me = me, length = length )
    end if

    me%val(pos1, newpos) = val

    if( present( pos2 ) ) pos2 = newpos

  end subroutine append_singlega2d_int
  ! ****************************************************************************


  ! ****************************************************************************
  !> append an array of values to the growing 2d array.
  !!
  subroutine append_arrayga2d_int(me, val, length, pos)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_int2darray_type) :: me
    !> array of values to append
    integer, intent(in) :: val(:)
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> the position in second dimension the elements were added to
    integer, intent(out), optional :: pos
    ! --------------------------------------------------------------------------
    integer :: newpos
    ! --------------------------------------------------------------------------

    if (me%nvals == huge(me%nvals)) then
       write(*,*) "reached end of integer range for growing array!"
       write(*,*) "aborting!!"
       stop
    end if

    newpos = me%nvals + 1
    if (newpos > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand( me = me, length = length )
    end if

    me%nvals = max(newpos, me%nvals)
    me%val( : , newpos ) = val(:)

    if( present( pos ) ) pos = newpos

  end subroutine append_arrayga2d_int
  ! ****************************************************************************


  ! ****************************************************************************
  !> append a single value to the growing 2d array.
  !!
  subroutine placeat_singlega2d_int(me, val, pos1, pos2, length)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_int2darray_type) :: me
    !> value to append
    integer, intent(in) :: val
    !> position in first dimension (cannot grow)
    integer, intent(in) :: pos1
    !> position in second dimension (can grow)
    integer, intent(in) :: pos2
    !> optional length to expand the array
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------
    integer :: expandlength
    ! --------------------------------------------------------------------------

    if (pos2 > me%containersize) then
      expandlength = pos2 - me%containersize
      if( present( length ) ) expandlength = max(expandlength, length)
      ! expand the array, if its boundary is reached
      call expand( me = me, length = expandlength )
    end if
    me%nvals = max(pos2, me%nvals)
    me%val(pos1,pos2) = val

  end subroutine placeat_singlega2d_int
  ! ****************************************************************************


  ! ****************************************************************************
  !> append an array of values to the growing 2d array.
  !!
  subroutine placeat_arrayga2d_int(me, val, pos2, length)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_int2darray_type) :: me
    !> array of values to append
    integer, intent(in) :: val(:)
    !> position in second dimension (can grow)
    integer, intent(in) :: pos2
    !> optional length to expand the array
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------
    integer :: expandlength
    ! --------------------------------------------------------------------------

    if (me%nvals == huge(me%nvals)) then
       write(*,*) "reached end of integer range for growing array!"
       write(*,*) "aborting!!"
       stop
    end if
    if (pos2 > me%containersize) then
      expandlength = pos2 - me%containersize
      if( present( length ) ) expandlength = max(expandlength, length)
      ! expand the array, if its boundary is reached
      call expand( me = me, length = expandlength )
    end if

    me%nvals = max(pos2, me%nvals)
    me%val( : , pos2 ) = val(:)

  end subroutine placeat_arrayga2d_int
  ! ****************************************************************************


  ! ****************************************************************************
  !> expand the growing 2d array
  subroutine expand_ga2d_int( me, length )
    ! --------------------------------------------------------------------------
    !> array to resize
    type(grw_int2darray_type) :: me
    !> optional length to expand the array
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------
    integer, allocatable :: swpval(:,:)
    integer :: explen
    ! --------------------------------------------------------------------------

    explen = 0
    ! increase the container by the requested length of double it
    if( present(length) ) then
      explen = length
    else
      ! set the global minimum length, if doubling would be smaller than that
      explen = max(me%containersize, minlength)
    end if


    ! if the current size plus explen exceeds the max container size,
    ! reduce the size to the max container size.
    if( (huge(me%containersize) - explen) <= me%containersize) then
      ! set max container size
      me%containersize = huge(me%containersize)
    else
      ! set the calculated size
      me%containersize = me%containersize + explen
    end if


    if ( me%nvals > 0 ) then
      allocate(swpval(me%containerwidth,me%containersize))
      !> first reset all entries
      swpval(:,:) = -1
      swpval(:,:me%nvals) = me%val(:,:me%nvals)
      call move_alloc( swpval, me%val )
    else ! me%nvals == 0
      if ( allocated(me%val) ) deallocate( me%val )
      allocate( me%val(me%containerwidth,me%containersize) )
    end if

  end subroutine expand_ga2d_int
  ! ****************************************************************************

! growing array 2d
! -----------------------------------------------------------------


  ! ****************************************************************************
  subroutine init_ga2d_real(me, width, length)
    ! --------------------------------------------------------------------------
    !> dynamic array to init
    type(grw_real2darray_type), intent(out) :: me
    !> width of the container
    integer, intent(in)           :: width
    !> initial length of the container
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------

    me%containerwidth = width
    if (present(length)) me%containersize = length

    allocate(me%val(me%containerwidth,me%containersize))
    ! reset all values (only positive values are valid)
    me%val(:,:) = -1.0_rk
    me%nvals = 0

  end subroutine init_ga2d_real
  ! ****************************************************************************


  ! ****************************************************************************
  !> destroy the 2d growing array
  !!
  subroutine destroy_ga2d_real(me)
    ! --------------------------------------------------------------------------
    !> dynamic array to init
    type(grw_real2darray_type), intent(inout) :: me
    ! --------------------------------------------------------------------------

    deallocate(me%val)
    me%nvals  = 0
    me%containersize = 0

  end subroutine destroy_ga2d_real
  ! ****************************************************************************


  ! ****************************************************************************
  !> append a single value to the growing 2d array.
  !!
  subroutine append_singlega2d_real(me, val, pos1, pos2, length)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_real2darray_type) :: me
    !> value to append
    real(kind=rk), intent(in) :: val
    !> position in first dimension (cannot grow)
    integer, intent(in) :: pos1
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> the position in second dimension the element were added to
    integer, intent(out), optional :: pos2
    ! --------------------------------------------------------------------------
    integer :: newpos
    ! --------------------------------------------------------------------------

    newpos = me%nvals + 1
    if (newpos > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand( me = me, length = length )
    end if

    me%val(pos1, newpos) = val

    if( present( pos2 ) ) pos2 = newpos

  end subroutine append_singlega2d_real
  ! ****************************************************************************


  ! ****************************************************************************
  !> append an array of values to the growing 2d array.
  !!
  subroutine append_arrayga2d_real(me, val, length, pos)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_real2darray_type) :: me
    !> array of values to append
    real(kind=rk), intent(in) :: val(:)
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> the position in second dimension the elements were added to
    integer, intent(out), optional :: pos
    ! --------------------------------------------------------------------------
    integer :: newpos
    ! --------------------------------------------------------------------------

    if (me%nvals == huge(me%nvals)) then
       write(*,*) "reached end of integer range for growing array!"
       write(*,*) "aborting!!"
       stop
    end if

    newpos = me%nvals + 1
    if (newpos > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand( me = me, length = length )
    end if

    me%nvals = max(newpos, me%nvals)
    me%val( : , newpos ) = val(:)

    if( present( pos ) ) pos = newpos

  end subroutine append_arrayga2d_real
  ! ****************************************************************************


  ! ****************************************************************************
  !> append a single value to the growing 2d array.
  !!
  subroutine placeat_singlega2d_real(me, val, pos1, pos2, length)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_real2darray_type) :: me
    !> value to append
    real(kind=rk), intent(in) :: val
    !> position in first dimension (cannot grow)
    integer, intent(in) :: pos1
    !> position in second dimension (can grow)
    integer, intent(in) :: pos2
    !> optional length to expand the array
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------
    integer :: expandlength
    ! --------------------------------------------------------------------------

    if (pos2 > me%containersize) then
      expandlength = pos2 - me%containersize
      if( present( length ) ) expandlength = max(expandlength, length)
      ! expand the array, if its boundary is reached
      call expand( me = me, length = expandlength )
    end if
    me%nvals = max(pos2, me%nvals)
    me%val(pos1,pos2) = val

  end subroutine placeat_singlega2d_real
  ! ****************************************************************************


  ! ****************************************************************************
  !> append an array of values to the growing 2d array.
  !!
  subroutine placeat_arrayga2d_real(me, val, pos2, length)
    ! --------------------------------------------------------------------------
    !> array to append the value to
    type(grw_real2darray_type) :: me
    !> array of values to append
    real(kind=rk), intent(in) :: val(:)
    !> position in second dimension (can grow)
    integer, intent(in) :: pos2
    !> optional length to expand the array
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------
    integer :: expandlength
    ! --------------------------------------------------------------------------

    if (me%nvals == huge(me%nvals)) then
       write(*,*) "reached end of integer range for growing array!"
       write(*,*) "aborting!!"
       stop
    end if
    if (pos2 > me%containersize) then
      expandlength = pos2 - me%containersize
      if( present( length ) ) expandlength = max(expandlength, length)
      ! expand the array, if its boundary is reached
      call expand( me = me, length = expandlength )
    end if

    me%nvals = max(pos2, me%nvals)
    me%val( : , pos2 ) = val(:)

  end subroutine placeat_arrayga2d_real
  ! ****************************************************************************


  ! ****************************************************************************
  !> expand the growing 2d array
  subroutine expand_ga2d_real( me, length )
    ! --------------------------------------------------------------------------
    !> array to resize
    type(grw_real2darray_type) :: me
    !> optional length to expand the array
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------
    real(kind=rk), allocatable :: swpval(:,:)
    integer :: explen
    ! --------------------------------------------------------------------------

    explen = 0
    ! increase the container by the requested length of double it
    if( present(length) ) then
      explen = length
    else
      ! set the global minimum length, if doubling would be smaller than that
      explen = max(me%containersize, minlength)
    end if


    ! if the current size plus explen exceeds the max container size,
    ! reduce the size to the max container size.
    if( (huge(me%containersize) - explen) <= me%containersize) then
      ! set max container size
      me%containersize = huge(me%containersize)
    else
      ! set the calculated size
      me%containersize = me%containersize + explen
    end if


    if ( me%nvals > 0 ) then
      allocate(swpval(me%containerwidth,me%containersize))
      !> first reset all entries
      swpval(:,:) = -1.0_rk
      swpval(:,:me%nvals) = me%val(:,:me%nvals)
      call move_alloc( swpval, me%val )
    else ! me%nvals == 0
      if ( allocated(me%val) ) deallocate( me%val )
      allocate( me%val(me%containerwidth,me%containersize) )
    end if

  end subroutine expand_ga2d_real
  ! ****************************************************************************

! growing array 2d
! -----------------------------------------------------------------



! Growing array 2D
! -----------------------------------------------------------------


end module tem_grow_array_module
! ****************************************************************************** !
