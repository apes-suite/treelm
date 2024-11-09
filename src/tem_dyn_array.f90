! Copyright (c) 2011-2013, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011-2012 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2014 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> Module to provide smart growing data structures.
!!
!! The dynamic arrays provided by this module are
!! capable of handling lists of values, which might
!! need to grow over time.
!! Removal of entries is not possible directly.
!! A ranking list for the indices is maintained to
!! fast lookups of given values by a binary search.
!! This allows for the efficient handling of lists
!! with unique entries.
!! The complete module might be put into a CoCo Text
!! template, to create new modules of this object
!! for different types. For now, two different
!! templates are used for the declaration part and
!! the implementation part.
!!
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
module tem_dyn_array_module

  ! include treelm modules
  use env_module, only: long_k, rk, minLength, labelLen, zeroLength

  implicit none

! -----------------------------------------------------------------
! In order to define this data structure easily for several data
! type, we use the Coco text copying feature here, to duplicate
! the necessary declarations.
! tname: indicates type of dynamic array (long, int, real, ...)
! tstring: is the actual string describing the type specification

! Here the actual declarations are put in by CoCo:
! \brief smart dynamic array (da) for integer(kind=long_k)
!
! this datatype implements a dynamic array, which is capable of
! growing and adding of unique elements. it is available for
! various types, here we deal with $tstring$.
! sorted array contains the pointers to val array, instead of
! the actual values in val array. for example:
! val:     8, 6, 7, 5
! sorted:  4, 2, 3, 1
  !> dynamic array (da) type for integer(kind=long_k)
  type dyn_longarray_type
    integer :: nvals = 0
    integer :: containersize = 0
    integer(kind=long_k), allocatable :: val(:)
    integer, allocatable :: sorted(:) !< pointers, not values
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_da_long
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_da_long
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_da_long
    module procedure append_da_veclong
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_da_long
  end interface

  !> empty the array, reset nvals to be 0
  interface empty
    module procedure empty_da_long
  end interface


  !> fix the dynamic array, meaning:
  !! store the array in the sorted order and cut
  !! off the trailing empty entries
  interface sorttruncate
    module procedure sorttruncate_da_long
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_da_long
  end interface

  !> return the position of a given value
  !! in the array val, which is what you usually want to know.
  !! it is the index of a given value
  interface positionofval
    module procedure posofval_long
  end interface

  !> return the position of a given value
  !! in the list 'sorted'. this is mainly for internal usage.
  !! the sorted list is only a pointer list to the actual values
  !! thus, in order to get the index of a given value, you
  !! need to look up the entry in the sorted list.
  !! this is done by the positionofval routine
  interface sortedposofval
    module procedure sortposofval_long
  end interface
! \brief smart dynamic array (da) for integer
!
! this datatype implements a dynamic array, which is capable of
! growing and adding of unique elements. it is available for
! various types, here we deal with $tstring$.
! sorted array contains the pointers to val array, instead of
! the actual values in val array. for example:
! val:     8, 6, 7, 5
! sorted:  4, 2, 3, 1
  !> dynamic array (da) type for integer
  type dyn_intarray_type
    integer :: nvals = 0
    integer :: containersize = 0
    integer, allocatable :: val(:)
    integer, allocatable :: sorted(:) !< pointers, not values
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_da_int
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_da_int
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_da_int
    module procedure append_da_vecint
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_da_int
  end interface

  !> empty the array, reset nvals to be 0
  interface empty
    module procedure empty_da_int
  end interface


  !> fix the dynamic array, meaning:
  !! store the array in the sorted order and cut
  !! off the trailing empty entries
  interface sorttruncate
    module procedure sorttruncate_da_int
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_da_int
  end interface

  !> return the position of a given value
  !! in the array val, which is what you usually want to know.
  !! it is the index of a given value
  interface positionofval
    module procedure posofval_int
  end interface

  !> return the position of a given value
  !! in the list 'sorted'. this is mainly for internal usage.
  !! the sorted list is only a pointer list to the actual values
  !! thus, in order to get the index of a given value, you
  !! need to look up the entry in the sorted list.
  !! this is done by the positionofval routine
  interface sortedposofval
    module procedure sortposofval_int
  end interface
! \brief smart dynamic array (da) for real(kind=rk)
!
! this datatype implements a dynamic array, which is capable of
! growing and adding of unique elements. it is available for
! various types, here we deal with $tstring$.
! sorted array contains the pointers to val array, instead of
! the actual values in val array. for example:
! val:     8, 6, 7, 5
! sorted:  4, 2, 3, 1
  !> dynamic array (da) type for real(kind=rk)
  type dyn_realarray_type
    integer :: nvals = 0
    integer :: containersize = 0
    real(kind=rk), allocatable :: val(:)
    integer, allocatable :: sorted(:) !< pointers, not values
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_da_real
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_da_real
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_da_real
    module procedure append_da_vecreal
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_da_real
  end interface

  !> empty the array, reset nvals to be 0
  interface empty
    module procedure empty_da_real
  end interface


  !> fix the dynamic array, meaning:
  !! store the array in the sorted order and cut
  !! off the trailing empty entries
  interface sorttruncate
    module procedure sorttruncate_da_real
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_da_real
  end interface

  !> return the position of a given value
  !! in the array val, which is what you usually want to know.
  !! it is the index of a given value
  interface positionofval
    module procedure posofval_real
  end interface

  !> return the position of a given value
  !! in the list 'sorted'. this is mainly for internal usage.
  !! the sorted list is only a pointer list to the actual values
  !! thus, in order to get the index of a given value, you
  !! need to look up the entry in the sorted list.
  !! this is done by the positionofval routine
  interface sortedposofval
    module procedure sortposofval_real
  end interface
! \brief smart dynamic array (da) for character(len=labellen)
!
! this datatype implements a dynamic array, which is capable of
! growing and adding of unique elements. it is available for
! various types, here we deal with $tstring$.
! sorted array contains the pointers to val array, instead of
! the actual values in val array. for example:
! val:     8, 6, 7, 5
! sorted:  4, 2, 3, 1
  !> dynamic array (da) type for character(len=labellen)
  type dyn_labelarray_type
    integer :: nvals = 0
    integer :: containersize = 0
    character(len=labellen), allocatable :: val(:)
    integer, allocatable :: sorted(:) !< pointers, not values
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_da_label
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_da_label
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_da_label
    module procedure append_da_veclabel
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_da_label
  end interface

  !> empty the array, reset nvals to be 0
  interface empty
    module procedure empty_da_label
  end interface


  !> fix the dynamic array, meaning:
  !! store the array in the sorted order and cut
  !! off the trailing empty entries
  interface sorttruncate
    module procedure sorttruncate_da_label
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_da_label
  end interface

  !> return the position of a given value
  !! in the array val, which is what you usually want to know.
  !! it is the index of a given value
  interface positionofval
    module procedure posofval_label
  end interface

  !> return the position of a given value
  !! in the list 'sorted'. this is mainly for internal usage.
  !! the sorted list is only a pointer list to the actual values
  !! thus, in order to get the index of a given value, you
  !! need to look up the entry in the sorted list.
  !! this is done by the positionofval routine
  interface sortedposofval
    module procedure sortposofval_label
  end interface
! -----------------------------------------------------------------

contains

! -----------------------------------------------------------------
! Also for the implementation, we use the copy feature, to provide
! the necessary duplications to deal with the various types.
! tname ... indicates type of dynamic array (long, int, real, ...)

! ****************************************************************************** !
! ***************************************************************************** !
  !> initialization of a dynamic array
  !!
  !! before a dynamic array can be used, it has to be initialized
  !! with this routine. the initial length provided here, can
  !! avoid reallocations and memory copying, if approximated
  !! correctly enough. if none is specified, the provided container
  !! initially will be of size 0.
  subroutine init_da_long(me, length)
    !-----------------------------------------------------------------
    type(dyn_longarray_type), intent(out) :: me !< dynamic array to init
    integer, intent(in), optional :: length !< initial length of the container
    !-----------------------------------------------------------------

    if (present(length)) then
      me%containersize = length
    else
      me%containersize = zerolength
    end if

    ! deallocate ...
    if( allocated( me%val ) ) deallocate(me%val)
    if( allocated( me%sorted ) ) deallocate(me%sorted)
    ! ... and reallocate
    allocate(me%val(me%containersize))
    allocate(me%sorted(me%containersize))
    me%nvals = 0

  end subroutine init_da_long

  !> destruction of a dynamic array
  !!
  !! this subroutine takes care of a proper destruction of a
  !! dynamic array, it frees the allocated memory and resets
  !! the internal counts to 0.
  subroutine destroy_da_long(me)
    type(dyn_longarray_type), intent(inout) :: me !< dynamic array to init

    me%containersize = 0
    me%nvals         = 0

    if( allocated( me%val ) ) deallocate(me%val)
    if( allocated( me%sorted ) ) deallocate(me%sorted)
  end subroutine destroy_da_long
! ***************************************************************************** !


! ***************************************************************************** !
  !> appending a value to the dynamic array
  !!
  !! with this subroutine, a given value can be added to the
  !! dynamic array. the actual position of this value in the
  !! dynamic array will be returned, so it can be found again
  !! easily later. with the wasadded flag, it is indicated,\n
  !! wasadded = true,  if this entry had to be added,\n
  !! wasadded = false, if this was already found in the array.
  subroutine append_da_long(me, val, length, pos, wasadded )
    !------------------------------------------------------------------------
    type(dyn_longarray_type) :: me   !< array to append the value to
    integer(kind=long_k), intent(in)     :: val  !< value to append
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> position in the array, if the value is found
    integer, intent(out), optional :: pos
    !> flag to indicate, if val was newly added
    logical, intent(out), optional :: wasadded
    !------------------------------------------------------------------------
    integer :: foundpos
    integer :: i
    !------------------------------------------------------------------------

    ! do a binary search on existing entries (returns closest entry next to
    ! it if not found).
    foundpos = sortedposofval(me, val, .true.)
    if( present( wasadded ) ) wasadded = .false.

    ! if it found the value, the position is smaller than nvals
    if (foundpos <= me%nvals) then

      ! the returned position might actually be the right entry already or
      ! not, check for this here.
      if ( me%val(me%sorted(foundpos)) == val ) then

        ! found the value in a list of unique values,
        ! nothing to do, just return its position.
        if( present( pos ) ) pos = me%sorted(foundpos)

      else

        ! need to append a new value!

        if (me%nvals == huge(me%nvals)) then
           write(*,*) "reached end of integer range for dynamic array!"
           write(*,*) "aborting!!"
           stop
        end if

        if( present( wasadded ) ) wasadded = .true.
        if (me%nvals == me%containersize) then

          ! container is full, need to expand it
          call expand(me = me, length = length)
        end if
        me%nvals = me%nvals + 1

        ! put the new value into the last position in the
        ! array.
        me%val(me%nvals) = val
        do while( foundpos < me%nvals )
          if(me%val(me%sorted(foundpos)) /= val) then
            exit
          end if
          ! in case of multiple entries with the same value
          ! move on to the first differing entry.
          foundpos = foundpos + 1
        end do
        ! shift the sorted list of indices, to create a
        ! whole for the value to be inserted, at position
        ! foundpos.
        do i=me%nvals-1,foundpos,-1
          me%sorted(i+1) = me%sorted(i)
        end do
        ! put the index of the new value into the
        ! sorted list at the now freed position.
        me%sorted(foundpos) = me%nvals
        if( present( pos ) ) pos = me%nvals

      end if

    else

      ! value to append is larger than all existing ones,
      ! just put it to the end of the list, this captures
      ! also the case of empty lists.
      ! in this case foundpos = me%nvals + 1 holds.
      if( present( wasadded ) ) wasadded = .true.
      if (foundpos > me%containersize) then
        ! expand the array, if its boundary is reached
        call expand(me = me, length = length)
      end if
      me%nvals = foundpos
      me%val(foundpos) = val
      me%sorted(foundpos) = foundpos
      if( present( pos ) ) pos = foundpos

    end if

  end subroutine append_da_long
! ***************************************************************************** !


! ***************************************************************************** !
  !> appending a sorted list of values to the dynamic array
  !!
  !! with this subroutine, a given list of sorted values can be added to the
  !! dynamic array. the actual positions of these values in the
  !! dynamic array will be returned, so it can be found again
  !! easily later. with the wasadded flag, it is indicated,\n
  !! wasadded = true,  if this entry had to be added,\n
  !! wasadded = false, if this was already found in the array.
  subroutine append_da_veclong(me, val, length, pos, wasadded )
    !------------------------------------------------------------------------
    type(dyn_longarray_type) :: me   !< array to append the value to
    integer(kind=long_k), intent(in)     :: val(:)  !< values to append
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> position in the array, the values are found at.
    integer, intent(out), optional :: pos(:)
    !> flag to indicate, if val was newly added
    logical, intent(out), optional :: wasadded(:)
    !------------------------------------------------------------------------
    integer(kind=long_k) :: lastval
    logical :: addedval(size(val))
    integer :: i
    integer :: veclen
    integer :: maxlen
    integer :: nappend
    integer :: rem_app
    integer :: curval, ival, iold, iadd
    integer, allocatable :: newsorted(:)
    !------------------------------------------------------------------------

    if (size(val) == 0) return

    veclen = size(val)
    maxlen = veclen + me%nvals

    allocate(newsorted(maxlen))

    addedval = .false.

    iold = 1
    iadd = 1

    nappend = 0
    curval = 0

    ! select the first entry before the loop unconditionally without checks
    ! for uniqueness (nothing to check against yet).
    if ( me%val(me%sorted(iold)) <= val(iadd) ) then
      curval = curval + 1
      newsorted(curval) = me%sorted(iold)
      lastval = me%val(me%sorted(iold))
      iold = iold + 1
    else
      curval = curval + 1
      nappend = nappend + 1
      newsorted(curval) = me%nvals + nappend
      lastval = val(iadd)
      if (present(pos))  pos(iadd) = newsorted(curval)
      addedval(iadd) = .true.
      iadd = iadd + 1
    end if

    do ival=2,maxlen

      if ( (iadd <= veclen) .and. (iold <= me%nvals) ) then

        if ( me%val(me%sorted(iold)) <= val(iadd) ) then

          ! the original list's values are appended to newsorted before
          ! the additional list is appended.
          curval = curval + 1
          newsorted(curval) = me%sorted(iold)
          lastval = me%val(me%sorted(iold))
          iold = iold + 1

        else

          ! only append the value to unique lists, if it is not yet in the list.
          ! (if it is already in the list, it has to be the previous (curval-1)
          !  entry.)
          if ( lastval < val(iadd) ) then
            nappend = nappend + 1
            curval = curval + 1
            newsorted(curval) = me%nvals + nappend
            lastval = val(iadd)
            addedval(iadd) = .true.
          end if
          if (present(pos)) pos(iadd) = newsorted(curval)
          iadd = iadd + 1

        end if

      else

        ! reached the end of one or both of the sorted lists.
        exit

      end if

    end do

    if (iold <= me%nvals) then
      ! still some values from the original list left.
      newsorted(curval+1:me%nvals+nappend) = me%sorted(iold:me%nvals)
    end if

    if (iadd <= veclen) then
      ! still some values from the list to append left.
      rem_app = iadd
      do i = rem_app,veclen
        if ( lastval < val(iadd) ) then
          nappend = nappend + 1
          curval = curval + 1
          newsorted(curval) = me%nvals + nappend
          lastval = val(iadd)
          addedval(iadd) = .true.
        end if
        if (present(pos)) pos(iadd) = newsorted(curval)
        iadd = iadd + 1
      end do
    end if

    if (me%nvals > huge(me%nvals)-nappend) then
       write(*,*) "reached end of integer range for dynamic array!"
       write(*,*) "aborting!!"
       stop
    end if

    if (me%nvals + nappend > me%containersize) then
      call expand( me        = me,      &
        &          increment = nappend, &
        &          length    = length   )
    end if
    me%sorted(:me%nvals+nappend) = newsorted(:me%nvals+nappend)
    curval = me%nvals
    do iadd=1,veclen
      if (addedval(iadd)) then
        curval = curval + 1
        me%val(curval) = val(iadd)
      end if
    end do
    me%nvals = me%nvals + nappend

    if( present( wasadded ) ) wasadded = addedval

  end subroutine append_da_veclong
! ***************************************************************************** !


! ***************************************************************************** !
  !> truncate the array after the last valid entry and hence cut off the empty
  !! trailing empty entries
  !!
  subroutine truncate_da_long(me)
    !------------------------------------------------------------------------
    type(dyn_longarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------
    integer(kind=long_k), allocatable :: swpval(:)
    integer, allocatable :: swpsort(:)
    !------------------------------------------------------------------------

    if (me%nvals < me%containersize) then
      allocate(swpval(me%nvals))
      allocate(swpsort(me%nvals))

      swpval = me%val(:me%nvals)
      swpsort = me%sorted(:me%nvals)

      call move_alloc(swpval, me%val)
      call move_alloc(swpsort, me%sorted)

      me%containersize = me%nvals
    end if

  end subroutine truncate_da_long
! ***************************************************************************** !


! ***************************************************************************** !
  !> empty all contents of the array without changing the size or status of any
  !! array
  !!
  subroutine empty_da_long(me)
    !------------------------------------------------------------------------
    type(dyn_longarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------
    ! reset the number of entries
    me%nvals = 0
  end subroutine empty_da_long
! ***************************************************************************** !


! ***************************************************************************** !
  !> fixing the dynamic array
  !!
  !! truncate the array after the last valid entry and hence cut off the empty
  !! trailing empty entries
  !! store the array in the sorted order according to the sorted( ) array
  !!
  subroutine sorttruncate_da_long(me)
    !------------------------------------------------------------------------
    type(dyn_longarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------
    type(dyn_longarray_type) :: tarray !< temporary array
    integer :: ival
    integer :: dpos
    !------------------------------------------------------------------------
    ! allocate the temporary array
    call init( me = tarray, length = me%nvals )
    ! copy the entries in a sorted fashion into the temporary array
    do ival = 1, me%nvals
      call append( me = tarray, val = me%val( me%sorted( ival )), &
           &       pos = dpos)
    enddo
    call destroy( me = me )

    me = tarray
    call destroy( me = tarray )

  end subroutine sorttruncate_da_long
! ***************************************************************************** !


! ***************************************************************************** !
  !> expanding the dynamic array
  !!
  !! this is a helping subroutine, which doubles the container
  !! of the given dynamic array. as the container might be
  !! initially 0-sized, a module variable minlength has been introduced, which
  !! is used here, to at least create a container of this size.
  subroutine expand_da_long(me, increment, length)
    !------------------------------------------------------------------------
    type(dyn_longarray_type) :: me !< array to resize
    integer, optional :: increment !< used for vector append
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !------------------------------------------------------------------------
    integer(kind=long_k), allocatable :: swpval(:)
    integer, allocatable :: swpsort(:)
    !------------------------------------------------------------------------
    integer :: addvals, explen
    !------------------------------------------------------------------------

    addvals = 1
    if (present(increment)) addvals = increment

    if (addvals > 0) then

      ! if length is present, use that, otherwise double the size
      if( present( length ) ) then
        explen = length
      else
        ! set the global minimum length, if doubling would be smaller than that
        explen = max(me%containersize, minlength)
      end if

      ! check whether all elements will fit
      if( addvals > explen ) then
        explen = addvals
      end if

      ! check whether the new size will exceed the max container size.
      if( (huge(me%containersize) - explen) <= me%containersize ) then
        ! if so, expand to the maximum size
        me%containersize = huge(me%containersize)
      else
        ! if not, expand to the calculated size
        me%containersize = me%containersize + explen
      end if

      ! only need to do something, if there are actually values to append.
      if (me%nvals > 0) then

        allocate(swpval(me%containersize))
        swpval(1:me%nvals) = me%val(1:me%nvals)
        call move_alloc( swpval, me%val )

        allocate(swpsort(me%containersize))
        swpsort(1:me%nvals) = me%sorted(1:me%nvals)
        call move_alloc( swpsort, me%sorted )

      else ! me%nvals == 0

        if( allocated(me%val) ) &
          deallocate(me%val)
        allocate(me%val(me%containersize))
        if( allocated(me%sorted) ) &
          deallocate(me%sorted)
        allocate(me%sorted(me%containersize))

      end if

    end if

  end subroutine expand_da_long
! ***************************************************************************** !


! ***************************************************************************** !
  !> return the sorted position of a value in the given dynamic array
  !!
  !! if the value was not found,
  !!  - return 0 if nextifnotfound = .false.
  !!  - return position at the end if nextifnotfound = .true.
  function sortposofval_long(me, val, nextifnotfound, lower, upper) result(pos)
    !------------------------------------------------------------------------
    type(dyn_longarray_type), intent(in) :: me !< dynamic array
    integer(kind=long_k), intent(in) :: val !< value to look for
    !> flag to indicate, if the next entry in the list should be returned,
    !! if the searched one is not found.
    logical, intent(in), optional :: nextifnotfound
    integer, intent(in), optional :: lower !< lower search limit
    integer, intent(in), optional :: upper !< upper search limit
    integer :: pos !< position of val in the sorted list, 0 if not found
    !------------------------------------------------------------------------
    logical :: retnext
    integer :: lb, ub
    integer :: mid
    integer(kind=long_k) :: lb_val, ub_val
    integer(kind=long_k) :: mid_val
    !------------------------------------------------------------------------

    retnext = .false.
    if (present(nextifnotfound)) retnext = nextifnotfound

    lb = 1
    ub = me%nvals

    if( present( lower ) ) lb = lower
    if( present( upper ) ) ub = upper

    pos = 0
    if (retnext) pos = lb

    !> binary search on sorted list
    do while(ub >= lb)
      lb_val = me%val(me%sorted(lb))

      if (val < lb_val) then
        if (retnext) pos = lb
        exit
      end if

      ub_val = me%val(me%sorted(ub))

      if (val > ub_val) then
        if (retnext) pos = ub+1
        exit
      end if

      ! safe guard against integer limit overflow
      mid = lb + (ub-lb) / 2
      mid_val = me%val(me%sorted(mid))
      if (val == mid_val) then
        pos = mid
        exit
      end if
      if (val > mid_val) then
        lb = mid + 1
      else
        ub = mid - 1
      end if
    end do
  end function sortposofval_long
! ***************************************************************************** !


! ***************************************************************************** !
  !> the actual position of a given value in the dynamic array
  !!
  !! most likely this is what you need in codes, using this
  !! data structure, it first does the binary search on the sorted
  !! values with sortposofval_long and then returns the looked
  !! up position in the original unsorted array, which corresponds
  !! to the position returned by the append routine.
  function posofval_long(me, val, nextifnotfound, lower, upper) result(pos)
    !------------------------------------------------------------------------
    type(dyn_longarray_type), intent(in) :: me !< dynamic array
    integer(kind=long_k), intent(in) :: val !< value to search for
    !> flag to indicate, if the position of the next entry in the sorted
    !! list should be returned instead, if val is not found.
    logical, intent(in), optional :: nextifnotfound
    integer, intent(in), optional :: lower !< lower search limit
    integer, intent(in), optional :: upper !< upper search limit
    integer :: pos !< position in the array of the searche value, 0 if not found
    !------------------------------------------------------------------------
    integer :: sortpos
    integer :: lb, ub
    !------------------------------------------------------------------------

    lb = 1
    ub = me%nvals

    if( present( lower ) ) lb = lower
    if( present( upper ) ) ub = upper

    sortpos = sortedposofval(me, val, nextifnotfound, lb, ub)

    ! if result (sorted pos)
    if ((sortpos <= me%nvals) .and. (sortpos > 0)) then
      pos = me%sorted(sortpos)
    else
      pos = sortpos
    end if
  end function posofval_long
! ***************************************************************************** !
! ***************************************************************************** !
  !> initialization of a dynamic array
  !!
  !! before a dynamic array can be used, it has to be initialized
  !! with this routine. the initial length provided here, can
  !! avoid reallocations and memory copying, if approximated
  !! correctly enough. if none is specified, the provided container
  !! initially will be of size 0.
  subroutine init_da_int(me, length)
    !-----------------------------------------------------------------
    type(dyn_intarray_type), intent(out) :: me !< dynamic array to init
    integer, intent(in), optional :: length !< initial length of the container
    !-----------------------------------------------------------------

    if (present(length)) then
      me%containersize = length
    else
      me%containersize = zerolength
    end if

    ! deallocate ...
    if( allocated( me%val ) ) deallocate(me%val)
    if( allocated( me%sorted ) ) deallocate(me%sorted)
    ! ... and reallocate
    allocate(me%val(me%containersize))
    allocate(me%sorted(me%containersize))
    me%nvals = 0

  end subroutine init_da_int

  !> destruction of a dynamic array
  !!
  !! this subroutine takes care of a proper destruction of a
  !! dynamic array, it frees the allocated memory and resets
  !! the internal counts to 0.
  subroutine destroy_da_int(me)
    type(dyn_intarray_type), intent(inout) :: me !< dynamic array to init

    me%containersize = 0
    me%nvals         = 0

    if( allocated( me%val ) ) deallocate(me%val)
    if( allocated( me%sorted ) ) deallocate(me%sorted)
  end subroutine destroy_da_int
! ***************************************************************************** !


! ***************************************************************************** !
  !> appending a value to the dynamic array
  !!
  !! with this subroutine, a given value can be added to the
  !! dynamic array. the actual position of this value in the
  !! dynamic array will be returned, so it can be found again
  !! easily later. with the wasadded flag, it is indicated,\n
  !! wasadded = true,  if this entry had to be added,\n
  !! wasadded = false, if this was already found in the array.
  subroutine append_da_int(me, val, length, pos, wasadded )
    !------------------------------------------------------------------------
    type(dyn_intarray_type) :: me   !< array to append the value to
    integer, intent(in)     :: val  !< value to append
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> position in the array, if the value is found
    integer, intent(out), optional :: pos
    !> flag to indicate, if val was newly added
    logical, intent(out), optional :: wasadded
    !------------------------------------------------------------------------
    integer :: foundpos
    integer :: i
    !------------------------------------------------------------------------

    ! do a binary search on existing entries (returns closest entry next to
    ! it if not found).
    foundpos = sortedposofval(me, val, .true.)
    if( present( wasadded ) ) wasadded = .false.

    ! if it found the value, the position is smaller than nvals
    if (foundpos <= me%nvals) then

      ! the returned position might actually be the right entry already or
      ! not, check for this here.
      if ( me%val(me%sorted(foundpos)) == val ) then

        ! found the value in a list of unique values,
        ! nothing to do, just return its position.
        if( present( pos ) ) pos = me%sorted(foundpos)

      else

        ! need to append a new value!

        if (me%nvals == huge(me%nvals)) then
           write(*,*) "reached end of integer range for dynamic array!"
           write(*,*) "aborting!!"
           stop
        end if

        if( present( wasadded ) ) wasadded = .true.
        if (me%nvals == me%containersize) then

          ! container is full, need to expand it
          call expand(me = me, length = length)
        end if
        me%nvals = me%nvals + 1

        ! put the new value into the last position in the
        ! array.
        me%val(me%nvals) = val
        do while( foundpos < me%nvals )
          if(me%val(me%sorted(foundpos)) /= val) then
            exit
          end if
          ! in case of multiple entries with the same value
          ! move on to the first differing entry.
          foundpos = foundpos + 1
        end do
        ! shift the sorted list of indices, to create a
        ! whole for the value to be inserted, at position
        ! foundpos.
        do i=me%nvals-1,foundpos,-1
          me%sorted(i+1) = me%sorted(i)
        end do
        ! put the index of the new value into the
        ! sorted list at the now freed position.
        me%sorted(foundpos) = me%nvals
        if( present( pos ) ) pos = me%nvals

      end if

    else

      ! value to append is larger than all existing ones,
      ! just put it to the end of the list, this captures
      ! also the case of empty lists.
      ! in this case foundpos = me%nvals + 1 holds.
      if( present( wasadded ) ) wasadded = .true.
      if (foundpos > me%containersize) then
        ! expand the array, if its boundary is reached
        call expand(me = me, length = length)
      end if
      me%nvals = foundpos
      me%val(foundpos) = val
      me%sorted(foundpos) = foundpos
      if( present( pos ) ) pos = foundpos

    end if

  end subroutine append_da_int
! ***************************************************************************** !


! ***************************************************************************** !
  !> appending a sorted list of values to the dynamic array
  !!
  !! with this subroutine, a given list of sorted values can be added to the
  !! dynamic array. the actual positions of these values in the
  !! dynamic array will be returned, so it can be found again
  !! easily later. with the wasadded flag, it is indicated,\n
  !! wasadded = true,  if this entry had to be added,\n
  !! wasadded = false, if this was already found in the array.
  subroutine append_da_vecint(me, val, length, pos, wasadded )
    !------------------------------------------------------------------------
    type(dyn_intarray_type) :: me   !< array to append the value to
    integer, intent(in)     :: val(:)  !< values to append
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> position in the array, the values are found at.
    integer, intent(out), optional :: pos(:)
    !> flag to indicate, if val was newly added
    logical, intent(out), optional :: wasadded(:)
    !------------------------------------------------------------------------
    integer :: lastval
    logical :: addedval(size(val))
    integer :: i
    integer :: veclen
    integer :: maxlen
    integer :: nappend
    integer :: rem_app
    integer :: curval, ival, iold, iadd
    integer, allocatable :: newsorted(:)
    !------------------------------------------------------------------------

    if (size(val) == 0) return

    veclen = size(val)
    maxlen = veclen + me%nvals

    allocate(newsorted(maxlen))

    addedval = .false.

    iold = 1
    iadd = 1

    nappend = 0
    curval = 0

    ! select the first entry before the loop unconditionally without checks
    ! for uniqueness (nothing to check against yet).
    if ( me%val(me%sorted(iold)) <= val(iadd) ) then
      curval = curval + 1
      newsorted(curval) = me%sorted(iold)
      lastval = me%val(me%sorted(iold))
      iold = iold + 1
    else
      curval = curval + 1
      nappend = nappend + 1
      newsorted(curval) = me%nvals + nappend
      lastval = val(iadd)
      if (present(pos))  pos(iadd) = newsorted(curval)
      addedval(iadd) = .true.
      iadd = iadd + 1
    end if

    do ival=2,maxlen

      if ( (iadd <= veclen) .and. (iold <= me%nvals) ) then

        if ( me%val(me%sorted(iold)) <= val(iadd) ) then

          ! the original list's values are appended to newsorted before
          ! the additional list is appended.
          curval = curval + 1
          newsorted(curval) = me%sorted(iold)
          lastval = me%val(me%sorted(iold))
          iold = iold + 1

        else

          ! only append the value to unique lists, if it is not yet in the list.
          ! (if it is already in the list, it has to be the previous (curval-1)
          !  entry.)
          if ( lastval < val(iadd) ) then
            nappend = nappend + 1
            curval = curval + 1
            newsorted(curval) = me%nvals + nappend
            lastval = val(iadd)
            addedval(iadd) = .true.
          end if
          if (present(pos)) pos(iadd) = newsorted(curval)
          iadd = iadd + 1

        end if

      else

        ! reached the end of one or both of the sorted lists.
        exit

      end if

    end do

    if (iold <= me%nvals) then
      ! still some values from the original list left.
      newsorted(curval+1:me%nvals+nappend) = me%sorted(iold:me%nvals)
    end if

    if (iadd <= veclen) then
      ! still some values from the list to append left.
      rem_app = iadd
      do i = rem_app,veclen
        if ( lastval < val(iadd) ) then
          nappend = nappend + 1
          curval = curval + 1
          newsorted(curval) = me%nvals + nappend
          lastval = val(iadd)
          addedval(iadd) = .true.
        end if
        if (present(pos)) pos(iadd) = newsorted(curval)
        iadd = iadd + 1
      end do
    end if

    if (me%nvals > huge(me%nvals)-nappend) then
       write(*,*) "reached end of integer range for dynamic array!"
       write(*,*) "aborting!!"
       stop
    end if

    if (me%nvals + nappend > me%containersize) then
      call expand( me        = me,      &
        &          increment = nappend, &
        &          length    = length   )
    end if
    me%sorted(:me%nvals+nappend) = newsorted(:me%nvals+nappend)
    curval = me%nvals
    do iadd=1,veclen
      if (addedval(iadd)) then
        curval = curval + 1
        me%val(curval) = val(iadd)
      end if
    end do
    me%nvals = me%nvals + nappend

    if( present( wasadded ) ) wasadded = addedval

  end subroutine append_da_vecint
! ***************************************************************************** !


! ***************************************************************************** !
  !> truncate the array after the last valid entry and hence cut off the empty
  !! trailing empty entries
  !!
  subroutine truncate_da_int(me)
    !------------------------------------------------------------------------
    type(dyn_intarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------
    integer, allocatable :: swpval(:)
    integer, allocatable :: swpsort(:)
    !------------------------------------------------------------------------

    if (me%nvals < me%containersize) then
      allocate(swpval(me%nvals))
      allocate(swpsort(me%nvals))

      swpval = me%val(:me%nvals)
      swpsort = me%sorted(:me%nvals)

      call move_alloc(swpval, me%val)
      call move_alloc(swpsort, me%sorted)

      me%containersize = me%nvals
    end if

  end subroutine truncate_da_int
! ***************************************************************************** !


! ***************************************************************************** !
  !> empty all contents of the array without changing the size or status of any
  !! array
  !!
  subroutine empty_da_int(me)
    !------------------------------------------------------------------------
    type(dyn_intarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------
    ! reset the number of entries
    me%nvals = 0
  end subroutine empty_da_int
! ***************************************************************************** !


! ***************************************************************************** !
  !> fixing the dynamic array
  !!
  !! truncate the array after the last valid entry and hence cut off the empty
  !! trailing empty entries
  !! store the array in the sorted order according to the sorted( ) array
  !!
  subroutine sorttruncate_da_int(me)
    !------------------------------------------------------------------------
    type(dyn_intarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------
    type(dyn_intarray_type) :: tarray !< temporary array
    integer :: ival
    integer :: dpos
    !------------------------------------------------------------------------
    ! allocate the temporary array
    call init( me = tarray, length = me%nvals )
    ! copy the entries in a sorted fashion into the temporary array
    do ival = 1, me%nvals
      call append( me = tarray, val = me%val( me%sorted( ival )), &
           &       pos = dpos)
    enddo
    call destroy( me = me )

    me = tarray
    call destroy( me = tarray )

  end subroutine sorttruncate_da_int
! ***************************************************************************** !


! ***************************************************************************** !
  !> expanding the dynamic array
  !!
  !! this is a helping subroutine, which doubles the container
  !! of the given dynamic array. as the container might be
  !! initially 0-sized, a module variable minlength has been introduced, which
  !! is used here, to at least create a container of this size.
  subroutine expand_da_int(me, increment, length)
    !------------------------------------------------------------------------
    type(dyn_intarray_type) :: me !< array to resize
    integer, optional :: increment !< used for vector append
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !------------------------------------------------------------------------
    integer, allocatable :: swpval(:)
    integer, allocatable :: swpsort(:)
    !------------------------------------------------------------------------
    integer :: addvals, explen
    !------------------------------------------------------------------------

    addvals = 1
    if (present(increment)) addvals = increment

    if (addvals > 0) then

      ! if length is present, use that, otherwise double the size
      if( present( length ) ) then
        explen = length
      else
        ! set the global minimum length, if doubling would be smaller than that
        explen = max(me%containersize, minlength)
      end if

      ! check whether all elements will fit
      if( addvals > explen ) then
        explen = addvals
      end if

      ! check whether the new size will exceed the max container size.
      if( (huge(me%containersize) - explen) <= me%containersize ) then
        ! if so, expand to the maximum size
        me%containersize = huge(me%containersize)
      else
        ! if not, expand to the calculated size
        me%containersize = me%containersize + explen
      end if

      ! only need to do something, if there are actually values to append.
      if (me%nvals > 0) then

        allocate(swpval(me%containersize))
        swpval(1:me%nvals) = me%val(1:me%nvals)
        call move_alloc( swpval, me%val )

        allocate(swpsort(me%containersize))
        swpsort(1:me%nvals) = me%sorted(1:me%nvals)
        call move_alloc( swpsort, me%sorted )

      else ! me%nvals == 0

        if( allocated(me%val) ) &
          deallocate(me%val)
        allocate(me%val(me%containersize))
        if( allocated(me%sorted) ) &
          deallocate(me%sorted)
        allocate(me%sorted(me%containersize))

      end if

    end if

  end subroutine expand_da_int
! ***************************************************************************** !


! ***************************************************************************** !
  !> return the sorted position of a value in the given dynamic array
  !!
  !! if the value was not found,
  !!  - return 0 if nextifnotfound = .false.
  !!  - return position at the end if nextifnotfound = .true.
  function sortposofval_int(me, val, nextifnotfound, lower, upper) result(pos)
    !------------------------------------------------------------------------
    type(dyn_intarray_type), intent(in) :: me !< dynamic array
    integer, intent(in) :: val !< value to look for
    !> flag to indicate, if the next entry in the list should be returned,
    !! if the searched one is not found.
    logical, intent(in), optional :: nextifnotfound
    integer, intent(in), optional :: lower !< lower search limit
    integer, intent(in), optional :: upper !< upper search limit
    integer :: pos !< position of val in the sorted list, 0 if not found
    !------------------------------------------------------------------------
    logical :: retnext
    integer :: lb, ub
    integer :: mid
    integer :: lb_val, ub_val
    integer :: mid_val
    !------------------------------------------------------------------------

    retnext = .false.
    if (present(nextifnotfound)) retnext = nextifnotfound

    lb = 1
    ub = me%nvals

    if( present( lower ) ) lb = lower
    if( present( upper ) ) ub = upper

    pos = 0
    if (retnext) pos = lb

    !> binary search on sorted list
    do while(ub >= lb)
      lb_val = me%val(me%sorted(lb))

      if (val < lb_val) then
        if (retnext) pos = lb
        exit
      end if

      ub_val = me%val(me%sorted(ub))

      if (val > ub_val) then
        if (retnext) pos = ub+1
        exit
      end if

      ! safe guard against integer limit overflow
      mid = lb + (ub-lb) / 2
      mid_val = me%val(me%sorted(mid))
      if (val == mid_val) then
        pos = mid
        exit
      end if
      if (val > mid_val) then
        lb = mid + 1
      else
        ub = mid - 1
      end if
    end do
  end function sortposofval_int
! ***************************************************************************** !


! ***************************************************************************** !
  !> the actual position of a given value in the dynamic array
  !!
  !! most likely this is what you need in codes, using this
  !! data structure, it first does the binary search on the sorted
  !! values with sortposofval_int and then returns the looked
  !! up position in the original unsorted array, which corresponds
  !! to the position returned by the append routine.
  function posofval_int(me, val, nextifnotfound, lower, upper) result(pos)
    !------------------------------------------------------------------------
    type(dyn_intarray_type), intent(in) :: me !< dynamic array
    integer, intent(in) :: val !< value to search for
    !> flag to indicate, if the position of the next entry in the sorted
    !! list should be returned instead, if val is not found.
    logical, intent(in), optional :: nextifnotfound
    integer, intent(in), optional :: lower !< lower search limit
    integer, intent(in), optional :: upper !< upper search limit
    integer :: pos !< position in the array of the searche value, 0 if not found
    !------------------------------------------------------------------------
    integer :: sortpos
    integer :: lb, ub
    !------------------------------------------------------------------------

    lb = 1
    ub = me%nvals

    if( present( lower ) ) lb = lower
    if( present( upper ) ) ub = upper

    sortpos = sortedposofval(me, val, nextifnotfound, lb, ub)

    ! if result (sorted pos)
    if ((sortpos <= me%nvals) .and. (sortpos > 0)) then
      pos = me%sorted(sortpos)
    else
      pos = sortpos
    end if
  end function posofval_int
! ***************************************************************************** !
! ***************************************************************************** !
  !> initialization of a dynamic array
  !!
  !! before a dynamic array can be used, it has to be initialized
  !! with this routine. the initial length provided here, can
  !! avoid reallocations and memory copying, if approximated
  !! correctly enough. if none is specified, the provided container
  !! initially will be of size 0.
  subroutine init_da_real(me, length)
    !-----------------------------------------------------------------
    type(dyn_realarray_type), intent(out) :: me !< dynamic array to init
    integer, intent(in), optional :: length !< initial length of the container
    !-----------------------------------------------------------------

    if (present(length)) then
      me%containersize = length
    else
      me%containersize = zerolength
    end if

    ! deallocate ...
    if( allocated( me%val ) ) deallocate(me%val)
    if( allocated( me%sorted ) ) deallocate(me%sorted)
    ! ... and reallocate
    allocate(me%val(me%containersize))
    allocate(me%sorted(me%containersize))
    me%nvals = 0

  end subroutine init_da_real

  !> destruction of a dynamic array
  !!
  !! this subroutine takes care of a proper destruction of a
  !! dynamic array, it frees the allocated memory and resets
  !! the internal counts to 0.
  subroutine destroy_da_real(me)
    type(dyn_realarray_type), intent(inout) :: me !< dynamic array to init

    me%containersize = 0
    me%nvals         = 0

    if( allocated( me%val ) ) deallocate(me%val)
    if( allocated( me%sorted ) ) deallocate(me%sorted)
  end subroutine destroy_da_real
! ***************************************************************************** !


! ***************************************************************************** !
  !> appending a value to the dynamic array
  !!
  !! with this subroutine, a given value can be added to the
  !! dynamic array. the actual position of this value in the
  !! dynamic array will be returned, so it can be found again
  !! easily later. with the wasadded flag, it is indicated,\n
  !! wasadded = true,  if this entry had to be added,\n
  !! wasadded = false, if this was already found in the array.
  subroutine append_da_real(me, val, length, pos, wasadded )
    !------------------------------------------------------------------------
    type(dyn_realarray_type) :: me   !< array to append the value to
    real(kind=rk), intent(in)     :: val  !< value to append
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> position in the array, if the value is found
    integer, intent(out), optional :: pos
    !> flag to indicate, if val was newly added
    logical, intent(out), optional :: wasadded
    !------------------------------------------------------------------------
    integer :: foundpos
    integer :: i
    !------------------------------------------------------------------------

    ! do a binary search on existing entries (returns closest entry next to
    ! it if not found).
    foundpos = sortedposofval(me, val, .true.)
    if( present( wasadded ) ) wasadded = .false.

    ! if it found the value, the position is smaller than nvals
    if (foundpos <= me%nvals) then

      ! the returned position might actually be the right entry already or
      ! not, check for this here.
      if ( me%val(me%sorted(foundpos)) == val ) then

        ! found the value in a list of unique values,
        ! nothing to do, just return its position.
        if( present( pos ) ) pos = me%sorted(foundpos)

      else

        ! need to append a new value!

        if (me%nvals == huge(me%nvals)) then
           write(*,*) "reached end of integer range for dynamic array!"
           write(*,*) "aborting!!"
           stop
        end if

        if( present( wasadded ) ) wasadded = .true.
        if (me%nvals == me%containersize) then

          ! container is full, need to expand it
          call expand(me = me, length = length)
        end if
        me%nvals = me%nvals + 1

        ! put the new value into the last position in the
        ! array.
        me%val(me%nvals) = val
        do while( foundpos < me%nvals )
          if(me%val(me%sorted(foundpos)) /= val) then
            exit
          end if
          ! in case of multiple entries with the same value
          ! move on to the first differing entry.
          foundpos = foundpos + 1
        end do
        ! shift the sorted list of indices, to create a
        ! whole for the value to be inserted, at position
        ! foundpos.
        do i=me%nvals-1,foundpos,-1
          me%sorted(i+1) = me%sorted(i)
        end do
        ! put the index of the new value into the
        ! sorted list at the now freed position.
        me%sorted(foundpos) = me%nvals
        if( present( pos ) ) pos = me%nvals

      end if

    else

      ! value to append is larger than all existing ones,
      ! just put it to the end of the list, this captures
      ! also the case of empty lists.
      ! in this case foundpos = me%nvals + 1 holds.
      if( present( wasadded ) ) wasadded = .true.
      if (foundpos > me%containersize) then
        ! expand the array, if its boundary is reached
        call expand(me = me, length = length)
      end if
      me%nvals = foundpos
      me%val(foundpos) = val
      me%sorted(foundpos) = foundpos
      if( present( pos ) ) pos = foundpos

    end if

  end subroutine append_da_real
! ***************************************************************************** !


! ***************************************************************************** !
  !> appending a sorted list of values to the dynamic array
  !!
  !! with this subroutine, a given list of sorted values can be added to the
  !! dynamic array. the actual positions of these values in the
  !! dynamic array will be returned, so it can be found again
  !! easily later. with the wasadded flag, it is indicated,\n
  !! wasadded = true,  if this entry had to be added,\n
  !! wasadded = false, if this was already found in the array.
  subroutine append_da_vecreal(me, val, length, pos, wasadded )
    !------------------------------------------------------------------------
    type(dyn_realarray_type) :: me   !< array to append the value to
    real(kind=rk), intent(in)     :: val(:)  !< values to append
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> position in the array, the values are found at.
    integer, intent(out), optional :: pos(:)
    !> flag to indicate, if val was newly added
    logical, intent(out), optional :: wasadded(:)
    !------------------------------------------------------------------------
    real(kind=rk) :: lastval
    logical :: addedval(size(val))
    integer :: i
    integer :: veclen
    integer :: maxlen
    integer :: nappend
    integer :: rem_app
    integer :: curval, ival, iold, iadd
    integer, allocatable :: newsorted(:)
    !------------------------------------------------------------------------

    if (size(val) == 0) return

    veclen = size(val)
    maxlen = veclen + me%nvals

    allocate(newsorted(maxlen))

    addedval = .false.

    iold = 1
    iadd = 1

    nappend = 0
    curval = 0

    ! select the first entry before the loop unconditionally without checks
    ! for uniqueness (nothing to check against yet).
    if ( me%val(me%sorted(iold)) <= val(iadd) ) then
      curval = curval + 1
      newsorted(curval) = me%sorted(iold)
      lastval = me%val(me%sorted(iold))
      iold = iold + 1
    else
      curval = curval + 1
      nappend = nappend + 1
      newsorted(curval) = me%nvals + nappend
      lastval = val(iadd)
      if (present(pos))  pos(iadd) = newsorted(curval)
      addedval(iadd) = .true.
      iadd = iadd + 1
    end if

    do ival=2,maxlen

      if ( (iadd <= veclen) .and. (iold <= me%nvals) ) then

        if ( me%val(me%sorted(iold)) <= val(iadd) ) then

          ! the original list's values are appended to newsorted before
          ! the additional list is appended.
          curval = curval + 1
          newsorted(curval) = me%sorted(iold)
          lastval = me%val(me%sorted(iold))
          iold = iold + 1

        else

          ! only append the value to unique lists, if it is not yet in the list.
          ! (if it is already in the list, it has to be the previous (curval-1)
          !  entry.)
          if ( lastval < val(iadd) ) then
            nappend = nappend + 1
            curval = curval + 1
            newsorted(curval) = me%nvals + nappend
            lastval = val(iadd)
            addedval(iadd) = .true.
          end if
          if (present(pos)) pos(iadd) = newsorted(curval)
          iadd = iadd + 1

        end if

      else

        ! reached the end of one or both of the sorted lists.
        exit

      end if

    end do

    if (iold <= me%nvals) then
      ! still some values from the original list left.
      newsorted(curval+1:me%nvals+nappend) = me%sorted(iold:me%nvals)
    end if

    if (iadd <= veclen) then
      ! still some values from the list to append left.
      rem_app = iadd
      do i = rem_app,veclen
        if ( lastval < val(iadd) ) then
          nappend = nappend + 1
          curval = curval + 1
          newsorted(curval) = me%nvals + nappend
          lastval = val(iadd)
          addedval(iadd) = .true.
        end if
        if (present(pos)) pos(iadd) = newsorted(curval)
        iadd = iadd + 1
      end do
    end if

    if (me%nvals > huge(me%nvals)-nappend) then
       write(*,*) "reached end of integer range for dynamic array!"
       write(*,*) "aborting!!"
       stop
    end if

    if (me%nvals + nappend > me%containersize) then
      call expand( me        = me,      &
        &          increment = nappend, &
        &          length    = length   )
    end if
    me%sorted(:me%nvals+nappend) = newsorted(:me%nvals+nappend)
    curval = me%nvals
    do iadd=1,veclen
      if (addedval(iadd)) then
        curval = curval + 1
        me%val(curval) = val(iadd)
      end if
    end do
    me%nvals = me%nvals + nappend

    if( present( wasadded ) ) wasadded = addedval

  end subroutine append_da_vecreal
! ***************************************************************************** !


! ***************************************************************************** !
  !> truncate the array after the last valid entry and hence cut off the empty
  !! trailing empty entries
  !!
  subroutine truncate_da_real(me)
    !------------------------------------------------------------------------
    type(dyn_realarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------
    real(kind=rk), allocatable :: swpval(:)
    integer, allocatable :: swpsort(:)
    !------------------------------------------------------------------------

    if (me%nvals < me%containersize) then
      allocate(swpval(me%nvals))
      allocate(swpsort(me%nvals))

      swpval = me%val(:me%nvals)
      swpsort = me%sorted(:me%nvals)

      call move_alloc(swpval, me%val)
      call move_alloc(swpsort, me%sorted)

      me%containersize = me%nvals
    end if

  end subroutine truncate_da_real
! ***************************************************************************** !


! ***************************************************************************** !
  !> empty all contents of the array without changing the size or status of any
  !! array
  !!
  subroutine empty_da_real(me)
    !------------------------------------------------------------------------
    type(dyn_realarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------
    ! reset the number of entries
    me%nvals = 0
  end subroutine empty_da_real
! ***************************************************************************** !


! ***************************************************************************** !
  !> fixing the dynamic array
  !!
  !! truncate the array after the last valid entry and hence cut off the empty
  !! trailing empty entries
  !! store the array in the sorted order according to the sorted( ) array
  !!
  subroutine sorttruncate_da_real(me)
    !------------------------------------------------------------------------
    type(dyn_realarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------
    type(dyn_realarray_type) :: tarray !< temporary array
    integer :: ival
    integer :: dpos
    !------------------------------------------------------------------------
    ! allocate the temporary array
    call init( me = tarray, length = me%nvals )
    ! copy the entries in a sorted fashion into the temporary array
    do ival = 1, me%nvals
      call append( me = tarray, val = me%val( me%sorted( ival )), &
           &       pos = dpos)
    enddo
    call destroy( me = me )

    me = tarray
    call destroy( me = tarray )

  end subroutine sorttruncate_da_real
! ***************************************************************************** !


! ***************************************************************************** !
  !> expanding the dynamic array
  !!
  !! this is a helping subroutine, which doubles the container
  !! of the given dynamic array. as the container might be
  !! initially 0-sized, a module variable minlength has been introduced, which
  !! is used here, to at least create a container of this size.
  subroutine expand_da_real(me, increment, length)
    !------------------------------------------------------------------------
    type(dyn_realarray_type) :: me !< array to resize
    integer, optional :: increment !< used for vector append
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !------------------------------------------------------------------------
    real(kind=rk), allocatable :: swpval(:)
    integer, allocatable :: swpsort(:)
    !------------------------------------------------------------------------
    integer :: addvals, explen
    !------------------------------------------------------------------------

    addvals = 1
    if (present(increment)) addvals = increment

    if (addvals > 0) then

      ! if length is present, use that, otherwise double the size
      if( present( length ) ) then
        explen = length
      else
        ! set the global minimum length, if doubling would be smaller than that
        explen = max(me%containersize, minlength)
      end if

      ! check whether all elements will fit
      if( addvals > explen ) then
        explen = addvals
      end if

      ! check whether the new size will exceed the max container size.
      if( (huge(me%containersize) - explen) <= me%containersize ) then
        ! if so, expand to the maximum size
        me%containersize = huge(me%containersize)
      else
        ! if not, expand to the calculated size
        me%containersize = me%containersize + explen
      end if

      ! only need to do something, if there are actually values to append.
      if (me%nvals > 0) then

        allocate(swpval(me%containersize))
        swpval(1:me%nvals) = me%val(1:me%nvals)
        call move_alloc( swpval, me%val )

        allocate(swpsort(me%containersize))
        swpsort(1:me%nvals) = me%sorted(1:me%nvals)
        call move_alloc( swpsort, me%sorted )

      else ! me%nvals == 0

        if( allocated(me%val) ) &
          deallocate(me%val)
        allocate(me%val(me%containersize))
        if( allocated(me%sorted) ) &
          deallocate(me%sorted)
        allocate(me%sorted(me%containersize))

      end if

    end if

  end subroutine expand_da_real
! ***************************************************************************** !


! ***************************************************************************** !
  !> return the sorted position of a value in the given dynamic array
  !!
  !! if the value was not found,
  !!  - return 0 if nextifnotfound = .false.
  !!  - return position at the end if nextifnotfound = .true.
  function sortposofval_real(me, val, nextifnotfound, lower, upper) result(pos)
    !------------------------------------------------------------------------
    type(dyn_realarray_type), intent(in) :: me !< dynamic array
    real(kind=rk), intent(in) :: val !< value to look for
    !> flag to indicate, if the next entry in the list should be returned,
    !! if the searched one is not found.
    logical, intent(in), optional :: nextifnotfound
    integer, intent(in), optional :: lower !< lower search limit
    integer, intent(in), optional :: upper !< upper search limit
    integer :: pos !< position of val in the sorted list, 0 if not found
    !------------------------------------------------------------------------
    logical :: retnext
    integer :: lb, ub
    integer :: mid
    real(kind=rk) :: lb_val, ub_val
    real(kind=rk) :: mid_val
    !------------------------------------------------------------------------

    retnext = .false.
    if (present(nextifnotfound)) retnext = nextifnotfound

    lb = 1
    ub = me%nvals

    if( present( lower ) ) lb = lower
    if( present( upper ) ) ub = upper

    pos = 0
    if (retnext) pos = lb

    !> binary search on sorted list
    do while(ub >= lb)
      lb_val = me%val(me%sorted(lb))

      if (val < lb_val) then
        if (retnext) pos = lb
        exit
      end if

      ub_val = me%val(me%sorted(ub))

      if (val > ub_val) then
        if (retnext) pos = ub+1
        exit
      end if

      ! safe guard against integer limit overflow
      mid = lb + (ub-lb) / 2
      mid_val = me%val(me%sorted(mid))
      if (val == mid_val) then
        pos = mid
        exit
      end if
      if (val > mid_val) then
        lb = mid + 1
      else
        ub = mid - 1
      end if
    end do
  end function sortposofval_real
! ***************************************************************************** !


! ***************************************************************************** !
  !> the actual position of a given value in the dynamic array
  !!
  !! most likely this is what you need in codes, using this
  !! data structure, it first does the binary search on the sorted
  !! values with sortposofval_real and then returns the looked
  !! up position in the original unsorted array, which corresponds
  !! to the position returned by the append routine.
  function posofval_real(me, val, nextifnotfound, lower, upper) result(pos)
    !------------------------------------------------------------------------
    type(dyn_realarray_type), intent(in) :: me !< dynamic array
    real(kind=rk), intent(in) :: val !< value to search for
    !> flag to indicate, if the position of the next entry in the sorted
    !! list should be returned instead, if val is not found.
    logical, intent(in), optional :: nextifnotfound
    integer, intent(in), optional :: lower !< lower search limit
    integer, intent(in), optional :: upper !< upper search limit
    integer :: pos !< position in the array of the searche value, 0 if not found
    !------------------------------------------------------------------------
    integer :: sortpos
    integer :: lb, ub
    !------------------------------------------------------------------------

    lb = 1
    ub = me%nvals

    if( present( lower ) ) lb = lower
    if( present( upper ) ) ub = upper

    sortpos = sortedposofval(me, val, nextifnotfound, lb, ub)

    ! if result (sorted pos)
    if ((sortpos <= me%nvals) .and. (sortpos > 0)) then
      pos = me%sorted(sortpos)
    else
      pos = sortpos
    end if
  end function posofval_real
! ***************************************************************************** !
! ***************************************************************************** !
  !> initialization of a dynamic array
  !!
  !! before a dynamic array can be used, it has to be initialized
  !! with this routine. the initial length provided here, can
  !! avoid reallocations and memory copying, if approximated
  !! correctly enough. if none is specified, the provided container
  !! initially will be of size 0.
  subroutine init_da_label(me, length)
    !-----------------------------------------------------------------
    type(dyn_labelarray_type), intent(out) :: me !< dynamic array to init
    integer, intent(in), optional :: length !< initial length of the container
    !-----------------------------------------------------------------

    if (present(length)) then
      me%containersize = length
    else
      me%containersize = zerolength
    end if

    ! deallocate ...
    if( allocated( me%val ) ) deallocate(me%val)
    if( allocated( me%sorted ) ) deallocate(me%sorted)
    ! ... and reallocate
    allocate(me%val(me%containersize))
    allocate(me%sorted(me%containersize))
    me%nvals = 0

  end subroutine init_da_label

  !> destruction of a dynamic array
  !!
  !! this subroutine takes care of a proper destruction of a
  !! dynamic array, it frees the allocated memory and resets
  !! the internal counts to 0.
  subroutine destroy_da_label(me)
    type(dyn_labelarray_type), intent(inout) :: me !< dynamic array to init

    me%containersize = 0
    me%nvals         = 0

    if( allocated( me%val ) ) deallocate(me%val)
    if( allocated( me%sorted ) ) deallocate(me%sorted)
  end subroutine destroy_da_label
! ***************************************************************************** !


! ***************************************************************************** !
  !> appending a value to the dynamic array
  !!
  !! with this subroutine, a given value can be added to the
  !! dynamic array. the actual position of this value in the
  !! dynamic array will be returned, so it can be found again
  !! easily later. with the wasadded flag, it is indicated,\n
  !! wasadded = true,  if this entry had to be added,\n
  !! wasadded = false, if this was already found in the array.
  subroutine append_da_label(me, val, length, pos, wasadded )
    !------------------------------------------------------------------------
    type(dyn_labelarray_type) :: me   !< array to append the value to
    character(len=*), intent(in)     :: val  !< value to append
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> position in the array, if the value is found
    integer, intent(out), optional :: pos
    !> flag to indicate, if val was newly added
    logical, intent(out), optional :: wasadded
    !------------------------------------------------------------------------
    integer :: foundpos
    integer :: i
    !------------------------------------------------------------------------

    ! do a binary search on existing entries (returns closest entry next to
    ! it if not found).
    foundpos = sortedposofval(me, val, .true.)
    if( present( wasadded ) ) wasadded = .false.

    ! if it found the value, the position is smaller than nvals
    if (foundpos <= me%nvals) then

      ! the returned position might actually be the right entry already or
      ! not, check for this here.
      if ( me%val(me%sorted(foundpos)) == val ) then

        ! found the value in a list of unique values,
        ! nothing to do, just return its position.
        if( present( pos ) ) pos = me%sorted(foundpos)

      else

        ! need to append a new value!

        if (me%nvals == huge(me%nvals)) then
           write(*,*) "reached end of integer range for dynamic array!"
           write(*,*) "aborting!!"
           stop
        end if

        if( present( wasadded ) ) wasadded = .true.
        if (me%nvals == me%containersize) then

          ! container is full, need to expand it
          call expand(me = me, length = length)
        end if
        me%nvals = me%nvals + 1

        ! put the new value into the last position in the
        ! array.
        me%val(me%nvals) = val
        do while( foundpos < me%nvals )
          if(me%val(me%sorted(foundpos)) /= val) then
            exit
          end if
          ! in case of multiple entries with the same value
          ! move on to the first differing entry.
          foundpos = foundpos + 1
        end do
        ! shift the sorted list of indices, to create a
        ! whole for the value to be inserted, at position
        ! foundpos.
        do i=me%nvals-1,foundpos,-1
          me%sorted(i+1) = me%sorted(i)
        end do
        ! put the index of the new value into the
        ! sorted list at the now freed position.
        me%sorted(foundpos) = me%nvals
        if( present( pos ) ) pos = me%nvals

      end if

    else

      ! value to append is larger than all existing ones,
      ! just put it to the end of the list, this captures
      ! also the case of empty lists.
      ! in this case foundpos = me%nvals + 1 holds.
      if( present( wasadded ) ) wasadded = .true.
      if (foundpos > me%containersize) then
        ! expand the array, if its boundary is reached
        call expand(me = me, length = length)
      end if
      me%nvals = foundpos
      me%val(foundpos) = val
      me%sorted(foundpos) = foundpos
      if( present( pos ) ) pos = foundpos

    end if

  end subroutine append_da_label
! ***************************************************************************** !


! ***************************************************************************** !
  !> appending a sorted list of values to the dynamic array
  !!
  !! with this subroutine, a given list of sorted values can be added to the
  !! dynamic array. the actual positions of these values in the
  !! dynamic array will be returned, so it can be found again
  !! easily later. with the wasadded flag, it is indicated,\n
  !! wasadded = true,  if this entry had to be added,\n
  !! wasadded = false, if this was already found in the array.
  subroutine append_da_veclabel(me, val, length, pos, wasadded )
    !------------------------------------------------------------------------
    type(dyn_labelarray_type) :: me   !< array to append the value to
    character(len=*), intent(in)     :: val(:)  !< values to append
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> position in the array, the values are found at.
    integer, intent(out), optional :: pos(:)
    !> flag to indicate, if val was newly added
    logical, intent(out), optional :: wasadded(:)
    !------------------------------------------------------------------------
    character(len=labellen) :: lastval
    logical :: addedval(size(val))
    integer :: i
    integer :: veclen
    integer :: maxlen
    integer :: nappend
    integer :: rem_app
    integer :: curval, ival, iold, iadd
    integer, allocatable :: newsorted(:)
    !------------------------------------------------------------------------

    if (size(val) == 0) return

    veclen = size(val)
    maxlen = veclen + me%nvals

    allocate(newsorted(maxlen))

    addedval = .false.

    iold = 1
    iadd = 1

    nappend = 0
    curval = 0

    ! select the first entry before the loop unconditionally without checks
    ! for uniqueness (nothing to check against yet).
    if ( me%val(me%sorted(iold)) <= val(iadd) ) then
      curval = curval + 1
      newsorted(curval) = me%sorted(iold)
      lastval = me%val(me%sorted(iold))
      iold = iold + 1
    else
      curval = curval + 1
      nappend = nappend + 1
      newsorted(curval) = me%nvals + nappend
      lastval = val(iadd)
      if (present(pos))  pos(iadd) = newsorted(curval)
      addedval(iadd) = .true.
      iadd = iadd + 1
    end if

    do ival=2,maxlen

      if ( (iadd <= veclen) .and. (iold <= me%nvals) ) then

        if ( me%val(me%sorted(iold)) <= val(iadd) ) then

          ! the original list's values are appended to newsorted before
          ! the additional list is appended.
          curval = curval + 1
          newsorted(curval) = me%sorted(iold)
          lastval = me%val(me%sorted(iold))
          iold = iold + 1

        else

          ! only append the value to unique lists, if it is not yet in the list.
          ! (if it is already in the list, it has to be the previous (curval-1)
          !  entry.)
          if ( lastval < val(iadd) ) then
            nappend = nappend + 1
            curval = curval + 1
            newsorted(curval) = me%nvals + nappend
            lastval = val(iadd)
            addedval(iadd) = .true.
          end if
          if (present(pos)) pos(iadd) = newsorted(curval)
          iadd = iadd + 1

        end if

      else

        ! reached the end of one or both of the sorted lists.
        exit

      end if

    end do

    if (iold <= me%nvals) then
      ! still some values from the original list left.
      newsorted(curval+1:me%nvals+nappend) = me%sorted(iold:me%nvals)
    end if

    if (iadd <= veclen) then
      ! still some values from the list to append left.
      rem_app = iadd
      do i = rem_app,veclen
        if ( lastval < val(iadd) ) then
          nappend = nappend + 1
          curval = curval + 1
          newsorted(curval) = me%nvals + nappend
          lastval = val(iadd)
          addedval(iadd) = .true.
        end if
        if (present(pos)) pos(iadd) = newsorted(curval)
        iadd = iadd + 1
      end do
    end if

    if (me%nvals > huge(me%nvals)-nappend) then
       write(*,*) "reached end of integer range for dynamic array!"
       write(*,*) "aborting!!"
       stop
    end if

    if (me%nvals + nappend > me%containersize) then
      call expand( me        = me,      &
        &          increment = nappend, &
        &          length    = length   )
    end if
    me%sorted(:me%nvals+nappend) = newsorted(:me%nvals+nappend)
    curval = me%nvals
    do iadd=1,veclen
      if (addedval(iadd)) then
        curval = curval + 1
        me%val(curval) = val(iadd)
      end if
    end do
    me%nvals = me%nvals + nappend

    if( present( wasadded ) ) wasadded = addedval

  end subroutine append_da_veclabel
! ***************************************************************************** !


! ***************************************************************************** !
  !> truncate the array after the last valid entry and hence cut off the empty
  !! trailing empty entries
  !!
  subroutine truncate_da_label(me)
    !------------------------------------------------------------------------
    type(dyn_labelarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------
    character(len=labellen), allocatable :: swpval(:)
    integer, allocatable :: swpsort(:)
    !------------------------------------------------------------------------

    if (me%nvals < me%containersize) then
      allocate(swpval(me%nvals))
      allocate(swpsort(me%nvals))

      swpval = me%val(:me%nvals)
      swpsort = me%sorted(:me%nvals)

      call move_alloc(swpval, me%val)
      call move_alloc(swpsort, me%sorted)

      me%containersize = me%nvals
    end if

  end subroutine truncate_da_label
! ***************************************************************************** !


! ***************************************************************************** !
  !> empty all contents of the array without changing the size or status of any
  !! array
  !!
  subroutine empty_da_label(me)
    !------------------------------------------------------------------------
    type(dyn_labelarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------
    ! reset the number of entries
    me%nvals = 0
  end subroutine empty_da_label
! ***************************************************************************** !


! ***************************************************************************** !
  !> fixing the dynamic array
  !!
  !! truncate the array after the last valid entry and hence cut off the empty
  !! trailing empty entries
  !! store the array in the sorted order according to the sorted( ) array
  !!
  subroutine sorttruncate_da_label(me)
    !------------------------------------------------------------------------
    type(dyn_labelarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------
    type(dyn_labelarray_type) :: tarray !< temporary array
    integer :: ival
    integer :: dpos
    !------------------------------------------------------------------------
    ! allocate the temporary array
    call init( me = tarray, length = me%nvals )
    ! copy the entries in a sorted fashion into the temporary array
    do ival = 1, me%nvals
      call append( me = tarray, val = me%val( me%sorted( ival )), &
           &       pos = dpos)
    enddo
    call destroy( me = me )

    me = tarray
    call destroy( me = tarray )

  end subroutine sorttruncate_da_label
! ***************************************************************************** !


! ***************************************************************************** !
  !> expanding the dynamic array
  !!
  !! this is a helping subroutine, which doubles the container
  !! of the given dynamic array. as the container might be
  !! initially 0-sized, a module variable minlength has been introduced, which
  !! is used here, to at least create a container of this size.
  subroutine expand_da_label(me, increment, length)
    !------------------------------------------------------------------------
    type(dyn_labelarray_type) :: me !< array to resize
    integer, optional :: increment !< used for vector append
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !------------------------------------------------------------------------
    character(len=labellen), allocatable :: swpval(:)
    integer, allocatable :: swpsort(:)
    !------------------------------------------------------------------------
    integer :: addvals, explen
    !------------------------------------------------------------------------

    addvals = 1
    if (present(increment)) addvals = increment

    if (addvals > 0) then

      ! if length is present, use that, otherwise double the size
      if( present( length ) ) then
        explen = length
      else
        ! set the global minimum length, if doubling would be smaller than that
        explen = max(me%containersize, minlength)
      end if

      ! check whether all elements will fit
      if( addvals > explen ) then
        explen = addvals
      end if

      ! check whether the new size will exceed the max container size.
      if( (huge(me%containersize) - explen) <= me%containersize ) then
        ! if so, expand to the maximum size
        me%containersize = huge(me%containersize)
      else
        ! if not, expand to the calculated size
        me%containersize = me%containersize + explen
      end if

      ! only need to do something, if there are actually values to append.
      if (me%nvals > 0) then

        allocate(swpval(me%containersize))
        swpval(1:me%nvals) = me%val(1:me%nvals)
        call move_alloc( swpval, me%val )

        allocate(swpsort(me%containersize))
        swpsort(1:me%nvals) = me%sorted(1:me%nvals)
        call move_alloc( swpsort, me%sorted )

      else ! me%nvals == 0

        if( allocated(me%val) ) &
          deallocate(me%val)
        allocate(me%val(me%containersize))
        if( allocated(me%sorted) ) &
          deallocate(me%sorted)
        allocate(me%sorted(me%containersize))

      end if

    end if

  end subroutine expand_da_label
! ***************************************************************************** !


! ***************************************************************************** !
  !> return the sorted position of a value in the given dynamic array
  !!
  !! if the value was not found,
  !!  - return 0 if nextifnotfound = .false.
  !!  - return position at the end if nextifnotfound = .true.
  function sortposofval_label(me, val, nextifnotfound, lower, upper) result(pos)
    !------------------------------------------------------------------------
    type(dyn_labelarray_type), intent(in) :: me !< dynamic array
    character(len=*), intent(in) :: val !< value to look for
    !> flag to indicate, if the next entry in the list should be returned,
    !! if the searched one is not found.
    logical, intent(in), optional :: nextifnotfound
    integer, intent(in), optional :: lower !< lower search limit
    integer, intent(in), optional :: upper !< upper search limit
    integer :: pos !< position of val in the sorted list, 0 if not found
    !------------------------------------------------------------------------
    logical :: retnext
    integer :: lb, ub
    integer :: mid
    character(len=labellen) :: lb_val, ub_val
    character(len=labellen) :: mid_val
    !------------------------------------------------------------------------

    retnext = .false.
    if (present(nextifnotfound)) retnext = nextifnotfound

    lb = 1
    ub = me%nvals

    if( present( lower ) ) lb = lower
    if( present( upper ) ) ub = upper

    pos = 0
    if (retnext) pos = lb

    !> binary search on sorted list
    do while(ub >= lb)
      lb_val = me%val(me%sorted(lb))

      if (val < lb_val) then
        if (retnext) pos = lb
        exit
      end if

      ub_val = me%val(me%sorted(ub))

      if (val > ub_val) then
        if (retnext) pos = ub+1
        exit
      end if

      ! safe guard against integer limit overflow
      mid = lb + (ub-lb) / 2
      mid_val = me%val(me%sorted(mid))
      if (val == mid_val) then
        pos = mid
        exit
      end if
      if (val > mid_val) then
        lb = mid + 1
      else
        ub = mid - 1
      end if
    end do
  end function sortposofval_label
! ***************************************************************************** !


! ***************************************************************************** !
  !> the actual position of a given value in the dynamic array
  !!
  !! most likely this is what you need in codes, using this
  !! data structure, it first does the binary search on the sorted
  !! values with sortposofval_label and then returns the looked
  !! up position in the original unsorted array, which corresponds
  !! to the position returned by the append routine.
  function posofval_label(me, val, nextifnotfound, lower, upper) result(pos)
    !------------------------------------------------------------------------
    type(dyn_labelarray_type), intent(in) :: me !< dynamic array
    character(len=*), intent(in) :: val !< value to search for
    !> flag to indicate, if the position of the next entry in the sorted
    !! list should be returned instead, if val is not found.
    logical, intent(in), optional :: nextifnotfound
    integer, intent(in), optional :: lower !< lower search limit
    integer, intent(in), optional :: upper !< upper search limit
    integer :: pos !< position in the array of the searche value, 0 if not found
    !------------------------------------------------------------------------
    integer :: sortpos
    integer :: lb, ub
    !------------------------------------------------------------------------

    lb = 1
    ub = me%nvals

    if( present( lower ) ) lb = lower
    if( present( upper ) ) ub = upper

    sortpos = sortedposofval(me, val, nextifnotfound, lb, ub)

    ! if result (sorted pos)
    if ((sortpos <= me%nvals) .and. (sortpos > 0)) then
      pos = me%sorted(sortpos)
    else
      pos = sortpos
    end if
  end function posofval_label
! ***************************************************************************** !
! ****************************************************************************** !


end module tem_dyn_array_module
! ****************************************************************************** !
!> Dynamic data structures
!!
!! During the processing of the mesh to create the solver data structures,
!! common tasks arise.
!! These include among others
!! a unique collection of elements in lists, i.e. without duplicates,
!! a localization of a stored element,
!! sorting of a list of records.
!! As Fortran does not intrinsically include dynamically growing data
!! structures, they first have to be introduced manually.
!!
!! \section tem_dyn_retrievePos Retrieving the position of a record
!! You might want to retrieve the position of a (unique) value within a list.
!! In order to achieve this efficiently, the list with \f$n\f$ entries must be
!! ordered so we can rely on the binary search, which performs with
!! \f$\mathcal{O}( \log n )\f$.
!!
!! \section tem_dyn_dataStructures Dynamic data structures
!! Hence we introduce a derived data type which allows an efficient and simple
!! element adding and retrieving process, as well as a sorted order of the
!! entries.
!! A [module] (@ref tem_dyn_array_module) provides the type definitions and all
!! the related procedures to act on this data type.
!! The data type includes the size of the container within, the amount of its
!! actual entries, a flag for uniqueness and the container array itself as
!! well as a conforming array for the sorted retrieval of the entries in the
!! value array.
!!
!! \todo What shall be listed in lst:tem_dyn_array ??
!!
!! The definition is depicted in listing \ref lst:tem_dyn_array where macros
!! are used for the label of a dynamic type `?tname?`and the used data type
!! `?tstring?`.
!! The introduction of macros allow for a generic
!! definition for all data types of the container,
!! whereas the compiled source code contains the definitions for each required
!! data type.
!!
!! In order to use such data structures they first have to be initialized with
!! options such as the uniqueness of the entries or the initial size of the
!! container.
!! The initialization simply allocates the container with the given options and
!! resets the number of records.
!! Correspondingly, the object can be destroyed by `destroy( me = valName )`.
!!
!!
!! If an object has been initialized, records can be added by simply calling
!! `append`.
!! Subsequently, checks are performed, for the existence of the entry in the
!! list and the sufficiency of the size of the container to incorporate another
!! entry.
!! In case the size does not suffice, the array is doubled in size through a
!! temporary array, and then the record is added.
!! Information about if and where the entry was stored is returned to the
!! caller.
