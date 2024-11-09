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
! ******************************************************************************
!> author: Kannan Masilamani
!! This module contains ellipsoid definition and routines related to ellipsoids

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

module tem_ellipsoid_module
  use env_module,                only: rk, minLength, labelLen, zeroLength
  use tem_aux_module,            only: tem_abort
  use tem_logging_module,        only: logunit

  use aotus_module,              only: flu_State, aot_get_val, &
    &                                  aoterr_Fatal, aoterr_NonExistent,       &
    &                                  aoterr_WrongType
  use aot_table_module,          only: aot_table_open, aot_table_close,        &
    &                                  aot_table_length

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

  public :: grw_ellipsoidArray_type
  public :: init, append, truncate, destroy, empty, placeAt
  public :: tem_ellipsoid_type, tem_load_ellipsoid
  public :: tem_ellipsoidCubeOverlap
  public :: tem_ellipsoid_out

  type tem_ellipsoid_type
    !> origin of the ellipsoid
    real(kind=rk) :: origin(3)

    !> radius of the ellipsoid
    real(kind=rk) :: radius(3)

    !> To choose what to do with intersection of this object
    !! if only_surface = true than the only the surface of the object
    !! is intersected
    !! if only_surface = false then the whole object is intersected
    !! default is set to false
    logical :: only_surface
  end type tem_ellipsoid_type

  !> growing array type for type(tem_ellipsoid_type)
  type grw_ellipsoidarray_type
    integer :: nvals = 0
    integer :: containersize = 0
    type(tem_ellipsoid_type), allocatable :: val(:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga_ellipsoid
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_ga_ellipsoid
  end interface

  !> empty the entries  without changing arrays
  interface empty
    module procedure empty_ga_ellipsoid
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga_ellipsoid
  end interface

  !> insert an element at a given position
  interface placeat
    module procedure placeat_ga_ellipsoid
    module procedure placeat_ga_ellipsoid_vec
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_ga_ellipsoid
    module procedure append_ga_ellipsoid_vec
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga_ellipsoid
  end interface


  !> interface to write out ellipsoids in lua format to a file
  interface tem_ellipsoid_out
    module procedure tem_ellipsoid_out_scal
    module procedure tem_ellipsoid_out_vec
  end interface tem_ellipsoid_out

  !> interface to load ellipsoids
  interface tem_load_ellipsoid
    module procedure tem_load_ellipsoid
    module procedure tem_load_ellipsoid_single
  end interface tem_load_ellipsoid

contains

  ! *****************************************************************************
  !> Load ellipsoid information from config file.
  subroutine tem_load_ellipsoid(me, transform, conf, thandle)
    ! --------------------------------------------------------------------------!
    !inferface variables
    !> array of ellipsoids
    type(tem_ellipsoid_type), allocatable, intent(out) :: me(:)
    !> transformation for spatial object
    type(tem_transformation_type), intent(in) :: transform
    !> lua state
    type(flu_state) :: conf
    integer, intent(in) :: thandle !< handle for canonical objects
    ! --------------------------------------------------------------------------!
    ! local varaibles
    integer :: sph_handle, sph_subHandle
    integer :: iObj, nObjects
    ! --------------------------------------------------------------------------!

    write(logunit(1),*) 'Loading ellipsoid: '

    call aot_table_open(L = conf, parent = thandle, thandle = sph_handle, &
      &                 key = 'object')
    call aot_table_open(L=conf, parent = sph_handle, thandle = sph_subHandle, &
      & pos = 1 )

    if ( sph_subHandle .eq. 0) then
      !object is a single table
      allocate(me(1))
      call aot_table_close(L=conf, thandle=sph_subHandle)
      call tem_load_ellipsoid_single( me(1), transform, conf, sph_handle )
    else
      !object is a multiple table
      call aot_table_close(L=conf, thandle=sph_subHandle)
      nObjects = aot_table_length(L=conf, thandle=sph_handle)
      allocate(me(nObjects))
      do iObj=1,nObjects
        call aot_table_open(L=conf, parent=sph_handle, thandle=sph_suBHandle,&
          & pos=iObj)
        call tem_load_ellipsoid_single( me(iObj), transform, conf, &
          &                             sph_Subhandle )
        call aot_table_close(L=conf, thandle=sph_subHandle)
      end do
    end if

    call aot_table_close(L=conf, thandle=sph_Handle)


  end subroutine tem_load_ellipsoid
  ! *****************************************************************************

  ! *****************************************************************************
  !> This routine single ellipsoid from object table
  subroutine tem_load_ellipsoid_single(me, transform, conf, thandle )
    ! --------------------------------------------------------------------------!
    !inferface variables
    !> single ellipsoid
    type(tem_ellipsoid_type), intent(out) :: me
    !> transformation for spatial object
    type(tem_transformation_type), intent(in) :: transform
    !> lua state
    type(flu_state) :: conf
    integer, intent(in) :: thandle !< handle for canonical objects
    ! --------------------------------------------------------------------------!
    integer :: ii, iError, vError(3), errFatal(3)
    ! --------------------------------------------------------------------------!
    errFatal = aoterr_fatal

    ! read origin of ellipsoid
    call aot_get_val(L=conf, thandle=thandle, val=me%origin, &
      &              ErrCode=vError, key='origin', pos = 1)
    if (any(btest(vError, errFatal))) then
      write(logunit(0),*) &
        &  ' Error in configuration: origin is not given to define a ellipsoid'
      call tem_abort()
    end if

    !read radius of  ellipsoid
    call aot_get_val(L=conf, thandle=thandle, val=me%radius, &
      &              ErrCode=vError, key='radius', pos=2 )
    if (any(btest(vError, errFatal))) then
      write(logunit(0),*) &
        &  ' Error in configuration: radius is not given to define a ellipsoid'
      call tem_abort()
    end if

    call aot_get_val(L=conf, thandle=thandle, val=me%only_surface, &
      &              ErrCode=iError, key='only_surface', &
      &              pos=3, default=.false.)

    if (btest(iError, aoterr_WrongType)) then
      write(logunit(0),*) 'Error occured, while retrieving ellipsoid only_surface'
      write(logunit(0),*) 'Variable has wrong type!'
      write(logunit(0),*) 'Should be a LOGICAL!'
      call tem_abort()
    endif

    write(logunit(1),"(A,3E12.5)") '        origin: ', me%origin
    write(logunit(1),"(A,3E12.5)") '        radius: ', me%radius
    write(logunit(1),"(A,L5)")      '  only surface: ', me%only_surface

    ! apply transformation to ellipsoid
    if(transform%active) then
      if(transform%deform%active) then
        write(logunit(5),"(A)")   '  apply deformation ...'
        do ii = 1, 3
          me%radius(ii) = me%radius(ii)     &
            &                      * transform%deform%matrix(ii,ii)
        end do
        me%origin = matmul(transform%deform%matrix, me%origin)
      endif
      if(transform%translate%active) then
        write(logunit(5),"(A)")   '  apply translation ...'
        me%origin = me%origin + transform%translate%vec
      endif
    endif

  end subroutine tem_load_ellipsoid_single
  ! ****************************************************************************

  ! ****************************************************************************
  !> This function checks intesection of solid cube and ellipsoid
  function tem_ellipsoidCubeOverlap(ellipsoid, cube) result(overlap)
    ! --------------------------------------------------------------------------!
    !inferface variables
    type(tem_ellipsoid_type), intent(in) :: ellipsoid !< spacer geometry data
    type(tem_cube_type), intent(in) :: cube
    logical :: overlap !< return value
    ! --------------------------------------------------------------------------!
    if(ellipsoid%only_surface) then
      overlap = hollowellipsoidCubeOverlap(ellipsoid, cube)
    else
      overlap = solidellipsoidCubeOverlap(ellipsoid, cube)
    endif

  end function tem_ellipsoidCubeOverlap
  ! ****************************************************************************

  ! ****************************************************************************
  !> This function checks intesection of solid cube and hollow ellipsoid
  !!
  !! This algorithm is taken from
  !! http://tog.acm.org/resources/GraphicsGems/gems/Boxellipsoid.c
  !!
  pure function hollowellipsoidCubeOverlap(me, cube) result(overlap)
    ! --------------------------------------------------------------------------!
    !> ellipsoid type
    type(tem_ellipsoid_type), intent(in) :: me
    !> cube type
    type(tem_cube_type), intent(in) :: cube
    logical :: overlap !< return value
    ! --------------------------------------------------------------------------!
    ! local variables
    real(kind=rk) :: rsqr,a, b
    integer :: i
    real(kind=rk) :: dmin, dmax
    ! --------------------------------------------------------------------------!

    !minimum distance
    dmin = 0.0_rk
    !maximum distance
    dmax = 0.0_rk

    rsqr = 1.0_rk

    do i=1,3
      ! a
      a = ( me%origin(i) - cube%origin(i) )**2 / me%radius(i)**2
      b = ( me%origin(i) - cube%endPnt(i) )**2 / me%radius(i)**2
      dmax = dmax + max(a,b)
      if ( me%origin(i) < cube%origin(i) ) then
        dmin = dmin + a
      else if ( me%origin(i) > cube%endPnt(i) ) then
        dmin = dmin + b
      end if
    end do

    overlap = ( (dmin <= rsqr) .and. (dmax >= rsqr) )

  end function hollowellipsoidCubeOverlap
  ! ****************************************************************************

  ! ****************************************************************************
  !> This function checks intesection of solid cube and solid ellipsoid
  !!
  !! This algorithm is taken from
  !! http://tog.acm.org/resources/GraphicsGems/gems/Boxellipsoid.c
  !!
  function solidellipsoidCubeOverlap(me, cube) result(overlap)
    ! --------------------------------------------------------------------------!
    !> ellipsoid object
    type(tem_ellipsoid_type), intent(in) :: me
    !> cube object
    type(tem_cube_type), intent(in) :: cube
    logical :: overlap !< return value
    ! --------------------------------------------------------------------------!
    ! local variables
    integer :: i
    real(kind=rk) :: dmin
    ! --------------------------------------------------------------------------!

    ! minimum distance
    dmin = 0.0_rk

    do i=1,3
      if ( me%origin(i) < cube%origin(i) ) then
        dmin = dmin + &
          &  ( me%origin(i) - cube%origin(i) )**2 / me%radius(i)**2
        ! dmin = dmin + ( me%origin(i) - cube%origin(i) )**2
      else if ( me%origin(i) > cube%endPnt(i) ) then
        dmin = dmin + &
          &  ( me%origin(i) - cube%endPnt(i) )**2 / me%radius(i)**2
        ! dmin = dmin +  ( me%origin(i) - cube%endPnt(i) )**2
      end if
    end do

    overlap = ( dmin <= 1.0_rk )

  end function solidellipsoidCubeOverlap
  ! ****************************************************************************

  ! ************************************************************************** !
  !> Write out an array of ellipsoids in lua format
  !!
  subroutine tem_ellipsoid_out_vec( me, conf )
    ! --------------------------------------------------------------------------
    !> ellipsoid types to write out
    type( tem_ellipsoid_type ), intent(in) :: me(:)
    !> Aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    ! --------------------------------------------------------------------------
    ! counter
    integer :: i
    ! --------------------------------------------------------------------------

    ! create a table with name ellipsoid
    call aot_out_open_table( put_conf = conf, tname = 'object' )

    do i = 1, size(me)
      call tem_ellipsoid_out_scal( me(i), conf )
    end do

    call aot_out_close_table( put_conf = conf )

  end subroutine tem_ellipsoid_out_vec
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Write out a ellipsoid shape in lua format
  !!
  subroutine tem_ellipsoid_out_scal( me, conf )
    ! --------------------------------------------------------------------------
    !> ellipsoid types to write out
    type( tem_ellipsoid_type ), intent(in) :: me
    !> Aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    ! --------------------------------------------------------------------------

    ! create a table with name ellipsoid if not exist
    if( conf%level == 0 ) then
      call aot_out_open_table( put_conf = conf, tname = 'object' )
    else
      call aot_out_open_table( put_conf = conf )
    end if

    call aot_out_val( put_conf = conf, vname = 'origin', val = me%origin )
    call aot_out_val( put_conf = conf, vname = 'radius', val = me%radius )
    call aot_out_val( put_conf = conf, vname = 'only_surface', &
      &               val = me%only_surface )

    call aot_out_close_table( put_conf = conf )

  end subroutine tem_ellipsoid_out_scal
  ! ************************************************************************** !


  subroutine init_ga_ellipsoid(me, length)
    type(grw_ellipsoidarray_type), intent(out) :: me !< dynamic array to init
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

  end subroutine init_ga_ellipsoid

  subroutine destroy_ga_ellipsoid(me)
    type(grw_ellipsoidarray_type), intent(inout) :: me !< dynamic array to destroy

    me%containersize = 0
    me%nvals = 0
    if( allocated( me%val ) ) deallocate(me%val)

  end subroutine destroy_ga_ellipsoid


  subroutine truncate_ga_ellipsoid(me)
    !------------------------------------------------------------------------
    type(grw_ellipsoidarray_type) :: me !< array to truncate
    !------------------------------------------------------------------------
    type(tem_ellipsoid_type), allocatable :: tarray(:)
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

  end subroutine truncate_ga_ellipsoid


  subroutine empty_ga_ellipsoid(me)
    !------------------------------------------------------------------------
    type(grw_ellipsoidarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------

    me%nvals = 0

  end subroutine empty_ga_ellipsoid

  !> adds the value to a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! element at the requested position will be replaced.
  subroutine placeat_ga_ellipsoid(me, val, pos, length)
    type(grw_ellipsoidarray_type) :: me !< array to place the value into
    type(tem_ellipsoid_type), intent(in) :: val !< value to place at the given position
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

  end subroutine placeat_ga_ellipsoid


  !> adds the values starting from a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! elements starting from the requested position will be replaced up to
  !! the element at position `pos + size(val) - 1`.
  subroutine placeat_ga_ellipsoid_vec(me, val, pos, length)
    type(grw_ellipsoidarray_type) :: me !< array to append the value to
    type(tem_ellipsoid_type), intent(in) :: val(:) !< values to append
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

  end subroutine placeat_ga_ellipsoid_vec


  subroutine append_ga_ellipsoid(me, val, length)
    type(grw_ellipsoidarray_type) :: me !< array to append the value to
    type(tem_ellipsoid_type), intent(in) :: val !< value to append
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

  end subroutine append_ga_ellipsoid

  subroutine append_ga_ellipsoid_vec(me, val, length)
    type(grw_ellipsoidarray_type) :: me !< array to append the value to
    type(tem_ellipsoid_type), intent(in) :: val(:) !< values to append
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

  end subroutine append_ga_ellipsoid_vec


  subroutine expand_ga_ellipsoid(me, pos, length)
    type(grw_ellipsoidarray_type) :: me !< array to resize
    integer, intent(in), optional :: pos !< optional predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    type(tem_ellipsoid_type), allocatable :: swpval(:)
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

  end subroutine expand_ga_ellipsoid


end module tem_ellipsoid_module

!> \page ellipsoid ellipsoid
!! ellipsoids are defined by an origin and radius.
!! ellipsoid is considered to be solid as default i.e. all the cubes inside the
!! ellipsoid are marked as intersected cubes.
!! It is possible to created hollow ellipsoids by setting only_surface = true,
!! it will mark only the cubes intersect with ellipsoid surface as intersected
!! cubes
!!
!! Valid definition:
!! \li Single ellipsoid
!! \verbatim
!! geometry={
!!   kind='ellipsoid',
!!     object={
!!       origin={0.0,0.0,0.0},
!!       radius=0.25,
!!       only_surface = true, -- If not defined default is set to false
!!     }
!! }
!! \endverbatim
!!
!! \li Multiple ellipsoid
!! \verbatim
!! geometry={
!!   kind='ellipsoid',
!!     object={
!!       {
!!       origin={0.0,0.0,0.0},
!!       radius=0.25
!!       },
!!       {
!!       origin={-2.0,0.0,0.0},
!!       radius=0.25
!!       }
!!     }
!! }
!! \endverbatim
!! \n\n
!! The following seeder file is to generate mesh with hollow ellipsoid (hollow => only_surface=true)
!!  inside:
!! \include testsuite/ellipsoid/seeder.lua
!! \n\n
!! Mesh generated with hollow ellipsoid by the seeder file:
!! \image html ellipsoid.png
!! \image html ellipsoid_withedges.png
!! \n\n
!! Cutview of mesh with hollow ellipsoid:
!! \image html ellipsoid_hollow.png
!! \n\n
!! Cutview of mesh with solid ellipsoid (solid => only_surface=false):
!! \image html ellipsoid_solid.png
!! \n\n
!! Example lua file is available at \link testsuite/ellipsoid/seeder.lua
!! \example testsuite/ellipsoid/seeder.lua

