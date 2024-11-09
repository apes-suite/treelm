! Copyright (c) 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!! This module contains cylinder definition and routines related to cylinders

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

module tem_cylinder_module
  use env_module,                only: rk, minLength, labelLen, zeroLength
  use tem_param_module,          only: PI
  use tem_aux_module,            only: tem_abort
  use tem_logging_module,        only: logunit
  use tem_cube_module,           only: tem_cube_type
  use tem_sphere_module,         only: tem_sphereCubeOverlap, tem_sphere_type
  use tem_line_module,           only: tem_line_type, tem_lineCubeOverlap
  use tem_transformation_module, only: tem_transformation_type
  use tem_float_module,          only: operator(.fle.), operator(.fge.)


  ! include aotus modules
  use aotus_module,     only: aot_get_val, aoterr_Fatal, aoterr_WrongType,     &
    &                         flu_State, aoterr_NonExistent
  use aot_table_module, only: aot_table_open, aot_table_close, aot_table_length
  use aot_out_module,   only: aot_out_type, aot_out_val,                       &
    &                         aot_out_open_table, aot_out_close_table

  implicit none
  private

  public :: grw_cylinderArray_type
  public :: init, append, truncate, destroy, empty, placeAt
  public :: tem_cylinder_type, tem_load_cylinder, tem_cylinderCubeOverlap
  public :: tem_cylinder_out

  !> This type provides information to
  !! create cylinder geometry
  type tem_cylinder_type
    !> vector defining length and axis of cylinder
    real(kind=rk) :: vec(3)
    real(kind=rk) :: radius !< radius of the cylinder
    real(kind=rk) :: origin(3) !< origin of the cylinder
    !> To choose what to do with intersection of this object
    !! if only_surface = true than the only the surface of the object
    !! is intersected
    !! if only_surface = false then the whole object is intersected
    !! default is set to false
    logical :: only_surface
  end type tem_cylinder_type

  !> growing array type for type(tem_cylinder_type)
  type grw_cylinderarray_type
    integer :: nvals = 0
    integer :: containersize = 0
    type(tem_cylinder_type), allocatable :: val(:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga_cylinder
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_ga_cylinder
  end interface

  !> empty the entries  without changing arrays
  interface empty
    module procedure empty_ga_cylinder
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga_cylinder
  end interface

  !> insert an element at a given position
  interface placeat
    module procedure placeat_ga_cylinder
    module procedure placeat_ga_cylinder_vec
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_ga_cylinder
    module procedure append_ga_cylinder_vec
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga_cylinder
  end interface


  !> interface to write out cylinders in lua format to a file
  interface tem_cylinder_out
    module procedure tem_cylinder_out_scal
    module procedure tem_cylinder_out_vec
  end interface tem_cylinder_out

  !> interface to load cylinders
  interface tem_load_cylinder
    module procedure tem_load_cylinder
    module procedure tem_load_cylinder_single
  end interface tem_load_cylinder

contains
  ! ****************************************************************************
  !> \brief Loading cylinder information from config file \n
  subroutine tem_load_cylinder(me, transform, conf, thandle)
    ! --------------------------------------------------------------------------!
    !inferface variables
    !> array of cylinders
    type(tem_cylinder_type), allocatable, intent(out) :: me(:)
    !> transformation for spatial object
    type(tem_transformation_type), intent(in) :: transform
    !> lua state
    type(flu_state) :: conf
    integer, intent(in) :: thandle !< handle for canonical objects
    ! --------------------------------------------------------------------------!
    ! local varaibles
    integer :: cyl_handle, cyl_subHandle
    integer :: iObj, nObjects
    ! --------------------------------------------------------------------------!

    write(logunit(1),*) 'Loading cylinder:'

    call aot_table_open(L = conf, parent = thandle, thandle = cyl_handle, &
      &                 key = 'object')
    call aot_table_open(L=conf, parent = cyl_handle, thandle = cyl_subHandle, &
      & pos = 1 )

    if ( cyl_subHandle .eq. 0) then
      !object is a single table
      allocate(me(1))
      call aot_table_close(L=conf, thandle=cyl_subHandle)
      call tem_load_cylinder_single( me(1), transform, conf, cyl_handle )
    else
      !object is a multiple table
      call aot_table_close(L=conf, thandle=cyl_subHandle)
      nObjects = aot_table_length(L=conf, thandle=cyl_handle)
      allocate(me(nObjects))
      do iObj=1,nObjects
        call aot_table_open(L=conf, parent=cyl_handle, thandle=cyl_suBHandle,&
          & pos=iObj)
        call tem_load_cylinder_single(me(iObj), transform, conf, cyl_Subhandle)
        call aot_table_close(L=conf, thandle=cyl_subHandle)
      end do
    end if

    call aot_table_close(L=conf, thandle=cyl_Handle)

  end subroutine tem_load_cylinder
  ! *****************************************************************************

  ! *****************************************************************************
  !> This routine single cylinder from object table
  subroutine tem_load_cylinder_single(me, transform, conf, thandle)
    ! --------------------------------------------------------------------------!
    !inferface variables
    type(tem_cylinder_type), intent(out) :: me !< cylinder data type
    type(flu_state) :: conf !< flu state
    integer, intent(in) :: thandle !< parent handle
    !> transformation for spatial object
    type(tem_transformation_type), intent(in) :: transform
    ! --------------------------------------------------------------------------!
    !local variable
    integer :: iError, vError(3), errFatal(3)
    ! --------------------------------------------------------------------------!

    errFatal = aoterr_fatal

    !read startPos distance of length from bounding box
    call aot_get_val(L=conf, thandle=thandle, &
    &              val=me%origin, ErrCode=vError, &
    &              key='origin', default=[0.0_rk,0.0_rk,0.0_rk] )
    if (any(btest(vError, errFatal))) then
      write(logunit(0),*) &
        &  'FATAL Error occured, while retrieving cylinder origin :'
      call tem_abort()
    end if

    !read radius of cylinder filament in length
    call aot_get_val(L=conf, thandle=thandle, &
      &              val=me%radius, ErrCode=iError, key='radius')
    if (btest(iError, aoterr_Fatal)) then
      write(logunit(0),*) &
        &  'FATAL Error occured, while retrieving cylinder radius'
      if (btest(iError, aoterr_NonExistent)) &
        &  write(logunit(0),*) 'Variable not existent!'
      if (btest(iError, aoterr_WrongType)) &
        &  write(logunit(0),*) 'Variable has wrong type!'
      call tem_abort()
    end if

    !read cylinder normal which defines the direction of the filament
    call aot_get_val(L=conf, thandle=thandle, &
      &              val=me%vec, ErrCode=vError, key='vec')
    if (any(btest(vError, errFatal))) then
      write(logunit(0),*) 'FATAL Error occured, while retrieving cylinder vec'
      call tem_abort()
    end if

    !cylinder type
    call aot_get_val(L=conf, thandle=thandle, val=me%only_surface, &
      &              ErrCode=iError, key='only_surface', &
      &              pos=4, default=.false.)

    if (btest(iError, aoterr_WrongType)) then
      write(logunit(0),*) &
        &  'Error occured, while retrieving cylinder only_surface'
      write(logunit(0),*) 'Variable has wrong type!'
      write(logunit(0),*) 'Should be a LOGICAL!'
      call tem_abort()
    endif

    write(logunit(1),*) '                  origin:', me%origin
    write(logunit(1),*) '                     vec:', me%vec
    write(logunit(1),*) '                  radius:', me%radius
    write(logunit(1),*) '            only_surface:', me%only_surface

    !apply transformation
    if(transform%active) then
      if(transform%deform%active) then
        me%vec = matmul(transform%deform%matrix, me%vec)
        me%origin = matmul(transform%deform%matrix, me%origin)
      endif
      if(transform%translate%active) then
        me%origin = transform%translate%vec + me%origin
      endif
    endif

  end subroutine tem_load_cylinder_single
  ! ****************************************************************************


  ! ****************************************************************************
  !> This function checks intesection of solid cube and cylinder.
  !!
  !! The test is done by projecting each cube vertices on the cylinder axis
  !! and check whether the projected point is within the cylinder length.
  !! If yes then use the projected point as the origin of sphere and
  !! do sphere-cube intersection.
  !!@todo HK: The algorithm used in here is not correct!
  !!@todo KM: BUG in defining cylinder with fixed length.
  !!          This implementation works only for cylinder with infinite length
  function tem_cylinderCubeOverlap(cylinder, cube) result(overlap)
    ! --------------------------------------------------------------------------!
    !inferface variables
    type(tem_cylinder_type), intent(in) :: cylinder !< cylinder geometry data
    type(tem_cube_type), intent(in) :: cube
    logical :: overlap !< return value
    ! --------------------------------------------------------------------------!
    ! local variables
    real(kind=rk) :: proj
    real(kind=rk) :: cubeVer(15,3)
    type(tem_sphere_type) :: sphere
    type(tem_line_type) :: line
    integer :: iVer
    real(kind=rk) :: pntIntersect(3)
    ! --------------------------------------------------------------------------!

    overlap = .false.

    !! we check for all vertices to find if any of the cube vertices
    !! intersect with cylinder
    cubeVer(1,:) = cube%center
    cubeVer(2,:) = cube%center + [ cube%halfwidth, 0.0_rk, 0.0_rk]
    cubeVer(3,:) = cube%center + [-cube%halfwidth, 0.0_rk, 0.0_rk]
    cubeVer(4,:) = cube%center + [0.0_rk,  cube%halfwidth, 0.0_rk]
    cubeVer(5,:) = cube%center + [0.0_rk, -cube%halfwidth, 0.0_rk]
    cubeVer(6,:) = cube%center + [0.0_rk, 0.0_rk,  cube%halfwidth]
    cubeVer(7,:) = cube%center + [0.0_rk, 0.0_rk, -cube%halfwidth]
    cubeVer(8,:) = cube%origin
    cubeVer(9,:) = [cube%endPnt(1), cube%origin(2), cube%origin(3)]
    cubeVer(10,:) = [cube%origin(1), cube%endPnt(2), cube%origin(3)]
    cubeVer(11,:) = [cube%origin(1), cube%origin(2), cube%endPnt(3)]
    cubeVer(12,:) = [cube%endPnt(1), cube%endPnt(2), cube%origin(3)]
    cubeVer(13,:) = [cube%endPnt(1), cube%origin(2), cube%endPnt(3)]
    cubeVer(14,:) = [cube%origin(1), cube%endPnt(2), cube%endPnt(3)]
    cubeVer(15,:) = cube%endPnt


    sphere%radius = cylinder%radius
    sphere%only_surface = cylinder%only_surface

    !! check whether cylinder axis intersect by cube
    !! if true then no need to check for intersection of each cube vertices
    !! on cylinder
    line%origin = cylinder%origin
    line%vec = cylinder%vec
    if(tem_lineCubeOverlap(line, cube, pntIntersect)) then
      ! project intersected point on cylinder to check for
      ! sphere cube overlap
      proj = dot_product((pntIntersect - cylinder%origin), &
        &                cylinder%vec) &
        &  / dot_product(cylinder%vec, cylinder%vec)
      ! compute the actual coordinate position of projected
      ! point on the line
      sphere%origin = cylinder%origin + proj*cylinder%vec
      overlap = tem_sphereCubeOverlap( sphere, cube )
      ! return if overlap is true
      if(overlap) return
    endif

    do iVer=1,15
      ! Find the projection of cubever on line and check
      ! whether projection point is between 0 and 1.
      ! If projection < 0 then the point is before the line
      ! if projection > 0 then the point is after the line
      proj = dot_product((cubeVer(iVer,:) - cylinder%origin), &
        &                cylinder%vec) &
        &  / dot_product(cylinder%vec, cylinder%vec)
      if( (proj .fge. 0.0_rk) .and. (proj .fle. 1.0_rk)) then
        ! compute the actual coordinate position of projected
        ! point on the line
        sphere%origin = cylinder%origin + proj*cylinder%vec
        overlap = tem_sphereCubeOverlap( sphere, cube )
        ! return if overlap is true
        if(overlap) return
      endif
    enddo

  end function tem_cylinderCubeOverlap
  ! ****************************************************************************

  ! ************************************************************************** !
  !> Write out an array of cylinders in lua format
  !!
  subroutine tem_cylinder_out_vec( me, conf )
    ! --------------------------------------------------------------------------
    !> cylinder types to write out
    type( tem_cylinder_type ), intent(in) :: me(:)
    !> Aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    ! --------------------------------------------------------------------------
    ! counter
    integer :: i
    ! --------------------------------------------------------------------------

    ! create a table with name cylinder
    call aot_out_open_table( put_conf = conf, tname = 'object' )

    do i = 1, size(me)
      call tem_cylinder_out_scal( me(i), conf )
    end do

    call aot_out_close_table( put_conf = conf )

  end subroutine tem_cylinder_out_vec
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Write out a cylinder shape in lua format
  !!
  subroutine tem_cylinder_out_scal( me, conf )
    ! --------------------------------------------------------------------------
    !> cylinder types to write out
    type( tem_cylinder_type ), intent(in) :: me
    !> Aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    ! --------------------------------------------------------------------------
    ! ---------------------------------------------------------------------------

    ! create a table with name cylinder if not exist
    if( conf%level == 0 ) then
      call aot_out_open_table( put_conf = conf, tname = 'object' )
    else
      call aot_out_open_table( put_conf = conf )
    end if

    call aot_out_val( put_conf = conf, vname = 'origin', val = me%origin )
    call aot_out_val( put_conf = conf, vname = 'vec', val = me%vec )
    call aot_out_val( put_conf = conf, vname = 'radius', val = me%radius )
    call aot_out_val( put_conf = conf, vname = 'only_surface', &
      &               val = me%only_surface )

    call aot_out_close_table( put_conf = conf )

  end subroutine tem_cylinder_out_scal
  ! ************************************************************************** !

  subroutine init_ga_cylinder(me, length)
    type(grw_cylinderarray_type), intent(out) :: me !< dynamic array to init
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

  end subroutine init_ga_cylinder

  subroutine destroy_ga_cylinder(me)
    type(grw_cylinderarray_type), intent(inout) :: me !< dynamic array to destroy

    me%containersize = 0
    me%nvals = 0
    if( allocated( me%val ) ) deallocate(me%val)

  end subroutine destroy_ga_cylinder


  subroutine truncate_ga_cylinder(me)
    !------------------------------------------------------------------------
    type(grw_cylinderarray_type) :: me !< array to truncate
    !------------------------------------------------------------------------
    type(tem_cylinder_type), allocatable :: tarray(:)
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

  end subroutine truncate_ga_cylinder


  subroutine empty_ga_cylinder(me)
    !------------------------------------------------------------------------
    type(grw_cylinderarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------

    me%nvals = 0

  end subroutine empty_ga_cylinder

  !> adds the value to a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! element at the requested position will be replaced.
  subroutine placeat_ga_cylinder(me, val, pos, length)
    type(grw_cylinderarray_type) :: me !< array to place the value into
    type(tem_cylinder_type), intent(in) :: val !< value to place at the given position
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

  end subroutine placeat_ga_cylinder


  !> adds the values starting from a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! elements starting from the requested position will be replaced up to
  !! the element at position `pos + size(val) - 1`.
  subroutine placeat_ga_cylinder_vec(me, val, pos, length)
    type(grw_cylinderarray_type) :: me !< array to append the value to
    type(tem_cylinder_type), intent(in) :: val(:) !< values to append
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

  end subroutine placeat_ga_cylinder_vec


  subroutine append_ga_cylinder(me, val, length)
    type(grw_cylinderarray_type) :: me !< array to append the value to
    type(tem_cylinder_type), intent(in) :: val !< value to append
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

  end subroutine append_ga_cylinder

  subroutine append_ga_cylinder_vec(me, val, length)
    type(grw_cylinderarray_type) :: me !< array to append the value to
    type(tem_cylinder_type), intent(in) :: val(:) !< values to append
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

  end subroutine append_ga_cylinder_vec


  subroutine expand_ga_cylinder(me, pos, length)
    type(grw_cylinderarray_type) :: me !< array to resize
    integer, intent(in), optional :: pos !< optional predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    type(tem_cylinder_type), allocatable :: swpval(:)
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

  end subroutine expand_ga_cylinder


end module tem_cylinder_module

!> \page cylinder Cylinders
!! Cylinders are defined by an origin, vector defining the length and the
!! axis and the radius.
!! Cylinder is considered to be solid as default i.e. all the cubes inside the
!! cylinder are marked as intersected cubes.
!! It is possible to created hollow cylinders by setting only_surface = true,
!! it will mark only the cubes intersect with cylinder surface as intersected
!! cubes
!!
!! Valid definition:
!! \li Single cylinder
!! \verbatim
!! geometry={
!!   kind='cylinder',
!!     object={
!!       origin={0.0,0.0,0.0},
!!       vec={1.0,0.0,0.0},
!!       radius=0.25,
!!       only_surface = true, -- If not defined default is set to false
!!     }
!! }
!! \endverbatim
!!
!! \li Multiple cylinder
!! \verbatim
!! geometry={
!!   kind='cylinder',
!!     object={
!!       {
!!       origin={0.0,0.0,0.0},
!!       vec={1.0,0.0,0.0},
!!       radius=0.25
!!       },
!!       {
!!       origin={0.0,0.0,0.0},
!!       vec={1.0,1.0,0.0},
!!       radius=0.25
!!       }
!!     }
!! }
!! \endverbatim
!! \n\n
!! Seeder file to generate mesh with single cylinder (only_surface=true) is below:
!! include testsuite/plane/seeder.lua
!! \n\n
!! Mesh with hollow cylinder (Hollow => only_surface = true)
!! \image html cylinder.png
!! \n\n
!! \image html cylinder_withedges.png
!! \n\n
!! Cutview of mesh with hollow cylinder
!! \image html cylinder_hollow.png
!! \n\n
!! As said earlier, cylinder can be created as solid one using 'only_surface=false'.
!! Cutview of Mesh generated with 'only_surface=false':
!! \image html cylinder_solid.png
!! \n\n
!! Example lua file is available at \link testsuite/cylinder/seeder.lua
!! \example testsuite/cylinder/seeder.lua

