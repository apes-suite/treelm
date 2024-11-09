! Copyright (c) 2015-2016, 2018-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2019 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
! Copyright (c) 2019, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!
! Copyright (c) 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
!
! Parts of this file were written by Peter Vitt for University of Siegen.
!
!> This module provides a spatial function to describe a 2D polygon.
!!
!! The spatial function will provide different values for points within the
!! polygon and those outside.
!!
!! A polygon is defined in the configuration by the following definition:
!!
!!```lua
!!   -- Vector to return inside the polygon
!!   inval = {1.0, 0.0, 2.0}
!!
!!   -- It also can be a scalar, if none is provided a scalar 1 is assumed.
!!   inval = 1.0 -- This is the value to return for points inside the polygon
!!               -- defaults to 1.
!!
!!   -- Vector to return outside the polygon
!!   outval = {0.0, 0.0, 0.0}
!!   -- Needs to conform to the definition of inval, that is, if inval is scalar
!!   -- or not given, outval needs to be a scalar.
!!   -- Defaults to all zero components for a vector of the length of inval.
!!   -- If inval is a scalar, outval has to be given as a scalar aswell:
!!   outval = 0.0 -- The value to return for points outside the polygon,
!!                -- defaults to 0.
!!
!!   -- List of 2D Points to be used as vertices for the polygon.
!!   -- The polygon will be closed by going from the last point back to the
!!   -- first one.
!!   vertex = { { 1.0,  0.0},
!!              { 0.0, -1.0},
!!              {-1.0,  0.0},
!!              { 0.0,  1.0},
!!            }
!!```
!!
module tem_polygon_material_module
  use env_module,         only: rk, labellen
  use aotus_module,       only: flu_state, aot_get_val, aoterr_Fatal
  use aot_table_module,   only: aot_table_open, aot_table_close, &
    &                           aot_table_length, aot_get_val
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logunit
  use tem_float_module,   only: operator(.feq.)

  implicit none

  private

  !> Type to store the vertices for the polygon
  type tem_polygon_vertex_type
    !> Number of vertices in the polygon.
    integer :: nVertices
    !> 2D Coordinates of the vertices.
    !! First index runs over number of vertices, second from 1 to 2.
    real(kind=rk), allocatable :: vertex(:,:)
  end type

  !> Type to store information regarding the movement of the
  !! polygon
  type tem_polygon_movement_type
    !> Linear movement of the polygon
    !! Include the values for the velocity, the first
    !! entry is the velocity in X direction, de second in Y
    !! The third one is the Z component of the velocity
    real(kind=rk), allocatable :: lin_parameter(:)
    !> Move the polygon with a sine fuction
    !! First and second entry belong to X direction and
    !! are the amplitude and the frequency. The third and
    !! forth entry are the devoted to the Y direction, for
    !! amplitude and the frequency respectivly.
    real(kind=rk), allocatable :: sin_parameter(:)
    !> Rotation of the polygon
    !! the first two entries belong to the directions, the first
    !! entry is the rot pointX and the second entry the rot pointY )
    !! the third one is the rot_speed omega
    real(kind=rk), allocatable :: rot_parameter(:)
    !> chaning the angle of attack for airfoild with a sinus
    !! the first entry is the phase shift (y), the second is
    !! the amplitude and the third one is the angular velocity (omega)
    real(kind=rk), allocatable :: angle_parameter(:)
    !> Kind of movement
    character(len=labellen) :: movement_kind
  end type tem_polygon_movement_type

  !> Description of a 2D closed polygon.
  type tem_polygon_material_type
    !> poly_list, we can have multiply of them
    type(tem_polygon_vertex_type),allocatable :: poly_list(:)
    !> Movement of each polygon
    type(tem_polygon_movement_type) :: moving
    !> Number of poly_list
    integer :: nPoly
    !> Extrude in z direction
    real(kind=rk) :: zmin
    real(kind=rk) :: zmax
    !> how many components inval/outval have,
    !> they are defined as vectors and might
    !> have more than 1 entries, defined by
    !> the user in the config file!
    integer :: nComponents
    !> Value of Material inside the polygon.
    real(kind=rk), allocatable :: inval(:)
    !> Value of Material outside the polygon.
    real(kind=rk), allocatable :: outval(:)
  end type tem_polygon_material_type

  public :: tem_polygon_material_type
  public :: tem_polygon_material_single_load
  public :: tem_polygon_material_movement_single
  public :: tem_polygon_material_movement_multi
  public :: tem_polygon_material_load
  public :: tem_polygon_material_multi_load
  public :: tem_polygon_material_value
  public :: tem_eval_polygon_material
  public :: tem_eval_polygon_material_3d
  public :: tem_eval_polygon_material_scal
  public :: tem_eval_polygon_material_scal_3d
  public :: tem_polygon_material_test_angle
  public :: tem_polygon_material_test_value

  ! Some private constants

  !> Definition of Pi
  real(kind=rk), parameter :: PI = 2*asin(1.0_rk)

  !> Overestimating tolerance factor for comparisons of reals
  real(kind=rk), parameter :: tolp = (1.0_rk + epsilon(1.0_rk))

  !> Underestimating tolerance factor for comparisons of reals
  real(kind=rk), parameter :: tolm = (1.0_rk - epsilon(1.0_rk))


contains

  subroutine tem_polygon_material_load(me, conf, thandle)
    ! ----------------------------------------------------------------------
    !> Polygon data structure to fill with information provided
    !! by the user in config.
    type(tem_polygon_material_type), intent(out) :: me

    !> Handle to the Lua script containing the polygon definition
    type(flu_state) :: conf

    !> Handle for the table containing the polygon definition.
    integer, intent(in), optional :: thandle
    ! ----------------------------------------------------------------------
    real(kind=rk), allocatable :: defout(:)
    integer :: vertex_table, valtable, vertices_table
    integer :: iVertex
    integer :: iError
    integer, allocatable :: vError(:)
    integer :: iError_v(2)
    integer :: iPoly
    ! ----------------------------------------------------------------------

    write(logUnit(1),*) 'Loading predefined function polygonal material:'
    valtable = 0
    call aot_get_val( L       = conf,    &
      &               thandle = thandle, &
      &               key     = 'zmin',  &
      &               val     = me%zmin, &
      &               default = 0.0_rk,  &
      &               ErrCode = iError   )

    if (btest(iError, aoterr_fatal)) then
      write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
      write(logunit(1),*) '      get a value for zmin!'
      write(logunit(1),*) '      This also applies if no inval is provided.'
      call tem_abort()
    end if
    call aot_get_val( L       = conf,    &
      &               thandle = thandle, &
      &               key     = 'zmax',  &
      &               val     = me%zmax, &
      &               default = 0.0_rk,  &
      &               ErrCode = iError   )

    if (btest(iError, aoterr_fatal)) then
      write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
      write(logunit(1),*) '      get a value for zmax!'
      write(logunit(1),*) '      This also applies if no inval is provided.'
      call tem_abort()
    end if
    call aot_table_open( L       = conf,    &
      &                  parent  = thandle, &
      &                  key     = 'inval', &
      &                  thandle = valtable )
    if (valtable  == 0) then
      ! inval not provided as a table, try to read it as a scalar.
      allocate(me%inval(1))
      call aot_get_val( L       = conf,        &
        &               thandle = thandle,     &
        &               key     = 'inval',     &
        &               val     = me%inval(1), &
        &               default = 1.0_rk,      &
        &               ErrCode = iError       )

      if (btest(iError, aoterr_Fatal)) then
        write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
        write(logunit(1),*) '      get value for inval!'
        call tem_abort()
      end if

      me%nComponents = 1
      allocate(me%outval(me%nComponents))
      allocate(vError(me%nComponents))

      ! Outval needs to be consistent with the inval definition, if inval was
      ! defined as a scalar, outval also has to be a scalar!
      ! We do not check for tables with single entries in this case.
      call aot_get_val( L       = conf,         &
        &               thandle = thandle,      &
        &               key     = 'outval',     &
        &               val     = me%outval(1), &
        &               default = 0.0_rk,       &
        &               ErrCode = iError        )

      if (btest(iError, aoterr_fatal)) then
        write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
        write(logunit(1),*) '      get a value for outval!'
        write(logunit(1),*) '      Note, that outval needs be a scalar, as'
        write(logunit(1),*) '      inval is provided as a scalar.'
        write(logunit(1),*) '      This also applies if no inval is provided.'
        call tem_abort()
      end if

    else

      ! Intable is a table, close it an read it into an array.
      call aot_table_close(L = conf, thandle = valtable)

      ! Value to use inside the polygon
      call aot_get_val( L         = conf,     &
        &               thandle   = thandle,  &
        &               key       = 'inval',  &
        &               val       = me%inval, &
        &               maxlength = 20,       &
        &               default   = [1.0_rk], &
        &               ErrCode   = vError    )

      if (any(btest(vError, aoterr_Fatal))) then
        write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
        write(logunit(1),*) '      get values for inval!'
        call tem_abort()
      end if
      me%nComponents = size(me%inval)
      deallocate(vError)

      ! Definition of outval needs to be consistent with inval, it has to have
      ! the same number of components, and also needs to be a vector.
      ! However, we define a default of all zeroes, so if outval is 0 for all
      ! components, this definition can be omitted in the user definition.
      allocate(me%outval(me%nComponents))
      allocate(vError(me%nComponents))
      allocate(defout(me%nComponents))

      defout = 0.0_rk

      ! Value to use outside the polygon
      call aot_get_val( L       = conf,      &
        &               thandle = thandle,   &
        &               key     = 'outval',  &
        &               val     = me%outval, &
        &               default = defout,    &
        &               ErrCode = vError     )

      if (any(btest(vError,aoterr_Fatal))) then
        write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
        write(logunit(1),*) '      get a value for outval!'
        write(logunit(1),*) '      Note, that outval needs to have the same'
        write(logunit(1),*) '      length as inval.'
        call tem_abort()
      end if

      deallocate(vError)
      deallocate(defout)

    end if

    !> read list of vertices
    call aot_table_open( L       = conf,          &
      &                  parent  = thandle,       &
      &                  key     = 'vertices',    &
      &                  thandle = vertices_table )

    if (vertices_table == 0) then
      write(logunit(1),*) 'ERROR in tem_polygon_material_load: No vertices'
      write(logunit(1),*) '      defined, unable to set up a polynomial!'
      write(logunit(1),*) '      Please define vertices tables with lists of'
      write(logunit(1),*) '      2D points.'
      call tem_abort()
    end if
    me%nPoly = aot_table_length( L       = conf,          &
      &                          thandle = vertices_table )
    
    allocate(me%poly_list(me%npoly))
   
   ! me%npoly = 1
    do iPoly = 1, me%npoly
      call aot_table_open( L       = conf,           &
        &                  parent  = vertices_table, &
        &                  pos     = ipoly,          &
        &                  thandle = vertex_table    )
      if (vertex_table == 0) then
        write(logunit(1),*) 'ERROR in tem_polygon_material_load:'
        write(logunit(1),*) '      Please define vertex tables with lists of'
        write(logunit(1),*) '      2D points.'
        call tem_abort()
      end if
      me%poly_list(ipoly)%nVertices =              &
        & aot_table_length( L       = conf,        &
        &                   thandle = vertex_table )

      allocate(me%poly_list(ipoly)%vertex(me%poly_list(ipoly)%nVertices, 2))
      do iVertex=1,me%poly_list(ipoly)%nVertices
        call aot_get_val( val     = me%poly_list(ipoly)%vertex(iVertex, :), &
          &               ErrCode = iError_v,                                  &
          &               L       = conf,                                      &
          &               thandle = vertex_table,                              &
          &               pos     = iVertex                                    )
        if (any(iError_v /= 0)) then
          write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
          write(logunit(1),*) '      obtain vertex ', iVertex, '!'
          write(logunit(1),*) '      Vertices have to be vectors of length 2,'
          write(logunit(1),*) '      with real numbers as entries.'
          call tem_abort()
        end if
      end do
      call aot_table_close( L       = conf,        &
        &                   thandle = vertex_table )

    end do
    call aot_table_close( L       = conf,          &
      &                   thandle = vertices_table )

  end subroutine tem_polygon_material_load

  subroutine tem_polygon_material_single_load(me, conf, thandle)
    ! ----------------------------------------------------------------------
    !> Polygon data structure to fill with information provided
    !! by the user in config.
    type(tem_polygon_material_type), intent(out) :: me

    !> Handle to the Lua script containing the polygon definition
    type(flu_state) :: conf

    !> Handle for the table containing the polygon definition.
    integer, intent(in), optional :: thandle
    ! ----------------------------------------------------------------------
    real(kind=rk), allocatable :: defout(:)
    integer :: vertex_table, valtable, vertices_table
    integer :: movement_table, parameter_table
    integer :: iVertex
    integer :: iError
    integer, allocatable :: vError(:)
    integer :: iError_v(2)
    ! ----------------------------------------------------------------------
    write(logUnit(1),*) 'Loading predefined function for single polygon:'
    valtable = 0
    !> get the z component
    call aot_get_val( L       = conf,    &
      &               thandle = thandle, &
      &               key     = 'zmin',  &
      &               val     = me%zmin, &
      &               default = 0.0_rk,  &
      &               ErrCode = iError   )

    if (btest(iError, aoterr_fatal)) then
      write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
      write(logunit(1),*) '      get a value for zmin!'
      write(logunit(1),*) '      This also applies if no inval is provided.'
      call tem_abort()
    end if
    call aot_get_val( L       = conf,    &
      &               thandle = thandle, &
      &               key     = 'zmax',  &
      &               val     = me%zmax, &
      &               default = 0.0_rk,  &
      &               ErrCode = iError   )
    if (btest(iError, aoterr_fatal)) then
      write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
      write(logunit(1),*) '      get a value for zmax!'
      write(logunit(1),*) '      This also applies if no inval is provided.'
      call tem_abort()
    end if
    !> read list of vertices
    call aot_table_open( L       = conf,          &
      &                  parent  = thandle,       &
      &                  key     = 'vertices',    &
      &                  thandle = vertices_table )
    if (vertices_table == 0) then
      write(logunit(1),*) 'ERROR in tem_polygon_material_load: No vertices'
      write(logunit(1),*) '      defined, unable to set up a polynomial!'
      write(logunit(1),*) '      Please define vertices tables with lists of'
      write(logunit(1),*) '      2D points.'
      call tem_abort()
    end if
    !> for a single polygon me%npoly = 1
    call aot_table_open( L       = conf,           &
      &                  parent  = vertices_table, &
      &                  pos     = 1,              &
      &                  thandle = vertex_table    )
    if (vertex_table == 0) then
      write(logunit(1),*) 'ERROR in tem_polygon_material_load: No vertices'
      write(logunit(1),*) '      Please define vertex tables with lists of'
      write(logunit(1),*) '      2D points.'
      call tem_abort()
    end if
    me%npoly = 1
    allocate(me%poly_list(me%npoly))

    me%poly_list(me%npoly)%nVertices =           &
      & aot_table_length( L       = conf,        &
      &                   thandle = vertex_table )

    allocate(me%poly_list(me%npoly)%vertex(me%poly_list(me%npoly)%nVertices, 2))
    do iVertex=1,me%poly_list(me%npoly)%nVertices
      call aot_get_val( val     = me%poly_list(me%npoly)%vertex(iVertex, :), &
        &               ErrCode = iError_v,                                  &
        &               L       = conf,                                      &
        &               thandle = vertex_table,                              &
        &               pos     = iVertex                                    )
      if (any(iError_v /= 0)) then
        write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
        write(logunit(1),*) '      obtain vertex ', iVertex, '!'
        write(logunit(1),*) '      Vertices have to be vectors of length 2,'
        write(logunit(1),*) '      with real numbers as entries.'
        call tem_abort()
      end if
    end do
    call aot_table_close( L       = conf,        &
      &                   thandle = vertex_table )


    call aot_table_close( L       = conf,          &
      &                   thandle = vertices_table )

    call aot_table_open( L       = conf,          &
      &                  parent  = thandle,       &
      &                  key     = 'movement',    &
      &                  thandle = movement_table )
    if (movement_table == 0) then
      write(logunit(1),*) 'ERROR in tem_polygon_material_load: No movement'
      write(logunit(1),*) '      defined, unable to continue!'
      write(logunit(1),*)                                         &
        & '      Please define movement table with a movement kind'
      call tem_abort()
    end if
    call aot_get_val( L         = conf,                    &
      &               thandle   = movement_table,          &
      &               key       = 'movement_kind',         &
      &               val       = me%moving%movement_kind, &
      &               default   = 'NO Movement',           &
      &               ErrCode   = iError                   )
    call aot_table_close( L       = conf,          &
      &                   thandle = movement_table )


    select case(me%moving%movement_kind)
      case( 'lin_movement_2d', 'lin_movement_3d')
        call aot_table_open( L       = conf,            &
          &                  parent  = thandle,         &
          &                  key     = 'lin_parameter', &
          &                  thandle = parameter_table  )
        if (parameter_table == 0) then
          write(logunit(1),*)                                   &
            & 'ERROR in tem_polygon_material_load: lin_parameter'
          write(logunit(1),*) '      not defined, unable to set movement!'
          write(logunit(1),*) '      Please define lin_parameter table'
          call tem_abort()
        end if
        me%nComponents = aot_table_length( L       = conf,           &
              &                            thandle = parameter_table )

        call aot_table_close( L       = conf,           &
          &                   thandle = parameter_table )

        allocate(vError(me%nComponents))
        allocate(me%moving%lin_parameter(me%nComponents))
        call aot_get_val( L         = conf,                    &
          &               thandle   = thandle,                 &
          &               key       = 'lin_parameter',         &
          &               val       = me%moving%lin_parameter, &
          &               maxlength = me%nComponents,          &
          &               ErrCode   = VError                   )

        if (any(btest(vError, aoterr_Fatal))) then
          write(logunit(1),*)                                 &
            & 'ERROR in tem_polygon_material_load: No movement'
          write(logunit(1),*)                                      &
            & '      values defined, unable to set up the movement!'
          write(logunit(1),*) '      Please define lin_parameter!'
          write(logunit(1),*)                                       &
            & '      For linear movement set velocity in X,Y for 2D!'
          call tem_abort()
        end if
        deallocate(vError)

      case( 'sin_movement_2d', 'sin_movement_3d')
        call aot_table_open( L       = conf,            &
          &                  parent  = thandle,         &
          &                  key     = 'sin_parameter', &
          &                  thandle = parameter_table  )
        if (parameter_table == 0) then
          write(logunit(1),*)                                   &
            & 'ERROR in tem_polygon_material_load: sin_parameter'
          write(logunit(1),*) '      not defined, unable to set movement!'
          write(logunit(1),*) '      Please define sin_parameter table'
          call tem_abort()
        end if
        me%nComponents = aot_table_length( L       = conf,           &
              &                            thandle = parameter_table )

        call aot_table_close( L       = conf,        &
          &                   thandle = parameter_table )

        allocate(vError(me%nComponents))
        allocate(me%moving%sin_parameter(me%nComponents))
        call aot_get_val( L         = conf,                    &
          &               thandle   = thandle,                 &
          &               key       = 'sin_parameter',         &
          &               val       = me%moving%sin_parameter, &
          &               maxlength = me%nComponents,          &
          &               ErrCode   = VError                   )

        if (any(btest(vError, aoterr_Fatal))) then
          write(logunit(1),*)                                 &
            & 'ERROR in tem_polygon_material_load: No movement'
          write(logunit(1),*)                                     &
            &'      values defined, unable to set up the movement!'
          write(logunit(1),*) '      Please define sin_parameter!'
          write(logunit(1),*)                                         &
            & '      For sinefuction set amplitude and frequency in X,'
          write(logunit(1),*) '      than the information for the Y direction!'
        end if
        deallocate(vError)

      case( 'angleofAttack_2d', 'angleofAttack_3d')
        call aot_table_open( L       = conf,            &
          &                  parent  = thandle,         &
          &                  key     = 'angle_parameter', &
          &                  thandle = parameter_table  )
        if (parameter_table == 0) then
          write(logunit(1),*)                                   &
            & 'ERROR in tem_polygon_material_load: angle_parameter'
          write(logunit(1),*) '      not defined, unable to set movement!'
          write(logunit(1),*) '      Please define sin_parameter table'
          call tem_abort()
        end if
        me%nComponents = aot_table_length( L       = conf,           &
              &                            thandle = parameter_table )

        call aot_table_close( L       = conf,        &
          &                   thandle = parameter_table )
        allocate(vError(me%nComponents))
        allocate(me%moving%sin_parameter(me%nComponents))
        call aot_get_val( L         = conf,                    &
          &               thandle   = thandle,                 &
          &               key       = 'angle_parameter',         &
          &               val       = me%moving%angle_parameter, &
          &               maxlength = me%nComponents,          &
          &               ErrCode   = VError                   )

        if (any(btest(vError, aoterr_Fatal))) then
          write(logunit(1),*) 'ERROR in tem_polygon_material_load: No' &
            & // ' movement for angle of attack'
          write(logunit(1),*)                                     &
            & '      values defined, unable to set up the movement!'
          write(logunit(1),*) '      Please define sin_parameter!'
          write(logunit(1),*)                                         &
            & '      For sinefunction set amplitude and frequency in X,'
          write(logunit(1),*) '      than the information for the Y direction!'
        end if
        deallocate(vError)
      case( 'rot_movement_2d', 'rot_movement_3d')
        call aot_table_open( L       = conf,            &
          &                  parent  = thandle,         &
          &                  key     = 'rot_parameter', &
          &                  thandle = parameter_table  )
        if (parameter_table == 0) then
          write(logunit(1),*)                                   &
            & 'ERROR in tem_polygon_material_load: rot_parameter'
          write(logunit(1),*) '      not defined, unable to set movement!'
          write(logunit(1),*) '      Please define rot_parameter table'
          call tem_abort()
        end if
        me%nComponents = aot_table_length( L       = conf,           &
          &                                thandle = parameter_table )

        call aot_table_close( L       = conf,        &
          &                   thandle = parameter_table )

        allocate(me%moving%sin_parameter(me%nComponents))
        allocate(vError(me%nComponents))
        call aot_get_val( L         = conf,                    &
          &               thandle   = thandle,                 &
          &               key       = 'rot_parameter',         &
          &               val       = me%moving%rot_parameter, &
          &               maxlength = me%nComponents,          &
          &               ErrCode   = VError                   )

        if (any(btest(vError, aoterr_Fatal))) then
          write(logunit(1),*) 'ERROR in tem_polygon_material_load: No movement'
          write(logunit(1),*)                                      &
            & '      values defined, unable to set up the movement!'
          write(logunit(1),*) '      Please define rot_parameter!'
          write(logunit(1),*)                                         &
            & '      For rotation provide the X and Y direction of the'
          write(logunit(1),*) '      rotation_point and the rotation speed! '
        end if
        deallocate(vError)
      case default
        write(logUnit(1),*) 'NO movement prescribed for the polygon'
    end select
    write(logUnit(1),*) 'The single geometry movement is defined by: ', &
      & me%moving%movement_kind

    call aot_table_open( L       = conf,    &
      &                  parent  = thandle, &
      &                  key     = 'inval', &
      &                  thandle = valtable )

    if (valtable == 0) then
      ! inval not provided as a table, try to read it as a scalar.
      allocate(me%inval(1))

      call aot_get_val( L         = conf,        &
        &               thandle   = thandle,     &
        &               key       = 'inval',     &
        &               val       = me%inval(1), &
        &               default   = 1.0_rk,      &
        &               ErrCode   = iError       )

      if (btest(iError, aoterr_Fatal)) then
        write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
        write(logunit(1),*) '      get value for inval!'
        call tem_abort()
      end if

      me%nComponents = 1
      allocate(me%outval(me%nComponents))
      allocate(vError(me%nComponents))

      ! Outval needs to be consistent with the inval definition, if inval was
      ! defined as a scalar, outval also has to be a scalar!
      ! We do not check for tables with single entries in this case.
      call aot_get_val( L       = conf,         &
        &               thandle = thandle,      &
        &               key     = 'outval',     &
        &               val     = me%outval(1), &
        &               default = 0.0_rk,       &
        &               ErrCode = iError        )

      if (btest(iError, aoterr_fatal)) then
        write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
        write(logunit(1),*) '      get a value for outval!'
        write(logunit(1),*) '      Note, that outval needs be a scalar, as'
        write(logunit(1),*) '      inval is provided as a scalar.'
        write(logunit(1),*) '      This also applies if no inval is provided.'
        call tem_abort()
      end if

    else

      ! Intable is a table, close it an read it into an array.
      call aot_table_close(L = conf, thandle = valtable)

      ! Value to use inside the polygon
      call aot_get_val( L         = conf,     &
        &               thandle   = thandle,  &
        &               key       = 'inval',  &
        &               val       = me%inval, &
        &               maxlength = 20,       &
        &               default   = [1.0_rk], &
        &               ErrCode   = vError    )

      if (any(btest(vError, aoterr_Fatal))) then
        write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
        write(logunit(1),*) '      get values for inval!'
        call tem_abort()
      end if
      me%nComponents = size(me%inval)
      deallocate(vError)

      ! Definition of outval needs to be consistent with inval, it has to have
      ! the same number of components, and also needs to be a vector.
      ! However, we define a default of all zeroes, so if outval is 0 for all
      ! components, this definition can be omitted in the user definition.
      allocate(me%outval(me%nComponents))
      allocate(vError(me%nComponents))
      allocate(defout(me%nComponents))

      defout = 0.0_rk

      ! Value to use outside the polygon
      call aot_get_val( L       = conf,      &
        &               thandle = thandle,   &
        &               key     = 'outval',  &
        &               val     = me%outval, &
        &               default = defout,    &
        &               ErrCode = vError     )

      if (any(btest(vError,aoterr_Fatal))) then
        write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
        write(logunit(1),*) '      get a value for outval!'
        write(logunit(1),*) '      Note, that outval needs to have the same'
        write(logunit(1),*) '      length as inval.'
        call tem_abort()
      end if

      deallocate(vError)
      deallocate(defout)

    end if
  end subroutine tem_polygon_material_single_load

  subroutine tem_polygon_material_multi_load(me, conf, thandle)
    ! ----------------------------------------------------------------------
    !> Polygon data structure to fill with information provided
    !! by the user in config.
    type(tem_polygon_material_type), intent(out) :: me

    !> Handle to the Lua script containing the polygon definition
    type(flu_state) :: conf

    !> Handle for the table containing the polygon definition.
    integer, intent(in), optional :: thandle
    ! ----------------------------------------------------------------------
    real(kind=rk), allocatable :: defout(:)
    integer :: vertex_table, valtable, vertices_table
    integer :: movement_table, parameter_table
    integer :: iVertex, ipoly
    integer :: iError
    integer, allocatable :: vError(:)
    integer :: iError_v(2)
    ! ----------------------------------------------------------------------

    write(logUnit(1),*) 'Loading predefined function for multi body polygon:'

    valtable =0
    !> get the z component
    call aot_get_val( L       = conf,    &
      &               thandle = thandle, &
      &               key     = 'zmin',  &
      &               val     = me%zmin, &
      &               default = 0.0_rk,  &
      &               ErrCode = iError   )

    if (btest(iError, aoterr_fatal)) then
      write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
      write(logunit(1),*) '      get a value for zmin!'
      write(logunit(1),*) '      This also applies if no inval is provided.'
      call tem_abort()
    end if
    call aot_get_val( L       = conf,    &
      &               thandle = thandle, &
      &               key     = 'zmax',  &
      &               val     = me%zmax, &
      &               default = 0.0_rk,  &
      &               ErrCode = iError   )
    if (btest(iError, aoterr_fatal)) then
      write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
      write(logunit(1),*) '      get a value for zmax!'
      write(logunit(1),*) '      This also applies if no inval is provided.'
      call tem_abort()
    end if
    !> read list of vertices
    call aot_table_open( L       = conf,          &
      &                  parent  = thandle,       &
      &                  key     = 'vertices',    &
      &                  thandle = vertices_table )

    if (vertices_table == 0) then
      write(logunit(1),*) 'ERROR in tem_polygon_material_load: No vertices'
      write(logunit(1),*) '      defined, unable to set up a polynomial!'
      write(logunit(1),*) '      Please define vertices tables with lists of'
      write(logunit(1),*) '      2D points.'
      call tem_abort()
    end if

     me%nPoly = aot_table_length( L       = conf,          &
       &                          thandle = vertices_table )
    !> how many lists provided, here we just have one
    allocate(me%poly_list(me%npoly))
    do iPoly = 1, me%npoly
      call aot_table_open( L       = conf,           &
        &                  parent  = vertices_table, &
        &                  pos     = ipoly,          &
        &                  thandle = vertex_table    )
      if (vertex_table == 0) then
        write(logunit(1),*) 'ERROR in tem_polygon_material_load: No vertices'
        write(logunit(1),*) '      Please define vertex tables with lists of'
        write(logunit(1),*) '      2D points.'
        call tem_abort()
      end if
      me%poly_list(ipoly)%nVertices = aot_table_length( L       = conf,        &
        &                                               thandle = vertex_table )

      allocate(me%poly_list(ipoly)%vertex(me%poly_list(ipoly)%nVertices, 2))
      do iVertex=1,me%poly_list(ipoly)%nVertices
        call aot_get_val( val     = me%poly_list(ipoly)%vertex(iVertex, :), &
          &               ErrCode = iError_v,                               &
          &               L       = conf,                                   &
          &               thandle = vertex_table,                           &
          &               pos     = iVertex                                 )
        if (any(iError_v /= 0)) then
          write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
          write(logunit(1),*) '      obtain vertex ', iVertex, '!'
          write(logunit(1),*) '      Vertices have to be vectors of length 2,'
          write(logunit(1),*) '      with real numbers as entries.'
          call tem_abort()
        end if
      end do
      call aot_table_close( L       = conf,        &
        &                   thandle = vertex_table )
    end do

    call aot_table_close( L       = conf,          &
      &                   thandle = vertices_table )

    call aot_table_open( L       = conf,          &
      &                  parent  = thandle,       &
      &                  key     = 'movement',    &
      &                  thandle = movement_table )
    if (movement_table == 0) then
      write(logunit(1),*) 'ERROR in tem_polygon_material_load: No movement'
      write(logunit(1),*) '      defined, unable to continue!'
      write(logunit(1),*)                                         &
        & '      Please define movement table with a movement kind'
      call tem_abort()
    end if

    call aot_get_val( L         = conf,                    &
      &               thandle   = movement_table,          &
      &               key       = 'movement_kind',         &
      &               val       = me%moving%movement_kind, &
      &               default   = 'No movement',           &
      &               ErrCode   = iError                   )

    call aot_table_close( L       = conf,          &
      &                   thandle = movement_table )

    select case(me%moving%movement_kind)
      case( 'lin_multi_body_2d', 'lin_multi_body_3d')
        call aot_table_open( L       = conf,            &
          &                  parent  = thandle,         &
          &                  key     = 'lin_parameter', &
          &                  thandle = parameter_table  )
        if (parameter_table == 0) then
          write(logunit(1),*)                                   &
            & 'ERROR in tem_polygon_material_load: lin_parameter'
          write(logunit(1),*) '      not defined, unable to set movement!'
          write(logunit(1),*) '      Define lin_parameter table for multi body'
          call tem_abort()
        end if
        me%nComponents = aot_table_length( L       = conf,           &
          &                                thandle = parameter_table )

        call aot_table_close( L       = conf,           &
          &                   thandle = parameter_table )

        allocate(me%moving%lin_parameter(me%nComponents))
        call aot_get_val( L         = conf,                    &
          &               thandle   = thandle,                 &
          &               key       = 'lin_parameter',         &
          &               val       = me%moving%lin_parameter, &
          &               maxlength = me%nComponents,          &
          &               ErrCode   = VError                   )

        if (any(btest(vError, aoterr_Fatal))) then
          write(logunit(1),*)                                 &
            & 'ERROR in tem_polygon_material_load: No movement'
          write(logunit(1),*)                                      &
            & '      values defined, unable to set up the movement!'
          write(logunit(1),*)                    &
            & '      Please define lin_parameter!'
          write(logunit(1),*)                                       &
            & '      For linear movement set velocity in X,Y for 2D!'
          call tem_abort()
        end if
      case( 'sin_multi_body_2d','sin_multi_body_3d')
        call aot_table_open( L       = conf,            &
          &                  parent  = thandle,         &
          &                  key     = 'sin_parameter', &
          &                  thandle = parameter_table  )
        if (parameter_table == 0) then
          write(logunit(1),*)                                   &
            & 'ERROR in tem_polygon_material_load: sin_parameter'
          write(logunit(1),*) '      not defined, unable to set movement!'
          write(logunit(1),*)                                        &
            & '      Please define sin_parameter table for multi body'
          call tem_abort()
        end if
        me%nComponents = aot_table_length( L       = conf,           &
          &                                thandle = parameter_table )

        call aot_table_close( L       = conf,        &
          &                   thandle = parameter_table )

        allocate(me%moving%sin_parameter(me%nComponents))
        call aot_get_val( L         = conf,                    &
          &               thandle   = thandle,                 &
          &               key       = 'sin_parameter',         &
          &               val       = me%moving%sin_parameter, &
          &               maxlength = me%nComponents,          &
          &               ErrCode   = VError                   )

        if (any(btest(vError, aoterr_Fatal))) then
          write(logunit(1),*) 'ERROR in tem_polygon_material_load: No movement'
          write(logunit(1),*)                                      &
            & '      values defined, unable to set up the movement!'
          write(logunit(1),*) '      Please define sin_parameter!'
          write(logunit(1),*)                                         &
            & '      For sinefuction set amplitude and frequency in X,'
          write(logunit(1),*) '      than the information for the Y direction!'
        end if

      case( 'rot_multi_body_2d','rot_multi_body_3d')
        call aot_table_open( L       = conf,            &
          &                  parent  = thandle,         &
          &                  key     = 'rot_parameter', &
          &                  thandle = parameter_table  )
        if (parameter_table == 0) then
          write(logunit(1),*)                                   &
            & 'ERROR in tem_polygon_material_load: rot_parameter'
          write(logunit(1),*) '      not defined, unable to set movement!'
          write(logunit(1),*) '      Please define rot_parameter table'
          call tem_abort()
        end if
        me%nComponents = aot_table_length( L       = conf,           &
              &                            thandle = parameter_table )

        call aot_table_close( L       = conf,        &
          &                   thandle = parameter_table )

        allocate(me%moving%rot_parameter(me%nComponents))
        call aot_get_val( L         = conf,                    &
          &               thandle   = thandle,                 &
          &               key       = 'rot_parameter',         &
          &               val       = me%moving%rot_parameter, &
          &               maxlength = me%nComponents,          &
          &               ErrCode   = VError                   )

        if (any(btest(vError, aoterr_Fatal))) then
          write(logunit(1),*) 'ERROR in tem_polygon_material_load: No movement'
          write(logunit(1),*)                                      &
            & '      values defined, unable to set up the movement!'
          write(logunit(1),*) '      Please define rot_parameter!'
          write(logunit(1),*)                                         &
            & '      For rotation provide the X and Y direction of the'
          write(logunit(1),*) '      rotation_point and the rotation speed! '
        end if
      case default
        write(logUnit(1),*) 'NO movement prescribed for the polygon'
      end select
      deallocate(vError)
      write(logUnit(1),*) 'The multi geometry movement is defined by: ', &
        & me%moving%movement_kind

    !! here we need to read out the kind and than point to the specific
    !! movement of the polygon
    call aot_table_open( L       = conf,    &
      &                  parent  = thandle, &
      &                  key     = 'inval', &
      &                  thandle = valtable )

    if (valtable  == 0) then
      ! inval not provided as a table, try to read it as a scalar.
      allocate(me%inval(1))
      call aot_get_val( L         = conf,        &
        &               thandle   = thandle,     &
        &               key       = 'inval',     &
        &               val       = me%inval(1), &
        &               default   = 1.0_rk,      &
        &               ErrCode   = iError       )

      if (btest(iError, aoterr_Fatal)) then
        write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
        write(logunit(1),*) '      get value for inval!'
        call tem_abort()
      end if

      me%nComponents = 1
      allocate(me%outval(me%nComponents))

      ! Outval needs to be consistent with the inval definition, if inval was
      ! defined as a scalar, outval also has to be a scalar!
      ! We do not check for tables with single entries in this case.
      call aot_get_val( L       = conf,         &
        &               thandle = thandle,      &
        &               key     = 'outval',     &
        &               val     = me%outval(1), &
        &               default = 0.0_rk,       &
        &               ErrCode = iError        )

      if (btest(iError, aoterr_fatal)) then
        write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
        write(logunit(1),*) '      get a value for outval!'
        write(logunit(1),*) '      Note, that outval needs be a scalar, as'
        write(logunit(1),*) '      inval is provided as a scalar.'
        write(logunit(1),*) '      This also applies if no inval is provided.'
        call tem_abort()
      end if

    else

      ! Intable is a table, close it an read it into an array.
      call aot_table_close(L = conf, thandle = valtable)

      ! Value to use inside the polygon
      call aot_get_val( L         = conf,     &
        &               thandle   = thandle,  &
        &               key       = 'inval',  &
        &               val       = me%inval, &
        &               maxlength = 20,       &
        &               default   = [1.0_rk], &
        &               ErrCode   = vError    )

      if (any(btest(vError, aoterr_Fatal))) then
        write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
        write(logunit(1),*) '      get values for inval!'
        call tem_abort()
      end if
      me%nComponents = size(me%inval)
      deallocate(vError)

      ! Definition of outval needs to be consistent with inval, it has to have
      ! the same number of components, and also needs to be a vector.
      ! However, we define a default of all zeroes, so if outval is 0 for all
      ! components, this definition can be omitted in the user definition.
      allocate(me%outval(me%nComponents))
      allocate(vError(me%nComponents))
      allocate(defout(me%nComponents))

      defout = 0.0_rk

      ! Value to use outside the polygon
      call aot_get_val( L       = conf,      &
        &               thandle = thandle,   &
        &               key     = 'outval',  &
        &               val     = me%outval, &
        &               default = defout,    &
        &               ErrCode = vError     )

      if (any(btest(vError,aoterr_Fatal))) then
        write(logunit(1),*) 'ERROR in tem_polygon_material_load: Not able to'
        write(logunit(1),*) '      get a value for outval!'
        write(logunit(1),*) '      Note, that outval needs to have the same'
        write(logunit(1),*) '      length as inval.'
        call tem_abort()
      end if

      deallocate(vError)
      deallocate(defout)

    end if
  end subroutine tem_polygon_material_multi_load


  !> Evaluate a list of points, and return inval for each that is within me
  !! and outval for all other points.
  function tem_eval_polygon_material(me, coord, n) result(res)
    ! ----------------------------------------------------------------------
    !> Description of the polygon to evaluate
    type(tem_polygon_material_type), intent(in) :: me

    !> Number of points to get a value for.
    integer, intent(in) :: n

    !> Coordinates for which the function should be evaluated.
    real(kind=rk), intent(in) :: coord(n,3)

    !> Resulting value at each point.
    real(kind=rk) :: res(n, me%nComponents)
    ! ----------------------------------------------------------------------
    integer :: iPoint, iPoly, icomp
    ! ----------------------------------------------------------------------
    do icomp = 1, me%nComponents
      res(:,iComp) = me%outval(iComp)
    end do

    do iPoly = 1, me%nPoly
      do iPoint=1,n
        if (res(iPoint,1) .feq. me%outval(1) ) then
          res(iPoint,:) = tem_polygon_material_value(  &
            & me          = me%poly_list(ipoly),       &
            & nComponents = me%nComponents,            &
            & inVal       = me%inVal,                  &
            & outVal      = me%outVal,                 &
            & point       = coord(iPoint,:2)           )
         endif   
      end do
    end do
  end function tem_eval_polygon_material


  !> Evaluate a list of points, and return inval for each that is within me
  !! and outval for all other points.
  function tem_eval_polygon_material_3d(me, coord, n) result(res)
    ! ----------------------------------------------------------------------
    !> Description of the polygon to evaluate
    type(tem_polygon_material_type), intent(in) :: me

    !> Number of points to get a value for.
    integer, intent(in) :: n

    !> Coordinates for which the function should be evaluated.
    real(kind=rk), intent(in) :: coord(n,3)

    !> Resulting value at each point.
    real(kind=rk) :: res(n, me%nComponents)
    ! ----------------------------------------------------------------------
    integer :: iPoint, iPoly
    ! ----------------------------------------------------------------------
    do iPoly = 1, me%nPoly
      do iPoint=1,n
        if (coord(iPoint,3) >= me%zmin .and. coord(iPoint,3) <= me%zmax ) then
          res(iPoint,:) = tem_polygon_material_value( &
            & me          = me%poly_list(iPoly),      &
            & nComponents = me%nComponents,           &
            & inVal       = me%inVal,                 &
            & outVal      = me%outVal,                &
            & point       = coord(iPoint,:2)          )
        else
          res(iPoint,:) = me%outval
        end if
      end do
    end do
  end function tem_eval_polygon_material_3d


  !> Evaluate a list of points, and return first component of inval for
  !! each that is within me and first component of outval for all
  !! other points.
  function tem_eval_polygon_material_scal(me, coord, n) result(res)
    ! ----------------------------------------------------------------------
    !> Description of the polygon to evaluate
    type(tem_polygon_material_type), intent(in) :: me

    !> Number of points to get a value for.
    integer, intent(in) :: n

    !> Coordinates for which the function should be evaluated.
    real(kind=rk), intent(in) :: coord(n,3)

    !> Resulting value at each point.
    real(kind=rk) :: res(n)
    ! ----------------------------------------------------------------------
    real(kind=rk) :: loc(me%nComponents)
    integer :: iPoint, iPoly
    ! ----------------------------------------------------------------------
    do iPoly = 1, me%nPoly
      do iPoint=1,n
        loc = tem_polygon_material_value(      &
          & me          = me%poly_list(ipoly), &
          & nComponents = me%nComponents,      &
          & inVal       = me%inVal,            &
          & outVal      = me%outVal,           &
          & point       = coord(iPoint,:2)     )
        res(iPoint) = loc(1)
      end do
    end do 
  end function tem_eval_polygon_material_scal


  !> Evaluate a list of points, and return first component of inval for
  !! each that is within me and first component of outval for all
  !! other points.
  function tem_eval_polygon_material_scal_3d(me, coord, n) result(res)
    ! ----------------------------------------------------------------------
    !> Description of the polygon to evaluate
    type(tem_polygon_material_type), intent(in) :: me

    !> Number of points to get a value for.
    integer, intent(in) :: n

    !> Coordinates for which the function should be evaluated.
    real(kind=rk), intent(in) :: coord(n,3)

    !> Resulting value at each point.
    real(kind=rk) :: res(n)
    ! ----------------------------------------------------------------------
    real(kind=rk) :: loc(me%nComponents)
    integer :: iPoint, iPoly
    ! ----------------------------------------------------------------------
    do iPoly = 1, me%nPoly
      do iPoint=1,n
        if (coord(iPoint,3) >= me%zmin .and. coord(iPoint,3) <= me%zmax ) then
          loc = tem_polygon_material_value(      & 
            & me          = me%poly_list(ipoly), &
            & nComponents = me%nComponents,      &
            & inVal       = me%inVal,            &
            & outVal      = me%outVal,           &
            & point       = coord(iPoint,:2)     )
        else
          loc = me%outval
        end if
        res(iPoint) = loc(1)
      end do
    end do
  end function tem_eval_polygon_material_scal_3d


  !> Return the material value for point, based on the position in relation to
  !! the polygon.
  !!
  !! The point containment is decided with the help of the Gauss-Bonnet
  !! theorem, which is highly robust, but requires the polygon vertices to
  !! follow an ordering along the polygon lines.
  !!
  !! The returned value will be polygon%inval, if the point is inside and
  !! polygon%outside otherwise.
  function tem_polygon_material_value(me, nComponents, inVal, outVal, point) &
    & result(res)
    ! ----------------------------------------------------------------------
    !> Polygon to describe the material shape.
    type(tem_polygon_vertex_type), intent(in) :: me
    integer, intent(in) :: nComponents
    real(kind=rk), intent(in) :: inVal(nComponents)
    real(kind=rk), intent(in) :: outVal(nComponents)

    !> Point to check against the polygon.
    real(kind=rk), intent(in) :: point(:)
    ! ----------------------------------------------------------------------
    !> Material value at point, as defined by the polygon.
    real(kind=rk) :: res(nComponents)
    real(kind=rk) :: diffa(me%nVertices,2)
    real(kind=rk) :: diffb(me%nVertices,2)
    real(kind=rk) :: anglesum
    integer :: iDir
    ! ----------------------------------------------------------------------

    ! Compute difference between all vertices of the polygon and the given
    ! point.
    diffa(:,1) = me%vertex(:,1) - point(1)
    diffa(:,2) = me%vertex(:,2) - point(2)
    ! Copy the differences to a shifted array, to allow the simultaneous
    ! computation of all the angles for polygon segments.
    ! After the last point, we jump back to the first index to close the
    ! polygon again.
    do iDir=1,2
      diffb(:me%nVertices-1,iDir) = diffa(2:,iDir)
      diffb(me%nVertices,iDir) = diffa(1,iDir)
    end do

    ! With the Gauss-Bonnet theorem, we can decide the containment of the
    ! point within the polygon, based on the sum of all angles under which
    ! the polygon segments are seen from the point.
    anglesum = sum( angle_between(va_x = diffa(:,1), va_y = diffa(:,2), &
      &                           vb_x = diffb(:,1), vb_y = diffb(:,2)) )

    if (abs(anglesum) > 3.0_rk) then
      ! The point is inside, if the sum of angles is 2 Pi actually, but to
      ! allow for numerical inaccuracies and points on the polygon segments
      ! we just check for greater 3 here.
      res = inval
    else
      ! Points with a smaller value are deemed to be outside the polygon.
      ! They might actually be located on a vertex, but the angle at the
      ! vertex would be greater for the outside part than the inside, so it
      ! should definitely be ok to assume this point as outside.
      res = outval
    end if

  end function tem_polygon_material_value


  function tem_polygon_material_movement_single(me, time, nPoint, coord) &
    & result(res)
    ! ----------------------------------------------------------------------
    type(tem_polygon_material_type), intent(in) :: me
    !>velocity value
    real(kind=rk), intent(in) :: time
    !> number of points to get value for
    integer, intent(in) :: nPoint
    !> points
    real(kind=rk), intent(in) :: coord(nPoint,3)
    !> List of values of each point
    real(kind=rk) :: res(nPoint*me%nComponents)
    ! ----------------------------------------------------------------------
    real(kind=rk) :: alpha
    integer :: iPoint
    type(tem_polygon_material_type) :: loc_polygon
    ! ----------------------------------------------------------------------
    loc_polygon = me
    select case(me%moving%movement_kind)
    case ('lin_movement_2d','lin_movement_3d')
      loc_polygon%poly_list(1)%vertex(:,1) = me%poly_list(1)%vertex(:,1) &
        & + me%moving%lin_parameter(1) * time
      loc_polygon%poly_list(1)%vertex(:,2) = me%poly_list(1)%vertex(:,2) &
        & + me%moving%lin_parameter(2) * time

    case('sin_movement_2d','sin_movement_3d')
      loc_polygon%poly_list(1)%vertex(:,1) = me%poly_list(1)%vertex(:,1) &
        & + me%moving%sin_parameter(1)                                   &
        & *sin(2*PI*me%moving%sin_parameter(2)*time)
      loc_polygon%poly_list(1)%vertex(:,2) = me%poly_list(1)%vertex(:,2) &
        & + me%moving%sin_parameter(3)                                   &
        & *sin(2*PI*me%moving%sin_parameter(4)*time)

    case('angleofAttack_2d','angleofAttack_3d' )
     ! rotation around y-axis
      alpha = me%moving%angle_parameter(1) + me%moving%angle_parameter(2) &
        & * sin( me%moving%angle_parameter(3) * time)

      loc_polygon%poly_list(1)%vertex(:,1) = me%poly_list(1)%vertex(:,1)  &
        & * cos(alpha) + me%poly_list(1)%vertex(:,2) * sin(alpha)
      loc_polygon%poly_list(1)%vertex(:,2) = me%poly_list(1)%vertex(:,2)  &
        & * cos(alpha) - me%poly_list(1)%vertex(:,1) * sin(alpha)

    case('rot_movement_2d','rot_movement_3d')
      loc_polygon%poly_list(1)%vertex(:,1) =                                  &
        & cos(me%moving%rot_parameter(3)*time)*(me%poly_list(1)%vertex(:,1)   &
        & - me%moving%rot_parameter(1))                                       &
        & - sin(me%moving%rot_parameter(3)*time)*(me%poly_list(1)%vertex(:,2) &
        & - me%moving%rot_parameter(2)) + me%moving%rot_parameter(1)

      loc_polygon%poly_list(1)%vertex(:,2) =                                  &
        & sin(me%moving%rot_parameter(3)*time)*(me%poly_list(1)%vertex(:,1)   &
        & - me%moving%rot_parameter(1))                                       &
        & + cos(me%moving%rot_parameter(3)*time)*(me%poly_list(1)%vertex(:,2) &
        & - me%moving%rot_parameter(2)) + me%moving%rot_parameter(2)
    case default
      call tem_abort( 'ERROR in tem_polygon_material_module: UNKNOWN movement' &
        & // ' for the polygon in 2D' )
    end select

    select case(me%moving%movement_kind)
    case('lin_movement_2d', 'sin_movement_2d', 'rot_movement_2d', &
        & 'angleofAttack_2d'                                      )
      do iPoint=1, nPoint
        res((iPoint-1)*me%nComponents+1:iPoint*me%nComponents) &
          & = tem_polygon_material_value(                      &
          &   me          = loc_polygon%poly_list(1),          &
          &   nComponents = loc_polygon%nComponents,           &
          &   inVal       = loc_polygon%inVal,                 &
          &   outVal      = loc_polygon%outVal,                &
          &   point       = coord(iPoint,:2)                   )
      end do
    case('lin_movement_3d', 'sin_movement_3d', 'rot_movement_3d', &
        & 'angleofAttack_3d'                                      )

      do iPoint=1, nPoint
        if (coord(iPoint,3) >= me%zmin .and. coord(iPoint,3) <= me%zmax ) then
          res((iPoint-1)*me%nComponents+1:iPoint*me%nComponents) &
            & = tem_polygon_material_value(                      &
            &   me          = loc_polygon%poly_list(1),          &
            &   nComponents = loc_polygon%nComponents,           &
            &   inVal       = loc_polygon%inVal,                 &
            &   outVal      = loc_polygon%outVal,                &
            &   point       = coord(iPoint,:2)                   )
        else
          res((iPoint-1)*me%nComponents+1:iPoint*me%nComponents) = me%outval
        end if
      end do
    case default
      call tem_abort( 'ERROR in tem_polygon_material_module: UNKNOWN movement' &
        & // 'for the polygon in 3D' )
    end select
  end function tem_polygon_material_movement_single


  function tem_polygon_material_movement_multi(me, time, nPoint, coord) &
    & result(res)
    ! ----------------------------------------------------------------------
    type(tem_polygon_material_type), intent(in) :: me
    !>velocity value
    real(kind=rk), intent(in) :: time
    !> number of points to get value for
    integer, intent(in) :: nPoint
    !> points
    real(kind=rk), intent(in) :: coord(nPoint,3)
    !> List of values of each point
    real(kind=rk) :: res(nPoint*me%nComponents)
    ! ----------------------------------------------------------------------
    integer :: iPoint, iPoly
    integer :: iComp
    type(tem_polygon_material_type) :: loc_polygon
    ! ----------------------------------------------------------------------

    loc_polygon = me
    select case(me%moving%movement_kind)
    case ('lin_multi_body_2d','lin_multi_body_3d')
      do iPoly = 1, me%nPoly
        loc_polygon%poly_list(iPoly)%vertex(:,1) =                            &
          & me%poly_list(iPoly)%vertex(:,1) + me%moving%lin_parameter(1) * time
        loc_polygon%poly_list(iPoly)%vertex(:,2) =                            &
          & me%poly_list(iPoly)%vertex(:,2) + me%moving%lin_parameter(2) * time
      end do
    case('sin_multi_body_2d', 'sin_multi_body_3d')
      do ipoly = 1, me%nPoly
        loc_polygon%poly_list(ipoly)%vertex(:,1) =                       &
          & me%poly_list(ipoly)%vertex(:,1) + me%moving%sin_parameter(1) &
          & *sin(2*PI*me%moving%sin_parameter(2)*time)
        loc_polygon%poly_list(ipoly)%vertex(:,2) =                       &
          & me%poly_list(ipoly)%vertex(:,2) + me%moving%sin_parameter(3) &
          & *sin(2*PI*me%moving%sin_parameter(4)*time)
      end do
    case('rot_multi_body_2d', 'rot_multi_body_3d')
      do ipoly = 1, me%nPoly
        loc_polygon%poly_list(ipoly)%vertex(:,1) =                   &
          & cos(me%moving%rot_parameter(3)*time)                     &
          & *(me%poly_list(ipoly)%vertex(:,1)                        &
          & - me%moving%rot_parameter(1))                            &
          & - sin(me%moving%rot_parameter(3)*time)                   &
          & *(me%poly_list(ipoly)%vertex(:,2)                        &
          & - me%moving%rot_parameter(2)) + me%moving%rot_parameter(1)

        loc_polygon%poly_list(ipoly)%vertex(:,2) =                   &
          & sin(me%moving%rot_parameter(3)*time)                     &
          & *(me%poly_list(ipoly)%vertex(:,1)                        &
          & - me%moving%rot_parameter(1))                            &
          & + cos(me%moving%rot_parameter(3)*time)                   &
          & *(me%poly_list(ipoly)%vertex(:,2)                        &
          & - me%moving%rot_parameter(2)) + me%moving%rot_parameter(2)
      end do
    case default
      call tem_abort( 'ERROR in tem_polygon_material_module: UNKNOWN movement' &
        & // 'for the polygon' )
    end select
    do icomp = 1, me%nComponents
      res(iComp::me%nComponents) = me%outval(iComp)
    end do
    select case(me%moving%movement_kind)
    case('lin_multi_body_2d', 'sin_multi_body_2d', 'rot_multi_body_2d')
      do ipoly = 1, me%nPoly
        do iPoint=1, nPoint
          if (res((iPoint-1)*me%nComponents+1) .feq. me%outval(1) ) then
            res((iPoint-1)*me%nComponents+1:iPoint*me%nComponents)   &
              & = tem_polygon_material_value(                        &
              &     me          = loc_polygon%poly_list(iPoly),      &
              &     nComponents = loc_polygon%nComponents,           &
              &     inVal       = loc_polygon%inVal,                 &
              &     outVal      = loc_polygon%outVal,                &
              &     point       = coord(iPoint,:2)                   )
          !else
          ! res = res((iPoint-1)*me%nComponents+1:iPoint*me%nComponents)=me%inval
          end if
        end do
      end do
    case('lin_multi_body_3d', 'sin_multi_body_3d', 'rot_multi_body_3d')
      do ipoly = 1, me%nPoly
        do iPoint=1, nPoint
          if (res((iPoint-1)*me%nComponents+1) .feq. me%outval(1) ) then
            if (coord(iPoint,3) >= me%zmin .and. coord(iPoint,3) <= me%zmax ) &
              & then
              res((iPoint-1)*me%nComponents+1:iPoint*me%nComponents)   &
                & = tem_polygon_material_value(                        &
                &     me          = loc_polygon%poly_list(iPoly),      &
                &     nComponents = loc_polygon%nComponents,           &
                &     inVal       = loc_polygon%inVal,                 &
                &     outVal      = loc_polygon%outVal,                &
                &     point       = coord(iPoint,:2)                   )
            else
              res((iPoint-1)*me%nComponents+1:iPoint*me%nComponents) = me%outval
            end if
          end if
        end do
      end do
    case default
      call tem_abort( 'ERROR in tem_polygon_material_module: UNKNOWN movement' &
        & // 'for the multi body polygon' )
    end select
  end function tem_polygon_material_movement_multi


  !> A subroutine to test the tem_polygon_material_value function in
  !! tem_polygon_material_test.
  !!
  !! The only argument will be true, if all calculations result in the
  !! expected values.
  subroutine tem_polygon_material_test_value(success)
    ! ----------------------------------------------------------------------
    !> Indicator if all tests were computed correctly.
    logical, intent(out) :: success
    ! ----------------------------------------------------------------------
    type(tem_polygon_material_type) :: polygon
    type(tem_polygon_vertex_type) :: me
    real(kind=rk) :: matval(1)
    ! ----------------------------------------------------------------------

    ! Use simple polygon with 5 vertices to check the containment function.
    ! It has roughly this shape:
    !
    !  x-------x
    !   \      |
    !    \     |
    !     x    |
    !    /     |
    !   /      |
    !  x-------x
    !
    ! With the 'x' indicating the 5 vertices.
    me%nVertices = 5
    allocate(me%Vertex(me%nVertices,2))

    me%Vertex(1,:) = [ 0.0_rk,  0.0_rk]
    me%Vertex(2,:) = [-1.0_rk, -1.0_rk]
    me%Vertex(3,:) = [ 1.0_rk, -1.0_rk]
    me%Vertex(4,:) = [ 1.0_rk,  1.0_rk]
    me%Vertex(5,:) = [-1.0_rk,  1.0_rk]

    ! Defining returned values for inside and outside the polygon:
    polygon%nComponents = 1
    allocate(polygon%inval(polygon%nComponents))
    allocate(polygon%outval(polygon%nComponents))

    polygon%inval = 1.0_rk
    polygon%outval = 0.0_rk

    success = .true.
    ! (If any subsequent test fails, we set this to false.)

    write(*,*) 'Checking point containment for this 5 point polygon:'
    write(*,*)
    write(*,*) '(-1,1)---(1,1)'
    write(*,*) '    \      |'
    write(*,*) '     \     |'
    write(*,*) '    (0,0)  |'
    write(*,*) '     /     |'
    write(*,*) '    /      |'
    write(*,*) '(-1,-1)--(1,-1)'
    write(*,*)

    ! Check a point inside the polygon:
    matval = tem_polygon_material_value( me          = me,                 &
      &                                  point       = [0.5_rk, 0.5_rk],   &
      &                                  inval       = polygon%inval,      &
      &                                  outval      = polygon%outval,     &
      &                                  nComponents = polygon%nComponents )

    if (matval(1) < 1.0_rk) success = .false.
    write(*,*) 'Point (0.5, 0.5) has value:', matval

    ! Check a point outside the polygon:
    matval = tem_polygon_material_value( me          = me,                 &
      &                                  point       = [-1.0_rk, 0.0_rk],  &
      &                                  inval       = polygon%inval,      &
      &                                  outval      = polygon%outval,     &
      &                                  nComponents = polygon%nComponents )

    if (matval(1) > epsilon(1.0_rk)) success = .false.
    write(*,*) 'Point (-1.0, 0.0) has value:', matval

    ! Check a point on a side of the polygon (should be considered inside):
    matval = tem_polygon_material_value( me          = me,                 &
      &                                  point       = [-0.5_rk, -0.5_rk], &
      &                                  inval       = polygon%inval,      &
      &                                  outval      = polygon%outval,     &
      &                                  nComponents = polygon%nComponents )

    if (matval(1) < 1.0_rk) success = .false.
    write(*,*) 'Point (-0.5, -0.5) has value:', matval

    ! Check an outward pointing corner (should be outside):
    matval = tem_polygon_material_value( me          = me,                 &
      &                                  point       = [1.0_rk, 1.0_rk],   &
      &                                  inval       = polygon%inval,      &
      &                                  outval      = polygon%outval,     &
      &                                  nComponents = polygon%nComponents )

    if (matval(1) > epsilon(1.0_rk)) success = .false.
    write(*,*) 'Point (1.0, 1.0) has value:', matval

    ! Check an inward pointing corner (should be inside):
    matval = tem_polygon_material_value( me          = me,                 &
      &                                  point       = [0.0_rk, 0.0_rk],   &
      &                                  inval       = polygon%inval,      &
      &                                  outval      = polygon%outval,     &
      &                                  nComponents = polygon%nComponents )

    if (matval(1) < 1.0_rk) success = .false.
    write(*,*) 'Point (0.0, 0.0) has value:', matval

  end subroutine tem_polygon_material_test_value



  !> Compute the angle between to vectors (they should not both be the 0
  !! vector).
  !!
  !! This function uses the crossproduct and arctan to find the angle between
  !! two 2D vectors. It takes the four components of the vectors as scalar
  !! arguments to allow vectorized evaluation of multiple vector pairs at once.
  !! If one of the vectors is zero, the returned angle is also zero.
  !! Also multiplicities of Pi are ignored and a 0 angle will be returned for
  !! them.
  elemental function angle_between(va_x, va_y, vb_x, vb_y) result(angle)
    ! ----------------------------------------------------------------------
    !> The first vector va
    real(kind=rk), intent(in) :: va_x, va_y

    !> The second vector vb
    real(kind=rk), intent(in) :: vb_x, vb_y

    !> The angle betweend va and vb
    real(kind=rk) :: angle
    ! ----------------------------------------------------------------------

    angle = 0.0_rk

    ! Only proceed if both vectors are non-zero.
    ! If one of the vectors is 0, we assume a zero angle.
    if (      (abs(va_x)+abs(va_y) > 2*tiny(va_x)) &
      & .and. (abs(vb_x)+abs(vb_y) > 2*tiny(vb_x)) ) then
      ! The angle is arctan( (va X vb)/(va*vb) ):
      ! Using atan2 here for numerical stability and correct signs.
      angle = atan2( (va_x*vb_y - va_y*vb_x), &
        &            (va_x*vb_x + va_y*vb_y)  )

      ! If the angle is close enough to Pi, consider this as a colinear
      ! 0 angle (avoid multiplicities of Pi in the summed angle).
      if (abs(angle) > tolm*Pi) angle = 0.0_rk
    end if

  end function angle_between


  !> A subroutine to test the angle_between function in
  !! tem_polygon_material_test.
  !!
  !! The only argument will be true, if all calculations result in the
  !! expected values.
  subroutine tem_polygon_material_test_angle(success)
    ! ----------------------------------------------------------------------
    !> Indicator if all tests were computed correctly.
    logical, intent(out) :: success
    ! ----------------------------------------------------------------------
    real(kind=rk) :: res
    ! ----------------------------------------------------------------------

    success = .true.

    write(*,*) 'Checking angle computation between two 2D Vectors:'
    ! 90 degree
    res = angle_between( va_x = 1.0_rk, va_y = 0.0_rk, &
      &                  vb_x = 0.0_rk, vb_y = 1.0_rk  )

    if (res < 0.5_rk*tolm*PI .or. res > 0.5_rk*tolp*PI) success = .false.
    write(*,*) '  90 deg:', res


    ! 45 degree
    res = angle_between( va_x = 10.0_rk, va_y = 0.0_rk, &
      &                  vb_x =  2.0_rk, vb_y = 2.0_rk  )

    if (res < 0.25_rk*tolm*PI .or. res > 0.25_rk*tolp*PI) success = .false.
    write(*,*) '  45 deg:', res


    ! -90 degree
    res = angle_between( va_x = 10.0_rk, va_y =  0.0_rk, &
      &                  vb_x =  0.0_rk, vb_y = -3.0_rk  )

    if (res < -0.5_rk*tolp*PI .or. res > -0.5_rk*tolm*PI) success = .false.
    write(*,*) ' -90 deg:', res


    ! -135 degree
    res = angle_between( va_x = 10.0_rk, va_y =  0.0_rk, &
      &                  vb_x = -3.0_rk, vb_y = -3.0_rk  )

    if (res < -0.75_rk*tolp*PI .or. res > -0.75_rk*tolm*PI) success = .false.
    write(*,*) '-135 deg:', res


    ! 120 degree
    res = angle_between( va_x =  0.0_rk,       va_y =  4.5_rk, &
      &                  vb_x = -sqrt(3.0_rk), vb_y = -1.0_rk  )

    if ( res < 2.0_rk/3.0_rk*tolm*PI &
      &  .or. res > 2.0_rk/3.0_rk*tolp*PI) success = .false.
    write(*,*) ' 120 deg:', res

    ! 0 degree
    res = angle_between( va_x = 3.0_rk, va_y = -1.0_rk, &
      &                  vb_x = 3.0_rk, vb_y = -1.0_rk  )

    if (res > epsilon(1.0_rk)) success = .false.
    write(*,*) '   0 deg:', res

    ! 180 degree
    res = angle_between( va_x = -1.0_rk, va_y = -1.0_rk, &
      &                  vb_x =  1.0_rk, vb_y =  1.0_rk  )

    if (res > epsilon(1.0_rk)) success = .false.
    write(*,*) ' 180 deg:', res

    ! 0 Vector
    res = angle_between( va_x = 1.0_rk, va_y = -1.0_rk, &
      &                  vb_x = 0.0_rk, vb_y =  0.0_rk  )

    if (res > epsilon(1.0_rk)) success = .false.
    write(*,*) 'Zero vec:', res

  end subroutine tem_polygon_material_test_angle

end module tem_polygon_material_module
