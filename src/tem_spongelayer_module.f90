! Copyright (c) 2013-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2013,2015-2016,2021-2022 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2018 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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
! ****************************************************************************** !
!> This module gathers the various spatial spongeLayer definitions
!!
module tem_spongeLayer_module

  ! include treelm modules
  use env_module,          only: rk, long_k, labelLen
  use tem_param_module,    only: PI
  use tem_aux_module,      only: tem_abort
  use tem_logging_module,  only: logUnit
  use treelmesh_module,    only: treelmesh_type
  use tem_geometry_module, only: tem_BaryOfId

  ! include aotus modules
  use aotus_module,       only: flu_State, aot_get_val,           &
    &                           aoterr_NonExistent, aoterr_Fatal
  use aot_table_module,   only: aot_table_open, aot_table_close,  &
   &                            aot_table_length, aot_get_val

  implicit none
  private

  public :: tem_spongeLayer_plane_type, tem_load_spongeLayer_plane
  public :: tem_spongeLayer_radial_type, tem_load_spongeLayer_radial
  public :: tem_spongeLayer_box_type, tem_load_spongeLayer_box
  public :: tem_spongelayer_plane_for
  public :: tem_spongelayer_box_for
  public :: tem_spongelayer_box2d_for
  public :: tem_spongeLayer_radial_for
  public :: tem_viscSpongeLayer_plane_for
  public :: tem_viscSpongeLayer_box_for
  public :: tem_viscSpongeLayer_box2d_for
  public :: tem_viscSpongeLayer_radial_for


  !> This type contains base data defined for all sponge layers
  type spongeLayer_base_type
    !> Thickness of the sponge layer.
    !! For planar sponge thickness is defined implicitly in place_normal
    real(kind=rk) :: thickness
    !> Damp factor or strength for the sponge Layer
    real(kind=rk) :: dampFactor
    !> damping exponent for the sponge layer
    real(kind=rk) :: dampExponent
    !> damping profile
    character(len=labelLen) :: dampProfile
    !> target states.
    !! For viscous sponge, viscosity is stored and multiplied with sponge
    !! strength
    real(kind=rk), allocatable :: targetState(:)
  end type spongeLayer_base_type

  !> This type contains data to define spongeLayer plane
  type, extends(spongeLayer_base_type) :: tem_spongeLayer_plane_type
    !> Sponge Plane origin
    real(kind=rk) :: origin(3)
    !> Sponge Plane normal
    real(kind=rk) :: normal(3)
  end type tem_spongeLayer_plane_type

  !> This type contains data to define spongeLayer box
  type, extends(spongeLayer_base_type) :: tem_spongeLayer_box_type
    !> Box origin, bottom left corner of sponge layer
    real(kind=rk) :: origin(3)
    !> Length of box in each dimension 
    real(kind=rk) :: extent(3)
    !> To create sponge box with rounded corners
    logical :: rounded_corner
    !> Corner radius for rounded box
    real(kind=rk) :: corner_radius
  end type tem_spongeLayer_box_type

  !> This type contains data to define spongeLayer radial
  type, extends(spongeLayer_base_type) :: tem_spongeLayer_radial_type
    !> Sponge radial origin
    real(kind=rk) :: origin(3)
    !> Sponge inner radius i.e. sponge start.
    !! Outer radius is computed by adding thickness to inner radius.
    real(kind=rk) :: radius
  end type tem_spongeLayer_radial_type

  !> Interface for sponge layer plane
  interface tem_spongeLayer_plane_for
    module procedure spongeLayer_plane_scalar_for_coord
    module procedure spongeLayer_plane_scalar_for_treeIDs
    module procedure spongeLayer_plane_vector_for_coord
    module procedure spongeLayer_plane_vector_for_treeIDs
  end interface tem_spongeLayer_plane_for

  !> Interface for sponge layer box
  interface tem_spongeLayer_box_for
    module procedure spongeLayer_box_scalar_for_coord
    module procedure spongeLayer_box_scalar_for_treeIDs
    module procedure spongeLayer_box_vector_for_coord
    module procedure spongeLayer_box_vector_for_treeIDs
  end interface tem_spongeLayer_box_for

  !> Interface for sponge layer box 2d
  interface tem_spongeLayer_box2d_for
    module procedure spongeLayer_box2d_scalar_for_coord
    module procedure spongeLayer_box2d_scalar_for_treeIDs
    module procedure spongeLayer_box2d_vector_for_coord
    module procedure spongeLayer_box2d_vector_for_treeIDs
  end interface tem_spongeLayer_box2d_for

  !> Interface for sponge layer radial
  interface tem_spongeLayer_radial_for
    module procedure spongeLayer_radial_scalar_for_coord
    module procedure spongeLayer_radial_scalar_for_treeIDs
    module procedure spongeLayer_radial_vector_for_coord
    module procedure spongeLayer_radial_vector_for_treeIDs
  end interface tem_spongeLayer_radial_for

  !> Interface for viscous sponge layer plane
  interface tem_viscSpongeLayer_plane_for
    module procedure viscSpongeLayer_plane_for_coord
    module procedure viscSpongeLayer_plane_for_treeIDs
  end interface tem_viscSpongeLayer_plane_for

  !> Interface for viscous sponge layer box
  interface tem_viscSpongeLayer_box_for
    module procedure viscSpongeLayer_box_for_coord
    module procedure viscSpongeLayer_box_for_treeIDs
  end interface tem_viscSpongeLayer_box_for

  !> Interface for viscous sponge layer box
  interface tem_viscSpongeLayer_box2d_for
    module procedure viscSpongeLayer_box2d_for_coord
    module procedure viscSpongeLayer_box2d_for_treeIDs
  end interface tem_viscSpongeLayer_box2d_for

  !> Interface for viscous sponge layer radial
  interface tem_viscSpongeLayer_radial_for
    module procedure viscSpongeLayer_radial_for_coord
    module procedure viscSpongeLayer_radial_for_treeIDs
  end interface tem_viscSpongeLayer_radial_for


contains


  ! -------------------------------------------------------------------------- !
  !> This subroutine load data for standard plane sponge layer
  !! Example:
  !!
  !!```lua
  !! spatial = {
  !!   -- supported options: 'spongelayer_plane', 'spongelayer_plane_1d', 
  !!   --                    'spongelayer_plane_2d', 'viscous_spongelayer_plane'
  !!   predefined = 'spongelayer',
  !!   origin = {0.0,0.0,0.0},
  !!   normal = {1.0, 0.0, 0.0},
  !!   thickness = 0.5,
  !!   damp_profile = 'linear', --'exponential', 'polynomial_n5', 'polynomial_n6'
  !!   damp_factor = 0.5,
  !!   damp_exponent = 1.0,
  !!   target_state = {
  !!     Default: density, velocityX, velocityY, velocityZ and pressure
  !!     density = 1.0,
  !!     pressure = 1.0,
  !!     velocityX = 0.0, velocityY = 0.0, velocityZ = 0.0
  !! }
  !!```
  subroutine tem_load_spongeLayer_plane(me, conf, thandle, ndim, nComp, &
    &                                   stateName)
    ! --------------------------------------------------------------------------
    !> Plane spongeLayer data type
    type(tem_spongeLayer_plane_type), intent(out) :: me
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: thandle
    !> number of Dimension for nonViscous sponges
    integer, intent(in) :: nDim
    !> Number of component of St-Fun variable under which this spatial function
    !! is defined
    integer, intent(in) :: nComp
    !> Load stateName from target_state table
    character(len=*), intent(in), optional :: stateName
    ! --------------------------------------------------------------------------
    integer :: vError(3), errfatal(3)
    real(kind=rk) :: thickness
    ! --------------------------------------------------------------------------
    errfatal = aotErr_Fatal
    ! Plane_origin
    call aot_get_val( L       = conf,      &
      &               thandle = thandle,   &
      &               key     = 'origin',  &
      &               val     = me%origin, &
      &               ErrCode = vError     )
    if (any(btest( vError, errFatal )) ) then
      write(logUnit(1),*) 'ERROR reading the plane_origin of sponge layer. ' &
        &              // 'It should have 3 entries for each coordinate.'
      call tem_abort()
    end if

    write(logUnit(1),*) ' * Origin =', me%origin

    ! Plane_normal
    call aot_get_val( L       = conf,      &
      &               thandle = thandle,   &
      &               key     = 'normal',  &
      &               val     = me%normal, &
      &               ErrCode = vError     )
    if (any(btest( vError, errFatal )) ) then
      write(logUnit(1),*) 'ERROR reading the plane_normal of sponge layer. ' &
        &              // 'It should have 3 entries for each coordinate.'
      call tem_abort()
    end if
    write(logUnit(1),*) ' * Normal =', me%normal

    ! Compute thickness from normal and use it only if thickness is not
    ! defined seperately.
    thickness = sqrt(me%normal(1)**2 + me%normal(2)**2 + me%normal(3)**2)

    ! Normalize the normal
    me%normal = me%normal/thickness
    write(logUnit(1),*) ' * Normalized normal =', me%normal

    ! Load base information required for sponge layer definition like
    ! damp_factor, damp_exponent and target_state
    call load_spongeLayer( conf      = conf,                     &
      &                    thandle   = thandle,                  &
      &                    me        = me%spongeLayer_base_type, &
      &                    nDim      = nDim,                     &
      &                    nComp     = nComp,                    &
      &                    stateName = stateName,                &
      &                    thickness = thickness                 )

  end subroutine tem_load_spongeLayer_plane
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This subroutine load data for radial sponge layer
  !! Example:
  !!
  !!```lua
  !! spatial = {
  !!   --supported options: 'spongelayer_radial','sponge_radial_2d',
  !!   --                   'viscous_spongelayer_radial',
  !!   --                   'viscous_spongelayer_radial_2d'
  !!   predefined = 'viscous_spongelayer_radial', 
  !!   origin = {0.0,0.0,0.0},
  !!   radius = 1.0, -- Sponge start
  !!   thickness = 0.3,
  !!   damp_profile = 'linear', --'exponential', 'polynomial_n5', 'polynomial_n6'
  !!   damp_factor = 0.5,
  !!   damp_exponent = 1.0,
  !!   target_state = {
  !!     Default: density, velocityX, velocityY, velocityZ and pressure
  !!     viscosity = 1e-3
  !! }
  !!```
  subroutine tem_load_spongeLayer_radial(me, conf, thandle, nDim, nComp, &
    &                                    stateName)
    ! --------------------------------------------------------------------------
    !> Radial spongeLayer data type
    type(tem_spongeLayer_radial_type), intent(out) :: me
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: thandle
    !> number of Dimension for nonViscous sponges
    integer, intent(in) :: nDim
    !> Number of component of St-Fun variable under which this spatial function
    !! is defined
    integer, intent(in) :: nComp
    !> Load stateName from target_state table
    character(len=*), intent(in), optional :: stateName
    ! --------------------------------------------------------------------------
    integer :: iError
    integer :: vError(3), errfatal(3)
    ! --------------------------------------------------------------------------

    errfatal = aotErr_Fatal
    ! Sponge origin
    call aot_get_val( L       = conf,      &
      &               thandle = thandle,   &
      &               key     = 'origin',  &
      &               val     = me%origin, &
      &               ErrCode = vError     )
    if (any(btest( vError, errFatal )) ) then
      write(logUnit(1),*) 'ERROR reading the sponge origin, ' &
        &              // 'origin is not well defined. '   &
        &              // 'It should have 3 entries for each coordinate.'
      call tem_abort()
    end if
    write(logUnit(1),*) ' * Origin =', me%origin

    ! Sponge inner radius
    call aot_get_val( L       = conf,      &
      &               thandle = thandle,   &
      &               key     = 'radius',  &
      &               val     = me%radius, &
      &               ErrCode = iError     )
    if (btest(iError, aotErr_Fatal)) then
      write(logUnit(1),*) 'ERROR reading the sponge inner radius. '
      call tem_abort()
    end if
    write(logUnit(1),*) ' * Inner radius =', me%radius

    ! Load base information required for sponge layer definition like
    ! damp_factor, damp_exponent and target_state
    call load_spongeLayer( conf      = conf,                     &
      &                    thandle   = thandle,                  &
      &                    me        = me%spongeLayer_base_type, &
      &                    nDim      = nDim,                     &
      &                    nComp     = nComp,                    &
      &                    stateName = stateName                 )

  end subroutine tem_load_spongeLayer_radial
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This subroutine load data for sponge layer box
  !! Example:
  !!
  !!```lua
  !! spatial = {
  !!   --supported options: 'spongelayer_box', 'spongelayer_box_2d'
  !!   --                   'viscous_spongelayer_box',
  !!   predefined = 'spongelayer_box',
  !!   origin = {0.0,0.0,0.0},
  !!   extent = {1.0, 2.0, 3.0},
  !!   thickness = 0.3,
  !!   damp_profile = 'linear', --'exponential', 'polynomial_n5', 'polynomial_n6'
  !!   damp_factor = 0.5,
  !!   damp_exponent = 1.0,
  !!   target_state = {
  !!     Default: density, velocityX, velocityY, velocityZ and pressure
  !!     density = 1.0,
  !!     pressure = 1.0,
  !!     velocityX = 0.0, velocityY = 0.0, velocityZ = 0.0
  !! }
  !!```
  subroutine tem_load_spongeLayer_box(me, conf, thandle, ndim, nComp, stateName)
    ! --------------------------------------------------------------------------
    !> Box spongeLayer data type
    type(tem_spongeLayer_box_type), intent(out) :: me
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: thandle
    !> number of Dimension for nonViscous sponges
    integer, intent(in) :: nDim
    !> Number of component of St-Fun variable under which this spatial function
    !! is defined
    integer, intent(in) :: nComp
    !> Load stateName from target_state table
    character(len=*), intent(in), optional :: stateName
    ! --------------------------------------------------------------------------
    integer :: vError(3), errfatal(3), iError
    real(kind=rk) :: min_halfextent
    ! --------------------------------------------------------------------------
    errfatal = aotErr_Fatal
    ! Box origin
    call aot_get_val( L       = conf,      &
      &               thandle = thandle,   &
      &               key     = 'origin',  &
      &               val     = me%origin, &
      &               ErrCode = vError     )
    if (any(btest( vError, errFatal )) ) then
      write(logUnit(1),*) 'ERROR reading the origin of box sponge layer. ' &
        &              // 'It should have 3 entries for each coordinate.'
      call tem_abort()
    end if

    write(logUnit(1),*) ' * Origin =', me%origin

    ! Box extent
    call aot_get_val( L       = conf,      &
      &               thandle = thandle,   &
      &               key     = 'extent',  &
      &               val     = me%extent, &
      &               ErrCode = vError     )
    if (any(btest( vError, errFatal )) ) then
      write(logUnit(1),*) 'ERROR reading the extent of box sponge layer. ' &
        &              // 'It should have 3 entries for each dimension.'
      call tem_abort()
    end if

    write(logUnit(1),*) ' * Extent =', me%extent

    ! Additional parameters for box with rounded corner
    call aot_get_val( L       = conf,              &
      &               thandle = thandle,           &
      &               key     = 'rounded_corner',  &
      &               val     = me%rounded_corner, &
      &               default = .false.,           &
      &               ErrCode = iError             )
    if (btest( iError, aotErr_Fatal )) then
      write(logUnit(1),*) 'WARNING: reading the rounded_corner of box '&
        &               //'sponge layer. Setting it to false'
      me%rounded_corner = .false.
    end if
    write(logUnit(1),*) ' * Rounded_corner =', me%rounded_corner

    if (me%rounded_corner) then
      call aot_get_val( L       = conf,              &
        &               thandle = thandle,           &
        &               key     = 'corner_radius',   &
        &               val     = me%corner_radius,  &
        &               ErrCode = iError             )
      if (btest( iError, aotErr_Fatal )) then
        write(logUnit(1),*) 'ERROR reading the corner_radius of box '&
          &               //'sponge layer with rounded corner.'
        call tem_abort()
      end if
      ! if corner radius is greater than minimum half extent of the box
      ! then set corner_radius to minimum of half extent
      min_halfextent = minval(me%extent(1:nDim))*0.5_rk
      if (me%corner_radius > min_halfextent) then
        write(logUnit(1),*) 'WARNING: corner_radius is greater than half of '
        write(logUnit(1),*) '  min of box extent. '
        write(logUnit(1),*) 'Setting it to half of min extent'
        me%corner_radius = min_halfextent
      end if
      write(logUnit(1),*) ' * Corner radius =', me%corner_radius
    end if

    ! Load base information required for sponge layer definition like
    ! damp_factor, damp_exponent and target_state
    call load_spongeLayer( conf      = conf,                     &
      &                    thandle   = thandle,                  &
      &                    me        = me%spongeLayer_base_type, &
      &                    nDim      = nDim,                     &
      &                    nComp     = nComp,                    &
      &                    stateName = stateName                 )

  end subroutine tem_load_spongeLayer_box
  ! -------------------------------------------------------------------------- !


  ! -------------------------------------------------------------------------- !
  !> This routine load base info for sponge layer
  subroutine load_spongeLayer(conf, thandle, me, ndim, nComp, stateName, &
    &                         thickness)
    ! --------------------------------------------------------------------------
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: thandle
    !> base spongeLayer data type
    type(spongeLayer_base_type), intent(out) :: me
    !> number of spatial dimensions
    integer, intent(in) :: ndim
    !> Number of component of St-Fun variable under which this spatial function
    !! is defined
    integer, intent(in) :: nComp
    !> Load stateName from target_state table
    character(len=*), intent(in), optional :: stateName
    !> Thickness computed from sponge layer plane normal. Use this thickness
    !! If thickness is not defined.
    real(kind=rk), intent(in), optional :: thickness
    ! --------------------------------------------------------------------------
    integer :: iError, ts_handle
    ! --------------------------------------------------------------------------
    ! Thickness
    call aot_get_val( L       = conf,         &
      &               thandle = thandle,      &
      &               key     = 'thickness',  &
      &               val     = me%thickness, &
      &               ErrCode = iError        )
    if (btest( iError, aotErr_Fatal )) then
      if (present(thickness)) then
        me%thickness = thickness
      else
        write(logUnit(1),*) 'ERROR reading the thickness of sponge layer.'
        write(logUnit(1),*) 'Thickness is required to calculate sponge end.'
        call tem_abort()
      end if
    end if

    write(logUnit(1),*) ' * Thickness =', me%thickness

    !damp_factor
    call aot_get_val( L       = conf,          &
      &               thandle = thandle,       &
      &               key     = 'damp_factor', &
      &               val     = me%dampFactor, &
      &               ErrCode = iError         )
    if (btest( iError, aotErr_Fatal )) then
      write(logUnit(1),*) 'ERROR reading the damp_factor of sponge layer.'
      call tem_abort()
    end if
    write(logUnit(1),*) ' * Damp_factor =', me%dampFactor

    !damp_profile
    call aot_get_val( L       = conf,           &
      &               thandle = thandle,        &
      &               key     = 'damp_profile', &
      &               val     = me%dampProfile, &
      &               ErrCode = iError          )

    ! Viscous sponge works only with exponential profile so no need load 
    ! damp_profile
    if (present(stateName)) then
      if (trim(stateName) == 'viscosity') then
        me%dampProfile = 'exponential'
        iError = 0
      end if
    end if

    if (btest( iError, aotErr_Fatal )) then
      write(logUnit(1),*) 'ERROR reading the damp_profile of sponge layer.'
      call tem_abort()
    end if
    write(logUnit(1),*) ' * Damp_profile =', me%dampProfile

    select case (trim(me%dampProfile))
    case ('linear')
      me%dampExponent = 1.0
    case ('exponential')
      !damp_exponent
      call aot_get_val( L       = conf,            &
        &               thandle = thandle,         &
        &               key     = 'damp_exponent', &
        &               val     = me%dampExponent, &
        &               ErrCode = iError,          &
        &               default = 1.0_rk           )

      if (btest( iError, aotErr_Fatal )) then
        write(logUnit(1),*) 'ERROR reading the damp_exponent of exponential ' &
          &                 //'sponge layer.'
        call tem_abort()
      end if
    case ('polynomial_n5', 'polynomial_n6')
      me%dampExponent = 1.0
    case default
      write(logUnit(1),*) 'ERROR unknown damp_profile for sponge layer.'
      write(logUnit(1),*) 'Supported options: '
      write(logUnit(1),*) '  linear, exponential, polynomial_n5, polynomial_n6'
      call tem_abort()
    end select

    write(logUnit(1),*) ' * Damp_exponent =', me%dampExponent


    ! Load stateName provided by caller function
    if ( present(stateName) .and. nComp == 1) then
      allocate(me%targetState(1))
      call aot_table_open( L       = conf,          &
        &                  parent  = thandle,       &
        &                  thandle = ts_handle,     &
        &                  key     = 'target_state' )

      call aot_get_val( L       = conf,               &
        &               thandle = ts_handle,          &
        &               key     = trim(stateName),    &
        &               val     = me%targetState(1),  &
        &               ErrCode = iError              )

      if (btest(iError, aotErr_Fatal)) then
        write(logUnit(1),*) 'ERROR reading the target state: '//trim(stateName)
        call tem_abort()
      end if
      call aot_table_close(conf, thandle)

      write(logUnit(1),*) ' * Target state:'
      write(logUnit(1),*) '   '//trim(stateName)//' =', me%targetState(1)
    else if (nComp > 1) then
      call load_defaultTargetState( conf        = conf,          &
        &                           parent      = thandle,       &
        &                           nDim        = nDim,          &
        &                           nComp       = nComp,         &
        &                           targetState = me%targetState )
    else
      write(logUnit(1),*) 'WARNING: nComp = 1 so no target states are loaded'
    end if

  end subroutine load_spongeLayer
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This routine loads state names from target_state table
  subroutine load_defaultTargetState(conf, parent, nDim, nComp, targetState)
    ! --------------------------------------------------------------------------
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: parent
    !> Number of dimension
    integer, intent(in) :: nDim
    !> Number of component of St-Fun variable under which this spatial function
    !! is defined
    integer, intent(in) :: nComp
    !> Target state value
    real(kind=rk), allocatable, intent(out) :: targetState(:)
    ! --------------------------------------------------------------------------
    integer :: thandle, iState, nState
    integer :: iError
    !> List of stateNames
    character(len=labelLen), allocatable :: stateName(:)
    ! --------------------------------------------------------------------------
    nState = 0
    select case (nDim)
    case (3)
      nState = 5
      allocate(stateName(nState))
      stateName = [ 'density  ', 'velocityX', &
        &           'velocityY', 'velocityZ', &
        &           'pressure '               ]
    case (2)
      nState = 4
      allocate(stateName(nState))
      stateName = [ 'density  ', 'velocityX', &
        &           'velocityY', 'pressure '  ]

    case (1)
      nState = 3
      allocate(stateName(nState))
      stateName = [ 'density  ', 'velocityX', &
        &           'pressure '               ]
    end select

    ! nState must be nComp - 1 to return all target states when evaluating
    ! the space-time function
    if (nComp /= nState+1) then
      write(logUnit(1),*) 'Error: Expected ncomponents = ', nState+1
      write(logUnit(1),*) '  Defined ncomponents in st-fun is ',  ncomp
      call tem_abort()
    end if

    ! target_state
    allocate(targetState(nState))
    call aot_table_open( L       = conf,          &
      &                  parent  = parent,        &
      &                  thandle = thandle,       &
      &                  key     = 'target_state' )

    write(logUnit(1),*) ' * Target state:'
    do iState = 1, nState
      call aot_get_val( L       = conf,                    &
           &            thandle = thandle,                 &
           &            key     = trim(stateName(iState)), &
           &            val     = targetState(iState),     &
           &            ErrCode = iError                   )

      if (btest(iError, aoterr_Fatal)) then
        write(*,*) 'FATAL Error occured, when loading target state: ' &
          &      // trim(stateName(iState))// '! Aborting'
        call tem_abort()
      end if

      write(logUnit(1),*) '   '//trim(stateName(iState))//' =', &
        &                 targetState(iState)
    end do

    call aot_table_close(conf, thandle)

  end subroutine load_defaultTargetState
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !
  !                            SPONGE LAYER PLANE
  !
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function returns the sigma for the planar shape spongelayer
  function spongelayer_plane_scalar_for_coord(me, coord, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_plane_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    select case (trim(me%dampProfile))
    case ('linear', 'exponential')
      res(:) = spongeLayer_plane_expon_for_coord(me, coord, n)
    case ('polynomial_n5')
      res(:) = spongeLayer_plane_polyn5_for_coord(me, coord, n)
    case ('polynomial_n6')
      res(:) = spongeLayer_plane_polyn6_for_coord(me, coord, n)
    end select

  end function spongelayer_plane_scalar_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function returns the sigma for the spongelayer computed from
  !! exponential function
  function spongelayer_plane_expon_for_coord(me, coord, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_plane_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(3), normal(3), vec(3), proj_len
    ! --------------------------------------------------------------------------
    origin(:) = me%origin
    ! Normal was normalized in load routine so multiply thickness to convert
    ! normal to user defined normal
    normal(:) = me%normal * me%thickness

    do i = 1,n
      vec(:) = coord(i,:) - origin(:)
      proj_len = (vec(1)*normal(1)+ vec(2)*normal(2)+vec(3)*normal(3))/   &
                 (normal(1)**2 + normal(2)**2 + normal(3)**2)
      sigma = me%dampFactor*((proj_len)**me%dampExponent)
      if (proj_len > 0 .and. proj_len < 1) then
        res(i) = sigma
      else if (proj_len > 1) then
        res(i) = me%dampFactor
      else
        res(i) = 0.0_rk
      end if

    end do

  end function spongelayer_plane_expon_for_coord
  ! -------------------------------------------------------------------------- !

 
  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the spongeLayer from coord for
  !! the polynomial order n5.
  !! Sponge profile:
  !! $$ \sigma(x) = \sigma_A*3125*(L+x0-x)*(x-x0)^4)/(256*(L)^5 $$
  !! where, \sigma_A - sponge strength, L - thickness, x0 - start of sponge.
  !!
  !! Profile is taken from:
  !! Xu, Hui; Sagaut, Pierre (2013): Analysis of the absorbing layers for the 
  !! weakly-compressible lattice Boltzmann methods. In Journal of Computational 
  !! Physics 245, pp. 14-42. DOI: 10.1016/j.jcp.2013.02.051.
  function spongeLayer_plane_polyn5_for_coord(me, coord, n)  result(res)
    ! --------------------------------------------------------------------------
    !> Spatial sponge layer to evaluate
    type(tem_spongeLayer_plane_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(3), normal(3), vec1(3), vec2(3)
    real(kind=rk) :: proj_len1, proj_len2, const_fac
    ! --------------------------------------------------------------------------
    origin(:) = me%origin
    normal(:) = me%normal
    const_fac = 3125_rk/(256_rk*me%thickness**5)

    do i = 1,n
      vec1(:) = coord(i,:) - origin(:)
      vec2(:) = me%thickness*normal(:) + origin(:) - coord(i,:)
      proj_len1 = vec1(1)*normal(1) + vec1(2)*normal(2) + vec1(3)*normal(3)
      proj_len2 = vec2(1)*normal(1) + vec2(2)*normal(2) + vec2(3)*normal(3)

      if (proj_len1 > 0 .and. proj_len2 > 0) then
        sigma = const_fac * proj_len2 * (proj_len1**4)
        res(i) = sigma*me%dampFactor
      else if (proj_len2 < 0) then ! If coord is beyond thickness
        res(i) = me%dampFactor
      else
        res(i) = 0.0_rk
      end if
    end do

  end function spongeLayer_plane_polyn5_for_coord
  ! -------------------------------------------------------------------------- !


  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the spongeLayer from coord for
  !! the polynomial order n6.
  !! Sponge profile:
  !! $$ \sigma(x) = \sigma_A*729*(L+x0-x)^2*(x-x0)^4)/(16*(L)^6 $$
  !! where, \sigma_A - sponge strength, L - thickness, x0 - start of sponge.
  function spongeLayer_plane_polyn6_for_coord(me, coord, n)  result(res)
    ! --------------------------------------------------------------------------
    !> Spatial sponge layer to evaluate
    type(tem_spongeLayer_plane_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(3), normal(3), vec1(3), vec2(3)
    real(kind=rk) :: proj_len1, proj_len2, const_fac
    ! --------------------------------------------------------------------------
    origin(:) = me%origin
    normal(:) = me%normal
    const_fac = 729_rk/(16_rk*me%thickness**6)

    do i = 1,n
      vec1(:) = coord(i,:) - origin(:)
      vec2(:) = me%thickness*normal(:) + origin(:) - coord(i,:)
      proj_len1 = vec1(1)*normal(1) + vec1(2)*normal(2) + vec1(3)*normal(3)
      proj_len2 = vec2(1)*normal(1) + vec2(2)*normal(2) + vec2(3)*normal(3)

      sigma = const_fac * proj_len2**2 * (proj_len1**4)
      if (proj_len1 > 0) then
        res(i) = min(1.0_rk, sigma) * me%dampFactor
      else
        res(i) = 0.0_rk
      end if
    end do

  end function spongeLayer_plane_polyn6_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the spongelayer and fills up
  !!  the res with the target state
  function spongelayer_plane_vector_for_coord(me, nComp, coord, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_plane_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> Number of entrys in each array
    integer, intent(in) :: ncomp
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n,ncomp)
    ! --------------------------------------------------------------------------
    integer :: i
    ! --------------------------------------------------------------------------
    res(:, 1) = spongeLayer_plane_scalar_for_coord(me, coord, n)

    if (ncomp > 1) then
      do i = 1,n
        res(i,2:) = me%targetState(:)
      end do
    end if

  end function spongelayer_plane_vector_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function returns the sigma for the spongelayer from treeids
  function spongelayer_plane_scalar_for_treeIDs(me, treeids, tree, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_plane_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    select case (trim(me%dampProfile))
    case ('linear', 'exponential')
      res(:) = spongeLayer_plane_expon_for_treeIDs(me, treeIDs, tree, n)
    case ('polynomial_n5')
      res(:) = spongeLayer_plane_polyn5_for_treeIDs(me, treeIDs, tree, n)
    case ('polynomial_n6')
      res(:) = spongeLayer_plane_polyn6_for_treeIDs(me, treeIDs, tree, n)
    end select

  end function spongeLayer_plane_scalar_for_treeIDs
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function returns the sigma for the exponential spongelayer from 
  !! treeids
  function spongelayer_plane_expon_for_treeIDs(me, treeids, tree, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_plane_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(3), normal(3), vec(3), coord(3), proj_len
    ! --------------------------------------------------------------------------
    origin(:) = me%origin
    ! Normal was normalized in load routine so multiply thickness to convert
    ! normal to user defined normal
    normal(:) = me%normal * me%thickness

    do i = 1,n
      !barycentric coordinate
      coord = tem_BaryOfId( tree, treeIds(i) )
      vec (:) = coord(:) - origin(:)
      proj_len = (vec(1)*normal(1)+ vec(2)*normal(2)+vec(3)*normal(3))/   &
                 (normal(1)**2 + normal(2)**2 + normal(3)**2)
      sigma = me%dampFactor*((proj_len)**me%dampExponent)
      if (proj_len > 0 .and. proj_len < 1) then
        res(i) = sigma
      else if (proj_len > 1) then
        res(i) = me%dampFactor
      else
        res(i) = 0.0_rk
      end if

    end do

  end function spongeLayer_plane_expon_for_treeIDs
  ! -------------------------------------------------------------------------- !


  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the spongelayer using
  !! polynomial order 5 from treeids.
  !! Sponge profile:
  !! $$ \sigma(x) = \sigma_A*3125*(L+x0-x)*(x-x0)^4)/(256*(L)^5 $$
  !! where, \sigma_A - sponge strength, L - thickness, x0 - start of sponge.
  !!
  !! Profile is taken from:
  !! Xu, Hui; Sagaut, Pierre (2013): Analysis of the absorbing layers for the 
  !! weakly-compressible lattice Boltzmann methods. In Journal of Computational 
  !! Physics 245, pp. 14-42. DOI: 10.1016/j.jcp.2013.02.051.
  function spongeLayer_plane_polyn5_for_treeids(me, treeids, tree, n) &
    & result(res)
    ! --------------------------------------------------------------------------
    !> Spatial absorb layer to evaluate
    type(tem_spongeLayer_plane_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(3), normal(3), vec1(3), vec2(3), coord(3)
    real(kind=rk) :: proj_len1, proj_len2, const_fac
    ! --------------------------------------------------------------------------
    origin(:) = me%origin
    normal(:) = me%normal
    const_fac = 3125_rk/(256_rk*me%thickness**5)

    do i = 1,n
      !barycentric coordinate
      coord = tem_BaryOfId( tree, treeIds(i) )
      vec1(:) = coord(:) - origin(:)
      vec2(:) = me%thickness*normal(:) + origin(:) - coord(:)
      proj_len1 = vec1(1)*normal(1)+ vec1(2)*normal(2)+vec1(3)*normal(3)
      proj_len2 = vec2(1)*normal(1)+ vec2(2)*normal(2)+vec2(3)*normal(3)

      if (proj_len1 > 0 .and. proj_len2 > 0) then
        sigma = const_fac * proj_len2 * (proj_len1**4)
        res(i) = sigma*me%dampFactor
      else if (proj_len2 < 0) then ! If coord is beyond thickness
        res(i) = me%dampFactor
      else
        res(i) = 0.0_rk
      end if

    end do

  end function spongeLayer_plane_polyn5_for_treeids
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the spongeLayer from treeid for
  !! the polynomial order n6.
  !! Sponge profile:
  !! $$ \sigma(x) = \sigma_A*729*(L+x0-x)^2*(x-x0)^4)/(16*(L)^6 $$
  !! where, \sigma_A - sponge strength, L - thickness, x0 - start of sponge.
  function spongeLayer_plane_polyn6_for_treeids(me, treeids, tree, n) &
    & result(res)
    ! --------------------------------------------------------------------------
    !> Spatial sponge layer to evaluate
    type(tem_spongeLayer_plane_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(3), normal(3), vec1(3), vec2(3), coord(3)
    real(kind=rk) :: proj_len1, proj_len2, const_fac
    ! --------------------------------------------------------------------------
    origin(:) = me%origin
    normal(:) = me%normal
    const_fac = 729_rk/(16_rk*me%thickness**6)

    do i = 1,n
      !barycentric coordinate
      coord = tem_BaryOfId( tree, treeIds(i) )
      vec1(:) = coord(:) - origin(:)
      vec2(:) = me%thickness*normal(:) + origin(:) - coord(:)
      proj_len1 = vec1(1)*normal(1)+ vec1(2)*normal(2)+vec1(3)*normal(3)
      proj_len2 = vec2(1)*normal(1)+ vec2(2)*normal(2)+vec2(3)*normal(3)

      if (proj_len1 > 0) then
        sigma = const_fac * proj_len2**2 * (proj_len1**4)
        res(i) = min(1.0_rk, sigma) * me%dampFactor
      else
        res(i) = 0.0_rk
      end if

    end do

  end function spongeLayer_plane_polyn6_for_treeids
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the spongelayer and fills up
  !!  the res with the target state
  function spongelayer_plane_vector_for_treeIDs(me, ncomp, treeids, tree, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_plane_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> Number of entrys in each array
    integer, intent(in) :: ncomp
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n,ncomp)
    ! --------------------------------------------------------------------------
    integer :: i
    ! --------------------------------------------------------------------------
    res(:, 1) = spongeLayer_plane_scalar_for_treeIDs(me, treeids, tree, n)

    if (ncomp > 1) then
      do i = 1,n
        res(i,2:) = me%targetState(:)
      end do
    end if

  end function spongelayer_plane_vector_for_treeIDs
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !
  !                            SPONGE LAYER BOX
  !
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function returns the sigma for the box shape spongelayer
  function spongelayer_box_scalar_for_coord(me, coord, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    select case (trim(me%dampProfile))
    case ('linear', 'exponential')
      res(:) = spongeLayer_box_expon_for_coord(me, coord, n)
    case ('polynomial_n5')
      if (me%rounded_corner) then
        res(:) = spongeLayer_box_roundCornerPolyn5_for_coord(me, coord, n)
      else
        res(:) = spongeLayer_box_sharpCornerPolyn5_for_coord(me, coord, n)
      end if
    case ('polynomial_n6')
      if (me%rounded_corner) then
        res(:) = spongeLayer_box_roundCornerPolyn6_for_coord(me, coord, n)
      else
        res(:) = spongeLayer_box_sharpCornerPolyn6_for_coord(me, coord, n)
      end if
    end select

  end function spongelayer_box_scalar_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function returns sigma for the box shape spongelayer for coord for
  !! exponential profile.
  function spongelayer_box_expon_for_coord(me, coord, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(3), extent(3), box_max(3), proj_len
    real(kind=rk) :: coordLoc(3), normal, vec_min(3), vec_max(3)
    real(kind=rk) :: rad, vec_minSqr(3), vec_maxSqr(3)
    ! --------------------------------------------------------------------------
    origin(:) = me%origin
    extent(:) = me%extent
    box_max(:) = origin(:) + extent(:)

    do i = 1,n
      coordLoc = coord(i,:)
      vec_min(:) = coordLoc(:) - origin(:)
      vec_max(:) = coordLoc(:) - box_max(:)
      vec_minSqr(:) = vec_min(:)**2
      vec_maxSqr(:) = vec_max(:)**2
      normal = 1.0_rk
      rad = 0.0_rk

      ! Bottom-South-West: -x,-y,-z
      if (all(coordLoc(:) < origin(:))) then
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(2), vec_minSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) < origin(2) &
        & .and. coordLoc(3) < origin(3)) then ! Bottom-South-East: x, -y, -z
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(2), vec_minSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) > box_max(2) &
        & .and. coordLoc(3) < origin(3)) then ! Bottom-North-West: -x, y, -z
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(2), vec_minSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) > box_max(2) &
        & .and. coordLoc(3) < origin(3)) then ! Bottom-North-East: x, y, -z
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(2), vec_minSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) < origin(2) &
        & .and. coordLoc(3) > box_max(3)) then ! Top-South-West: -x, -y, z
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) < origin(2) &
        & .and. coordLoc(3) > box_max(3)) then ! Top-South-East: x, -y, z
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) > box_max(2) &
        & .and. coordLoc(3) > box_max(3)) then ! Top-North-West: -x, y, z
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(2), vec_maxSqr(3)) )

      else if (all(coordLoc(:) > box_max(:))) then ! Bottom-North-East: x, y, z
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(2) < origin(2) .and. coordLoc(3) < origin(3)) then
        ! Botton-South: -y,-z
        rad = sqrt( max(vec_minSqr(2), vec_minSqr(3)) )

      else if (coordLoc(2) < origin(2) .and. coordLoc(3) > box_max(3)) then
        ! Top-South: -y, z
        rad = sqrt( max(vec_minSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(2) > box_max(2) .and. coordLoc(3) < origin(3)) then
        ! Botom-North: -y, z
        rad = sqrt( max(vec_maxSqr(2), vec_minSqr(3)) )

      else if (coordLoc(2) > box_max(2) .and. coordLoc(3) > box_max(3)) then
        ! Top-North: y, z
        rad = sqrt( max(vec_maxSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(3) < origin(3)) then
        ! Botton-West: -x,-z
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(3) < origin(3)) then
        ! Bottom-East: x,-z
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(3) > box_max(3)) then
        ! Top-West: -x, z
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(3) > box_max(3)) then
        ! Top-East: x, z
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) < origin(2)) then
        ! South-West: -x,-y
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(2)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) > box_max(2)) then
        ! North-West: -x, y
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(2)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) < origin(2)) then
        ! South-East: x, -y
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(2)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) > box_max(2)) then
        ! North-East: x, y
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(2)) )

      else if (coordLoc(1) < origin(1)) then ! West: -x
        normal = -1_rk
        rad = vec_min(1)

      else if (coordLoc(2) < origin(2)) then ! South: -y
        normal = -1_rk
        rad = vec_min(2)

      else if (coordLoc(3) < origin(3)) then ! Bottom: -z
        normal = -1_rk
        rad = vec_min(3)

      else if (coordLoc(1) > box_max(1)) then ! East: x
        normal = 1_rk
        rad = vec_max(1)

      else if (coordLoc(2) > box_max(2)) then ! North: y
        normal = 1_rk
        rad = vec_max(2)

      else if (coordLoc(3) > box_max(3)) then ! Top: z
        normal = 1_rk
        rad = vec_max(3)
      end if

      proj_len = rad*normal/me%thickness
      sigma = me%dampFactor*((proj_len)**me%dampExponent)
      if (proj_len > 0 .and. proj_len < 1) then
        res(i) = sigma
      else if (proj_len > 1) then
        res(i) = me%dampFactor
      else
        res(i) = 0.0_rk
      end if

    end do

  end function spongelayer_box_expon_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the box shape sponge layer from 
  !! coord for polynomial n5.
  !! Sponge profile:
  !! $$ \sigma(x) = \sigma_A*3125*(L+x0-x)*(x-x0)^4)/(256*(L)^5 $$
  !! where, \sigma_A - sponge strength, L - thickness, x0 - start of sponge.
  !!
  !! Profile is taken from:
  !! Xu, Hui; Sagaut, Pierre (2013): Analysis of the absorbing layers for the 
  !! weakly-compressible lattice Boltzmann methods. In Journal of Computational 
  !! Physics 245, pp. 14-42. DOI: 10.1016/j.jcp.2013.02.051.
  function spongeLayer_box_sharpCornerPolyn5_for_coord(me, coord, n) &
    & result(res)
    ! --------------------------------------------------------------------------
    !> Spatial sponge layer to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(3), extent(3), coordLoc(3), normal
    real(kind=rk) :: proj_len1, proj_len2
    real(kind=rk) :: vec_min(3), vec_max(3), vec_minSqr(3), vec_maxSqr(3) 
    real(kind=rk) :: rad, const_fac, box_max(3)
    ! --------------------------------------------------------------------------
    origin(:) = me%origin
    extent(:) = me%extent
    box_max(:) = origin(:) + extent(:)

    const_fac = 3125_rk/(256_rk*me%thickness**5)

    do i = 1,n
      coordLoc = coord(i,:)
      vec_min(:) = coordLoc(:) - origin(:)
      vec_max(:) = coordLoc(:) - box_max(:)
      vec_minSqr(:) = vec_min(:)**2
      vec_maxSqr(:) = vec_max(:)**2
      proj_len1 = 0_rk
      proj_len2 = 0_rk
      normal = 1_rk
      rad = 0_rk

      ! Bottom-South-West: -x,-y,-z
      if (all(coordLoc(:) < origin(:))) then
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(2), vec_minSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) < origin(2) &
        & .and. coordLoc(3) < origin(3)) then ! Bottom-South-East: x, -y, -z
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(2), vec_minSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) > box_max(2) &
        & .and. coordLoc(3) < origin(3)) then ! Bottom-North-West: -x, y, -z
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(2), vec_minSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) > box_max(2) &
        & .and. coordLoc(3) < origin(3)) then ! Bottom-North-East: x, y, -z
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(2), vec_minSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) < origin(2) &
        & .and. coordLoc(3) > box_max(3)) then ! Top-South-West: -x, -y, z
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) < origin(2) &
        & .and. coordLoc(3) > box_max(3)) then ! Top-South-East: x, -y, z
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) > box_max(2) &
        & .and. coordLoc(3) > box_max(3)) then ! Top-North-West: -x, y, z
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(2), vec_maxSqr(3)) )

      else if (all(coordLoc(:) > box_max(:))) then ! Top-North-East: x, y, z
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(2) < origin(2) .and. coordLoc(3) < origin(3)) then
        ! Botton-South: -y,-z
        rad = sqrt( max(vec_minSqr(2), vec_minSqr(3)) )

      else if (coordLoc(2) < origin(2) .and. coordLoc(3) > box_max(3)) then
        ! Top-South: -y, z
        rad = sqrt( max(vec_minSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(2) > box_max(2) .and. coordLoc(3) < origin(3)) then
        ! Botom-North: -y, z
        rad = sqrt( max(vec_maxSqr(2), vec_minSqr(3)) )

      else if (coordLoc(2) > box_max(2) .and. coordLoc(3) > box_max(3)) then
        ! Top-North: y, z
        rad = sqrt( max(vec_maxSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(3) < origin(3)) then
        ! Botton-West: -x,-z
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(3) < origin(3)) then
        ! Bottom-East: x,-z
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(3) > box_max(3)) then
        ! Top-West: -x, z
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(3) > box_max(3)) then
        ! Top-East: x, z
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) < origin(2)) then
        ! South-West: -x,-y
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(2)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) > box_max(2)) then
        ! North-West: -x, y
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(2)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) < origin(2)) then
        ! South-East: x, -y
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(2)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) > box_max(2)) then
        ! North-East: x, y
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(2)) )

      else if (coordLoc(1) < origin(1)) then ! West: -x
        normal = -1_rk
        rad = vec_min(1)

      else if (coordLoc(2) < origin(2)) then ! South: -y
        normal = -1_rk
        rad = vec_min(2)

      else if (coordLoc(3) < origin(3)) then ! Bottom: -z
        normal = -1_rk
        rad = vec_min(3)

      else if (coordLoc(1) > box_max(1)) then ! East: x
        normal = 1_rk
        rad = vec_max(1)

      else if (coordLoc(2) > box_max(2)) then ! North: y
        normal = 1_rk
        rad = vec_max(2)

      else if (coordLoc(3) > box_max(3)) then ! Top: z
        normal = 1_rk
        rad = vec_max(3)
      end if

      proj_len1 = rad*normal
      proj_len2 = (me%thickness*normal - rad)*normal

      if (proj_len1 > 0 .and. proj_len2 > 0) then
        sigma = const_fac * proj_len2 * (proj_len1**4)
        res(i) = sigma*me%dampFactor
      else if (proj_len2 < 0) then ! If coord is beyond thickness
        res(i) = me%dampFactor
      else
        res(i) = 0.0_rk
      end if
    end do

  end function spongeLayer_box_sharpCornerPolyn5_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the box shape sponge layer from 
  !! coord for polynomial n5.
  !! Sponge profile:
  !! $$ \sigma(x) = \sigma_A*3125*(L+x0-x)*(x-x0)^4)/(256*(L)^5 $$
  !! where, \sigma_A - sponge strength, L - thickness, x0 - start of sponge.
  !!
  !! Profile is taken from:
  !! Xu, Hui; Sagaut, Pierre (2013): Analysis of the absorbing layers for the 
  !! weakly-compressible lattice Boltzmann methods. In Journal of Computational 
  !! Physics 245, pp. 14-42. DOI: 10.1016/j.jcp.2013.02.051.
  function spongeLayer_box_roundCornerPolyn5_for_coord(me, coord, n) &
    & result(res)
    ! --------------------------------------------------------------------------
    !> Spatial sponge layer to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(3), extent(3), coordLoc(3), normal
    real(kind=rk) :: proj_len1, proj_len2
    real(kind=rk) :: vec_corMin(3), vec_corMax(3)
    real(kind=rk) :: vec_corMinSqr(3), vec_corMaxSqr(3) 
    real(kind=rk) :: rad, const_fac, box_max(3)
    real(kind=rk) :: corInRad, corOutRad, corMin(3), corMax(3)
    logical :: isCorner
    ! --------------------------------------------------------------------------
    origin(:) = me%origin
    extent(:) = me%extent
    box_max(:) = origin(:) + extent(:)
    corInRad = me%corner_radius
    corOutRad = corInRad + me%thickness
    corMin(:) = origin(:) + corInRad
    corMax(:) = box_max(:) - corInRad

    const_fac = 3125_rk/(256_rk*me%thickness**5)

    do i = 1,n
      coordLoc = coord(i,:)
      vec_corMin(:) = coordLoc(:) - corMin(:)
      vec_corMax(:) = coordLoc(:) - corMax(:)
      vec_corMinSqr(:) = vec_corMin(:)**2
      vec_corMaxSqr(:) = vec_corMax(:)**2
      proj_len1 = 0_rk
      proj_len2 = 0_rk
      normal = 1_rk
      rad = 0_rk
      isCorner = .true.

      ! Bottom-South-West: -x,-y,-z
      if (all(coordLoc(:) < corMin(:))) then
        rad = sqrt( vec_corMinSqr(1) + vec_corMinSqr(2) + vec_corMinSqr(3))

      else if (coordLoc(1) > corMax(1) .and. coordLoc(2) < corMin(2) &
        & .and. coordLoc(3) < corMin(3)) then ! Bottom-South-East: x, -y, -z
        rad = sqrt( vec_corMaxSqr(1) + vec_corMinSqr(2) + vec_corMinSqr(3))

      else if (coordLoc(1) < corMin(1) .and. coordLoc(2) > corMax(2) &
        & .and. coordLoc(3) < corMin(3)) then ! Bottom-North-West: -x, y, -z
        rad = sqrt( vec_corMinSqr(1) + vec_corMaxSqr(2) + vec_corMinSqr(3))

      else if (coordLoc(1) > corMax(1) .and. coordLoc(2) > corMax(2) &
        & .and. coordLoc(3) < corMin(3)) then ! Bottom-North-East: x, y, -z
        rad = sqrt( vec_corMaxSqr(1) + vec_corMaxSqr(2) + vec_corMinSqr(3))

      else if (coordLoc(1) < corMin(1) .and. coordLoc(2) < corMin(2) &
        & .and. coordLoc(3) > corMax(3)) then ! Top-South-West: -x, -y, z
        rad = sqrt( vec_corMinSqr(1) + vec_corMinSqr(2) + vec_corMaxSqr(3))

      else if (coordLoc(1) > corMax(1) .and. coordLoc(2) < corMin(2) &
        & .and. coordLoc(3) > corMax(3)) then ! Top-South-East: x, -y, z
        rad = sqrt( vec_corMaxSqr(1) + vec_corMinSqr(2) + vec_corMaxSqr(3))

      else if (coordLoc(1) < corMin(1) .and. coordLoc(2) > corMax(2) &
        & .and. coordLoc(3) > corMax(3)) then ! Top-North-West: -x, y, z
        rad = sqrt( vec_corMinSqr(1) + vec_corMaxSqr(2) + vec_corMaxSqr(3))

      else if (all(coordLoc(:) > corMax(:))) then ! Top-North-East: x, y, z
        rad = sqrt( vec_corMaxSqr(1) + vec_corMaxSqr(2) + vec_corMaxSqr(3))

      else if (coordLoc(2) < corMin(2) .and. coordLoc(3) < corMin(3)) then
        ! Botton-South: -y,-z
        rad = sqrt( vec_corMinSqr(2) + vec_corMinSqr(3))

      else if (coordLoc(2) < corMin(2) .and. coordLoc(3) > corMax(3)) then
        ! Top-South: -y, z
        rad = sqrt( vec_corMinSqr(2) + vec_corMaxSqr(3))

      else if (coordLoc(2) > corMax(2) .and. coordLoc(3) < corMin(3)) then
        ! Bottom-North: y, -z
        rad = sqrt( vec_corMaxSqr(2) + vec_corMinSqr(3))

      else if (coordLoc(2) > corMax(2) .and. coordLoc(3) > corMax(3)) then
        ! Top-North: y, z
        rad = sqrt( vec_corMaxSqr(2) + vec_corMaxSqr(3))

      else if (coordLoc(1) < corMin(1) .and. coordLoc(3) < corMin(3)) then
        ! Botton-West: -x,-z
        rad = sqrt( vec_corMinSqr(1) + vec_corMinSqr(3))

      else if (coordLoc(1) > corMax(1) .and. coordLoc(3) < corMin(3)) then
        ! Bottom-East: x,-z
        rad = sqrt( vec_corMaxSqr(1) + vec_corMinSqr(3))

      else if (coordLoc(1) < corMin(1) .and. coordLoc(3) > corMax(3)) then
        ! Top-West: -x, z
        rad = sqrt( vec_corMinSqr(1) + vec_corMaxSqr(3))

      else if (coordLoc(1) > corMax(1) .and. coordLoc(3) > corMax(3)) then
        ! Top-East: x, z
        rad = sqrt( vec_corMaxSqr(1) + vec_corMaxSqr(3))

      else if (coordLoc(1) < corMin(1) .and. coordLoc(2) < corMin(2)) then
        ! South-West: -x,-y
        rad = sqrt(vec_corMinSqr(1) + vec_corMinSqr(2))

      else if (coordLoc(1) < corMin(1) .and. coordLoc(2) > corMax(2)) then
        ! North-West: -x, y
        rad = sqrt(vec_corMinSqr(1) + vec_corMaxSqr(2))

      else if (coordLoc(1) > corMax(1) .and. coordLoc(2) < corMin(2)) then
        ! South-East: x, -y
        rad = sqrt(vec_corMaxSqr(1) + vec_corMinSqr(2))

      else if (coordLoc(1) > corMax(1) .and. coordLoc(2) > corMax(2)) then
        ! North-East: x, y
        rad = sqrt(vec_corMaxSqr(1) + vec_corMaxSqr(2))

      else
        isCorner = .false.
        if (coordLoc(1) < origin(1)) then ! West: -x
          normal = -1_rk
          rad = coordLoc(1) - origin(1)

        else if (coordLoc(2) < origin(2)) then ! South: -y
          normal = -1_rk
          rad = coordLoc(2) - origin(2)

        else if (coordLoc(3) < origin(3)) then ! Bottom: -z
          normal = -1_rk
          rad = coordLoc(3) - origin(3)

        else if (coordLoc(1) > box_max(1)) then ! East: x
          normal = 1_rk
          rad = coordLoc(1) - box_max(1)

        else if (coordLoc(2) > box_max(2)) then ! North: y
          normal = 1_rk
          rad = coordLoc(2) - box_max(2)

        else if (coordLoc(3) > box_max(3)) then ! Top: z
          normal = 1_rk
          rad = coordLoc(3) - box_max(3)
        end if
      end if

      if (isCorner) then
        proj_len1 = rad - corInRad
        proj_len2 = corOutRad - rad
      else
        proj_len1 = rad*normal
        proj_len2 = (me%thickness*normal - rad)*normal
      end if

      if (proj_len1 > 0 .and. proj_len2 > 0) then
        sigma = const_fac * proj_len2 * (proj_len1**4)
        res(i) = sigma*me%dampFactor
      else if (proj_len2 < 0) then ! If coord is beyond thickness
        res(i) = me%dampFactor
      else
        res(i) = 0.0_rk
      end if
    end do

  end function spongeLayer_box_roundCornerPolyn5_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the box shape sponge layer from 
  !! coord for polynomial n6.
  !! Sponge profile:
  !! $$ \sigma(x) = \sigma_A*729*(L+x0-x)^2*(x-x0)^4)/(16*(L)^6 $$
  !! where, \sigma_A - sponge strength, L - thickness, x0 - start of sponge.
  function spongeLayer_box_sharpCornerPolyn6_for_coord(me, coord, n) &
    & result(res)
    ! --------------------------------------------------------------------------
    !> Spatial sponge layer to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(3), extent(3), coordLoc(3), normal
    real(kind=rk) :: proj_len1, proj_len2
    real(kind=rk) :: vec_min(3), vec_max(3), vec_minSqr(3), vec_maxSqr(3) 
    real(kind=rk) :: rad, const_fac, box_max(3)
    ! --------------------------------------------------------------------------
    origin(:) = me%origin
    extent(:) = me%extent
    box_max(:) = origin(:) + extent(:)

    const_fac = 729_rk/(16_rk*me%thickness**6)

    do i = 1,n
      coordLoc = coord(i,:)
      vec_min(:) = coordLoc(:) - origin(:)
      vec_max(:) = coordLoc(:) - box_max(:)
      vec_minSqr(:) = vec_min(:)**2
      vec_maxSqr(:) = vec_max(:)**2
      proj_len1 = 0_rk
      proj_len2 = 0_rk
      normal = 1_rk
      rad = 0_rk

      ! Bottom-South-West: -x,-y,-z
      if (all(coordLoc(:) < origin(:))) then
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(2), vec_minSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) < origin(2) &
        & .and. coordLoc(3) < origin(3)) then ! Bottom-South-East: x, -y, -z
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(2), vec_minSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) > box_max(2) &
        & .and. coordLoc(3) < origin(3)) then ! Bottom-North-West: -x, y, -z
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(2), vec_minSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) > box_max(2) &
        & .and. coordLoc(3) < origin(3)) then ! Bottom-North-East: x, y, -z
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(2), vec_minSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) < origin(2) &
        & .and. coordLoc(3) > box_max(3)) then ! Top-South-West: -x, -y, z
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) < origin(2) &
        & .and. coordLoc(3) > box_max(3)) then ! Top-South-East: x, -y, z
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) > box_max(2) &
        & .and. coordLoc(3) > box_max(3)) then ! Top-North-West: -x, y, z
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(2), vec_maxSqr(3)) )

      else if (all(coordLoc(:) > box_max(:))) then ! Top-North-East: x, y, z
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(2) < origin(2) .and. coordLoc(3) < origin(3)) then
        ! Botton-South: -y,-z
        rad = sqrt( max(vec_minSqr(2), vec_minSqr(3)) )

      else if (coordLoc(2) < origin(2) .and. coordLoc(3) > box_max(3)) then
        ! Top-South: -y, z
        rad = sqrt( max(vec_minSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(2) > box_max(2) .and. coordLoc(3) < origin(3)) then
        ! Botom-North: -y, z
        rad = sqrt( max(vec_maxSqr(2), vec_minSqr(3)) )

      else if (coordLoc(2) > box_max(2) .and. coordLoc(3) > box_max(3)) then
        ! Top-North: y, z
        rad = sqrt( max(vec_maxSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(3) < origin(3)) then
        ! Botton-West: -x,-z
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(3) < origin(3)) then
        ! Bottom-East: x,-z
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(3) > box_max(3)) then
        ! Top-West: -x, z
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(3) > box_max(3)) then
        ! Top-East: x, z
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) < origin(2)) then
        ! South-West: -x,-y
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(2)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) > box_max(2)) then
        ! North-West: -x, y
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(2)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) < origin(2)) then
        ! South-East: x, -y
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(2)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) > box_max(2)) then
        ! North-East: x, y
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(2)) )

      else if (coordLoc(1) < origin(1)) then ! West: -x
        normal = -1_rk
        rad = vec_min(1)

      else if (coordLoc(2) < origin(2)) then ! South: -y
        normal = -1_rk
        rad = vec_min(2)

      else if (coordLoc(3) < origin(3)) then ! Bottom: -z
        normal = -1_rk
        rad = vec_min(3)

      else if (coordLoc(1) > box_max(1)) then ! East: x
        normal = 1_rk
        rad = vec_max(1)

      else if (coordLoc(2) > box_max(2)) then ! North: y
        normal = 1_rk
        rad = vec_max(2)

      else if (coordLoc(3) > box_max(3)) then ! Top: z
        normal = 1_rk
        rad = vec_max(3)
      end if

      proj_len1 = rad*normal
      proj_len2 = (me%thickness*normal - rad)*normal

      if (proj_len1 > 0) then
        sigma = const_fac * proj_len2**2 * (proj_len1**4)
        res(i) = min(1.0_rk, sigma) * me%dampFactor
      else
        res(i) = 0.0_rk
      end if
    end do

  end function spongeLayer_box_sharpCornerPolyn6_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the box shape sponge layer from 
  !! coord for polynomial n6.
  !! Sponge profile:
  !! $$ \sigma(x) = \sigma_A*729*(L+x0-x)^2*(x-x0)^4)/(16*(L)^6 $$
  !! where, \sigma_A - sponge strength, L - thickness, x0 - start of sponge.
  function spongeLayer_box_roundCornerPolyn6_for_coord(me, coord, n) &
    & result(res)
    ! --------------------------------------------------------------------------
    !> Spatial sponge layer to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(3), extent(3), coordLoc(3), normal
    real(kind=rk) :: proj_len1, proj_len2
    real(kind=rk) :: vec_corMin(3), vec_corMax(3)
    real(kind=rk) :: vec_corMinSqr(3), vec_corMaxSqr(3) 
    real(kind=rk) :: rad, const_fac, box_max(3)
    real(kind=rk) :: corInRad, corOutRad, corMin(3), corMax(3)
    logical :: isCorner
    ! --------------------------------------------------------------------------
    origin(:) = me%origin
    extent(:) = me%extent
    box_max(:) = origin(:) + extent(:)
    corInRad = me%corner_radius
    corOutRad = corInRad + me%thickness
    corMin(:) = origin(:) + corInRad
    corMax(:) = box_max(:) - corInRad

    const_fac = 729_rk/(16_rk*me%thickness**6)

    do i = 1,n
      coordLoc = coord(i,:)
      vec_corMin(:) = coordLoc(:) - corMin(:)
      vec_corMax(:) = coordLoc(:) - corMax(:)
      vec_corMinSqr(:) = vec_corMin(:)**2
      vec_corMaxSqr(:) = vec_corMax(:)**2
      proj_len1 = 0_rk
      proj_len2 = 0_rk
      normal = 1_rk
      rad = 0_rk
      isCorner = .true.

      ! Bottom-South-West: -x,-y,-z
      if (all(coordLoc(:) < corMin(:))) then
        rad = sqrt( vec_corMinSqr(1) + vec_corMinSqr(2) + vec_corMinSqr(3))

      else if (coordLoc(1) > corMax(1) .and. coordLoc(2) < corMin(2) &
        & .and. coordLoc(3) < corMin(3)) then ! Bottom-South-East: x, -y, -z
        rad = sqrt( vec_corMaxSqr(1) + vec_corMinSqr(2) + vec_corMinSqr(3))

      else if (coordLoc(1) < corMin(1) .and. coordLoc(2) > corMax(2) &
        & .and. coordLoc(3) < corMin(3)) then ! Bottom-North-West: -x, y, -z
        rad = sqrt( vec_corMinSqr(1) + vec_corMaxSqr(2) + vec_corMinSqr(3))

      else if (coordLoc(1) > corMax(1) .and. coordLoc(2) > corMax(2) &
        & .and. coordLoc(3) < corMin(3)) then ! Bottom-North-East: x, y, -z
        rad = sqrt( vec_corMaxSqr(1) + vec_corMaxSqr(2) + vec_corMinSqr(3))

      else if (coordLoc(1) < corMin(1) .and. coordLoc(2) < corMin(2) &
        & .and. coordLoc(3) > corMax(3)) then ! Top-South-West: -x, -y, z
        rad = sqrt( vec_corMinSqr(1) + vec_corMinSqr(2) + vec_corMaxSqr(3))

      else if (coordLoc(1) > corMax(1) .and. coordLoc(2) < corMin(2) &
        & .and. coordLoc(3) > corMax(3)) then ! Top-South-East: x, -y, z
        rad = sqrt( vec_corMaxSqr(1) + vec_corMinSqr(2) + vec_corMaxSqr(3))

      else if (coordLoc(1) < corMin(1) .and. coordLoc(2) > corMax(2) &
        & .and. coordLoc(3) > corMax(3)) then ! Top-North-West: -x, y, z
        rad = sqrt( vec_corMinSqr(1) + vec_corMaxSqr(2) + vec_corMaxSqr(3))

      else if (all(coordLoc(:) > corMax(:))) then ! Top-North-East: x, y, z
        rad = sqrt( vec_corMaxSqr(1) + vec_corMaxSqr(2) + vec_corMaxSqr(3))

      else if (coordLoc(2) < corMin(2) .and. coordLoc(3) < corMin(3)) then
        ! Botton-South: -y,-z
        rad = sqrt( vec_corMinSqr(2) + vec_corMinSqr(3))

      else if (coordLoc(2) < corMin(2) .and. coordLoc(3) > corMax(3)) then
        ! Top-South: -y, z
        rad = sqrt( vec_corMinSqr(2) + vec_corMaxSqr(3))

      else if (coordLoc(2) > corMax(2) .and. coordLoc(3) < corMin(3)) then
        ! Bottom-North: y, -z
        rad = sqrt( vec_corMaxSqr(2) + vec_corMinSqr(3))

      else if (coordLoc(2) > corMax(2) .and. coordLoc(3) > corMax(3)) then
        ! Top-North: y, z
        rad = sqrt( vec_corMaxSqr(2) + vec_corMaxSqr(3))

      else if (coordLoc(1) < corMin(1) .and. coordLoc(3) < corMin(3)) then
        ! Botton-West: -x,-z
        rad = sqrt( vec_corMinSqr(1) + vec_corMinSqr(3))

      else if (coordLoc(1) > corMax(1) .and. coordLoc(3) < corMin(3)) then
        ! Bottom-East: x,-z
        rad = sqrt( vec_corMaxSqr(1) + vec_corMinSqr(3))

      else if (coordLoc(1) < corMin(1) .and. coordLoc(3) > corMax(3)) then
        ! Top-West: -x, z
        rad = sqrt( vec_corMinSqr(1) + vec_corMaxSqr(3))

      else if (coordLoc(1) > corMax(1) .and. coordLoc(3) > corMax(3)) then
        ! Top-East: x, z
        rad = sqrt( vec_corMaxSqr(1) + vec_corMaxSqr(3))

      else if (coordLoc(1) < corMin(1) .and. coordLoc(2) < corMin(2)) then
        ! South-West: -x,-y
        rad = sqrt(vec_corMinSqr(1) + vec_corMinSqr(2))

      else if (coordLoc(1) < corMin(1) .and. coordLoc(2) > corMax(2)) then
        ! North-West: -x, y
        rad = sqrt(vec_corMinSqr(1) + vec_corMaxSqr(2))

      else if (coordLoc(1) > corMax(1) .and. coordLoc(2) < corMin(2)) then
        ! South-East: x, -y
        rad = sqrt(vec_corMaxSqr(1) + vec_corMinSqr(2))

      else if (coordLoc(1) > corMax(1) .and. coordLoc(2) > corMax(2)) then
        ! North-East: x, y
        rad = sqrt(vec_corMaxSqr(1) + vec_corMaxSqr(2))

      else
        isCorner = .false.
        if (coordLoc(1) < origin(1)) then ! West: -x
          normal = -1_rk
          rad = coordLoc(1) - origin(1)

        else if (coordLoc(2) < origin(2)) then ! South: -y
          normal = -1_rk
          rad = coordLoc(2) - origin(2)

        else if (coordLoc(3) < origin(3)) then ! Bottom: -z
          normal = -1_rk
          rad = coordLoc(3) - origin(3)

        else if (coordLoc(1) > box_max(1)) then ! East: x
          normal = 1_rk
          rad = coordLoc(1) - box_max(1)

        else if (coordLoc(2) > box_max(2)) then ! North: y
          normal = 1_rk
          rad = coordLoc(2) - box_max(2)

        else if (coordLoc(3) > box_max(3)) then ! Top: z
          normal = 1_rk
          rad = coordLoc(3) - box_max(3)
        end if
      end if

      if (isCorner) then
        proj_len1 = rad - corInRad
        proj_len2 = corOutRad - rad
      else
        proj_len1 = rad*normal
        proj_len2 = (me%thickness*normal - rad)*normal
      end if

      if (proj_len1 > 0) then
        sigma = const_fac * proj_len2**2 * (proj_len1**4)
        res(i) = min(1.0_rk, sigma) * me%dampFactor
      else
        res(i) = 0.0_rk
      end if
    end do

  end function spongeLayer_box_roundCornerPolyn6_for_coord
  ! -------------------------------------------------------------------------- !


  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the spongelayer and fills up
  !!  the res with the target state
  function spongelayer_box_vector_for_coord(me, nComp, coord, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> Number of entrys in each array
    integer, intent(in) :: ncomp
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n,ncomp)
    ! --------------------------------------------------------------------------
    integer :: i
    ! --------------------------------------------------------------------------
    res(:, 1) = spongeLayer_box_scalar_for_coord(me, coord, n)

    if (ncomp > 1) then
      do i = 1,n
        res(i,2:) = me%targetState(:)
      end do
    end if

  end function spongelayer_box_vector_for_coord
  ! -------------------------------------------------------------------------- !


  ! -------------------------------------------------------------------------- !
  !> This function returns the sigma for the spongelayer from treeids
  function spongelayer_box_scalar_for_treeIDs(me, treeids, tree, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    select case (trim(me%dampProfile))
    case ('linear', 'exponential')
      res(:) = spongeLayer_box_expon_for_treeIDs(me, treeIDs, tree, n)
    case ('polynomial_n5')
      if (me%rounded_corner) then
        res(:) = spongeLayer_box_roundCornerPolyn5_for_treeIDs(me, treeIDs, &
          &                                                    tree, n)
      else
        res(:) = spongeLayer_box_sharpCornerPolyn5_for_treeIDs(me, treeIDs, &
          &                                                    tree, n)
      end if
    case ('polynomial_n6')
      if (me%rounded_corner) then
        res(:) = spongeLayer_box_roundCornerPolyn6_for_treeIDs(me, treeIDs, &
          &                                                    tree, n)
      else
        res(:) = spongeLayer_box_sharpCornerPolyn6_for_treeIDs(me, treeIDs, &
          &                                                    tree, n)
      end if
    end select

  end function spongeLayer_box_scalar_for_treeIDs
  ! -------------------------------------------------------------------------- !


  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the box shape spongelayer fom treeid
  !! for exponential profile
  function spongelayer_box_expon_for_treeIDs(me, treeIDs, tree, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: res_i(1), coord(3,1)
    ! --------------------------------------------------------------------------
    do i = 1,n
      !barycentric coordinate
      coord(:,1) = tem_BaryOfId( tree, treeIds(i) )
      res_i = spongeLayer_box_expon_for_coord(me, coord, 1)
      res(i) = res_i(1)
    end do

  end function spongelayer_box_expon_for_treeIDs
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the box shape spongelayer fom treeid
  !! for polynomial n5 profile
  function spongelayer_box_sharpCornerPolyn5_for_treeIDs(me, treeIDs, tree, n) &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: res_i(1), coord(3,1)
    ! --------------------------------------------------------------------------
    do i = 1,n
      !barycentric coordinate
      coord(:,1) = tem_BaryOfId( tree, treeIds(i) )
      res_i = spongeLayer_box_sharpCornerPolyn5_for_coord(me, coord, 1)
      res(i) = res_i(1)
    end do

  end function spongelayer_box_sharpCornerPolyn5_for_treeIDs
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the box shape spongelayer fom treeid
  !! for polynomial n5 profile
  function spongelayer_box_roundCornerPolyn5_for_treeIDs(me, treeIDs, tree, n) &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: res_i(1), coord(3,1)
    ! --------------------------------------------------------------------------
    do i = 1,n
      !barycentric coordinate
      coord(:,1) = tem_BaryOfId( tree, treeIds(i) )
      res_i = spongeLayer_box_roundCornerPolyn5_for_coord(me, coord, 1)
      res(i) = res_i(1)
    end do

  end function spongelayer_box_roundCornerPolyn5_for_treeIDs
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the box shape spongelayer fom treeid
  !! for polynomial n5 profile
  function spongelayer_box_sharpCornerPolyn6_for_treeIDs(me, treeIDs, tree, n) &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: res_i(1), coord(3,1)
    ! --------------------------------------------------------------------------
    do i = 1,n
      !barycentric coordinate
      coord(:,1) = tem_BaryOfId( tree, treeIds(i) )
      res_i = spongeLayer_box_sharpCornerPolyn6_for_coord(me, coord, 1)
      res(i) = res_i(1)
    end do

  end function spongelayer_box_sharpCornerPolyn6_for_treeIDs
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the box shape spongelayer fom treeid
  !! for polynomial n5 profile
  function spongelayer_box_roundCornerPolyn6_for_treeIDs(me, treeIDs, tree, n) &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: res_i(1), coord(3,1)
    ! --------------------------------------------------------------------------
    do i = 1,n
      !barycentric coordinate
      coord(:,1) = tem_BaryOfId( tree, treeIds(i) )
      res_i = spongeLayer_box_roundCornerPolyn6_for_coord(me, coord, 1)
      res(i) = res_i(1)
    end do

  end function spongelayer_box_roundCornerPolyn6_for_treeIDs
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the spongelayer and fills up
  !!  the res with the target state
  function spongelayer_box_vector_for_treeIDs(me, ncomp, treeids, tree, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> Number of entrys in each array
    integer, intent(in) :: ncomp
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n,ncomp)
    ! --------------------------------------------------------------------------
    integer :: i
    ! --------------------------------------------------------------------------
    res(:, 1) = spongeLayer_box_scalar_for_treeIDs(me, treeids, tree, n)

    if (ncomp > 1) then
      do i = 1,n
        res(i,2:) = me%targetState(:)
      end do
    end if

  end function spongelayer_box_vector_for_treeIDs
  ! -------------------------------------------------------------------------- !


  ! -------------------------------------------------------------------------- !
  !
  !                            SPONGE LAYER BOX 2D
  !
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function returns the sigma for the 2d box shape spongelayer
  function spongelayer_box2d_scalar_for_coord(me, coord, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    select case (trim(me%dampProfile))
    case ('linear', 'exponential')
      res(:) = spongeLayer_box2d_expon_for_coord(me, coord, n)
    case ('polynomial_n5')
      if (me%rounded_corner) then
        res(:) = spongeLayer_box2d_roundCornerPolyn5_for_coord(me, coord, n)
      else
        res(:) = spongeLayer_box2d_sharpCornerPolyn5_for_coord(me, coord, n)
      end if
    case ('polynomial_n6')
      if (me%rounded_corner) then
        res(:) = spongeLayer_box2d_roundCornerPolyn6_for_coord(me, coord, n)
      else
        res(:) = spongeLayer_box2d_sharpCornerPolyn6_for_coord(me, coord, n)
      end if
    end select

  end function spongelayer_box2d_scalar_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function returns sigma for the box shape spongelayer for coord for
  !! exponential profile.
  function spongelayer_box2d_expon_for_coord(me, coord, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(2), extent(2), box_max(2), proj_len
    real(kind=rk) :: coordLoc(2), normal, vec_min(2), vec_max(2)
    real(kind=rk) :: rad, vec_minSqr(2), vec_maxSqr(2)
    ! --------------------------------------------------------------------------
    origin(:) = me%origin(1:2)
    extent(:) = me%extent(1:2)
    box_max(:) = origin(:) + extent(:)

    do i = 1,n
      coordLoc = coord(i,1:2)
      vec_min(:) = coordLoc(:) - origin(:)
      vec_max(:) = coordLoc(:) - box_max(:)
      vec_minSqr(:) = vec_min(:)**2
      vec_maxSqr(:) = vec_max(:)**2
      normal = 1.0_rk
      rad = 0.0_rk

      if (coordLoc(1) < origin(1) .and. coordLoc(2) < origin(2)) then
        ! South-West: -x,-y
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(2)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) > box_max(2)) then
        ! North-West: -x, y
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(2)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) < origin(2)) then
        ! South-East: x, -y
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(2)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) > box_max(2)) then
        ! North-East: x, y
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(2)) )

      else if (coordLoc(1) < origin(1)) then ! West: -x
        normal = -1_rk
        rad = vec_min(1)

      else if (coordLoc(2) < origin(2)) then ! South: -y
        normal = -1_rk
        rad = vec_min(2)

      else if (coordLoc(1) > box_max(1)) then ! East: x
        normal = 1_rk
        rad = vec_max(1)

      else if (coordLoc(2) > box_max(2)) then ! North: y
        normal = 1_rk
        rad = vec_max(2)
      end if

      proj_len = rad*normal/me%thickness
      sigma = me%dampFactor*((proj_len)**me%dampExponent)
      if (proj_len > 0 .and. proj_len < 1) then
        res(i) = sigma
      else if (proj_len > 1) then
        res(i) = me%dampFactor
      else
        res(i) = 0.0_rk
      end if

    end do

  end function spongelayer_box2d_expon_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the 2d box shape sponge layer from 
  !! coord for polynomial n5.
  !! Sponge profile:
  !! $$ \sigma(x) = \sigma_A*3125*(L+x0-x)*(x-x0)^4)/(256*(L)^5 $$
  !! where, \sigma_A - sponge strength, L - thickness, x0 - start of sponge.
  !!
  !! Profile is taken from:
  !! Xu, Hui; Sagaut, Pierre (2013): Analysis of the absorbing layers for the 
  !! weakly-compressible lattice Boltzmann methods. In Journal of Computational 
  !! Physics 245, pp. 14-42. DOI: 10.1016/j.jcp.2013.02.051.
  function spongeLayer_box2d_sharpCornerPolyn5_for_coord(me, coord, n) &
    & result(res)
    ! --------------------------------------------------------------------------
    !> Spatial sponge layer to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(2), extent(2), coordLoc(2), normal
    real(kind=rk) :: proj_len1, proj_len2
    real(kind=rk) :: vec_min(2), vec_max(2), vec_minSqr(2), vec_maxSqr(2) 
    real(kind=rk) :: rad, const_fac, box_max(2)
    ! --------------------------------------------------------------------------
    origin(:) = me%origin(1:2)
    extent(:) = me%extent(1:2)
    box_max(:) = origin(:) + extent(:)

    const_fac = 3125_rk/(256_rk*me%thickness**5)

    do i = 1,n
      coordLoc = coord(i,1:2)
      vec_min(:) = coordLoc(:) - origin(:)
      vec_max(:) = coordLoc(:) - box_max(:)
      vec_minSqr(:) = vec_min(:)**2
      vec_maxSqr(:) = vec_max(:)**2
      proj_len1 = 0_rk
      proj_len2 = 0_rk
      normal = 1_rk
      rad = 0_rk

      if (coordLoc(1) < origin(1) .and. coordLoc(2) < origin(2)) then
        ! South-West: -x,-y
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(2)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) > box_max(2)) then
        ! North-West: -x, y
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(2)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) < origin(2)) then
        ! South-East: x, -y
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(2)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) > box_max(2)) then
        ! North-East: x, y
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(2)) )

      else if (coordLoc(1) < origin(1)) then ! West: -x
        normal = -1_rk
        rad = vec_min(1)

      else if (coordLoc(2) < origin(2)) then ! South: -y
        normal = -1_rk
        rad = vec_min(2)

      else if (coordLoc(1) > box_max(1)) then ! East: x
        normal = 1_rk
        rad = vec_max(1)

      else if (coordLoc(2) > box_max(2)) then ! North: y
        normal = 1_rk
        rad = vec_max(2)
      end if

      proj_len1 = rad*normal
      proj_len2 = (me%thickness*normal - rad)*normal

      if (proj_len1 > 0 .and. proj_len2 > 0) then
        sigma = const_fac * proj_len2 * (proj_len1**4)
        res(i) = sigma*me%dampFactor
      else if (proj_len2 < 0) then ! If coord is beyond thickness
        res(i) = me%dampFactor
      else
        res(i) = 0.0_rk
      end if
    end do

  end function spongeLayer_box2d_sharpCornerPolyn5_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the 2d box shape sponge layer from 
  !! coord for polynomial n5.
  !! Sponge profile:
  !! $$ \sigma(x) = \sigma_A*3125*(L+x0-x)*(x-x0)^4)/(256*(L)^5 $$
  !! where, \sigma_A - sponge strength, L - thickness, x0 - start of sponge.
  !!
  !! Profile is taken from:
  !! Xu, Hui; Sagaut, Pierre (2013): Analysis of the absorbing layers for the 
  !! weakly-compressible lattice Boltzmann methods. In Journal of Computational 
  !! Physics 245, pp. 14-42. DOI: 10.1016/j.jcp.2013.02.051.
  function spongeLayer_box2d_roundCornerPolyn5_for_coord(me, coord, n) &
    & result(res)
    ! --------------------------------------------------------------------------
    !> Spatial sponge layer to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(2), extent(2), coordLoc(2), normal
    real(kind=rk) :: proj_len1, proj_len2
    real(kind=rk) :: vec_corMin(2), vec_corMax(2)
    real(kind=rk) :: vec_corMinSqr(2), vec_corMaxSqr(2) 
    real(kind=rk) :: rad, const_fac, box_max(2)
    real(kind=rk) :: corInRad, corOutRad, corMin(2), corMax(2)
    logical :: isCorner
    ! --------------------------------------------------------------------------
    origin(:) = me%origin(1:2)
    extent(:) = me%extent(1:2)
    box_max(:) = origin(:) + extent(:)
    corInRad = me%corner_radius
    corOutRad = corInRad + me%thickness
    corMin(:) = origin(:) + corInRad
    corMax(:) = box_max(:) - corInRad

    const_fac = 3125_rk/(256_rk*me%thickness**5)

    do i = 1,n
      coordLoc = coord(i,1:2)
      vec_corMin(:) = coordLoc(:) - corMin(:)
      vec_corMax(:) = coordLoc(:) - corMax(:)
      vec_corMinSqr(:) = vec_corMin(:)**2
      vec_corMaxSqr(:) = vec_corMax(:)**2
      proj_len1 = 0_rk
      proj_len2 = 0_rk
      normal = 1_rk
      rad = 0_rk
      isCorner = .true.

      if (coordLoc(1) < corMin(1) .and. coordLoc(2) < corMin(2)) then
        ! South-West: -x,-y
        rad = sqrt(vec_corMinSqr(1) + vec_corMinSqr(2))

      else if (coordLoc(1) < corMin(1) .and. coordLoc(2) > corMax(2)) then
        ! North-West: -x, y
        rad = sqrt(vec_corMinSqr(1) + vec_corMaxSqr(2))

      else if (coordLoc(1) > corMax(1) .and. coordLoc(2) < corMin(2)) then
        ! South-East: x, -y
        rad = sqrt(vec_corMaxSqr(1) + vec_corMinSqr(2))

      else if (coordLoc(1) > corMax(1) .and. coordLoc(2) > corMax(2)) then
        ! North-East: x, y
        rad = sqrt(vec_corMaxSqr(1) + vec_corMaxSqr(2))

      else
        isCorner = .false.
        if (coordLoc(1) < origin(1)) then ! West: -x
          normal = -1_rk
          rad = coordLoc(1) - origin(1)

        else if (coordLoc(2) < origin(2)) then ! South: -y
          normal = -1_rk
          rad = coordLoc(2) - origin(2)

        else if (coordLoc(1) > box_max(1)) then ! East: x
          normal = 1_rk
          rad = coordLoc(1) - box_max(1)

        else if (coordLoc(2) > box_max(2)) then ! North: y
          normal = 1_rk
          rad = coordLoc(2) - box_max(2)
        end if
      end if

      if (isCorner) then
        proj_len1 = rad - corInRad
        proj_len2 = corOutRad - rad
      else
        proj_len1 = rad*normal
        proj_len2 = (me%thickness*normal - rad)*normal
      end if

      if (proj_len1 > 0 .and. proj_len2 > 0) then
        sigma = const_fac * proj_len2 * (proj_len1**4)
        res(i) = sigma*me%dampFactor
      else if (proj_len2 < 0) then ! If coord is beyond thickness
        res(i) = me%dampFactor
      else
        res(i) = 0.0_rk
      end if
    end do

  end function spongeLayer_box2d_roundCornerPolyn5_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the 2d box shape sponge layer from 
  !! coord for polynomial n6.
  !! Sponge profile:
  !! $$ \sigma(x) = \sigma_A*729*(L+x0-x)^2*(x-x0)^4)/(16*(L)^6 $$
  !! where, \sigma_A - sponge strength, L - thickness, x0 - start of sponge.
  function spongeLayer_box2d_sharpCornerPolyn6_for_coord(me, coord, n) &
    & result(res)
    ! --------------------------------------------------------------------------
    !> Spatial sponge layer to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(2), extent(2), coordLoc(2), normal
    real(kind=rk) :: proj_len1, proj_len2
    real(kind=rk) :: vec_min(2), vec_max(2), vec_minSqr(2), vec_maxSqr(2) 
    real(kind=rk) :: rad, const_fac, box_max(2)
    ! --------------------------------------------------------------------------
    origin(:) = me%origin(1:2)
    extent(:) = me%extent(1:2)
    box_max(:) = origin(:) + extent(:)

    const_fac = 729_rk/(16_rk*me%thickness**6)

    do i = 1,n
      coordLoc = coord(i,1:2)
      vec_min(:) = coordLoc(:) - origin(:)
      vec_max(:) = coordLoc(:) - box_max(:)
      vec_minSqr(:) = vec_min(:)**2
      vec_maxSqr(:) = vec_max(:)**2
      proj_len1 = 0_rk
      proj_len2 = 0_rk
      normal = 1_rk
      rad = 0_rk

      if (coordLoc(1) < origin(1) .and. coordLoc(2) < origin(2)) then
        ! South-West: -x,-y
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(2)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) > box_max(2)) then
        ! North-West: -x, y
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(2)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) < origin(2)) then
        ! South-East: x, -y
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(2)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) > box_max(2)) then
        ! North-East: x, y
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(2)) )

      else if (coordLoc(1) < origin(1)) then ! West: -x
        normal = -1_rk
        rad = vec_min(1)

      else if (coordLoc(2) < origin(2)) then ! South: -y
        normal = -1_rk
        rad = vec_min(2)

      else if (coordLoc(1) > box_max(1)) then ! East: x
        normal = 1_rk
        rad = vec_max(1)

      else if (coordLoc(2) > box_max(2)) then ! North: y
        normal = 1_rk
        rad = vec_max(2)
      end if

      proj_len1 = rad*normal
      proj_len2 = (me%thickness*normal - rad)*normal

      if (proj_len1 > 0) then
        sigma = const_fac * proj_len2**2 * (proj_len1**4)
        res(i) = min(1.0_rk, sigma) * me%dampFactor
      else
        res(i) = 0.0_rk
      end if
    end do

  end function spongeLayer_box2d_sharpCornerPolyn6_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the 2d box shape sponge layer from 
  !! coord for polynomial n6.
  !! Sponge profile:
  !! $$ \sigma(x) = \sigma_A*729*(L+x0-x)^2*(x-x0)^4)/(16*(L)^6 $$
  !! where, \sigma_A - sponge strength, L - thickness, x0 - start of sponge.
  function spongeLayer_box2d_roundCornerPolyn6_for_coord(me, coord, n) &
    & result(res)
    ! --------------------------------------------------------------------------
    !> Spatial sponge layer to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(2), extent(2), coordLoc(2), normal
    real(kind=rk) :: proj_len1, proj_len2
    real(kind=rk) :: vec_corMin(2), vec_corMax(2)
    real(kind=rk) :: vec_corMinSqr(2), vec_corMaxSqr(2) 
    real(kind=rk) :: rad, const_fac, box_max(2)
    real(kind=rk) :: corInRad, corOutRad, corMin(2), corMax(2)
    logical :: isCorner
    ! --------------------------------------------------------------------------
    origin(:) = me%origin(1:2)
    extent(:) = me%extent(1:2)
    box_max(:) = origin(:) + extent(:)
    corInRad = me%corner_radius
    corOutRad = corInRad + me%thickness
    corMin(:) = origin(:) + corInRad
    corMax(:) = box_max(:) - corInRad

    const_fac = 729_rk/(16_rk*me%thickness**6)

    do i = 1,n
      coordLoc = coord(i,1:2)
      vec_corMin(:) = coordLoc(:) - corMin(:)
      vec_corMax(:) = coordLoc(:) - corMax(:)
      vec_corMinSqr(:) = vec_corMin(:)**2
      vec_corMaxSqr(:) = vec_corMax(:)**2
      proj_len1 = 0_rk
      proj_len2 = 0_rk
      normal = 1_rk
      rad = 0_rk
      isCorner = .true.

      if (coordLoc(1) < corMin(1) .and. coordLoc(2) < corMin(2)) then
        ! South-West: -x,-y
        rad = sqrt(vec_corMinSqr(1) + vec_corMinSqr(2))

      else if (coordLoc(1) < corMin(1) .and. coordLoc(2) > corMax(2)) then
        ! North-West: -x, y
        rad = sqrt(vec_corMinSqr(1) + vec_corMaxSqr(2))

      else if (coordLoc(1) > corMax(1) .and. coordLoc(2) < corMin(2)) then
        ! South-East: x, -y
        rad = sqrt(vec_corMaxSqr(1) + vec_corMinSqr(2))

      else if (coordLoc(1) > corMax(1) .and. coordLoc(2) > corMax(2)) then
        ! North-East: x, y
        rad = sqrt(vec_corMaxSqr(1) + vec_corMaxSqr(2))

      else
        isCorner = .false.
        if (coordLoc(1) < origin(1)) then ! West: -x
          normal = -1_rk
          rad = coordLoc(1) - origin(1)

        else if (coordLoc(2) < origin(2)) then ! South: -y
          normal = -1_rk
          rad = coordLoc(2) - origin(2)

        else if (coordLoc(1) > box_max(1)) then ! East: x
          normal = 1_rk
          rad = coordLoc(1) - box_max(1)

        else if (coordLoc(2) > box_max(2)) then ! North: y
          normal = 1_rk
          rad = coordLoc(2) - box_max(2)
        end if
      end if

      if (isCorner) then
        proj_len1 = rad - corInRad
        proj_len2 = corOutRad - rad
      else
        proj_len1 = rad*normal
        proj_len2 = (me%thickness*normal - rad)*normal
      end if

      if (proj_len1 > 0) then
        sigma = const_fac * proj_len2**2 * (proj_len1**4)
        res(i) = min(1.0_rk, sigma) * me%dampFactor
      else
        res(i) = 0.0_rk
      end if
    end do

  end function spongeLayer_box2d_roundCornerPolyn6_for_coord
  ! -------------------------------------------------------------------------- !


  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the spongelayer and fills up
  !!  the res with the target state
  function spongelayer_box2d_vector_for_coord(me, nComp, coord, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> Number of entrys in each array
    integer, intent(in) :: ncomp
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n,ncomp)
    ! --------------------------------------------------------------------------
    integer :: i
    ! --------------------------------------------------------------------------
    res(:, 1) = spongeLayer_box2d_scalar_for_coord(me, coord, n)

    if (ncomp > 1) then
      do i = 1,n
        res(i,2:) = me%targetState(:)
      end do
    end if

  end function spongelayer_box2d_vector_for_coord
  ! -------------------------------------------------------------------------- !


  ! -------------------------------------------------------------------------- !
  !> This function returns the sigma for the spongelayer from treeids
  function spongelayer_box2d_scalar_for_treeIDs(me, treeids, tree, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    select case (trim(me%dampProfile))
    case ('linear', 'exponential')
      res(:) = spongeLayer_box2d_expon_for_treeIDs(me, treeIDs, tree, n)
    case ('polynomial_n5')
      if (me%rounded_corner) then
        res(:) = spongeLayer_box2d_roundCornerPolyn5_for_treeIDs(me, treeIDs, &
          &                                                      tree, n)
      else
        res(:) = spongeLayer_box2d_sharpCornerPolyn5_for_treeIDs(me, treeIDs, &
          &                                                      tree, n)
      end if
    case ('polynomial_n6')
      if (me%rounded_corner) then
        res(:) = spongeLayer_box2d_roundCornerPolyn6_for_treeIDs(me, treeIDs, &
          &                                                      tree, n)
      else
        res(:) = spongeLayer_box2d_sharpCornerPolyn6_for_treeIDs(me, treeIDs, &
          &                                                      tree, n)
      end if
    end select

  end function spongeLayer_box2d_scalar_for_treeIDs
  ! -------------------------------------------------------------------------- !


  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the box shape spongelayer fom treeid
  !! for exponential profile
  function spongelayer_box2d_expon_for_treeIDs(me, treeIDs, tree, n)  &
    &                           result(res)
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: res_i(1), coord(3,1)
    ! --------------------------------------------------------------------------
    do i = 1,n
      !barycentric coordinate
      coord(:,1) = tem_BaryOfId( tree, treeIds(i) )
      res_i = spongeLayer_box2d_expon_for_coord(me, coord, 1)
      res(i) = res_i(1)
    end do

  end function spongelayer_box2d_expon_for_treeIDs
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the box shape spongelayer fom treeid
  !! for polynomial n5 profile
  function spongelayer_box2d_sharpCornerPolyn5_for_treeIDs(me, treeIDs, tree, &
    &                                                      n) result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: res_i(1), coord(3,1)
    ! --------------------------------------------------------------------------
    do i = 1,n
      !barycentric coordinate
      coord(:,1) = tem_BaryOfId( tree, treeIds(i) )
      res_i = spongeLayer_box2d_sharpCornerPolyn5_for_coord(me, coord, 1)
      res(i) = res_i(1)
    end do

  end function spongelayer_box2d_sharpCornerPolyn5_for_treeIDs
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the box shape spongelayer fom treeid
  !! for polynomial n5 profile
  function spongelayer_box2d_roundCornerPolyn5_for_treeIDs(me, treeIDs, tree, &
    &                                                      n) result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: res_i(1), coord(3,1)
    ! --------------------------------------------------------------------------
    do i = 1,n
      !barycentric coordinate
      coord(:,1) = tem_BaryOfId( tree, treeIds(i) )
      res_i = spongeLayer_box2d_roundCornerPolyn5_for_coord(me, coord, 1)
      res(i) = res_i(1)
    end do

  end function spongelayer_box2d_roundCornerPolyn5_for_treeIDs
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the box shape spongelayer fom treeid
  !! for polynomial n5 profile
  function spongelayer_box2d_sharpCornerPolyn6_for_treeIDs(me, treeIDs, tree, &
    &                                                      n) result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: res_i(1), coord(3,1)
    ! --------------------------------------------------------------------------
    do i = 1,n
      !barycentric coordinate
      coord(:,1) = tem_BaryOfId( tree, treeIds(i) )
      res_i = spongeLayer_box2d_sharpCornerPolyn6_for_coord(me, coord, 1)
      res(i) = res_i(1)
    end do

  end function spongelayer_box2d_sharpCornerPolyn6_for_treeIDs
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the box shape spongelayer fom treeid
  !! for polynomial n5 profile
  function spongelayer_box2d_roundCornerPolyn6_for_treeIDs(me, treeIDs, tree, &
    &                                                      n) result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: res_i(1), coord(3,1)
    ! --------------------------------------------------------------------------
    do i = 1,n
      !barycentric coordinate
      coord(:,1) = tem_BaryOfId( tree, treeIds(i) )
      res_i = spongeLayer_box2d_roundCornerPolyn6_for_coord(me, coord, 1)
      res(i) = res_i(1)
    end do

  end function spongelayer_box2d_roundCornerPolyn6_for_treeIDs
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the spongelayer and fills up
  !!  the res with the target state
  function spongelayer_box2d_vector_for_treeIDs(me, ncomp, treeids, tree, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> Number of entrys in each array
    integer, intent(in) :: ncomp
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n,ncomp)
    ! --------------------------------------------------------------------------
    integer :: i
    ! --------------------------------------------------------------------------
    res(:, 1) = spongeLayer_box2d_scalar_for_treeIDs(me, treeids, tree, n)

    if (ncomp > 1) then
      do i = 1,n
        res(i,2:) = me%targetState(:)
      end do
    end if

  end function spongelayer_box2d_vector_for_treeIDs
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !
  !                            SPONGE LAYER RADIAL
  !
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function returns the sigma for the radial viscosity spongelayer
  !! for 2D and 3D
  function spongelayer_radial_scalar_for_coord(me, coord, nDim, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_radial_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> Dimension
    integer, intent(in) :: nDim
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    select case (trim(me%dampProfile))
    case ('linear', 'exponential')
      res(:) = spongeLayer_radial_expon_for_coord(me, coord, nDim, n)
    case ('polynomial_n5')
      res(:) = spongeLayer_radial_polyn5_for_coord(me, coord, nDim, n)
    case ('polynomial_n6')
      res(:) = spongeLayer_radial_polyn6_for_coord(me, coord, nDim, n)
    end select

  end function spongelayer_radial_scalar_for_coord
  ! -------------------------------------------------------------------------- !


  ! -------------------------------------------------------------------------- !
  !> This function returns the sigma for the radial viscosity spongelayer
  !! for 2D and 3D
  function spongelayer_radial_expon_for_coord(me, coord, nDim, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_radial_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> Dimension
    integer, intent(in) :: nDim
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, rad, origin(3), outerRadius, vec(3)
    ! --------------------------------------------------------------------------
    origin(:) = me%origin
    outerRadius = me%radius + me%thickness
    vec = 0.0_rk

    do i = 1, n
      vec(1:nDim) = coord(i, 1:nDim) - origin(1:nDim)
      rad = sqrt( vec(1)**2 + vec(2)**2 + vec(3)**2 )

      if ( rad > me%radius .and. rad < outerRadius ) then
        sigma = me%dampFactor                                        &
          &   * ( (rad-me%radius)/me%thickness )**me%dampExponent
      else if (rad > outerRadius) then
        sigma = me%dampFactor
      else
        sigma = 0.0_rk
      end if

      res(i) = sigma

    end do

  end function spongelayer_radial_expon_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the spongelayer radial from coord.
  !! Sponge profile:
  !! $$ \sigma(x) = \sigma_A*3125*(L+x0-x)*(x-x0)^4)/(256*(L)^5 $$
  !! where, \sigma_A - sponge strength, L - thickness, x0 - start of sponge.
  !!
  !! Profile is taken from:
  !! Xu, Hui; Sagaut, Pierre (2013): Analysis of the absorbing layers for the 
  !! weakly-compressible lattice Boltzmann methods. In Journal of Computational 
  !! Physics 245, pp. 14-42. DOI: 10.1016/j.jcp.2013.02.051.
  function spongeLayer_radial_polyn5_for_coord(me, coord, nDim, n) result(res)
    ! --------------------------------------------------------------------------
    !> Spatial sponge layer to evaluate
    type(tem_spongeLayer_radial_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> Dimension
    integer, intent(in) :: nDim
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: rad, outer_radius, sigma, vec(3)
    real(kind=rk) :: const_fac
    ! --------------------------------------------------------------------------
    outer_radius = me%radius + me%thickness
    vec = 0.0_rk
    const_fac = 3125_rk/(256_rk*me%thickness**5)

    do i = 1,n
      vec(1:nDim) = coord(i,1:nDim) - me%origin(1:nDim)
      rad = sqrt( vec(1)**2 + vec(2)**2 + vec(3)**2 )
      if (rad > me%radius .and. rad < outer_radius) then
        sigma = const_fac * (outer_radius - rad) * (rad - me%radius)**4
        res(i) = sigma * me%dampFactor
      else if (rad > outer_radius) then
        res(i) = me%dampFactor
      else
        res(i) = 0.0_rk
      end if
    end do

  end function spongeLayer_radial_polyn5_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the spongelayer radial from coord 
  !! for the polynomial order n6.
  !! Sponge profile:
  !! $$ \sigma(x) = \sigma_A*729*(L+x0-x)^2*(x-x0)^4)/(16*(L)^6 $$
  !! where, \sigma_A - sponge strength, L - thickness, x0 - start of sponge.
  function spongeLayer_radial_polyn6_for_coord(me, coord, nDim, n) result(res)
    ! --------------------------------------------------------------------------
    !> Spatial sponge layer to evaluate
    type(tem_spongeLayer_radial_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> Dimension
    integer, intent(in) :: nDim
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: rad, outer_radius, sigma, vec(3)
    real(kind=rk) :: const_fac
    ! --------------------------------------------------------------------------
    outer_radius = me%radius + me%thickness
    vec = 0.0_rk
    const_fac = 729_rk/(16_rk*me%thickness**6)

    do i = 1,n
      vec(1:nDim) = coord(i,1:nDim) - me%origin(1:nDim)
      rad = sqrt( vec(1)**2 + vec(2)**2 + vec(3)**2 )
      if (rad > me%radius) then
        sigma = const_fac * (outer_radius - rad)**2 * (rad - me%radius)**4
        res(i) = min(1.0_rk, sigma) * me%dampFactor
      else
        res(i) = 0.0_rk
      end if
    end do

  end function spongeLayer_radial_polyn6_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the radial viscosity spongelayer
  !! for 2D and 3D and fills up rest with target_state.
  !! This function is currectly used to define viscosity sponge in musubi.
  function spongelayer_radial_vector_for_coord(me, nComp, coord, nDim, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_radial_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> Number of entrys in each array
    integer, intent(in) :: nComp
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> Dimension
    integer, intent(in) :: nDim
    !> return value
    real(kind=rk) :: res(n, nComp)
    ! --------------------------------------------------------------------------
    integer :: i
    ! --------------------------------------------------------------------------
    res(:, 1) = spongeLayer_radial_scalar_for_coord(me, coord, nDim, n)  

    if (ncomp > 1) then
      do i = 1,n
        res(i,2:) = me%targetState(:)
      end do
    end if

  end function spongelayer_radial_vector_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function returns the sigma for the radial viscosity spongelayer
  !! for 2D and 3D from treeIDs
  function spongelayer_radial_scalar_for_treeIDs(me, treeids, tree, nDim, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_radial_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> Dimension
    integer, intent(in) :: nDim
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    select case (trim(me%dampProfile))
    case ('linear', 'exponential')
      res(:) = spongeLayer_radial_expon_for_treeIDs(me, treeIDs, tree, nDim, n)
    case ('polynomial_n5')
      res(:) = spongeLayer_radial_polyn5_for_treeIDs(me, treeIDs, tree, nDim, n)
    case ('polynomial_n6')
      res(:) = spongeLayer_radial_polyn6_for_treeIDs(me, treeIDs, tree, nDim, n)
    end select

  end function spongelayer_radial_scalar_for_treeIDs
  ! -------------------------------------------------------------------------- !


  ! -------------------------------------------------------------------------- !
  !> This function returns the sigma for the radial viscosity spongelayer
  !! for 2D and 3D from treeIDs
  function spongelayer_radial_expon_for_treeIDs(me, treeids, tree, nDim, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_radial_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> Dimension
    integer, intent(in) :: nDim
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, rad, origin(3), outerRadius, vec(3), coord(3)
    ! --------------------------------------------------------------------------
    origin(:) = me%origin
    outerRadius = me%radius + me%thickness
    vec = 0.0_rk

    do i = 1, n
      !barycentric coordinate
      coord = tem_BaryOfId( tree, treeIds(i) )
      vec(1:nDim) = coord(1:nDim) - origin(1:nDim)
      rad = sqrt( vec(1)**2 + vec(2)**2 + vec(3)**2 )

      if ( rad > me%radius .and. rad < outerRadius ) then
        sigma = me%dampFactor                                        & 
          &   * ( (rad-me%radius)/me%thickness )**me%dampExponent
      else
        sigma = 0.0_rk
      end if

      res(i) = sigma

    end do

  end function spongelayer_radial_expon_for_treeIDs
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the radial spongelayer
  !! polynomial order 5 from treeids.
  !! Sponge profile:
  !! $$ \sigma(x) = \sigma_A*3125*(L+x0-x)*(x-x0)^4)/(256*(L)^5 $$
  !! where, \sigma_A - sponge strength, L - thickness, x0 - start of sponge.
  !!
  !! Profile is taken from:
  !! Xu, Hui; Sagaut, Pierre (2013): Analysis of the absorbing layers for the 
  !! weakly-compressible lattice Boltzmann methods. In Journal of Computational 
  !! Physics 245, pp. 14-42. DOI: 10.1016/j.jcp.2013.02.051.
  function spongeLayer_radial_polyn5_for_treeids(me, treeids, tree, nDim, n) &
    &                                     result(res)
    ! --------------------------------------------------------------------------
    !> Spatial sponge layer to evaluate
    type(tem_spongeLayer_radial_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> Dimension
    integer, intent(in) :: nDim
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: rad, outer_radius, sigma, vec(3), coord(3)
    real(kind=rk) :: const_fac
    ! --------------------------------------------------------------------------

    outer_radius = me%radius + me%thickness
    vec = 0.0_rk
    const_fac = 3125_rk/(256_rk*me%thickness**5)

    do i = 1,n
      !barycentric coordinate
      coord = tem_BaryOfId( tree, treeIds(i) )
      vec(1:nDim) = coord(1:nDim) - me%origin(1:nDim)
      rad = sqrt( vec(1)**2 + vec(2)**2 + vec(3)**2 )
      if (rad > me%radius .and. rad < outer_radius) then
        sigma = const_fac * (outer_radius - rad) * (rad - me%radius)**4
        res(i) = sigma * me%dampFactor
      else if (rad > outer_radius) then
        res(i) = me%dampFactor
      else
        res(i) = 0.0_rk
      end if
    end do

  end function spongeLayer_radial_polyn5_for_treeids
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the spongeLayer from treeid for
  !! the polynomial order n6.
  !! Sponge profile:
  !! $$ \sigma(x) = \sigma_A*729*(L+x0-x)^2*(x-x0)^4)/(16*(L)^6 $$
  !! where, \sigma_A - sponge strength, L - thickness, x0 - start of sponge.
  function spongeLayer_radial_polyn6_for_treeids(me, treeids, tree, nDim, n) &
    &                                     result(res)
    ! --------------------------------------------------------------------------
    !> Spatial sponge layer to evaluate
    type(tem_spongeLayer_radial_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> Dimension
    integer, intent(in) :: nDim
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: rad, outer_radius, sigma, vec(3), coord(3)
    real(kind=rk) :: const_fac
    ! --------------------------------------------------------------------------

    outer_radius = me%radius + me%thickness
    vec = 0.0_rk
    const_fac = 729_rk/(16_rk*me%thickness**6)

    do i = 1,n
      !barycentric coordinate
      coord = tem_BaryOfId( tree, treeIds(i) )
      vec(1:nDim) = coord(1:nDim) - me%origin(1:nDim)
      rad = sqrt( vec(1)**2 + vec(2)**2 + vec(3)**2 )
      if (rad > me%radius) then
        sigma = const_fac * (outer_radius - rad)**2 * (rad - me%radius)**4
        res(i) = min(1.0_rk, sigma) * me%dampFactor
      else
        res(i) = 0.0_rk
      end if
    end do

  end function spongeLayer_radial_polyn6_for_treeids
  ! -------------------------------------------------------------------------- !


  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the radial viscosity spongelayer
  !! for 2D and 3D from treeIDs and fills up rest with target_state.
  !! This function is currectly used to define viscosity sponge in musubi.
  function spongelayer_radial_vector_for_treeIDs(me, nComp, treeids, tree, &
    &                                            nDim, n) result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_radial_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> Number of entrys in each array
    integer, intent(in) :: nComp
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> Dimension
    integer, intent(in) :: nDim
    !> return value
    real(kind=rk) :: res(n, nComp)
    ! --------------------------------------------------------------------------
    integer :: i
    ! --------------------------------------------------------------------------
    res(:, 1) = spongeLayer_radial_scalar_for_treeIDs(me, treeids, tree, &
      &                                               nDim, n)  

    if (ncomp > 1) then
      do i = 1,n
        res(i,2:) = me%targetState(:)
      end do
    end if

  end function spongelayer_radial_vector_for_treeIDs
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !
  !                            VISCOUS SPONGE LAYER
  !
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the planar viscosity spongelayer 
  !! and multiply with targetState 'viscosity'.
  !! This function is currectly used to define viscosity sponge in musubi.
  function viscSpongelayer_plane_for_coord(me, coord, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_plane_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(3), normal(3), vec(3), proj_len
    ! --------------------------------------------------------------------------
    origin(:) = me%origin
    normal(:) = me%normal

    do i = 1,n
      vec(:) = coord(i,:) - origin(:)
      proj_len = (vec(1)*normal(1)+ vec(2)*normal(2)+vec(3)*normal(3))/   &
                 (normal(1)**2 + normal(2)**2 + normal(3)**2)
      sigma = 1.0_rk + (me%dampFactor-1.0_rk)*((proj_len)**me%dampExponent)
      if (proj_len > 0 .and. proj_len < 1) then
        res(i) = sigma * me%targetState(1)
      else if (proj_len > 1) then
        res(i) = me%dampFactor * me%targetState(1)
      else
        res(i) = me%targetState(1)
      end if

    end do

  end function viscSpongelayer_plane_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the planar viscosity spongelayer 
  !! and multiply with targetState 'viscosity' from treeid.
  !! This function is currectly used to define viscosity sponge in musubi.
  function viscSpongelayer_plane_for_treeIDs(me, treeIDs, tree, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_plane_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(3), normal(3), vec(3), proj_len, coord(3)
    ! --------------------------------------------------------------------------
    origin(:) = me%origin
    normal(:) = me%normal

    do i = 1,n
      !barycentric coordinate
      coord = tem_BaryOfId( tree, treeIds(i) )
      vec(:) = coord(:) - origin(:)
      proj_len = (vec(1)*normal(1)+ vec(2)*normal(2)+vec(3)*normal(3))/   &
                 (normal(1)**2 + normal(2)**2 + normal(3)**2)
      sigma = 1.0_rk + (me%dampFactor-1.0_rk)*((proj_len)**me%dampExponent)
      if (proj_len > 0 .and. proj_len < 1) then
        res(i) = sigma * me%targetState(1)
      else if (proj_len > 1) then
        res(i) = me%dampFactor * me%targetState(1)
      else
        res(i) = me%targetState(1)
      end if

    end do

  end function viscSpongelayer_plane_for_treeIDs
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the box viscosity spongelayer 
  !! and multiply with targetState 'viscosity'.
  !! This function is currectly used to define viscosity sponge in musubi.
  function viscSpongelayer_box_for_coord(me, coord, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(3), extent(3), box_max(3), proj_len
    real(kind=rk) :: coordLoc(3), normal, vec_min(3), vec_max(3)
    real(kind=rk) :: rad, vec_minSqr(3), vec_maxSqr(3)
    ! --------------------------------------------------------------------------
    origin(:) = me%origin
    extent(:) = me%extent
    box_max(:) = origin(:) + extent(:)

    do i = 1,n
      coordLoc = coord(i,:)
      vec_min(:) = coordLoc(:) - origin(:)
      vec_max(:) = coordLoc(:) - box_max(:)
      vec_minSqr(:) = vec_min(:)**2
      vec_maxSqr(:) = vec_max(:)**2
      normal = 1.0_rk
      rad = 0.0_rk


      ! Bottom-South-West: -x,-y,-z
      if (all(coordLoc(:) < origin(:))) then
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(2), vec_minSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) < origin(2) &
        & .and. coordLoc(3) < origin(3)) then ! Bottom-South-East: x, -y, -z
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(2), vec_minSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) > box_max(2) &
        & .and. coordLoc(3) < origin(3)) then ! Bottom-North-West: -x, y, -z
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(2), vec_minSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) > box_max(2) &
        & .and. coordLoc(3) < origin(3)) then ! Bottom-North-East: x, y, -z
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(2), vec_minSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) < origin(2) &
        & .and. coordLoc(3) > box_max(3)) then ! Top-South-West: -x, -y, z
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) < origin(2) &
        & .and. coordLoc(3) > box_max(3)) then ! Top-South-East: x, -y, z
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) > box_max(2) &
        & .and. coordLoc(3) > box_max(3)) then ! Top-North-West: -x, y, z
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(2), vec_maxSqr(3)) )

      else if (all(coordLoc(:) > box_max(:))) then ! Bottom-North-East: x, y, z
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(2) < origin(2) .and. coordLoc(3) < origin(3)) then
        ! Botton-South: -y,-z
        rad = sqrt( max(vec_minSqr(2), vec_minSqr(3)) )

      else if (coordLoc(2) < origin(2) .and. coordLoc(3) > box_max(3)) then
        ! Top-South: -y, z
        rad = sqrt( max(vec_minSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(2) > box_max(2) .and. coordLoc(3) < origin(3)) then
        ! Botom-North: -y, z
        rad = sqrt( max(vec_maxSqr(2), vec_minSqr(3)) )

      else if (coordLoc(2) > box_max(2) .and. coordLoc(3) > box_max(3)) then
        ! Top-North: y, z
        rad = sqrt( max(vec_maxSqr(2), vec_maxSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(3) < origin(3)) then
        ! Botton-West: -x,-z
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(3) < origin(3)) then
        ! Bottom-East: x,-z
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(3) > box_max(3)) then
        ! Top-West: -x, z
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(3)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(3) > box_max(3)) then
        ! Top-East: x, z
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(3)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) < origin(2)) then
        ! South-West: -x,-y
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(2)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) > box_max(2)) then
        ! North-West: -x, y
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(2)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) < origin(2)) then
        ! South-East: x, -y
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(2)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) > box_max(2)) then
        ! North-East: x, y
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(2)) )

      else if (coordLoc(1) < origin(1)) then ! West: -x
        normal = -1_rk
        rad = vec_min(1)

      else if (coordLoc(2) < origin(2)) then ! South: -y
        normal = -1_rk
        rad = vec_min(2)

      else if (coordLoc(3) < origin(3)) then ! Bottom: -z
        normal = -1_rk
        rad = vec_min(3)

      else if (coordLoc(1) > box_max(1)) then ! East: x
        normal = 1_rk
        rad = vec_max(1)

      else if (coordLoc(2) > box_max(2)) then ! North: y
        normal = 1_rk
        rad = vec_max(2)

      else if (coordLoc(3) > box_max(3)) then ! Top: z
        normal = 1_rk
        rad = vec_max(3)
      end if

      proj_len = rad*normal/me%thickness
      sigma = 1.0_rk + (me%dampFactor-1.0_rk)*((proj_len)**me%dampExponent)
      if (proj_len > 0 .and. proj_len < 1) then
        res(i) = sigma * me%targetState(1)
      else if (proj_len > 1) then
        res(i) = me%dampFactor * me%targetState(1)
      else
        res(i) = me%targetState(1)
      end if

    end do

  end function viscSpongelayer_box_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the box viscosity spongelayer 
  !! and multiply with targetState 'viscosity' from treeid.
  !! This function is currectly used to define viscosity sponge in musubi.
  function viscSpongelayer_box_for_treeIDs(me, treeIDs, tree, n)  &
    &                           result(res)
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: res_i(1), coord(3,1)
    ! --------------------------------------------------------------------------
    do i = 1,n
      !barycentric coordinate
      coord(:,1) = tem_BaryOfId( tree, treeIds(i) )
      res_i = viscSpongeLayer_box_for_coord(me, coord, 1)
      res(i) = res_i(1)
    end do

  end function viscSpongelayer_box_for_treeIDs
  ! -------------------------------------------------------------------------- !


  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the box viscosity spongelayer box2d 
  !! and multiply with targetState 'viscosity'.
  !! This function is currectly used to define viscosity sponge in musubi.
  function viscSpongelayer_box2d_for_coord(me, coord, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, origin(2), extent(2), box_max(2), proj_len
    real(kind=rk) :: coordLoc(2), normal, vec_min(2), vec_max(2)
    real(kind=rk) :: rad, vec_minSqr(2), vec_maxSqr(2)
    ! --------------------------------------------------------------------------
    origin(:) = me%origin(1:2)
    extent(:) = me%extent(1:2)
    box_max(:) = origin(:) + extent(:)

    do i = 1,n
      coordLoc = coord(i,1:2)
      vec_min(:) = coordLoc(:) - origin(:)
      vec_max(:) = coordLoc(:) - box_max(:)
      vec_minSqr(:) = vec_min(:)**2
      vec_maxSqr(:) = vec_max(:)**2
      normal = 1.0_rk
      rad = 0.0_rk


      if (coordLoc(1) < origin(1) .and. coordLoc(2) < origin(2)) then
        ! South-West: -x,-y
        rad = sqrt( max(vec_minSqr(1), vec_minSqr(2)) )

      else if (coordLoc(1) < origin(1) .and. coordLoc(2) > box_max(2)) then
        ! North-West: -x, y
        rad = sqrt( max(vec_minSqr(1), vec_maxSqr(2)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) < origin(2)) then
        ! South-East: x, -y
        rad = sqrt( max(vec_maxSqr(1), vec_minSqr(2)) )

      else if (coordLoc(1) > box_max(1) .and. coordLoc(2) > box_max(2)) then
        ! North-East: x, y
        rad = sqrt( max(vec_maxSqr(1), vec_maxSqr(2)) )

      else if (coordLoc(1) < origin(1)) then ! West: -x
        normal = -1_rk
        rad = vec_min(1)

      else if (coordLoc(2) < origin(2)) then ! South: -y
        normal = -1_rk
        rad = vec_min(2)

      else if (coordLoc(1) > box_max(1)) then ! East: x
        normal = 1_rk
        rad = vec_max(1)

      else if (coordLoc(2) > box_max(2)) then ! North: y
        normal = 1_rk
        rad = vec_max(2)
      end if

      proj_len = rad*normal/me%thickness
      sigma = 1.0_rk + (me%dampFactor-1.0_rk)*((proj_len)**me%dampExponent)
      if (proj_len > 0 .and. proj_len < 1) then
        res(i) = sigma * me%targetState(1)
      else if (proj_len > 1) then
        res(i) = me%dampFactor * me%targetState(1)
      else
        res(i) = me%targetState(1)
      end if

    end do

  end function viscSpongelayer_box2d_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the box viscosity spongelayer box2d 
  !! and multiply with targetState 'viscosity' from treeid.
  !! This function is currectly used to define viscosity sponge in musubi.
  function viscSpongelayer_box2d_for_treeIDs(me, treeIDs, tree, n)  &
    &                           result(res)
    !> Spacetime function to evaluate
    type(tem_spongeLayer_box_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: res_i(1), coord(3,1)
    ! --------------------------------------------------------------------------
    do i = 1,n
      !barycentric coordinate
      coord(:,1) = tem_BaryOfId( tree, treeIds(i) )
      res_i = viscSpongeLayer_box2d_for_coord(me, coord, 1)
      res(i) = res_i(1)
    end do

  end function viscSpongelayer_box2d_for_treeIDs
  ! -------------------------------------------------------------------------- !


  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the radial viscosity spongelayer
  !! for 2D and 3D, and multiply with targetState.
  !! This function is currectly used to define viscosity sponge in musubi.
  function viscSpongelayer_radial_for_coord(me, coord, nDim, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_radial_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> Dimension
    integer, intent(in) :: nDim
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, rad, origin(3), vec(3), outerRadius
    ! --------------------------------------------------------------------------
    origin(:) = me%origin
    outerRadius = me%radius + me%thickness
    vec = 0.0_rk

    do i = 1, n
      vec(1:nDim) = coord(i, 1:nDim) - origin(1:nDim)
      rad = sqrt( vec(1)**2 + vec(2)**2 + vec(3)**2 )

      if ( rad > me%radius .and. rad < outerRadius ) then
        sigma = 1.0_rk + (me%dampFactor-1.0_rk)                               &
          &            * ( (rad-me%radius)/me%thickness )**me%dampExponent
      else if (rad > outerRadius) then
        sigma = me%dampFactor
      else
        sigma = 1.0_rk
      end if

      res(i) = sigma * me%targetState(1)

    enddo

  end function viscSpongelayer_radial_for_coord
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This function calculates the sigma for the radial viscosity spongelayer
  !! for 2D and 3D, and multiply with targetState using treeid.
  !! This function is currectly used to define viscosity sponge in musubi.
  function viscSpongelayer_radial_for_treeIDs(me, treeIDs, tree, nDim, n)  &
    &                           result(res)
    ! --------------------------------------------------------------------------
    !> Spacetime function to evaluate
    type(tem_spongeLayer_radial_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> Dimension
    integer, intent(in) :: nDim
    !> return value
    real(kind=rk) :: res(n)
    ! --------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: sigma, rad, origin(3), vec(3), coord(3)
    real(kind=rk) :: outerRadius
    ! --------------------------------------------------------------------------
    origin(:) = me%origin
    outerRadius = me%radius + me%thickness
    vec = 0.0_rk

    do i = 1, n
      !barycentric coordinate
      coord = tem_BaryOfId( tree, treeIds(i) )
      vec(1:nDim) = coord(1:nDim) - origin(1:nDim)
      rad = sqrt( vec(1)**2 + vec(2)**2 + vec(3)**2 )

      if ( rad > me%radius .and. rad < outerRadius ) then
        sigma = 1.0_rk + (me%dampFactor-1.0_rk)                               &
          &            * ( (rad-me%radius)/me%thickness )**me%dampExponent
      else if (rad > outerRadius) then
        sigma = me%dampFactor
      else
        sigma = 1.0_rk
      end if

      res(i) = sigma * me%targetState(1)

    enddo

  end function viscSpongelayer_radial_for_treeIDs
  ! -------------------------------------------------------------------------- !

end module tem_spongeLayer_module
