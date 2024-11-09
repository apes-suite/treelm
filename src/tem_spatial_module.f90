! Copyright (c) 2011-2016,2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2016,2021-2022 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012, 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2013-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014 Julia Moos <julia.moos@student.uni-siegen.de>
! Copyright (c) 2015-2016, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2019 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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

!> Spatial functions let you describe functions that may vary in space.
!! The are needed for initial conditions and volumetric information like
!! source terms.
!! A space-time function that depends also on time can be constructed by
!! the space-time function of kind `combined`, see the
!! [[tem_spacetime_fun_module]].
!!
!! The spatial functions may be constants, various predefined functions, or
!! arbitrary Lua functions that depend on the three space coordinates.
!! In the Lua script, the spatial function is either a constant scalar number,
!! a Lua function, or a table.
!!
!! If it is a table, it either is a constant definition for a vector value,
!! for example a velocity definition, or it is a table that contains a
!! definition of either a `fun`, `const` or `predefined` key.
!! Here, `fun` would in turn need to be a Lua function, and const the definition
!! of constant values (either a scalar, or a vector in a table).
!! The `predefined` key signals the use of one of the available predefined
!! functions and further parameters for their definition depend on the
!! respective functions.
!!
!! Available predefined functions are:
!!
!! - 'crvpx' [[tem_ic_predefs_module:load_ic_2dcrvp]]
!! - 'crvpy' [[tem_ic_predefs_module:load_ic_2dcrvp]]
!! - 'crvppressure' [[tem_ic_predefs_module:load_ic_2dcrvp]]
!! - 'cylindricalwave' [[tem_cylindricalWave_module]]
!! - 'gausspulse' [[tem_ic_predefs_module:load_ic_gausspulse]]
!! - 'heaviside_gibbs'
!!   [[tem_heaviside_gibbs_fun_module:tem_load_heaviside_gibbs]]
!! - 'miescatter_magneticfieldx' [[tem_miescatter_module:tem_load_miescatter]]
!! - 'miescatter_magneticfieldy' [[tem_miescatter_module:tem_load_miescatter]]
!! - 'parabol' [[tem_spatial_module:load_spatial_parabol]]
!! - 'pml' [[tem_pmlLayer_module]]
!! - 'polygon_material' [[tem_polygon_material_module]]
!! - 'polygon_material_3d' [[tem_polygon_material_module]]
!! - 'random' [[tem_spatial_module:load_spatial_random]]
!! - 'rectangular'
!! - 'spongelayer_plane' [[tem_spongeLayer_module]]
!! - 'spongelayer_plane_1d' [[tem_spongeLayer_module]]
!! - 'viscous_spongelayer_radial' [[tem_spongeLayer_module:tem_load_spongeLayer_radial]]
!! - 'viscous_spongelayer_radial_2d' [[tem_spongeLayer_module:tem_load_spongeLayer_radial]]
!! - 'viscous_spongelayer_box' [[tem_spongeLayer_module:tem_load_spongeLayer_box]]
!! - 'viscous_spongelayer_box_2d' [[tem_spongeLayer_module:tem_load_spongeLayer_box]]
!! - 'spongelayer_box' [[tem_spongeLayer_module:tem_load_spongeLayer_box]]
!! - 'spongelayer_box' [[tem_spongeLayer_module:tem_load_spongeLayer_box]]
!! - 'spongelayer_box_2d' [[tem_spongeLayer_module:tem_load_spongeLayer_box]]
!! - 'spongelayer_radial' [[tem_spongeLayer_module:tem_load_spongeLayer_radial]]
!! - 'tgv_p' [[tem_ic_predefs_module:load_ic_tgv]]
!! - 'tgv_ux' [[tem_ic_predefs_module:load_ic_tgv]]
!! - 'tgv_uy' [[tem_ic_predefs_module:load_ic_tgv]]
!! - 'tgv_sxx' [[tem_ic_predefs_module:load_ic_tgv]]
!! - 'tgv_sxz' [[tem_ic_predefs_module:load_ic_tgv]]
!! - 'tgv_syy' [[tem_ic_predefs_module:load_ic_tgv]]
!! - 'tgv_syz' [[tem_ic_predefs_module:load_ic_tgv]]
!!
!! The simplest spatial function definition is just a constant scalar number:
!!
!!```lua
!!  spatial = 1.0
!!```
!!
!! Similarly the definition of a spatial function by a Lua function can simply
!! be written as
!!
!!```lua
!!  spatial = fun(x,y,z)
!!    return x*y*z
!!  end fun
!!```
!!
!! The Lua function needs to take 3 arguments, representing the three spatial
!! coordinates.
!!
!! An example for the definition of a predefined function can take the following
!! form:
!!
!!```lua
!!  spatial = {
!!    predefined = 'random',
!!    min = 0,
!!    max = 1
!!  }
!!```
!!
module tem_spatial_module

  ! include treelm modules
  use env_module,             only: rk, long_k, LabelLen, globalMaxLevels
  use tem_aux_module,         only: tem_abort
  use tem_tools_module,       only: upper_to_lower
  use treelmesh_module,       only: treelmesh_type
  use tem_geometry_module,    only: tem_BaryOfId
  use tem_canonicalND_module, only: tem_canonicalND_type
  use tem_shape_module,       only: tem_shape_type, tem_load_shape
  use tem_pointData_module,   only: tem_grwPoints_type
  use tem_ic_predefs_module,  only: ic_gausspulse_type, load_ic_gausspulse, &
    &                               ic_gausspulse_for,                      &
    &                               ic_2dcrvp_type, load_ic_2dcrvp,         &
    &                               ic_2dcrvpX_for, ic_2dcrvpY_for,         &
    &                               ic_2dcrvpPressure_for,                  &
    &                               ic_tgv_type, load_ic_tgv,               &
    &                               ic_tgv_pressure_for,                    &
    &                               ic_tgv_ux_for, ic_tgv_uy_for,           &
    &                               ic_tgv_Sxx_for, ic_tgv_Syy_for,         &
    &                               ic_tgv_Sxz_for, ic_tgv_Syz_for
  use tem_logging_module,     only: logUnit, tem_toStr
  use tem_grow_array_module,  only: grw_realArray_type, append, truncate
  use tem_miescatter_module,  only: tem_miescatter_field_type,              &
    &                               tem_load_miescatter_displacementfieldz, &
    &                               tem_load_miescatter_magneticfieldx,     &
    &                               tem_load_miescatter_magneticfieldy,     &
    &                               tem_eval_miescatter_displz,             &
    &                               tem_eval_miescatter_magnx,              &
    &                               tem_eval_miescatter_magny
  use tem_heaviside_gibbs_fun_module, only: tem_heaviside_gibbs_type, &
    &                                       tem_load_heaviside_gibbs, &
    &                                       tem_eval_heaviside_gibbs
  use tem_pmlLayer_module,            only: tem_pmlLayer_type, &
    &                                       tem_evaluate_pml,  &
    &                                       tem_load_pmlLayer
  use tem_cylindricalWave_module,     only: tem_cylindricalWave_type, &
    &                                       tem_load_cylindricalWave, &
    &                                       tem_eval_cylindricalWave
  use tem_polygon_material_module,    only: tem_polygon_material_type,      &
    &                                       tem_polygon_material_load,      &
    &                                       tem_eval_polygon_material,      &
    &                                       tem_eval_polygon_material_3d,   &
    &                                       tem_eval_polygon_material_scal, &
    &                                       tem_eval_polygon_material_scal_3d
  use tem_spongeLayer_module,         only: tem_spongeLayer_plane_type,      &
    &                                       tem_spongeLayer_radial_type,     &
    &                                       tem_spongeLayer_box_type,        &
    &                                       tem_load_spongeLayer_plane,      &
    &                                       tem_load_spongeLayer_radial,     &
    &                                       tem_load_spongeLayer_box,        &
    &                                       tem_spongeLayer_plane_for,       &
    &                                       tem_spongeLayer_box_for,         &
    &                                       tem_spongeLayer_box2d_for,       &
    &                                       tem_spongeLayer_radial_for,      &
    &                                       tem_viscSpongeLayer_plane_for,   &
    &                                       tem_viscSpongeLayer_radial_for,  &
    &                                       tem_viscSpongelayer_box_for,     &
    &                                       tem_viscSpongelayer_box2d_for

  ! include aotus modules
  use aotus_module,     only: flu_State, aoterr_Fatal, aoterr_NonExistent, &
    &                         aoterr_WrongType, aot_top_get_val, aot_get_val
  use aot_table_module, only: aot_table_open, aot_table_close
  use aot_fun_module,   only: aot_fun_type, aot_fun_open, aot_fun_close, &
    &                         aot_fun_put, aot_fun_do
  use aot_references_module, only: aot_reference_for

  implicit none

  private

  real(kind=rk) :: ref_amp = 1.0_rk

  public :: tem_spatial_type
  public :: tem_load_spatial
  public :: tem_spatial_lua_for
  public :: tem_spatial_parabol3d_for
  public :: tem_spatial_parabol2d_for
  public :: tem_spatial_for
  public :: tem_spatial_storeVal

  !> stores values of spatial term during initialization to reduce
  !! computations during main loop.
  type spatial_value_type
    !> Evaluated variable value on each point.
    !! For vectorial variable, the values are stored as
    !! (iVal-1)*nComp + iComp
    type(grw_realArray_type) :: evalVal
  end type spatial_value_type

  !> defines parabola functions
  !! shape kind defines 2d or 3d parabola
  type spatial_parabol_type
    !> define point, line or plane spatial profile
    type( tem_shape_type ) :: geometry
    !> magnitude of the spatial value
    real(kind=rk), allocatable :: amplitude(:)
  end type spatial_parabol_type

  !> Defines a random spatial distribution within a given interval.
  type spatial_random_type
    !> Minimal value to produce
    real(kind=rk) :: val_min
    !> Maximal value to produce
    real(kind=rk) :: val_max
    !> Length of the value interval
    real(kind=rk) :: val_range
  end type spatial_random_type

  !> Defines a stationary solution of the Euler equation by the Hopf Fibration.
  type spatial_hopf_type
    !> Minimal value to produce
    real(kind=rk) :: radius
    !> Maximal value to produce
    real(kind=rk) :: p0
  end type spatial_hopf_type

  !> contains spatial state information
  type tem_spatial_type
    !> kind: how is the spatial defined?
    !!
    !! - 'none' = no spatial modifier defined
    !! - 'const' = a constant factor
    !! - 'lua_fun' = defined as a lua function
    !! - 'parabol' = parabolic function
    !!               shape = {object = line} defines 2d parabola
    !!               shape = {object = plane} defines 3d parabola
    !! Further predefined functions might be added here
    character(len=LabelLen) :: kind

    !> constant spatial value for nComponents
    real(kind=rk), allocatable :: const(:)

    !> Reference to the Lua function if the spatial function is defined
    !! as a Lua function.
    integer :: lua_fun_ref = 0

    type(flu_state) :: conf

    !> defines gausspulse
    type( ic_gausspulse_type ) :: gausspulse

    !> 2d co-rotating vortex pair
    type( ic_2dcrvp_type ) :: crvp

    !> Taylor-Green vortex type
    type( ic_tgv_type ) :: tgv

    !> Parabol type
    type( spatial_parabol_type ) :: parabol

    !> Random type
    type( spatial_random_type ) :: random

    !> Hopf Fibration type
    type( spatial_hopf_type ) :: hopf

    !> Spatial function for Mie-series solution
    !! of electrodynamic scattering at dielectric sphere.
    type(tem_miescatter_field_type) :: mie_fun

    !> Spatial function for Heaviside function including Gibbs
    !! oscillations.
    type(tem_heaviside_gibbs_type) :: heaviside_gibbs_fun

    !> store spatial value on each level
    type(spatial_value_type) :: valOnLvl(globalMaxLevels)

    !> type for the plane sponge layer
    type(tem_spongeLayer_plane_type) :: spongePlane

    !> type for the box sponge layer
    type(tem_spongeLayer_box_type) :: spongeBox

    !> type for the radial sponge layer
    type(tem_spongeLayer_radial_type) :: spongeRadial

    !> type for the pml damping medium
    type(tem_pmlLayer_type) :: pml

    !> type for a scalar cylindrical wave.
    type(tem_cylindricalWave_type) :: cylindricalWave

    !> Description of a material definition by a polygon.
    type(tem_polygon_material_type) :: polygon_material

    !> range of x and y dimention for rectangular function
    real(kind=rk) :: rect_ly
    real(kind=rk) :: rect_lz

    !> to store spatial values during initialization
    logical :: isStored = .false.
  end type tem_spatial_type

  interface tem_spatial_lua_for
    module procedure tem_spatial_lua_for_treeIDs
    module procedure tem_spatial_lua_for_coord
    module procedure tem_spatial_lua_for_index
    module procedure tem_spatial_lua_vector_for_treeIDs
    module procedure tem_spatial_lua_vector_for_coord
    module procedure tem_spatial_lua_vector_for_index
  end interface tem_spatial_lua_for

  interface tem_spatial_parabol2d_for
    module procedure tem_spatial_parabol2d_for_treeIDs
    module procedure tem_spatial_parabol2d_for_coord
  end interface tem_spatial_parabol2d_for

  interface tem_spatial_parabol3d_for
    module procedure tem_spatial_parabol3d_for_treeIDs
    module procedure tem_spatial_parabol3d_for_coord
  end interface tem_spatial_parabol3d_for

  interface tem_spatial_for
    module procedure tem_spatial_for_treeIDs
    module procedure tem_spatial_for_coord
    module procedure tem_spatial_vector_for_treeIDs
    module procedure tem_spatial_vector_for_coord
    module procedure tem_spatial_scalar_for_index
    module procedure tem_spatial_vector_for_index
  end interface tem_spatial_for

  interface tem_spatial_storeVal
    module procedure tem_spatial_scalar_storeVal
    module procedure tem_spatial_vector_storeVal
  end interface tem_spatial_storeVal


contains


  ! ------------------------------------------------------------------------ !
  !>  This subroutine load spatial boundary state variable.
  !!
  !! If spatial is defined as block than read block for predefined Fortran
  !! function variables else it is defined as constant.
  !! Valid definitions:
  !!
  !! - Constant
  !!```lua
  !! spatial = 1.0
  !!```
  !! - lua_function
  !!```lua
  !! spatial = lua_fun_name
  !!```
  !! - define lua function inside a table
  !!```lua
  !! spatial  = {fun=lua_fun_name, store=<>}
  !!```
  !! Note. Lua function take 3 input arguments (x,y,z) i.e barycentric
  !! coordinates of an element
  !! - Predefined Fortran function
  !!```lua
  !! spatial  = {predefined = "fun_name", fun_parameters}
  !!```
  !! example given in subroutine [[load_spatial_parabol]] definition
  !!
  subroutine tem_load_spatial( me, conf, parent, key, defaultValue, nComp, &
    &                          errCode                                     )
    ! -------------------------------------------------------------------- !
    !> spatial boundary state type
    type(tem_spatial_type), intent(out) :: me
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: parent
    !> state variable key string defined in lua
    character(len=*), intent(in), optional :: key
    !> What should be set s a default value for the quantities if no
    !! quantity was given in the lua file
    real(kind=rk), intent(in), optional :: defaultValue
    !> number of components of the variable
    integer, intent(in), optional :: nComp
    !> Error code from lua loading
    integer, intent(out), optional :: errCode
    ! -------------------------------------------------------------------- !
    integer :: thandle, iError
    integer :: loc_nComp
    type(aot_fun_type) :: fun
    character(len=labelLen) :: local_key
    real(kind=rk) :: local_default
    logical :: loadasConst
    ! -------------------------------------------------------------------- !

    if( present( nComp ))then
      loc_nComp = nComp
    else
      loc_nComp = 1
    end if

    if(present(defaultValue)) then
      local_default = defaultValue
    else
      local_default = 1._rk
    endif

    if(present(key)) then
      local_key = trim(key)
    else
      local_key = 'spatial'
    endif
    write(logUnit(1),"(A)") 'Loading spatial for '//trim(local_key)

    ! default values
    loadasConst = .false.

    ! First test for a lua function
    call aot_fun_open(L=conf, parent=parent, fun=fun, key=trim(local_key))

    me%conf = conf

    if_fun: if (fun%handle /= 0) then

      ! There was a function found for the spatial entry, go on using it
      me%kind = 'lua_fun'
      ! Create a reference for this function and store it.
      me%lua_fun_ref = aot_reference_for(conf)
      call aot_fun_close( L=conf, fun=fun )
      iError = 0

    else ! Not a LUA function

      call aot_fun_close( L=conf, fun=fun )

      ! Not a function, try to interpret it as a table'
      call aot_table_open( L       = conf,           &
        &                  thandle = thandle,        &
        &                  parent  = parent,         &
        &                  key     = trim(local_key) )

      ! table is defined
      if_table: if (thandle /= 0) then

        ! check whether to store spatial values during initialization or not
        call aot_get_val( L       = conf,        &
          &               thandle = thandle,     &
          &               key     = 'store',     &
          &               val     = me%isStored, &
          &               default = .false.,     &
          &               ErrCode = iError       )

        ! inside table try to open lua function with key 'fun'
        call aot_fun_open( L      = conf,    &
          &                parent = thandle, &
          &                fun    = fun,     &
          &                key    = 'fun'    )

        ! lua function defined with key 'fun'
        if (fun%handle /= 0) then
          me%kind = 'lua_fun'
          ! Create a reference for this function and store it.
          me%lua_fun_ref = aot_reference_for(conf)

          call aot_fun_close( L=conf, fun=fun )
          iError = 0
        else
          call aot_fun_close( L=conf, fun=fun )
          ! Not a lua function with key 'fun'.
          ! Try to interpret it as a predefined'
          call aot_get_val( L       = conf,         &
            &               thandle = thandle,      &
            &               key     = 'predefined', &
            &               val     = me%kind,      &
            &               default = 'unknown',    &
            &               ErrCode = iError        )
          if( btest(iError, aoterr_WrongType) ) then
            call tem_abort('Error in retrieving the function kind.')
          end if

          ! predefined key word is not defined.
          ! try to interpret as constant with key "const"
          if ( btest(iError, aoterr_NonExistent) ) then
            call load_spatial_asConst( const   = me%const, &
              &                        conf    = conf,     &
              &                        errCode = iError,   &
              &                        parent  = thandle,  &
              &                        key     = 'const',  &
              &                        nComp   = loc_nComp )
            if (iError == 0) then
              me%kind = 'const'
            else
              call aot_table_close( L=conf, thandle=thandle )
              ! Constant is not defined with key "const"
              ! Try to interpret directly from local_key
              loadasConst = .true.
            end if
          else ! predefined exist
            call load_spatial_predefined( me      = me,       &
              &                           conf    = conf,     &
              &                           thandle = thandle,  &
              &                           nComp   = loc_nComp )
            iError = 0
          end if ! not predefined
        end if ! not lua_fun with key "fun"

      else

        ! not a table try to interprect as constant
        loadasConst = .true.

      end if if_table

      call aot_table_close( L=conf, thandle=thandle )

    end if if_fun

    if (loadasConst) then
      ! Not a table with specific key word like "fun", "const", "predefined"
      ! Try to interpret it as a constant with key
      write(logUnit(7),"(A)") 'Try to load as constant with key ' &
        &                     // trim(local_key)
      call load_spatial_asConst( const   = me%const,        &
        &                        conf    = conf,            &
        &                        errCode = iError,          &
        &                        parent  = parent,          &
        &                        key     = trim(local_key), &
        &                        nComp   = loc_nComp        )

      if (iError == 0) then
        me%kind = 'const'
      else
        write(logUnit(1),"(A)") ' WARNING: variable is non-existent!'
        write(logUnit(1),"(A)") ' Setting kind to be "none" and value = ' &
          &                     // trim( tem_toStr(local_default) )       &
          &                     // ' for all ' //                         &
          &                     trim( tem_toStr(loc_nComp) ) // ' components'
        if (allocated(me%const)) deallocate(me%const)
        allocate(me%const(loc_nComp))
        me%const = local_default
        me%kind = 'none'
      end if
    end if

    write(logUnit(3),"(A)") '   Spatial for ' // trim(local_key) &
      & // ' is a defined as ' // trim(me%kind)
    if (trim(me%kind) == 'const') then
      if ( loc_nComp == 1 ) then
        write(logUnit(3),"(A)") 'value = ' // trim(tem_toStr(me%const(1)))
      else
        write(logUnit(3),"(A)") 'value = '              &
          & // trim(tem_toStr(me%const(1:loc_nComp),','))
      end if
    end if

    if (present(errCode)) errCode = iError

  end subroutine tem_load_spatial
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Load spatial as constant
  subroutine load_spatial_asConst( const, conf, errCode, parent, key, nComp )
    ! -------------------------------------------------------------------- !
    !> value to be filled
    real(kind=rk), allocatable :: const(:)
    !> lua state type
    type(flu_State) :: conf
    !> errCode = -1, if spatial function is not defined as constant
    integer, intent(out) :: errCode
    !> aotus parent handle
    integer, intent(in) :: parent
    !> key is either "local_key" or "const"
    character(len=*), intent(in) :: key
    !> number of components of the variable
    integer, intent(in) :: nComp
    ! -------------------------------------------------------------------- !
    integer, allocatable :: vError(:), vErr_NonExistent(:), vErr_WrongType(:)
    ! -------------------------------------------------------------------- !

    errCode = 0

    if (nComp == 1) then
      write(logUnit(5),*) 'Load spatial as single constant'
      allocate(const(nComp), vError(nComp))
      call aot_get_val( L       = conf,     &
        &               thandle = parent,   &
        &               val     = const(1), &
        &               key     = key,      &
        &               ErrCode = vError(1) )

      if (btest(vError(1), aoterr_NonExistent)) then
        ! if previous checks fails return errCode =-1
        ! and set default kind = none and value = 1
        write(logUnit(3),*) 'ERROR: Fatal error occured in loading ' &
          & // 'scalar constants'
        errCode = -1
        return
      end if
    else !nComp > 1
      write(logUnit(3),"(A,I0)") &
        & 'Load spatial as array of constant with nComp ', nComp
      !Try to load constant of nComp
      call aot_get_val( L         = conf,   &
        &               thandle   = parent, &
        &               val       = const,  &
        &               maxLength = nComp,  &
        &               key       = key,    &
        &               ErrCode   = vError  )
      if (size(const) == 0) then
        write(logUnit(3),*) 'Error loading array constant from key "'&
          & // trim(key) // '"'
        errCode = -1
        return
      else
        ! if any error index is fatal then return errCode - 1
        allocate(vErr_NonExistent(size(vError)))
        vErr_NonExistent = aoterr_NonExistent
        if(any(btest(vError, vErr_NonExistent))) then
          write(logUnit(3),*) 'ERROR: Fatal error occured in loading ' &
            & // 'array of constants with nComp:', nComp
          errCode = -1
          return
        else if (size(const) /= nComp) then
          ! no fatal error but number of constant defined /= nComp
          write(logUnit(1),*) 'Spatial as constant. table size /= nComp', &
            &                 nComp
          call tem_abort()
        end if
      end if
    end if

    allocate(vErr_WrongType(size(vError)))
    vErr_WrongType = aoterr_WrongType
    if (any(btest(vError, vErr_WrongType))) then
      write(logUnit(1),*)'FATAL Error occured in definition of ' &
        & // 'the spatial'
      write(logUnit(1),*)'while retrieving spatial constant:'
      write(logUnit(1),*)'Variable has wrong type (should be a ' &
        & // 'real number)!'
      call tem_abort()
    end if

  end subroutine load_spatial_asConst
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Load predefined spatial function
  subroutine load_spatial_predefined( me, conf, thandle, nComp )
    ! -------------------------------------------------------------------- !
    !> spatial fun information
    type(tem_spatial_type), intent(inout) :: me
    !> lua state type
    type(flu_State) :: conf
    !> spatial function handle
    integer, intent(in) :: thandle
    !> Number of components
    integer, intent(in) :: nComp
    ! -------------------------------------------------------------------- !
    integer :: nDim, iError
    character(len=labelLen) :: funkind
    ! -------------------------------------------------------------------- !

    write(logUnit(1),*) 'Spatial is defined as predefined:'

    funkind = adjustl(me%kind)
    funkind = upper_to_lower(funkind)

    select case(trim(funkind))
    case('parabol')
      write(logUnit(3),*) '  is a parabolic profile'
      call load_spatial_parabol( me      = me%parabol, &
        &                        conf    = conf,       &
        &                        thandle = thandle,    &
        &                        nComp   = nComp       )
    case('random')
      write(logUnit(3),*) '  is a random distribution'
      call load_spatial_random( me      = me%random, &
        &                       conf    = conf,      &
        &                       thandle = thandle    )
    case('gausspulse')
      write(logUnit(3),*) '  is a gaussian pulse'
      call load_ic_gausspulse( conf    = conf,         &
        &                      thandle = thandle,      &
        &                      me      = me%gausspulse )
    case('tgv_p')
      write(logUnit(5),*) '  is a Taylor-Green vortex pressure'
      call load_ic_tgv( conf    = conf,         &
        &               thandle = thandle,      &
        &               me      = me%tgv        )
    case('tgv_ux')
      write(logUnit(5),*) '  is a Taylor-Green vortex Ux'
      call load_ic_tgv( conf    = conf,         &
        &               thandle = thandle,      &
        &               me      = me%tgv        )
    case('tgv_uy')
      write(logUnit(5),*) '  is a Taylor-Green vortex Uy'
      call load_ic_tgv( conf    = conf,         &
        &               thandle = thandle,      &
        &               me      = me%tgv        )
    case('tgv_sxx')
      write(logUnit(5),*) '  is a Taylor-Green vortex Sxx'
      call load_ic_tgv( conf    = conf,         &
        &               thandle = thandle,      &
        &               me      = me%tgv        )
    case('tgv_syy')
      write(logUnit(5),*) '  is a Taylor-Green vortex Syy'
      call load_ic_tgv( conf    = conf,         &
        &               thandle = thandle,      &
        &               me      = me%tgv        )
    case('tgv_syz')
      write(logUnit(5),*) '  is a Taylor-Green vortex Syz'
      call load_ic_tgv( conf    = conf,         &
        &               thandle = thandle,      &
        &               me      = me%tgv        )
    case('tgv_sxz')
      write(logUnit(5),*) '  is a Taylor-Green vortex Sxz'
      call load_ic_tgv( conf    = conf,         &
        &               thandle = thandle,      &
        &               me      = me%tgv        )
!HK: Not implemented yet.
!HK!        case('hopf')
!HK!          write(logUnit(1),*) '  is a Hopf fibration solution'
!HK!          call load_spatial_hopf( me = me%hopf, conf = conf, &
!HK!            &                     thandle = thandle )
    case('crvpx')
      write(logUnit(3),*) '  is a Spinning vortex pair X component'
      call load_ic_2dcrvp( conf    = conf,    &
        &                  thandle = thandle, &
        &                  me      = me%crvp  )
    case('crvpy')
      write(logUnit(3),*) '  is a Spinning vortex pair Y component'
      call load_ic_2dcrvp( conf    = conf,    &
        &                  thandle = thandle, &
        &                  me      = me%crvp  )
    case('crvppressure')
      write(logUnit(3),*) '  is a Spinning vortex pair density'
      call load_ic_2dcrvp( conf    = conf,    &
        &                  thandle = thandle, &
        &                  me      = me%crvp  )
    case('miescatter_displacementfieldz')
      call tem_load_miescatter_displacementfieldz( conf    = conf,      &
        &                                          thandle = thandle,   &
        &                                          me      = me%mie_fun )
    case('miescatter_magneticfieldx')
      call tem_load_miescatter_magneticfieldx( conf    = conf,      &
        &                                      thandle = thandle,   &
        &                                      me      = me%mie_fun )
    case('miescatter_magneticfieldy')
      call tem_load_miescatter_magneticfieldy( conf    = conf,      &
        &                                      thandle = thandle,   &
        &                                      me      = me%mie_fun )
    case('heaviside_gibbs')
      call tem_load_heaviside_gibbs( conf    = conf,                  &
        &                            thandle = thandle,               &
        &                            me      = me%heaviside_gibbs_fun )
    case('spongelayer_plane')
      ndim = 3
      call tem_load_spongeLayer_plane( conf    = conf,           &
        &                              thandle = thandle,        &
        &                              me      = me%spongePlane, &
        &                              ndim    = ndim,           &
        &                              nComp   = nComp           )
    case('spongelayer_plane_2d')
      ndim = 2
      call tem_load_spongeLayer_plane( conf    = conf,           &
        &                              thandle = thandle,        &
        &                              me      = me%spongePlane, &
        &                              ndim    = ndim,           &
        &                              nComp   = nComp           )
    case('spongelayer_plane_1d')
      ndim = 1
      call tem_load_spongeLayer_plane( conf    = conf,           &
        &                              thandle = thandle,        &
        &                              me      = me%spongePlane, &
        &                              ndim    = ndim,           &
        &                              nComp   = nComp           )
    case('spongelayer_radial')
      write(logUnit(3),*) '  is a radial sponge layer'
      ndim = 3
      call tem_load_spongelayer_radial( me       = me%spongeRadial, &
        &                               conf     = conf,            &
        &                               thandle  = thandle,         &
        &                               nDim     = nDim,            &
        &                               nComp    = nComp            )
    case('spongelayer_radial_2d')
      write(logUnit(3),*) '  is a 2d radial sponge layer'
      ndim = 2
      call tem_load_spongelayer_radial( me       = me%spongeRadial, &
        &                               conf     = conf,            &
        &                               thandle  = thandle,         &
        &                               nDim     = nDim,            &
        &                               nComp    = nComp            )
    case('spongelayer_box')
      call tem_load_spongeLayer_box( conf    = conf,         &
        &                            thandle = thandle,      &
        &                            me      = me%spongeBox, &
        &                            ndim    = 3,            &
        &                            nComp   = nComp         )
    case('spongelayer_box_2d')
      call tem_load_spongeLayer_box( conf    = conf,         &
        &                            thandle = thandle,      &
        &                            me      = me%spongeBox, &
        &                            ndim    = 2,            &
        &                            nComp   = nComp         )
    case('viscous_spongelayer_plane')
      ndim = 3
      call tem_load_spongeLayer_plane( conf      = conf,           &
        &                              thandle   = thandle,        &
        &                              me        = me%spongePlane, &
        &                              ndim      = nDim,           &
        &                              nComp     = nComp,          &
        &                              stateName = 'viscosity'     )
    case('viscous_spongelayer_box')
      ndim = 3
      call tem_load_spongeLayer_box( conf      = conf,         &
        &                            thandle   = thandle,      &
        &                            me        = me%spongeBox, &
        &                            ndim      = nDim,         &
        &                            nComp     = nComp,        &
        &                            stateName = 'viscosity'   )
    case('viscous_spongelayer_box_2d')
      ndim = 2
      call tem_load_spongeLayer_box( conf      = conf,         &
        &                            thandle   = thandle,      &
        &                            me        = me%spongeBox, &
        &                            ndim      = nDim,         &
        &                            nComp     = nComp,        &
        &                            stateName = 'viscosity'   )
    case('viscous_spongelayer_radial')
      write(logUnit(3),*) '  is a radial viscous sponge layer'
      nDim = 3
      call tem_load_spongelayer_radial( me        = me%spongeRadial, &
        &                               conf      = conf,            &
        &                               thandle   = thandle,         &
        &                               nDim      = nDim,            &
        &                               nComp     = nComp,           &
        &                               stateName = 'viscosity'      )
    case('viscous_spongelayer_radial_2d')
      write(logUnit(3),*) '  is a 2d radial viscous sponge layer'
      nDim = 2
      call tem_load_spongelayer_radial( me        = me%spongeRadial, &
        &                               conf      = conf,            &
        &                               thandle   = thandle,         &
        &                               nDim      = nDim,            &
        &                               nComp     = nComp,           &
        &                               stateName = 'viscosity'      )
    case('pml')
      call tem_load_pmlLayer( conf    = conf,    &
        &                     thandle = thandle, &
        &                     me      = me%pml   )
    case('cylindricalwave')
      call tem_load_cylindricalWave( conf    = conf,              &
        &                            thandle = thandle,           &
        &                            me      = me%cylindricalWave )
    case('polygon_material','polygon_material_3d')
      call tem_polygon_material_load( conf    = conf,               &
        &                             thandle = thandle,            &
        &                             me      = me%polygon_material )
    case('rectangular', 'gate')
      call aot_get_val( L       = conf,       &
        &               thandle = thandle,    &
        &               key     = 'ly',       &
        &               val     = me%rect_ly, &
        &               ErrCode = iError,     &
        &               default = 0.5_rk      )
      call aot_get_val( L       = conf,       &
        &               thandle = thandle,    &
        &               key     = 'lz',       &
        &               val     = me%rect_lz, &
        &               ErrCode = iError,     &
        &               default = 0.5_rk      )
      write(logUnit(3),"(A)") '   this is a predefined rectangular shape'
      write(logUnit(3),"(A)") '     ly: '//trim(tem_toStr(me%rect_ly))
      write(logUnit(3),"(A)") '     lz: '//trim(tem_toStr(me%rect_lz))
    case default
      write(logUnit(1),*)'ERROR in definition of the spatial'//&
        &            ' boundary Conditions:'
      write(logUnit(1),*)'Selected an unknown spatial function: '//&
        &            trim(me%kind)
      call tem_abort()
    end select

  end subroutine load_spatial_predefined
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This subroutine load 3d parabola type variables from LUA file.
  !!
  !! Specify `amplitude` to size of nComp in spatial block.
  !! If amplitude is not defined, it set to 1.0_rk for all components
  !! Valid definition:
  !!
  !! 1. 2D parabola
  !!```lua
  !! spatial = {predefined='parabol', shape = { center/origin={-5.0, 0.0, 0.0},
  !!                                            vec={0.0, 1.0, 0.0}},
  !!                                            amplitude = 0.01}
  !! where, center - line center.
  !!        origin - line origin. (define either center or origin)
  !!        vec    - length of the line.
  !!```
  !! Parabolic 2D profile at channel inlet.
  !!
  !! ![parabol2d](|media|/parabol2d.png)
  !!
  !! 2. 3D parabola
  !!```lua
  !! spatial = {
  !!   predefined = 'parabol',
  !!   shape = {
  !!     center/origin = {-5.0, 0.0, 0.0},
  !!     vec = {
  !!       {0.0, 1.0, 0.0},
  !!       {1.0, 0.0, 0.0}
  !!     },
  !!     amplitude = 0.01
  !!  }
  !! where, center - plane center.
  !!        origin - plane origin. (define either center or origin)
  !!        vec    - length of the plane in 2 direction.
  !!        amplitude - maximum value
  !!```
  !! Parabolic 3D profile at channel inlet.
  !! ![parabol3d](|media|/parabolic3d.png)
  !!
  !! Example:
  !! Inlet velocity at inlet plane with velocity along x-dir with a maximum
  !! velocity 0.08.
  !! Inlet plane is normal to x-dir and located at offset of {-5,0,0} from
  !! origin. The plane spans -1 to 1 in both y and z-dir.
  !!
  !!```lua
  !! boundary_condition = {
  !!   {
  !!     label = 'inlet',
  !!     kind = 'inlet_ubb',
  !!     velocityX = {
  !!       kind = 'combined',
  !!       spatial = {
  !!         predefined='parabol',
  !!         shape = {
  !!           origin = {-5.0, 0.0, 0.0},
  !!           vec={ {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0} }
  !!         },
  !!         amplitude = 0.01
  !!       }
  !!     },
  !!     velocityY = 0.0
  !!     velocityZ = 0.0
  !!   }
  !! }
  !!```
  !!```
  !!  _______vecB____
  !! |       |      |
  !! |       |      |
  !! |       |______|vecA
  !! |     center   |
  !! |              |
  !! |______________|_
  !!```
  !!
  subroutine load_spatial_parabol( me, conf, thandle, nComp )
    ! -------------------------------------------------------------------- !
    !> parabola3d spatial datas
    type(spatial_parabol_type),intent(out) :: me
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: thandle
    !> Number of components
    integer, intent(in) :: nComp
    ! -------------------------------------------------------------------- !
    integer :: iError
    ! -------------------------------------------------------------------- !

    ! load geometry
    call tem_load_shape( conf    = conf,       &
      &                  parent  = thandle,    &
      &                  me      = me%geometry )

    if (size(me%geometry%canoND) /= 1) then
      write(logUnit(1),*) 'Error: Requires single shape for spatial "parabol"'
      call tem_abort()
    end if

    call load_spatial_asConst( const   = me%amplitude, &
      &                        conf    = conf,         &
      &                        errCode = iError,       &
      &                        parent  = thandle,      &
      &                        key     = 'amplitude',  &
      &                        nComp   = nComp         )

    if ( iError /= 0 ) then
      write(logUnit(1),*) 'Warning: amplitude is not defined.'
      write(logUnit(1),*) '         Set to default value 1.0 for all components'
      me%amplitude = ref_amp
    end if

    if (nComp == 1) then
      write(logUnit(5),"(A)") ' amplitude = ' &
        & // trim(tem_toStr(me%amplitude(1)))
    else
      write(logUnit(5),"(A)") ' amplitude = ' &
        & // trim(tem_toStr(me%amplitude(1:nComp),','))
    end if

  end subroutine load_spatial_parabol
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Reading the definition for a random distribution function.
  !!
  !! This spatial function will return random numbers within a given interval.
  !! The interval has to be specified by 'min' and 'max' in the Lua
  !! configuration.
  !! Resulting values will be given by: min + (max-min)*random_number
  !! The range defaults to min=0.0 - max=1.0.
  !!
  subroutine load_spatial_random( me, conf, thandle )
    ! -------------------------------------------------------------------- !
    !> random spatial specification
    type(spatial_random_type),intent(out) :: me
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: thandle
    ! -------------------------------------------------------------------- !
    integer :: iError
    ! -------------------------------------------------------------------- !

    ! val_min
    call aot_get_val( L       = conf,       &
      &               thandle = thandle,    &
      &               key     = 'min',      &
      &               val     = me%val_min, &
      &               ErrCode = iError,     &
      &               default = 0.0_rk      )
    if ( btest( iError, aoterr_Fatal ) ) then
      write(logUnit(1),*) 'FATAL Error occured in definition of the spatial' &
        &                 // ' function'
      write(logUnit(1),*) 'while retrieving the min for random numbers!'
      call tem_abort()
    end if
    if( btest( iError, aoterr_NonExistent ))then
      write(logUnit(1),*) ' WARNING: min is non-existent.' &
        &                 // ' Using default value 0.0.'
    end if
    write(logUnit(1),*) ' * val_min =', me%val_min

    ! val_max
    call aot_get_val( L       = conf,       &
      &               thandle = thandle,    &
      &               key     = 'max',      &
      &               val     = me%val_max, &
      &               ErrCode = iError,     &
      &               default = 1.0_rk      )
    if( btest( iError, aoterr_Fatal ))then
      write(logUnit(1),*) 'FATAL Error occured in definition of the spatial' &
        & // ' function'
      write(logUnit(1),*) 'while retrieving the max for random numbers!'
      call tem_abort()
    end if
    if( btest( iError, aoterr_NonExistent ))then
      write(logUnit(1),*) ' WARNING: max is non-existent.' &
        & // ' Using default value 1.0.'
    end if
    write(logUnit(1),*) ' * val_max =', me%val_max

    me%val_range = me%val_max - me%val_min

  end subroutine load_spatial_random
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !>  This function computes 3d parabola profile from treeIDs of an element.
  !!
  !! This profile is defined
  !! by element barycentric coordinate and 3d parabola parameters.
  !! 3D parabola profile at given plane is computed in the following way:
  !!
  !! - Project barycentric coordinate vector in a given plane via
  !! \( \alpha = ((baryCoord - center) \cdot vecA)/\sqrt{ vecA\cdot vecA }\)
  !! - Compute spatial value using
  !! \( res = (1+\alpha)*(1-\alpha)*(1+\beta)*(1-\beta) \)
  !! Example:
  !! Parabolic 3D profile at channel inlet.
  !!
  !! ![parabolic3d](|media|/parabolic3D.png)
  !!
  function tem_spatial_parabol3d_for_treeIds( me, treeIds, tree, n ) result(res)
    ! -------------------------------------------------------------------- !
    !> contains plane parameter for 3d parabola
    type( tem_canonicalND_type ) :: me
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> number of return values
    integer, intent(in) :: n
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value of a function
    real(kind=rk) :: res(n)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: coord(3), alpha, beta, diff(3)
    real(kind=rk) :: vecAsqr, vecBsqr
    real(kind=rk) :: center(3), halfvec(2,3)
    !loop variables
    integer :: iDir, jDir
    ! -------------------------------------------------------------------- !

    jDir = 0
    do iDir = 1, 3
      if( me%active( iDir )) then
        jDir = jDir + 1
        halfvec( jDir, : ) = me%vec( :, iDir ) / 2._rk
      end if
    end do


    center = me%origin + halfvec(1, :) + halfvec(2, :)

    vecAsqr = dot_product( halfvec(1, :), halfvec(1, :) )
    vecBsqr = dot_product( halfvec(2, :), halfvec(2, :) )

    !loop over number of return values
    do iDir = 1, n

      !barycentric coordinate
      coord = tem_BaryOfId( tree, treeIds(iDir) )

      !distance between parabola center and barycentric coordinates
      diff = coord - center
      !projection of diff in a plane on vecA
      alpha =  dot_product( diff, halfvec(1, :) ) / vecAsqr

      !projection of diff in a plane on vecB
      beta =  dot_product( diff, halfvec(2, :) ) / vecBsqr

      res( iDir )  = ( 1.0_rk-alpha ) * ( 1.0_rk+alpha ) &
        &          * ( 1.0_rk-beta ) * ( 1.0_rk+beta )

      if ( abs(alpha) .gt. 1.0_rk .or. abs(beta) .gt. 1.0_rk ) then
        res( iDir ) = 0.0_rk
      end if
    end do

  end function tem_spatial_parabol3d_for_treeIds
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !>  This function computes 3d parabola profile from coord of an element.
  !!
  !! This profile is defined
  !! by element barycentric coordinate and 3d parabola parameters.
  !! 3D parabola profile at given plane is computed in the following way:
  !!
  !! - Project barycentric coordinate vector in a given plane via
  !! \( \alpha = ((baryCoord - center) \cdot vecA)/\sqrt{ vecA\cdot vecA }\)
  !! - Compute spatial value using
  !! \( res = (1+\alpha)*(1-\alpha)*(1+\beta)*(1-\beta) \)
  !! Example:
  !! Parabolic 3D profile at channel inlet.
  !! ![parabolic 3D](|media|/parabolic3D.png)
  !!
  !! @todo use projection of point on line and point on plane
  !!
  function tem_spatial_parabol3d_for_coord( me, coord, n ) result(res)
    ! -------------------------------------------------------------------- !
    !> contains parameter for 3d parabola
    type( tem_canonicalND_type ) :: me
    !> the number of return values
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent(in) :: coord(n,3)
    !> return value of a function
    real(kind=rk) :: res(n)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: alpha, beta, diff(3)
    real(kind=rk) :: vecAsqr, vecBsqr
    real(kind=rk) :: center(3), halfvec(2,3)
    integer :: iDir, jDir
    ! -------------------------------------------------------------------- !

    jDir = 0
    do iDir = 1, 3
      if( me%active( iDir )) then
        jDir = jDir + 1
        halfvec(jDir, :) = me%vec(:, iDir) / 2._rk
      end if
    end do

    center = me%origin + halfvec(1, :) + halfvec(2, :)

    vecAsqr = dot_product( halfvec(1, :), halfvec(1, :) )
    vecBsqr = dot_product( halfvec(2, :), halfvec(2, :) )

    !loop over number of return values
    do iDir = 1, n

      !distance between parabola center and barycentric coordinates
      diff = coord(iDir,:) - center
      !projection of diff in a plane on vecA
      alpha =  dot_product( diff, halfvec(1, :) ) / vecAsqr

      ! projection of diff in a plane on vecB
      beta =  dot_product( diff, halfvec(2, :) ) / vecBsqr

      res(iDir)  = (1.0_rk - alpha) * (1.0_rk + alpha) &
        &         * (1.0_rk - beta) * (1.0_rk + beta)

      if ( (abs(alpha) .gt. 1.d0) .or. (abs(beta) .gt. 1.d0) ) then
        res(iDir) = 0.0_rk
      end if

    end do

  end function tem_spatial_parabol3d_for_coord
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !>  This function computes 2d parabola profile from treeIds of elements
  !!
  !! This profile is defined
  !! by element barycentric coordinate and 2d parabola parameters.
  !! 2D parabola profile at given plane is computed in the following way:
  !!
  !! - Project barycentric coordinate vector in a given plane via
  !! \( alpha = ((baryCoord - center) \cdot vecA)/\sqrt(vecA\cdot vecA)\)
  !! - Compute spatial value using
  !! \( res = (1+alpha)*(1-alpha) \)
  !! Actual parabola 2d formula: \f$ (y-y_0) = a(x-x_0)^2 \)
  !! parabola open towards opposite x and with vertex at (1,0) is given by
  !! \( y = (1-x)^2 \)
  !!
  function tem_spatial_parabol2d_for_treeIds( me, treeIds, tree, n ) result(res)
    ! -------------------------------------------------------------------- !
    !> contains parameters for 2d parabola
    type( tem_canonicalND_type ) :: me
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> number of return values
    integer, intent(in) :: n
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value of a function
    real(kind=rk) :: res(n)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: coord(3), alpha, diff(3)
    real(kind=rk) :: vecAsqr, center(3), halfdir(3)
    !loop variables
    integer :: iDir
    ! -------------------------------------------------------------------- !

    ! since this is a line only one of the three vectors in vec
    ! is active
    halfdir = 0.0_rk
    do iDir = 1, 3
      if( me%active( iDir )) halfdir = me%vec( :, iDir )/2.0_rk
    end do
    center = me%origin + halfdir
    vecAsqr = dot_product( halfdir, halfdir )
    !loop over number of return values
    do iDir = 1, n

      !barycentric coordinate
      coord = tem_BaryOfId( tree, treeIds(iDir) )

      !distance between parabola center and barycentric coordinates
      diff = coord - center

      !projection of diff in a plane on vecA
      alpha = dot_product( diff, halfdir ) / vecAsqr

      res(iDir)  = ( 1.0_rk-alpha ) * ( 1.0_rk+alpha )

      if ( abs(alpha) .gt. 1.0d0 ) then
        res(iDir) = 0.0_rk
      end if

    end do

  end function tem_spatial_parabol2d_for_treeIds
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !>  This function computes 2d parabola profile from coord of elements
  !!
  !! This profile is defined
  !! by element barycentric coordinate and 2d parabola parameters.
  !! 2D parabola profile at given plane is computed in the following way:
  !!
  !! - Project barycentric coordinate vector in a given plane via
  !! \( alpha = ((baryCoord - center) \cdot vecA)/\sqrt(vecA\cdot vecA)\)
  !! - Compute spatial value using
  !! \( res = (1+alpha)*(1-alpha) \)
  !!
  function tem_spatial_parabol2d_for_coord( me, coord, n ) result(res)
    ! -------------------------------------------------------------------- !
    !> contains line parameters for 2d parabola
    type( tem_canonicalND_type ) :: me
    !> number of return values
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent(in) :: coord(n,3)
    !> return value of a function
    real(kind=rk) :: res(n)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: alpha, diff(3)
    real(kind=rk) :: vecAsqr, center(3), halfdir(3)
    ! loop variable
    integer :: iDir
    ! -------------------------------------------------------------------- !

    halfdir = 0.0_rk
    do iDir = 1, 3
      if( me%active(iDir) ) halfdir = me%vec(:, iDir)/2.0_rk
    end do

    vecAsqr = dot_product(halfdir, halfdir)

    center = me%origin + halfdir
    !loop over number of return values
    do iDir = 1, n

      !distance between parabola center and barycentric coordinates
      diff = coord(iDir,:) - center

      !projection of diff in a plane on vecA
      alpha =  dot_product( diff, halfdir ) / vecAsqr

      res(iDir) = (1.0_rk - alpha) * (1.0_rk + alpha)

      if ( abs(alpha) .gt. 1.0d0 ) then
        res(iDir) = 0.0_rk
      end if

    end do

  end function tem_spatial_parabol2d_for_coord
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Produce random numbers
  !!
  function tem_spatial_random_for( me, n ) result(res)
    ! -------------------------------------------------------------------- !
    !> Interval definition to get the random values from.
    type(spatial_random_type), intent(in) :: me
    !> number of return values
    integer, intent(in) :: n
    !> return value of a function
    real(kind=rk) :: res(n)
    ! -------------------------------------------------------------------- !

    call random_number(res)
    res = me%val_min + me%val_range*res

  end function tem_spatial_random_for
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function invokes the Lua function,
  !! in which the barycentric coordinates are computed from treeIds of an
  !! element.
  !!
  !! Lua function defined in the script is
  !! connected to the conf handle and return the result of the function.
  !! The Lua function takes barycentric coordinate as input argument
  !! i.e fun_name(x,y,z)
  function tem_spatial_lua_for_treeIds( fun_ref, conf, treeIds, tree, n ) &
    &                                 result(res)
    ! -------------------------------------------------------------------- !
    !> Lua reference to the function to evaluate.
    integer, intent(in) :: fun_ref
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> number of return values
    integer, intent(in) :: n
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> lua state
    type(flu_State) :: conf
    !> return value
    real(kind=rk) :: res(n)
    ! -------------------------------------------------------------------- !
    type(aot_fun_type) :: fun
    integer :: iError
    integer :: iDir, jDir
    real(kind=rk) :: coord(3)
    ! -------------------------------------------------------------------- !

    call aot_fun_open(L=conf, fun=fun, ref=fun_ref)

    do iDir=1,n
      coord = tem_BaryOfId( tree, treeIds(iDir) )
      do jDir=1,3
        call aot_fun_put(L=conf, fun=fun, arg=coord(jDir))
      end do
      call aot_fun_do(L=conf, fun=fun, nresults=1)
      call aot_top_get_val(L=conf, val=res(iDir), ErrCode=iError)
      if ( btest(iError,aoterr_Fatal) ) then
        write(logunit(0),*) "ERROR Obtaining a spatial function"
        write(logunit(0),*) "Probably wrong number of components returned"
        write(logunit(0),*) "or a scalar was return as a lua table"
        write(logunit(0),*) 'Expected nComp: 1'
        write(logunit(0),*) 'ErrorCodes: ', iError
        write(logunit(0),*) "Check return values of your Lua functions!"
        call tem_abort()
      end if
    end do

    call aot_fun_close(L=conf, fun=fun)

  end function tem_spatial_lua_for_treeIds
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !>  This function invokes the vectorial Lua function,
  !! in which the barycentric coordinates are computed from treeIds of an
  !! element.
  !!
  !! Lua function defined in the script is
  !! connected to the conf handle and return the result of the function.
  !! The Lua function takes barycentric coordinate as input argument
  !! i.e fun_name(x,y,z)
  function tem_spatial_lua_vector_for_treeIds( fun_ref, conf, treeIds, tree, &
    &                                          n, ncomp ) result(res)
    ! -------------------------------------------------------------------- !
    !> Lua reference to the function to evaluate.
    integer, intent(in) :: fun_ref
    !> global treelm mesh
    type(treelmesh_type), intent(in) ::tree
    !> number of return values
    integer, intent(in) :: n
    !> Number of components in each returned value
    integer, intent(in) :: ncomp
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> lua state
    type(flu_State) :: conf
    !> return value
    real(kind=rk) :: res(n,ncomp)
    ! -------------------------------------------------------------------- !
    type(aot_fun_type) :: fun
    integer :: iError(nComp)
    integer :: iDir, jDir
    real(kind=rk) :: coord(3)
    ! -------------------------------------------------------------------- !

    call aot_fun_open(L=conf, fun=fun, ref=fun_ref)

    do iDir=1,n
      coord = tem_BaryOfId( tree, treeIds(iDir) )
      do jDir=1,3
        call aot_fun_put(L=conf, fun=fun, arg=coord(jDir))
      end do
      call aot_fun_do(L=conf, fun=fun, nresults=1)
      call aot_top_get_val(L=conf, val=res(iDir,:), ErrCode=iError)
      if ( any(btest(iError,aoterr_Fatal)) ) then
        write(logunit(0),*) "ERROR Obtaining a spatial function"
        write(logunit(0),*) "Probably wrong number of components returned"
        write(logunit(0),*) "or a scalar was return as a lua table"
        write(logunit(0),*) 'Expected nComp: ', nComp
        write(logunit(0),*) 'ErrorCodes: ', iError
        write(logunit(0),*) "Check return values of your Lua functions!"
        call tem_abort()
      end if
    end do

    call aot_fun_close(L=conf, fun=fun)

  end function tem_spatial_lua_vector_for_treeIds
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function invokes the Lua function,
  !! which takes barycentric coordinates of an element.
  !!
  !! Lua function defined in the script is
  !! connected to the conf handle and return the result of the function.
  !! The Lua function takes barycentric coordinate as input argument
  !! i.e fun_name(x,y,z)
  !!
  function tem_spatial_lua_for_coord( fun_ref, conf, coord, n ) result(res)
    ! -------------------------------------------------------------------- !
    !> Lua reference to the function to evaluate.
    integer, intent(in) :: fun_ref
    !> number of return values
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent(in) :: coord(n,3)
    !> lua state
    type(flu_State) :: conf
    !> return value
    real(kind=rk) :: res(n)
    ! -------------------------------------------------------------------- !
    type(aot_fun_type) :: fun
    integer :: iError
    integer :: iDir, jDir
    ! -------------------------------------------------------------------- !

    call aot_fun_open(L=conf, fun=fun, ref=fun_ref)

    do iDir = 1, n
      do jDir = 1, 3
        call aot_fun_put( L=conf, fun=fun, arg=coord(iDir,jDir) )
      end do
      call aot_fun_do(L=conf, fun=fun, nresults=1)
      call aot_top_get_val(L=conf, val=res(iDir), ErrCode=iError)
      if ( btest(iError,aoterr_Fatal) ) then
        write(logunit(0),*) "ERROR Obtaining a spatial function"
        write(logunit(0),*) "Probably wrong number of components returned"
        write(logunit(0),*) "or a scalar was return as a lua table"
        write(logunit(0),*) 'Expected nComp: 1'
        write(logunit(0),*) 'ErrorCodes: ', iError
        write(logunit(0),*) "Check return values of your Lua functions!"
        call tem_abort()
      end if
    end do

    call aot_fun_close(L=conf, fun=fun)

  end function tem_spatial_lua_for_coord
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function invokes the vectorial Lua function,
  !! which takes barycentric coordinates of an element.
  !!
  !! Lua function defined in the script is
  !! connected to the conf handle and return the result of the function.
  !! The Lua function takes barycentric coordinate as input argument
  !! i.e fun_name(x,y,z)
  !!
  function tem_spatial_lua_vector_for_coord( fun_ref, conf, coord, n, ncomp ) &
    &      result(res)
    ! -------------------------------------------------------------------- !
    !> Lua reference to the function to evaluate.
    integer, intent(in) :: fun_ref
    !> number of return values
    integer, intent(in) :: n
    !> number of components in the resulting vector
    integer, intent(in) :: ncomp
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent(in) :: coord(n,3)
    !> lua state
    type(flu_State), optional :: conf
    !> return value
    real(kind=rk) :: res(n,ncomp)
    ! -------------------------------------------------------------------- !
    type(aot_fun_type) :: fun
    integer :: iError(ncomp)
    integer :: iDir, jDir
    ! -------------------------------------------------------------------- !

    call aot_fun_open(L=conf, fun=fun, ref=fun_ref)

    do iDir = 1, n
      do jDir = 1, 3
        call aot_fun_put( L=conf, fun=fun, arg=coord(iDir,jDir) )
      end do
      call aot_fun_do(L=conf, fun=fun, nresults=1)
      call aot_top_get_val(L=conf, val=res(iDir,:), ErrCode=iError)
      if ( any(btest(iError,aoterr_Fatal)) ) then
        write(logunit(0),*) "ERROR Obtaining a spactial function"
        write(logunit(0),*) "Probably wrong number of components returned"
        write(logunit(0),*) "or a scalar was return as a lua table"
        write(logunit(0),*) 'Expected nComp: ', nComp
        write(logunit(0),*) 'ErrorCodes: ', iError
        write(logunit(0),*) "Check return values of your Lua functions!"
        call tem_abort()
      end if
    end do

    call aot_fun_close(L=conf, fun=fun)

  end function tem_spatial_lua_vector_for_coord
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function invokes the Lua function,
  !! which takes tem_grwPoints_type and evaluate a function at a point of
  !! given idx in grwPnt.
  !!
  !! Lua function defined in the script is
  !! connected to the conf handle and return the result of the function.
  !! The Lua function takes barycentric coordinate as input argument
  !! i.e fun_name(x,y,z)
  !!
  function tem_spatial_lua_for_index( fun_ref, conf, grwPnt, idx, nVals ) &
    &                                 result(res)
    ! -------------------------------------------------------------------- !
    !> Lua reference to the function to evaluate.
    integer, intent(in) :: fun_ref
    !> number of return values
    integer, intent(in) :: nVals
    !> growing array of all spatial point of a variable
    type(tem_grwPoints_type), intent(in) :: grwPnt
    !> Index position to return a pre-store value or to compute
    integer, intent(in) :: idx(nVals)
    !> lua state
    type(flu_State) :: conf
    !> return value
    real(kind=rk) :: res(nVals)
    ! -------------------------------------------------------------------- !
    type(aot_fun_type) :: fun
    integer :: iError
    integer :: iDir, jDir
    real(kind=rk) :: coord(3)
    ! -------------------------------------------------------------------- !

    call aot_fun_open(L=conf, fun=fun, ref=fun_ref)

    do iDir = 1, nVals
      coord(:) =  (/ grwPnt%coordX%val( idx(iDir) ), &
        &            grwPnt%coordY%val( idx(iDir) ), &
        &            grwPnt%coordZ%val( idx(iDir) ) /)
      do jDir = 1, 3
        call aot_fun_put( L=conf, fun=fun, arg=coord(jDir) )
      end do
      call aot_fun_do(L=conf, fun=fun, nresults=1)
      call aot_top_get_val(L=conf, val=res(iDir), ErrCode=iError)
      if ( btest(iError,aoterr_Fatal) ) then
        write(logunit(0),*) "ERROR Obtaining a spatial function"
        write(logunit(0),*) "Probably wrong number of components returned"
        write(logunit(0),*) "or a scalar was return as a lua table"
        write(logunit(0),*) 'Expected nComp: 1'
        write(logunit(0),*) 'ErrorCodes: ', iError
        write(logunit(0),*) "Check return values of your Lua functions!"
        call tem_abort()
      end if
    end do

    call aot_fun_close(L=conf, fun=fun)

  end function tem_spatial_lua_for_index
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function invokes the vectorial Lua function,
  !! which takes tem_grwPoints_type and evaluate a function at a point of
  !! given idx in grwPnt.
  !!
  !! Lua function defined in the script is
  !! connected to the conf handle and return the result of the function.
  !! The Lua function takes barycentric coordinate as input argument
  !! i.e fun_name(x,y,z)
  !!
  function tem_spatial_lua_vector_for_index( fun_ref, conf, grwPnt, idx, &
    &                                        nVals, nComps ) result(res)
    ! -------------------------------------------------------------------- !
    !> Lua reference to the function to evaluate.
    integer, intent(in) :: fun_ref
    !> number of return values
    integer, intent(in) :: nVals
    !> number of components per returned value
    integer, intent(in) :: nComps
    !> growing array of all spatial point of a variable
    type(tem_grwPoints_type), intent(in) :: grwPnt
    !> Index position to return a pre-store value or to compute
    integer, intent(in) :: idx(nVals)
    !> lua state
    type(flu_State) :: conf
    !> return value
    real(kind=rk) :: res(nVals,nComps)
    ! -------------------------------------------------------------------- !
    type(aot_fun_type) :: fun
    integer :: iError(ncomps)
    integer :: iDir, jDir
    real(kind=rk) :: coord(3)
    ! -------------------------------------------------------------------- !

    call aot_fun_open(L=conf, fun=fun, ref=fun_ref)

    do iDir = 1, nVals
      coord(:) =  (/ grwPnt%coordX%val( idx(iDir) ), &
        &            grwPnt%coordY%val( idx(iDir) ), &
        &            grwPnt%coordZ%val( idx(iDir) ) /)
      do jDir = 1, 3
        call aot_fun_put( L=conf, fun=fun, arg=coord(jDir) )
      end do
      call aot_fun_do(L=conf, fun=fun, nresults=1)
      call aot_top_get_val(L=conf, val=res(iDir,:), ErrCode=iError)
      if ( any(btest(iError,aoterr_Fatal)) ) then
        write(logunit(0),*) "ERROR Obtaining a spactial function"
        write(logunit(0),*) "Probably wrong number of components returned"
        write(logunit(0),*) "or a scalar was return as a lua table"
        write(logunit(0),*) 'Expected nComp: ', nComps
        write(logunit(0),*) 'ErrorCodes: ', iError
        write(logunit(0),*) "Check return values of your Lua functions!"
        call tem_abort()
      end if
    end do

    call aot_fun_close(L=conf, fun=fun)

  end function tem_spatial_lua_vector_for_index
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function invokes different spatial boundary kinds like constant, lua
  !! function or predefined Fortran function for given treeIDs
  !!
  !! If a spatial block is not defined, default value is set to 1.0
  !! or default value passed while loading tem_load_spatial.
  !! If both spatial and temporal block are not defined in the lua file, the
  !! return value = 1.0_rk.
  !! based spatial_kind(kind).
  !!
  !! 1. const - set constant value
  !! 2. lua_fun - lua function
  !! 3. gausspulse - fortran gauss pulse function
  !! 4. 2dcrvp     - fortran spinning vortex function
  !! 5. parabol    - fotran parabolic function
  !!
  function tem_spatial_for_treeIDs( me, treeIds, tree, n) result( res )
    ! -------------------------------------------------------------------- !
    !> spatial type for given boundary state
    type( tem_spatial_type ) :: me
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> number of return values
    integer, intent(in) :: n
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value of a function
    real( kind=rk ) :: res(n)
    ! -------------------------------------------------------------------- !

    select case( trim(adjustl(me%kind)) )
    case( 'none', 'const' )
      res = me%const(1)

    case( 'lua_fun' )
      res = tem_spatial_lua_for( me%lua_fun_ref, me%conf, treeIds, tree, n)

    case( 'parabol' )
      select case( trim(adjustl(me%parabol%geometry%canoND(1)%kind)) )
      case('line')
        res = tem_spatial_parabol2d_for( me%parabol%geometry%canoND(1), &
          &                              treeIds,                       &
          &                              tree,                          &
          &                              n                              )
      case('plane')
        res = tem_spatial_parabol3d_for( me%parabol%geometry%canoND(1), &
          &                              treeIds,                       &
          &                              tree,                          &
          &                              n                              )
      end select
      res = res*me%parabol%amplitude(1)

    case ('random')
      res = tem_spatial_random_for(me%random, n)

    case( 'spongelayer_plane', 'spongelayer_plane_2d', 'spongelayer_plane_1d' )
      res = tem_spongeLayer_plane_for(me%spongePlane, treeids, tree, n)

    case( 'spongelayer_box' )
      res = tem_spongeLayer_box_for(me%spongeBox, treeids, tree, n)

    case( 'spongelayer_box_2d' )
      res = tem_spongeLayer_box2d_for(me%spongeBox, treeids, tree, n)

    case( 'spongelayer_radial_2d')
      res = tem_spongeLayer_radial_for( &
        & me      = me%spongeRadial,    &
        & treeids = treeids,            &
        & tree    = tree,               &
        & n       = n,                  &
        & nDim    = 2                   )

    case( 'spongelayer_radial')
      res = tem_spongeLayer_radial_for( &
        & me      = me%spongeRadial,    &
        & treeids = treeids,            &
        & tree    = tree,               &
        & n       = n,                  &
        & nDim    = 3                   )

    case( 'viscous_spongelayer_plane' )
      res = tem_viscSpongeLayer_plane_for(me%spongePlane, treeids, tree, n)

    case( 'viscous_spongelayer_box' )
      res = tem_viscSpongeLayer_box_for(me%spongeBox, treeids, tree, n)

    case( 'viscous_spongelayer_box_2d' )
      res = tem_viscSpongeLayer_box2d_for(me%spongeBox, treeids, tree, n)

    case( 'viscous_spongelayer_radial_2d')
      res = tem_viscSpongeLayer_radial_for( &
        & me      = me%spongeRadial,        &
        & treeids = treeids,                &
        & tree    = tree,                   &
        & n       = n,                      &
        & nDim    = 2                       )

    case( 'viscous_spongelayer_radial')
      res = tem_viscSpongeLayer_radial_for( &
        & me      = me%spongeRadial,        &
        & treeids = treeids,                &
        & tree    = tree,                   &
        & n       = n,                      &
        & nDim    = 3                       )

    case default
      call tem_abort(                                                   &
        & 'ERROR: Unknown spatial function in tem_spatial_for_treeIDs.' )

    end select

   end function tem_spatial_for_treeIDs
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function tem_spatial_vector_for_treeIDs( me, treeIds, tree, n, ncomp) &
    &      result(res)
    ! -------------------------------------------------------------------- !
    !> spatial type for given boundary state
    type(tem_spatial_type) :: me
    !> global treelm mesh
    type(treelmesh_type), intent(in) ::tree
    !> number of return values
    integer, intent(in) :: n
    !> number of components in the resulting vector
    integer, intent(in) :: ncomp
    !> treeIds of elements in given level
    integer(kind=long_k), intent(in) :: treeIds(n)
    !> return value of a function
    real(kind=rk) :: res(n,ncomp)
    ! -------------------------------------------------------------------- !
    integer :: i
    ! -------------------------------------------------------------------- !

    select case( trim(adjustl(me%kind)) )
    case( 'none', 'const' )
      do i = 1, nComp
        res(:,i) = me%const(i)
      end do

    case( 'lua_fun' )
      res = tem_spatial_lua_for( fun_ref = me%lua_fun_ref, &
        &                        conf    = me%conf,        &
        &                        treeIds = treeIds,        &
        &                        tree    = tree,           &
        &                        n       = n,              &
        &                        ncomp   = ncomp           )

    case( 'spongelayer_plane', 'spongelayer_plane_2d', 'spongelayer_plane_1d' )
      res = tem_spongeLayer_plane_for(me%spongePlane, nComp, treeids, tree, n)

    case( 'spongelayer_box' )
      res = tem_spongeLayer_box_for(me%spongeBox, nComp, treeids, tree, n)

    case( 'spongelayer_box_2d' )
      res = tem_spongeLayer_box2d_for(me%spongeBox, nComp, treeids, tree, n)

    case( 'spongelayer_radial_2d')
      res = tem_spongeLayer_radial_for( &
        & me      = me%spongeRadial,    &
        & nComp   = nComp,              &
        & treeids = treeids,            &
        & tree    = tree,               &
        & n       = n,                  &
        & nDim    = 2                   )

    case( 'spongelayer_radial')
      res = tem_spongeLayer_radial_for( &
        & me      = me%spongeRadial,    &
        & nComp   = nComp,              &
        & treeids = treeids,            &
        & tree    = tree,               &
        & n       = n,                  &
        & nDim    = 3                   )

    case( 'parabol' )
      select case( trim(adjustl(me%parabol%geometry%canoND(1)%kind)) )
      case('line')
        res(:,1) = tem_spatial_parabol2d_for(        &
          & me      = me%parabol%geometry%canoND(1), &
          & treeIds = treeIds,                       &
          & tree    = tree,                          &
          & n       = n                              )

      case('plane')
        res(:,1) = tem_spatial_parabol3d_for(        &
          & me      = me%parabol%geometry%canoND(1), &
          & treeIds = treeIds,                       &
          & tree    = tree,                          &
          & n       = n                              )
      end select

      do i = 2, nComp
        res(:,i) = res(:,1)*me%parabol%amplitude(i)
      end do

      res(:,1) = res(:,1)*me%parabol%amplitude(1)

    case default
      write(logUnit(1),*)'ERROR: No vectorial routine for spatial functions of'
      write(logUnit(1),*)'       kind "' // trim(adjustl(me%kind)) // '"'
      write(logUnit(1),*)'       available! Have a look in the ' &
        & // 'tem_spatial_module'
      write(logUnit(1),*)'       for implemented vectorial functions.'
      call tem_abort()
    end select

   end function tem_spatial_vector_for_treeIDs
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function invokes different spatial boundary kinds like constant, lua
  !! function or predefined Fortran function for given coord
  !!
  !! If a spatial block is not defined and a temporal block is defined in the
  !! lua file, the return value is either 1.0 or default value provided
  !! in tem_load_spatial.
  !! If both spatial and temporal block are not defined in the lua file, the
  !! return value = 1.0_rk.
  !! based spatial_kind(kind).
  !!
  !! 1. const - set constant value
  !! 2. lua_fun - lua function
  !! 3. gausspulse - fortran gauss pulse function
  !! 4. 2dcrvp     - fortran spinning vortex function
  !! 5. parabol    - fotran parabolic function
  !!
  function tem_spatial_for_coord( me, coord, n ) result( res )
    ! -------------------------------------------------------------------- !
    !> spatial type for given boundary state
    type(tem_spatial_type), intent(in) :: me
    !> number of return values
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent(in) :: coord(n,3)
    !> return value of a function
    real( kind=rk ) :: res(n)
    ! -------------------------------------------------------------------- !

    select case( trim(adjustl(me%kind)) )
    case( 'none', 'const' )
      res = me%const(1)

    case( 'lua_fun' )
      res = tem_spatial_lua_for(me%lua_fun_ref, me%conf, coord, n)

    case ('random')
      res = tem_spatial_random_for(me%random, n)

    case( 'parabol' )
      select case( trim(adjustl(me%parabol%geometry%canoND(1)%kind)) )
      case('line')
        res = tem_spatial_parabol2d_for(           &
          & me    = me%parabol%geometry%canoND(1), &
          & coord = coord,                         &
          & n     = n                              )

      case('plane')
        res = tem_spatial_parabol3d_for(           &
          & me    = me%parabol%geometry%canoND(1), &
          & coord = coord,                         &
          & n     = n                              )

      end select

      res = res*me%parabol%amplitude(1)

    case( 'viscous_spongelayer_plane' )
      res = tem_viscSpongelayer_plane_for( &
        & me    = me%spongePlane,          &
        & coord = coord,                   &
        & n     = n                        )

    case( 'viscous_spongelayer_box' )
      res = tem_viscSpongelayer_box_for( &
        & me    = me%spongeBox,          &
        & coord = coord,                 &
        & n     = n                      )

    case( 'viscous_spongelayer_box_2d' )
      res = tem_viscSpongelayer_box2d_for( &
        & me    = me%spongeBox,            &
        & coord = coord,                   &
        & n     = n                        )

    case( 'viscous_spongelayer_radial_2d' )
      res = tem_viscSpongelayer_radial_for( &
        & me    = me%spongeRadial,          &
        & coord = coord,                    &
        & nDim  = 2,                        &
        & n     = n                         )

    case( 'viscous_spongelayer_radial' )
      res = tem_viscSpongelayer_radial_for( &
        & me    = me%spongeRadial,          &
        & coord = coord,                    &
        & nDim  = 3,                        &
        & n     = n                         )

    case( 'gausspulse' )
      res = ic_gausspulse_for(me%gausspulse, coord, n)

    case( 'crvpX' )
      res = ic_2dcrvpX_for(me%crvp, coord, n)

    case( 'crvpY' )
      res = ic_2dcrvpY_for(me%crvp, coord, n)

    case( 'crvpPressure' )
      res = ic_2dcrvpPressure_for(me%crvp, coord, n)

    case('miescatter_displacementfieldz')
      res = tem_eval_miescatter_displz( me    = me%mie_fun, &
        &                               coord = coord,      &
        &                               time  = 0.0_rk,     &
        &                               n     = n           )

    case('miescatter_magneticfieldx')
      res = tem_eval_miescatter_magnx( me    = me%mie_fun, &
        &                              coord = coord,      &
        &                              time  = 0.0_rk,     &
        &                              n     = n           )

    case('miescatter_magneticfieldy')
      res = tem_eval_miescatter_magny( me    = me%mie_fun, &
        &                              coord = coord,      &
        &                              time  = 0.0_rk,     &
        &                              n     = n           )

    case('heaviside_gibbs')
      res = tem_eval_heaviside_gibbs( me    = me%heaviside_gibbs_fun, &
        &                             coord = coord,                  &
        &                             n     = n                       )

    case('cylindricalwave')
      res = tem_eval_cylindricalWave( me    = me%cylindricalWave, &
        &                             coord = coord,              &
        &                             time  = 0.0_rk,             &
        &                             n     = n                   )

    case('tgv_p')
      res = ic_tgv_pressure_for( me%tgv, coord, n )

    case('tgv_ux')
      res = ic_tgv_ux_for( me%tgv, coord, n )

    case('tgv_uy')
      res = ic_tgv_uy_for( me%tgv, coord, n )

    case('tgv_sxx')
      res = ic_tgv_sxx_for( me%tgv, coord, n )

    case('tgv_syy')
      res = ic_tgv_syy_for( me%tgv, coord, n )

    case('tgv_sxz')
      res = ic_tgv_sxz_for( me%tgv, coord, n )

    case('tgv_syz')
      res = ic_tgv_syz_for( me%tgv, coord, n )

    case('polygon_material')
      res = tem_eval_polygon_material_scal( me    = me%polygon_material, &
        &                                   coord = coord,               &
        &                                   n     = n                    )

    case('polygon_material_3d')
      res = tem_eval_polygon_material_scal_3d( me    = me%polygon_material, &
        &                                      coord = coord,               &
        &                                      n     = n                    )

    case( 'spongelayer_plane', 'spongelayer_plane_2d', 'spongelayer_plane_1d' )
      res = tem_spongeLayer_plane_for(me%spongePlane, coord, n)

    case( 'spongelayer_box' )
      res = tem_spongeLayer_box_for(me%spongeBox, coord, n)

    case( 'spongelayer_box_2d' )
      res = tem_spongeLayer_box2d_for(me%spongeBox, coord, n)

    case( 'spongelayer_radial_2d' )
      res = tem_spongeLayer_radial_for( &
        & me    = me%spongeRadial,      &
        & coord = coord,                &
        & n     = n,                    &
        & nDim  = 2                     )

    case( 'spongelayer_radial' )
      res = tem_spongeLayer_radial_for( &
        & me    = me%spongeRadial,      &
        & coord = coord,                &
        & n     = n,                    &
        & nDim  = 3                     )

    case default
      call tem_abort(                                                 &
        & 'ERROR: Unknown spatial function in tem_spatial_for_coord.' )

    end select

   end function tem_spatial_for_coord
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function invokes different spatial boundary kinds like constant, lua
  !! function or predefined Fortran function for given coord
  !!
  !! If a spatial block is not defined and a temporal block is defined in the
  !! lua file, the return value = ref_value.
  !! If both spatial and temporal block are not defined in the lua file, the
  !! return value = 1.0_rk.
  !! based spatial_kind(kind).
  !!
  !! 1. const - set constant value
  !! 2. lua_fun - lua function
  !! 3. gausspulse - fortran gauss pulse function
  !! 4. 2dcrvp     - fortran spinning vortex function
  !! 5. parabol    - fotran parabolic function
  function tem_spatial_vector_for_coord( me, coord, n, ncomp ) result( res )
    ! -------------------------------------------------------------------- !
    !> spatial type for given boundary state
    type(tem_spatial_type) :: me
    !> number of return values
    integer, intent(in) :: n
    !> number of components per returned value
    integer, intent(in) :: ncomp
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent(in) :: coord(n,3)
    !> return value of a function
    real( kind=rk ) :: res(n,ncomp)
    ! -------------------------------------------------------------------- !
    integer :: i
    ! -------------------------------------------------------------------- !

    select case( trim(adjustl(me%kind)) )
    case( 'none', 'const' )
      do i = 1, nComp
        res(:,i) = me%const(i)
      end do

    case( 'lua_fun' )
      res = tem_spatial_lua_for(me%lua_fun_ref, me%conf, coord, n, ncomp)

    case( 'spongelayer_plane', 'spongelayer_plane_2d', 'spongelayer_plane_1d' )
      res = tem_spongeLayer_plane_for(me%spongePlane, nComp, coord, n)

    case( 'spongeLayer_box' )
      res = tem_spongeLayer_box_for(me%spongeBox, nComp, coord, n)

    case( 'spongeLayer_box_2d' )
      res = tem_spongeLayer_box2d_for(me%spongeBox, nComp, coord, n)

    case( 'spongelayer_radial_2d' )
      res = tem_spongeLayer_radial_for( &
        & me    = me%spongeRadial,      &
        & nComp = nComp,                &
        & coord = coord,                &
        & n     = n,                    &
        & nDim  = 2                     )

    case( 'spongelayer_radial' )
      res = tem_spongeLayer_radial_for( &
        & me    = me%spongeRadial,      &
        & nComp = nComp,                &
        & coord = coord,                &
        & n     = n,                    &
        & nDim  = 3                     )

    case( 'radial_spongelayer_2d' )
      res = tem_spongelayer_radial_for( &
        & me    = me%spongeRadial,      &
        & nComp = nComp,                &
        & coord = coord,                &
        & nDim  = 2,                    &
        & n     = n                     )

    case( 'radial_spongelayer' )
      res = tem_spongelayer_radial_for( &
        & me    = me%spongeRadial,      &
        & nComp = nComp,                &
        & coord = coord,                &
        & nDim  = 3,                    &
        & n     = n                     )

    case( 'pml' )
      res = tem_evaluate_pml(me%pml, nComp, coord, n)

    case( 'polygon_material' )
      res = tem_eval_polygon_material( me    = me%polygon_material, &
        &                              coord = coord,               &
        &                              n     = n                    )

    case( 'polygon_material_3d' )
      res = tem_eval_polygon_material_3d( me    = me%polygon_material, &
        &                                 coord = coord,               &
        &                                 n     = n                    )

    case( 'parabol' )
      select case( trim(adjustl(me%parabol%geometry%canoND(1)%kind)) )
      case( 'line' )
        res(:,1) = tem_spatial_parabol2d_for(      &
          & me    = me%parabol%geometry%canoND(1), &
          & coord = coord,                         &
          & n     = n                              )

      case( 'plane' )
        res(:,1) = tem_spatial_parabol3d_for(      &
          & me    = me%parabol%geometry%canoND(1), &
          & coord = coord,                         &
          & n     = n                              )

      end select

      do i = 2, nComp
        res(:,i) = res(:,1)*me%parabol%amplitude(i)
      end do

      res(:,1) = res(:,1)*me%parabol%amplitude(1)

    case('rectangular', 'gate')
      do i = 1, n
        if( (abs(coord(i,2)) < me%rect_ly)       &
          & .and. (abs(coord(i,3)) < me%rect_lz) ) then
          res(i,1) = 1.0_rk
        else
          res(i,1) = 0.0_rk
        end if
      end do

      res(:, 2:nComp) = 0.0_rk

    case default
      write(logUnit(1),*)'ERROR: No vectorial routine for spatial functions of'
      write(logUnit(1),*)'       kind "' // trim(adjustl(me%kind)) // '"'
      write(logUnit(1),*)'       available! Have a look in the ' &
        & // 'tem_spatial_module'
      write(logUnit(1),*)'       for implemented vectorial functions.'
      call tem_abort()

    end select

   end function tem_spatial_vector_for_coord
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function returns pre-stored value at given idx or evaluate a spatial
  !! function for a point at given idx in growing array of points.
  !! Return value is a scalar.
  function tem_spatial_scalar_for_index( me, grwPnt, idx, nVals, iLevel ) &
    &                                    result (res)
    ! -------------------------------------------------------------------- !
    !> spatial type
    type(tem_spatial_type), intent(in) :: me
    !> number of return values
    integer, intent(in) :: nVals
    !> growing array of all spatial point of a variable
    type(tem_grwPoints_type), intent(in) :: grwPnt
    !> Index position to return a pre-store value or to compute
    integer, intent(in) :: idx(nVals)
    !> return value of a function
    real( kind=rk ) :: res(nVals)
    !> Level to which the evaluated values to be returned
    integer, intent(in) :: iLevel
    ! -------------------------------------------------------------------- !
    integer :: iVal
    real(kind=rk) :: coord(1,3), res_tmp(1)
    ! -------------------------------------------------------------------- !

    if (me%isStored) then
      res(1:nVals) = me%valOnLvl(iLevel)%evalVal%val( idx(1:nVals) )
    else
      select case (trim(me%kind))
      case ('none', 'const')
        res = me%const(1)

      case ('lua_fun')
        res = tem_spatial_lua_for( fun_ref = me%lua_fun_ref, &
          &                        conf    = me%conf,        &
          &                        grwPnt  = grwPnt,         &
          &                        idx     = idx,            &
          &                        nVals   = nVals           )

      case default
        do iVal = 1, nVals
          coord(1,:) =  (/ grwPnt%coordX%val( idx(iVal) ), &
            &              grwPnt%coordY%val( idx(iVal) ), &
            &              grwPnt%coordZ%val( idx(iVal) ) /)

          res_tmp = tem_spatial_for_coord( me    = me,    &
            &                              coord = coord, &
            &                              n     = 1      )
          res(iVal) = res_tmp(1)
        end do

      end select
    end if

  end function tem_spatial_scalar_for_index
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function returns pre-stored value at given idx or evaluate a spatial
  !! function for a point at given idx.
  !! Return value is a vector.
  function tem_spatial_vector_for_index( me, grwPnt, idx, nVals, iLevel, &
    &                                    nComps ) result (res)
    ! -------------------------------------------------------------------- !
    !> spatial type
    type(tem_spatial_type), intent(in) :: me
    !> number of return values
    integer, intent(in) :: nVals
    !> number of components per returned value
    integer, intent(in) :: nComps
    !> growing array of all spatial point of a variable
    type(tem_grwPoints_type), intent(in) :: grwPnt
    !> Index position to return a pre-store value or to compute
    integer, intent(in) :: idx(nVals)
    !> return value of a function
    real( kind=rk ) :: res(nVals, nComps)
    !> Level to which the evaluated values to be returned
    integer, intent(in) :: iLevel
    ! -------------------------------------------------------------------- !
    integer :: iVal, iComp, offset
    real(kind=rk) :: coord(1,3), res_tmp(1, nComps)
    ! -------------------------------------------------------------------- !

    if (me%isStored) then

      do iVal = 1, nVals
        offset = (idx(iVal)-1)*nComps
        res(iVal, :) = me%valOnLvl(iLevel)%evalVal%val( offset+1       &
          &                                             :offset+nComps )
      end do

    else

      select case (trim(me%kind))
      case ('none', 'const')
        do iComp = 1, nComps
          res(:, iComp) = me%const(iComp)
        end do

      case ('lua_fun')
        res = tem_spatial_lua_for( fun_ref = me%lua_fun_ref, &
          &                        conf    = me%conf,        &
          &                        grwPnt  = grwPnt,         &
          &                        idx     = idx,            &
          &                        nVals   = nVals,          &
          &                        nComps  = nComps          )

      case default
        do iVal = 1, nVals
          coord(1,:) =  (/ grwPnt%coordX%val( idx(iVal) ), &
            &              grwPnt%coordY%val( idx(iVal) ), &
            &              grwPnt%coordZ%val( idx(iVal) ) /)

          res_tmp = tem_spatial_vector_for_coord( me    = me,    &
            &                                     coord = coord, &
            &                                     n     = 1,     &
            &                                     nComp = nComps )
          res(iVal,:) = res_tmp(1,:)
        end do

      end select

    end if

  end function tem_spatial_vector_for_index
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This routine evaluate scalar spatial function and store its value in
  !! growing array
  subroutine tem_spatial_scalar_storeVal( me, coord, nVals, iLevel )
    ! -------------------------------------------------------------------- !
    !> spatial type for given boundary state
    type(tem_spatial_type), intent(inout) :: me
    !> number of return values
    integer, intent(in) :: nVals
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent(in) :: coord(nVals,3)
    !> Level to which the evaluated values to be stored
    integer, intent(in) :: iLevel
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: evalVal(nVals)
    ! -------------------------------------------------------------------- !

    ! Do not store value for const or none
    if (trim(me%kind) /= 'const' .or. trim(me%kind) /= 'none') then

      me%isStored = .true.

      evalVal = tem_spatial_for( me    = me,    &
        &                        coord = coord, &
        &                        n     = nVals  )

      call append( me     = me%valOnLvl(iLevel)%evalVal, &
        &          val    = evalVal,                     &
        &          length = nVals                        )

      call truncate(me = me%valOnLvl(iLevel)%evalVal)

    else

      me%isStored = .false.

    end if

  end subroutine tem_spatial_scalar_storeVal
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This routine evaluate vector spatial function and store its value in
  !! growing array with access Array Of Structure pattern
  !! (iVal-1)*nComps + iComp
  subroutine tem_spatial_vector_storeVal( me, coord, nVals, iLevel, nComps )
    ! -------------------------------------------------------------------- !
    !> spatial type for given boundary state
    type(tem_spatial_type), intent(inout) :: me
    !> number of return values
    integer, intent(in) :: nVals
    !> number of components per returned value
    integer, intent(in) :: nComps
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent(in) :: coord(nVals,3)
    !> Level to which the evaluated values to be stored
    integer, intent(in) :: iLevel
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: evalVal(nVals, nComps)
    integer :: iVal
    ! -------------------------------------------------------------------- !

    ! Do not store value for const or none
    if (trim(me%kind) /= 'const' .or. trim(me%kind) /= 'none') then
      me%isStored = .true.

      evalVal = tem_spatial_for( me    = me,    &
        &                        coord = coord, &
        &                        n     = nVals, &
        &                        nComp = nComps )

      ! Append evaluated values in 1D growing array with AOS acces
      do iVal = 1, nVals
        call append( me     = me%valOnLvl(iLevel)%evalVal, &
          &          val    = evalVal(iVal,:),             &
          &          length = nVals * nComps               )
      end do

      call truncate(me = me%valOnLvl(iLevel)%evalVal)
    else
      me%isStored = .false.
    end if


  end subroutine tem_spatial_vector_storeVal
  ! ------------------------------------------------------------------------ !

end module tem_spatial_module
