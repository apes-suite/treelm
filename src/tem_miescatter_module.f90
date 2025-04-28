! Copyright (c) 2013 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
! **************************************************************************** !
!> This module gathers the definition of the analytical solution by
!! scattering of an electromagnetic wave at a dielectric sphere/cylinder. The
!! solution is given by means of Mie-series.
!! The solution is given in:
!! Cai, W., & Deng, S. (2003). An upwinding embedded boundary method
!! for Maxwell's equations in media with material interfaces: 2D case.
!! Journal of Computational Physics, 190(1), 159-183.
!! doi:10.1016/S0021-9991(03)00269-9
!!
module tem_miescatter_module

  ! include treelm modules
  use env_module,         only: rk
  use tem_aux_module,     only: tem_abort
  use tem_param_module,   only: PI
  use tem_logging_module, only: logUnit

  ! include aotus modules
  use flu_binding,      only: flu_state
  use aotus_module,     only: aoterr_NonExistent
  use aot_table_module, only: aot_table_open, aot_table_close,                 &
   &                          aot_table_length, aot_get_val

  implicit none

  private

  !> Parameters of the solution for Mie-Scatter at
  !! dielectric cylinder (infinite height in z direction)
  type tem_miescatter_type
    !> The center of the cylinder (in the x-y plane)
    real(kind=rk) :: center(2) = [ 0.0_rk, 0.0_rk ]

    !> The radius of the cylinder
    real(kind=rk) :: radius = 1.0_rk

    !> Permeability (mu) of the background.
    real(kind=rk) :: permeability_background = 1.0_rk

    !> Permitivity (epsilon) of the background.
    real(kind=rk) :: permitivity_background = 1.0_rk

    !> Permeability (mu) of the cylinder.
    real(kind=rk) :: permeability_cylinder = 1.0_rk

    !> Permitivity (epsilon) of the cylinder.
    real(kind=rk) :: permitivity_cylinder = 1.0_rk

    !> Wave number for the background
    real(kind=rk) :: wavenumber_background

    !> Wave number for the cylinder
    real(kind=rk) :: wavenumber_cylinder

    !> Angular frequency ( = 2 * pi * f , where
    !! f is the original frequency in Hertz)
    real(kind=rk) :: omega = 2.0_rk * PI * 1.0_rk

  end type tem_miescatter_type


  !> Expansion data for the Mier scatter solution.
  type tem_mieexpansion_type

    !> The number of expansion coefficients
    integer :: nCoeffs

    !> Expansion coefficients inside the cylinder
    complex(kind=rk), allocatable :: c_tot(:)

    !> Expansion coefficients outside the cylinder (for the scattered field)
    complex(kind=rk), allocatable :: c_scat(:)

  end type tem_mieexpansion_type


  !> Parameters of the solution for Mie-Scatter at
  !! dielectric cylinder (infinite height in z direction).
  type tem_miescatter_field_type

    !> Parameter of the geometrical and material setup
    type(tem_miescatter_type) :: miescatter

    !> The expansion parameters for the scattered solution
    type(tem_mieexpansion_type) :: mieexpansion

  end type tem_miescatter_field_type

  interface tem_load_miescatter_displacementfieldz
    module procedure tem_load_miescatter
  end interface tem_load_miescatter_displacementfieldz

  interface tem_load_miescatter_magneticfieldx
    module procedure tem_load_miescatter
  end interface tem_load_miescatter_magneticfieldx

  interface tem_load_miescatter_magneticfieldy
    module procedure tem_load_miescatter
  end interface tem_load_miescatter_magneticfieldy

  public :: tem_miescatter_field_type,              &
    &       tem_load_miescatter_displacementfieldz, &
    &       tem_load_miescatter_magneticfieldx,     &
    &       tem_load_miescatter_magneticfieldy,     &
    &       tem_eval_miescatter_displz,             &
    &       tem_eval_miescatter_magnx,              &
    &       tem_eval_miescatter_magny

contains

! ****************************************************************************** !
  !> load gauss pulse variables to set initial condition
  !!
  subroutine tem_load_miescatter(conf, thandle, me)
    ! ---------------------------------------------------------------------------
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: thandle
    !> Global gauss pulse data type
    type(tem_miescatter_field_type), intent(out) :: me
    ! ---------------------------------------------------------------------------
    integer :: iError, iPos, cent_handle, back_handle, circ_handle
    real(kind=rk) :: frequency
    ! ---------------------------------------------------------------------------

    ! Read center cooridnates
    call aot_table_open( L       = conf,                                       &
      &                  parent  = thandle,                                    &
      &                  thandle = cent_handle,                                &
      &                  key     = 'center' )
    if ( cent_handle.eq.0 ) then
      write(logUnit(1),*)'ERROR in tem_load_miescatter_displacementfieldz: '
      write(logUnit(1),*)'Not able to read center-table for Mie-Scatter '//    &
        &                'solution from lua, stopping ...'
      call tem_abort()
    end if
    do iPos = 1, 2
      call aot_get_val( L       = conf,                                        &
        &               thandle = cent_handle,                                 &
        &               pos     = iPos,                                        &
        &               val     = me%miescatter%center(iPos),                  &
        &               ErrCode = iError )
      if ( iError.ne.0 ) then
        write(logUnit(1),*)'ERROR in tem_load_miescatter_displacementfieldz: '
        write(logUnit(1),*)'Not able to read center for Mie-Scatter '//        &
          &                'solution from lua, stopping ...'
        call tem_abort()
      end if
    end do
    call aot_table_close(conf, cent_handle)


    ! Read radius of the cylinder
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'radius',                                      &
      &               val     = me%miescatter%radius,                          &
      &               ErrCode = iError )
    if ( iError.ne.0 ) then
      write(logUnit(1),*)'ERROR in tem_load_miescatter_displacementfieldz: '
      write(logUnit(1),*)'Not able to read radius for Mie-Scatter '//          &
         &                'solution from lua, stopping ...'
      call tem_abort()
    end if


    ! Read background permeability and permitivity
    call aot_table_open( L       = conf,                                       &
      &                  parent  = thandle,                                    &
      &                  thandle = back_handle,                                &
      &                  key     = 'permeaPermit_background' )
    if ( back_handle.eq.0 ) then
      write(logUnit(1),*)'ERROR in tem_load_miescatter_displacementfieldz: '
      write(logUnit(1),*) 'Not able to open permeaPermit_background-table for '&
                    & // 'Mie-Scatter solution from lua, stopping ...'
      call tem_abort()
    end if
    call aot_get_val( L       = conf,                                          &
      &               thandle = back_handle,                                   &
      &               pos     = 1,                                             &
      &               val     = me%miescatter%permeability_background,         &
      &               ErrCode = iError )
    if ( iError.ne.0 ) then
      write(logUnit(1),*)'ERROR in tem_load_miescatter_displacementfieldz: '
      write(logUnit(1),*) 'Not able to read background permeability for ' //   &
                    & 'Mie-Scatter solution from lua, stopping ...'
      call tem_abort()
    end if
    call aot_get_val( L       = conf,                                          &
      &               thandle = back_handle,                                   &
      &               pos     = 2,                                             &
      &               val     = me%miescatter%permitivity_background,          &
      &               ErrCode = iError )
    if ( iError.ne.0 ) then
      write(logUnit(1),*)'ERROR in tem_load_miescatter_displacementfieldz: '
      write(logUnit(1),*) 'Not able to read background permitivity for ' //    &
                    & 'Mie-Scatter solution from lua, stopping ...'
      call tem_abort()
    end if
    call aot_table_close(conf, back_handle)


    ! Read background permeability and permitivity
    call aot_table_open( L       = conf,                                       &
      &                  parent  = thandle,                                    &
      &                  thandle = circ_handle,                                &
      &                  key     = 'permeaPermit_cylinder' )
    if ( back_handle.eq.0 ) then
      write(logUnit(1),*)'ERROR in tem_load_miescatter_displacementfieldz: '
      write(logUnit(1),*) 'Not able to open permeaPermit_cylinder-table for '//&
                    & 'Mie-Scatter solution from lua, stopping ...'
      call tem_abort()
    end if
    call aot_get_val( L       = conf,                                          &
      &               thandle = circ_handle,                                   &
      &               pos     = 1,                                             &
      &               val     = me%miescatter%permeability_cylinder,           &
      &               ErrCode = iError )
    if ( iError.ne.0 ) then
      write(logUnit(1),*)'ERROR in tem_load_miescatter_displacementfieldz: '
      write(logUnit(1),*) 'Not able to read cylinder permeability for ' // &
                    & 'Mie-Scatter solution from lua, stopping ...'
      call tem_abort()
    end if
    call aot_get_val( L       = conf,                                          &
      &               thandle = circ_handle,                                   &
      &               pos     = 2,                                             &
      &               val     = me%miescatter%permitivity_cylinder,            &
      &               ErrCode = iError )
    if ( iError.ne.0 ) then
      write(logUnit(1),*)'ERROR in tem_load_miescatter_displacementfieldz: '
      write(logUnit(1),*) 'Not able to read cylinder permitivity for ' // &
                    & 'Mie-Scatter solution from lua, stopping ...'
      call tem_abort()
    end if
    call aot_table_close(conf, circ_handle)


    ! Read the frequency
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'frequency',                                   &
      &               val     = frequency,                                     &
      &               ErrCode = iError )
    if ( iError.ne.0 ) then
      write(logUnit(1),*)'ERROR in tem_load_miescatter_displacementfieldz: '
      write(logUnit(1),*)'Not able to read frequency for Mie-Scatter '//       &
        &                'solution from lua, stopping ...'
      call tem_abort()
    end if
    me%miescatter%omega = 2.0_rk * PI * frequency

    ! Calculate the wave numbers inside and outside of the cylinder
    me%miescatter%wavenumber_cylinder = me%miescatter%omega &
    & * sqrt(me%miescatter%permeability_cylinder *          &
    &        me%miescatter%permitivity_cylinder             )
    me%miescatter%wavenumber_background = me%miescatter%omega &
    & * sqrt(me%miescatter%permeability_background *          &
    &        me%miescatter%permitivity_background             )

    ! Setup the parameters for the solution inside and outside the
    ! dielectric cylinder...
    call tem_init_data( conf, thandle, me%miescatter, me%mieexpansion )

  end subroutine tem_load_miescatter

  !> Init the expansion coefficients for the Mie-Scattering.
  subroutine tem_init_data( conf, thandle, miescatter, mieexpansion )
    ! --------------------------------------------------------------------------
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: thandle
    !> The description of the Mie-Scattering setup.
    type(tem_miescatter_type), intent(in) :: miescatter
    !> The expansion coefficients to be filled
    type(tem_mieexpansion_type), intent(out) :: mieexpansion
    ! --------------------------------------------------------------------------
    integer :: nCoeffs, n
    ! --------------------------------------------------------------------------

    ! Get the number of expansion coefficients to be used ....
    nCoeffs = tem_get_ncoeffs_miescat(conf,thandle)
    mieexpansion%nCoeffs = nCoeffs

    ! Calculate the expansion coefficients inside the cylinder
    allocate(mieexpansion%c_tot(-nCoeffs:nCoeffs))
    innerLoop: do n = -nCoeffs, nCoeffs

      ! C_n^{tot}
      mieexpansion%c_tot(n)                                                    &
        & = (0.0_rk,1.0_rk)**(-n) *                                            &
        & ( (miescatter%wavenumber_background /                                &
        &               miescatter%permeability_background)                    &
        & * bessel_jn_derivative(n,                                            &
        &               miescatter%wavenumber_background * miescatter%radius)  &
        & * hankel2_n(n,miescatter%wavenumber_background * miescatter%radius)  &
        & - (miescatter%wavenumber_background /                                &
        &               miescatter%permeability_background)                    &
        & * hankel2_n_derivative(n,                                            &
        &               miescatter%wavenumber_background * miescatter%radius)  &
        & * bessel_jn(n,miescatter%wavenumber_background * miescatter%radius)  &
        & )                                                                    &
        & /                                                                    &
        & (                                                                    &
        &   (miescatter%wavenumber_cylinder /                                  &
        &               miescatter%permeability_cylinder)                      &
        & * bessel_jn_derivative(n,                                            &
        &               miescatter%wavenumber_cylinder * miescatter%radius)    &
        & * hankel2_n(n,miescatter%wavenumber_background * miescatter%radius)  &
        & - (miescatter%wavenumber_background /                                &
        &               miescatter%permeability_background)&
        & * hankel2_n_derivative(n,                                            &
        &               miescatter%wavenumber_background * miescatter%radius)  &
        & * bessel_jn(n,miescatter%wavenumber_cylinder * miescatter%radius)    &
        & )

    end do innerLoop

    ! Calculate the expansion coefficients outside of the cylinder, for the
    ! scattered field
    allocate(mieexpansion%c_scat(-nCoeffs:nCoeffs))
    outerLoop: do n = -nCoeffs, nCoeffs

      ! C_n^{scat}
      mieexpansion%c_scat(n) &
      & = (0.0_rk,1.0_rk)**(-n) *                                              &
      & ( (miescatter%wavenumber_background /                                  &
      &               miescatter%permeability_background)                      &
      & * bessel_jn_derivative(n,                                              &
      &               miescatter%wavenumber_background*miescatter%radius)      &
      & * bessel_jn(n,miescatter%wavenumber_cylinder*miescatter%radius)        &
      & - (miescatter%wavenumber_cylinder/miescatter%permeability_cylinder)    &
      & * bessel_jn_derivative(n,                                              &
      &               miescatter%wavenumber_cylinder*miescatter%radius)        &
      & * bessel_jn(n,miescatter%wavenumber_background*miescatter%radius)      &
      & )                                                                      &
      & /                                                                      &
      & (                                                                      &
      &   (miescatter%wavenumber_cylinder/miescatter%permeability_cylinder)    &
      & * bessel_jn_derivative(n,                                              &
      &               miescatter%wavenumber_cylinder*miescatter%radius)        &
      & * hankel2_n(n,miescatter%wavenumber_background*miescatter%radius)      &
      & - (miescatter%wavenumber_background /                                  &
      &               miescatter%permeability_background)                      &
      & * hankel2_n_derivative(n,                                              &
      &               miescatter%wavenumber_background*miescatter%radius)      &
      & * bessel_jn(n,miescatter%wavenumber_cylinder*miescatter%radius)        &
      & )

    end do outerLoop

  end subroutine tem_init_data


  !> Evaluate displacement field (z component) for Mie-Scattering of
  !! electromagnetic wave at dielectric cylinder.
  function tem_eval_miescatter_displz(me, coord, time, n) result(res)
    ! ---------------------------------------------------------------------------
    !> The function to evaluate
    type(tem_miescatter_field_type), intent(in) :: me
    !> Number of points to evaluate the function for.
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> The time to evaluate the function at.
    real(kind=rk), intent( in )  :: time
    !> return value of the function
    real(kind=rk) :: res(n)
    ! ---------------------------------------------------------------------------
    integer :: iPoint, iCoeff
    ! Polar coordinate vector in x-y plane.
    ! First entry: radius \n
    ! Second entry: angle
    real(kind=rk) :: polar(2)
    complex(kind=rk) :: tmp
    ! ---------------------------------------------------------------------------

    do iPoint = 1, n

      ! Convert to polar coordinates (relative to the center of the
      ! cylinder.
      polar = cart2polar( coord(iPoint,1)-me%miescatter%center(1), &
                        & coord(iPoint,2)-me%miescatter%center(2)  )

      ! Inside the cylinder
      if(polar(1).le.me%miescatter%radius) then

        tmp = (0.0_rk, 0.0_rk)
        do iCoeff = -me%mieexpansion%nCoeffs , me%mieexpansion%nCoeffs
          tmp = tmp                                                           &
            & + me%mieexpansion%c_tot(iCoeff)                                 &
            & * bessel_jn(iCoeff, me%miescatter%wavenumber_cylinder*polar(1)) &
            & * exp( (0.0_rk,1.0_rk) * iCoeff * polar(2)  )
        end do
        ! Twiddle by phase and get the real part
        res(iPoint) = real( tmp * exp( (0.0_rk,1.0_rk) * me%miescatter%omega  &
          &                * time ) )

        ! Multiply with permitivity to get the displacement field
        res(iPoint) = me%miescatter%permitivity_cylinder * res(iPoint)

      ! Outside the cylinder
      else

        ! Add up incoming and scattered field
        tmp = (0.0_rk, 0.0_rk)
        do iCoeff = -me%mieexpansion%nCoeffs , me%mieexpansion%nCoeffs
          tmp = tmp +                                                     &
            & ( (0.0_rk,1.0_rk)**(-iCoeff)                                &
            &   * bessel_jn(iCoeff,                                       &
            &               me%miescatter%wavenumber_background*polar(1)) &
            &   + me%mieexpansion%c_scat(iCoeff)                          &
            &   * hankel2_n(iCoeff,                                       &
            &               me%miescatter%wavenumber_background*polar(1)) &
            & ) * exp( (0.0_rk,1.0_rk) * iCoeff * polar(2)  )
        end do
        ! Twiddle by phase and get the real part
        res(iPoint) = real(tmp * exp( (0.0_rk,1.0_rk) * &
          &                me%miescatter%omega * time ) )

        ! Multiply with permitivity to get the displacement field
        res(iPoint) = me%miescatter%permitivity_background * res(iPoint)

      end if

    end do

  end function tem_eval_miescatter_displz

  !> Evaluate magnetic field (x-component) for Mie-Scattering of
  !! electromagnetic wave at dielectric cylinder.
  function tem_eval_miescatter_magnx(me, coord, time, n) result(res)
    ! ---------------------------------------------------------------------------
    !> The function to evaluate
    type(tem_miescatter_field_type), intent(in) :: me
    !> Number of points to evaluate the function for.
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> The time to evaluate the function at.
    real(kind=rk), intent( in )  :: time
    !> return value of the function
    real(kind=rk) :: res(n)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: H_r(n), H_theta(n), H_x, polar(2)
    integer :: iPoint
    ! ---------------------------------------------------------------------------

    ! Calculate magnetizing field in polar coordinate system
    ! ... radial component
    H_r = tem_eval_miescatter_magnradial(me, coord, time, n)
    !... angular component
    H_theta = tem_eval_miescatter_magnangular(me, coord, time, n)

    ! Convert from polar vector field to cartesian vector field
    ! and convert from magnetizing field (i.e. H) to magnetic field (i.e. B)
    do iPoint = 1, n

      ! Convert to polar coordinates (relative to the center of the
      ! cylinder.
      polar = cart2polar( coord(iPoint,1)-me%miescatter%center(1), &
                        & coord(iPoint,2)-me%miescatter%center(2)  )

      ! Vector field in cartesian coordinates
      H_x = H_r(iPoint) * cos(polar(2)) - H_theta(iPoint) * sin(polar(2))

      ! Inside the cylinder
      if(polar(1) .le. me%miescatter%radius) then
        res(iPoint) = me%miescatter%permeability_cylinder * H_x
      ! Outside the cylinder
      else
        res(iPoint) = me%miescatter%permeability_background * H_x
      end if

    end do

  end function tem_eval_miescatter_magnx


  !> Evaluate magnetic field (y-component) for Mie-Scattering of electromagnetic
  !! wave at dielectric cylinder.
  function tem_eval_miescatter_magny(me, coord, time, n) result(res)
    ! ---------------------------------------------------------------------------
    !> The function to evaluate
    type(tem_miescatter_field_type), intent(in) :: me
    !> Number of points to evaluate the function for.
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> The time to evaluate the function at.
    real(kind=rk), intent( in )  :: time
    !> return value of the function
    real(kind=rk) :: res(n)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: H_r(n), H_theta(n), H_y, polar(2)
    integer :: iPoint
    ! ---------------------------------------------------------------------------

    ! Calculate magnetizing field in polar coordinate system
    ! ... radial component
    H_r = tem_eval_miescatter_magnradial(me, coord, time, n)
    !... angular component
    H_theta = tem_eval_miescatter_magnangular(me, coord, time, n)

    ! Convert from polar vector field to cartesian vector field
    ! and convert from magnetizing field (i.e. H) to magnetic field (i.e. B)
    do iPoint = 1, n

      ! Convert to polar coordinates (relative to the center of the
      ! cylinder.
      polar = cart2polar( coord(iPoint,1)-me%miescatter%center(1), &
                        & coord(iPoint,2)-me%miescatter%center(2)  )

      ! Vector field in cartesian coordinates
      H_y = H_r(iPoint) * sin(polar(2)) + H_theta(iPoint) * cos(polar(2))

      ! Inside the cylinder
      if(polar(1) .le. me%miescatter%radius) then
        res(iPoint) = me%miescatter%permeability_cylinder * H_y
      ! Outside the cylinder
      else
        res(iPoint) = me%miescatter%permeability_background * H_y
      end if

    end do

  end function tem_eval_miescatter_magny


  !> Evaluate magnetizing field (angular-component) for Mie-Scattering of
  !! electromagnetic wave at dielectric cylinder.
  function tem_eval_miescatter_magnangular(me, coord, time, n) result(res)
    ! ---------------------------------------------------------------------------
    !> The function to evaluate
    type(tem_miescatter_field_type), intent(in) :: me
    !> Number of points to evaluate the function for.
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> The time to evaluate the function at.
    real(kind=rk), intent( in )  :: time
    !> return value of the function
    real(kind=rk) :: res(n)
    ! ---------------------------------------------------------------------------
    integer :: iPoint, iCoeff
    ! Polar coordinate vector in x-y plane.
    ! First entry: radius \n
    ! Second entry: angle
    real(kind=rk) :: polar(2)
    complex(kind=rk) :: tmp
    ! ---------------------------------------------------------------------------

    do iPoint = 1, n

      ! Convert to polar coordinates (relative to the center of the
      ! cylinder.
      polar = cart2polar( coord(iPoint,1)-me%miescatter%center(1), &
        &                 coord(iPoint,2)-me%miescatter%center(2)  )

      ! Inside the cylinder
      if(polar(1).le.me%miescatter%radius) then

        tmp = (0.0_rk, 0.0_rk)
        do iCoeff = -me%mieexpansion%nCoeffs , me%mieexpansion%nCoeffs
          tmp = tmp                                                   &
            & + me%mieexpansion%c_tot(iCoeff)                         &
            & * bessel_jn_derivative(iCoeff,                          &
            &             me%miescatter%wavenumber_cylinder*polar(1)) &
            & * exp( (0.0_rk,1.0_rk) * iCoeff * polar(2) )
        end do
        ! Twiddle by phase and get the real part
        res(iPoint) = real( tmp                                         &
          ! JZ: the paper has an additional factor of -1 here. However, I think
          ! that this factor is wrong at this place. So, I commented out
          ! the factor.
          ! & * (-1.0_rk) &
          & * (exp( (0.0_rk,1.0_rk)*me%miescatter%omega*time))          &
          & * (0.0_rk,-1.0_rk)*me%miescatter%wavenumber_cylinder        &
          & / (me%miescatter%omega*me%miescatter%permeability_cylinder) )

      ! Outside the cylinder
      else

        ! Add up incoming and scattered field
        tmp = (0.0_rk, 0.0_rk)
        do iCoeff = -me%mieexpansion%nCoeffs , me%mieexpansion%nCoeffs
          tmp = tmp +                                                     &
            & ( (0.0_rk,1.0_rk)**(-iCoeff)                                &
            &   * bessel_jn_derivative(iCoeff,                            &
            &               me%miescatter%wavenumber_background*polar(1)) &
            &   + me%mieexpansion%c_scat(iCoeff)                          &
            &   * hankel2_n_derivative(iCoeff,                            &
            &               me%miescatter%wavenumber_background*polar(1)) &
            & ) * exp( (0.0_rk,1.0_rk) * iCoeff * polar(2)  )
        end do
        ! Twiddle by phase and get the real part
        res(iPoint) = real( tmp                                            &
           ! JZ: the paper has an additional factor of -1 here. However, I think
           ! that this factor is wrong at this place. So, I commented out
           ! the factor.
           !      & * (-1.0_rk) &
           & * (exp( (0.0_rk,1.0_rk) * me%miescatter%omega * time ))       &
           & * (0.0_rk,-1.0_rk)*me%miescatter%wavenumber_background        &
           & / (me%miescatter%omega*me%miescatter%permeability_background) )

      end if

    end do

  end function tem_eval_miescatter_magnangular


  !> Evaluate magnetizing field (radial-component) for Mie-Scattering of
  !! electromagnetic wave at dielectric cylinder.
  function tem_eval_miescatter_magnradial(me, coord, time, n) result(res)
    ! ---------------------------------------------------------------------------
    !> The function to evaluate
    type(tem_miescatter_field_type), intent(in) :: me
    !> Number of points to evaluate the function for.
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> The time to evaluate the function at.
    real(kind=rk), intent( in )  :: time
    !> return value of the function
    real(kind=rk) :: res(n)
    ! ---------------------------------------------------------------------------
    integer :: iPoint, iCoeff
    ! Polar coordinate vector in x-y plane.
    ! First entry: radius \n
    ! Second entry: angle
    real(kind=rk) :: polar(2)
    complex(kind=rk) :: tmp
    ! ---------------------------------------------------------------------------

    do iPoint = 1, n

      ! Convert to polar coordinates (relative to the center of the
      ! cylinder.
      polar = cart2polar( coord(iPoint,1)-me%miescatter%center(1), &
        &                 coord(iPoint,2)-me%miescatter%center(2)  )

      ! Inside the cylinder
      if(polar(1).le.me%miescatter%radius) then

        tmp = (0.0_rk, 0.0_rk)
        do iCoeff = -me%mieexpansion%nCoeffs , me%mieexpansion%nCoeffs
          tmp = tmp                                                           &
            & + me%mieexpansion%c_tot(iCoeff)                                 &
            & * bessel_jn(iCoeff, me%miescatter%wavenumber_cylinder*polar(1)) &
            & * (0.0_rk,1.0_rk) * iCoeff                                      &
            & * exp( (0.0_rk,1.0_rk) * iCoeff * polar(2) )
        end do
        ! Twiddle by phase and get the real part, first we check if distance from
        ! center is above a certain threshold (avoid division by zero).
        if(abs(polar(1)) > epsilon(1.0_rk) ) then
          res(iPoint) = real( tmp                                        &
            ! JZ: the paper has an additional factor of -1 here. However, I think
            ! that this factor is wrong at this place. So, I commented out
            ! the factor.
            !     & * (-1.0_rk) &
            & * (exp( (0.0_rk,1.0_rk)*me%miescatter%omega*time))         &
            & * (0.0_rk,1.0_rk)                                          &
            & / (me%miescatter%omega*me%miescatter%permeability_cylinder &
            &    *polar(1))                                              &
            & )
        else
          res(iPoint) = 0.0_rk
        end if

      ! Outside the cylinder
      else

        ! Add up incoming and scattered field
        tmp = (0.0_rk, 0.0_rk)
        do iCoeff = -me%mieexpansion%nCoeffs , me%mieexpansion%nCoeffs
          tmp = tmp +                                                     &
            & ( (0.0_rk,1.0_rk)**(-iCoeff)                                &
            &   * bessel_jn(iCoeff,                                       &
            &               me%miescatter%wavenumber_background*polar(1)) &
            &   + me%mieexpansion%c_scat(iCoeff)                          &
            &   * hankel2_n(iCoeff,                                       &
            &               me%miescatter%wavenumber_background*polar(1)) &
            & )                                                           &
            & * (0.0_rk,1.0_rk) * iCoeff                                  &
            & * exp( (0.0_rk,1.0_rk) * iCoeff * polar(2) )
        end do
        ! Twiddle by phase and get the real part, first we check if distance from
        ! center is above a certain threshold (avoid division by zero).
        if(abs(polar(1)) > epsilon(1.0_rk) ) then
          res(iPoint) = real( tmp                                          &
            ! JZ: the paper has an additional factor of -1 here. However, I think
            ! that this factor is wrong at this place. So, I commented out
            ! the factor.
            !     & * (-1.0_rk) &
            & * (exp( (0.0_rk,1.0_rk) * me%miescatter%omega * time ))      &
            & * (0.0_rk,1.0_rk)                                            &
            & / (me%miescatter%omega*me%miescatter%permeability_background &
            &    *polar(1))                                                &
            & )
        else
          res(iPoint) = 0.0_rk
        end if

      end if

    end do

  end function tem_eval_miescatter_magnradial

  !> Convert from cartesian coordinates (in the x-y plane) to
  !! polar coordinates (radius,angle)
  function cart2polar(x,y) result(polar)
    ! --------------------------------------------------------------------------
    !> X coordinate
    real(kind=rk) :: x
    !> Y coordinate
    real(kind=rk) :: y
    !> Polar coordinates, radius (first entry) and angle (second entry)
    real(kind=rk) :: polar(2)
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    polar(1) = sqrt(x*x + y*y)

    ! Atan2 is not defined when both coordinates are zero. To cover this
    ! situation correctly, we define the angle to be 0.
    if(polar(1) > epsilon(1.0_rk) ) then
      polar(2) = atan2(y,x)
    else
      polar(2) = 0.0_rk
    end if

  end function cart2polar


  !> Compute derivative for Bessel function of first kind of order n.
  function bessel_jn_derivative(n,x) result(der)
    ! --------------------------------------------------------------------------
    !> The order of the function.
    integer :: n
    !> The evaluation point.
    real(kind=rk), intent(in) :: x
    !> The function value.
    real(kind=rk) :: der
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    ! Calculate the derivative by recursive relation
    der = 0.5_rk*(bessel_jn(n-1,x)-bessel_jn(n+1,x))

  end function bessel_jn_derivative



  !> Compute derivative for Hankel function of second kind of order n.
  function hankel2_n_derivative(n,x) result(der)
    ! --------------------------------------------------------------------------
    !> The order of the function.
    integer :: n
    !> The evaluation point.
    real(kind=rk), intent(in) :: x
    !> The function value.
    complex(kind=rk) :: der
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    ! Calculate the derivative by recursive relation
    der = 0.5_rk*( hankel2_n(n-1,x)-hankel2_n(n+1,x) )

  end function hankel2_n_derivative



  !> Compute Hankel function of second kind of order n.
  function hankel2_n(n,x) result(val)
    ! --------------------------------------------------------------------------
    !> The order of the function.
    integer :: n
    !> The evaluation point.
    real(kind=rk), intent(in) :: x
    !> The function value.
    complex(kind=rk) :: val
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    val = bessel_jn(n, x) - (0.0_rk, 1.0_rk)*bessel_yn(n, x)

  end function hankel2_n


  function tem_get_ncoeffs_miescat(conf, thandle ) result(nCoeffs)
    ! --------------------------------------------------------------------------
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: thandle
    !> The number of coefficients necessary for the Fourier series
    integer :: nCoeffs
    ! --------------------------------------------------------------------------
    integer :: iError
    ! --------------------------------------------------------------------------

    !> @todo JZ: Currently, we read the number of necessary coefficients
    !! from the lua script. However, there are also some automatic ways
    !! to estimate the number of coefficients necessary to obtain a certain
    !! numerical accuaracy, e.g. in [Wiscombe1980] .

    ! Read the frequency
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'nCoeffs',                                     &
      &               val     = nCoeffs,                                       &
      &               ErrCode = iError )
    if ( iError.ne.0 ) then
      call tem_abort( 'ERROR in tem_get_ncoeffs_miescat: Not able to read' &
        & // ' nCoeffs for Mie-Scatter solution from lua'                  )
    end if

  end function tem_get_ncoeffs_miescat

end module tem_miescatter_module
! ****************************************************************************** !

