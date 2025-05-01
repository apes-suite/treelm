! Copyright (c) 2013, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
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
  use tem_logging_module, only: logUnit
  use tem_aux_module,     only: tem_abort
  use tem_param_module,   only: PI

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

    ! The number of expansion coefficients
    integer :: nCoeffs

    ! Expansion coefficients inside the cylinder
    complex(kind=rk), allocatable :: c_tot(:)

    ! Expansion coefficients outside the cylinder (for the scattered field)
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

  public :: tem_miescatter_field_type, &
          & tem_load_miescatter_displacementfieldz, &
          & tem_load_miescatter_magneticfieldx, &
          & tem_load_miescatter_magneticfieldy, &
          & tem_eval_miescatter_displz, &
          & tem_eval_miescatter_magnx, &
          & tem_eval_miescatter_magny

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

    write(logUnit(0),*)'ERROR: the Mie-Scatter function is not available in '//&
      &                'this'
    write(logUnit(0),*)'       executable, due to the compiler lacking '//     &
      &                'support of'
    write(logUnit(0),*)'       Bessel functions!'
    call tem_abort()

  end subroutine tem_load_miescatter
! ****************************************************************************** !


! ****************************************************************************** !
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

    res = 0.0_rk

  end function tem_eval_miescatter_displz
! ****************************************************************************** !


! ****************************************************************************** !
  !> Evaluate magnetic field (x-component) for Mie-Scattering of electromagnetic
  !! wave at dielectric cylinder.
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

    res = 0.0_rk

  end function tem_eval_miescatter_magnx
! ****************************************************************************** !


! ****************************************************************************** !
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

    res = 0.0_rk

  end function tem_eval_miescatter_magny
! ****************************************************************************** !


! ****************************************************************************** !
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

    res = 0.0_rk

  end function tem_eval_miescatter_magnangular
! ****************************************************************************** !


! ****************************************************************************** !
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

    res = 0.0_rk

  end function tem_eval_miescatter_magnradial
! ****************************************************************************** !

end module tem_miescatter_module
! ****************************************************************************** !

