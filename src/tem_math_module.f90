! Copyright (c) 2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012, 2014 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2013, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013, 2018-2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
!> Some generic matrix and vector function
!!
module tem_math_module

  ! include treelm modules
  use env_module,         only: rk
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit
  use tem_matrix_module,  only: tem_matrix_type, invert_matrix

  implicit none

  private

  public :: invert_matrix
  public :: tem_matrix_type
  public :: tem_intp_trilinearReduced
  public :: cross_product3D
  public :: inamuroDelta3D

  interface tem_intp_trilinearReduced
    module procedure tem_intp_trilinearReduced_scal
    module procedure tem_intp_trilinearReduced_vect
  end interface

  contains

! ****************************************************************************** !
  !> This function returns the tri-linearly interpolated values from the seven
  !! source points to the target position located at targetCoord.
  !! The source points are arranged in a square from (0,0,0)x(1,1,1)
  !! The order of the source points are according to the morton curve
  !! \[
  !! \begin{matrix}
  !!  1 & 2 & 3 & 4 \\
  !! (0,0,0); & (1,0,0); & (0,1,0); & (1,1,0) \\
  !!  5 & 6 & 7 \\
  !! (0,0,1); & (1,0,1); & (0,1,1) \\
  !! \end{matrix}
  !! \]
  !!
  function tem_intp_trilinearReduced_scal( srcVal, targetCoord ) result( phi )
    ! ---------------------------------------------------------------------------
    !> source values of the square corners
    real(kind=rk), intent(in) :: srcVal(7)
    !> interpolation location within the square
    real(kind=rk), intent(in) :: targetCoord(3)
    !> interpolated value
    real(kind=rk) :: phi
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: phi_northFront, phi_southFront
    real(kind=rk) :: phi_northBack, phi_southBack
    real(kind=rk) :: phi_front, phi_back
    ! ---------------------------------------------------------------------------
    ! Linear interpolation on the cube front side (z = 0 )
    phi_northFront = (1._rk - targetCoord(1))* srcVal(3)                       &
      &                     + targetCoord(1) * srcVal(4)
    phi_southFront = (1._rk - targetCoord(1))* srcVal(1)                       &
      &                     + targetCoord(1) * srcVal(2)
    ! Linear interpolation on the cube back side (z = 1 )
    phi_northBack  = srcVal(7) !(1._rk - targetCoord(1))* srcVal(7) ! + targetCoord(1)*srcVal(8)
    phi_southBack  = (1._rk - targetCoord(1))* srcVal(5)                       &
      &                     + targetCoord(1) * srcVal(6)
    ! Linear interpolation on the cube front side (z = 0 )
    phi_front = (1._rk - targetCoord(2))* phi_southFront                       &
      &                + targetCoord(2) * phi_northFront
    ! Linear interpolation on the cube back side (z = 1 )
    phi_back  = (1._rk - targetCoord(2))* phi_southBack                        &
      &                + targetCoord(2) * phi_northBack
    phi       = (1._rk - targetCoord(3))* phi_front                            &
      &                + targetCoord(3) * phi_back

  end function tem_intp_trilinearReduced_scal
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function returns the tri-linearly interpolated values from the seven
  !! source points to the target position located at targetCoord.
  !! The source points are arranged in a square from (0,0,0)x(1,1,1)
  !! The order of the source points are according to the morton curve
  !!
  !! \[
  !! \begin{matrix}
  !!  1 & 2 & 3 & 4 \\
  !! (0,0,0); & (1,0,0); & (0,1,0); & (1,1,0) \\
  !!  5 & 6 & 7 \\
  !! (0,0,1); & (1,0,1); & (0,1,1) \\
  !! \end{matrix}
  !! \]
  !!
  function tem_intp_trilinearReduced_vect( srcVal, targetCoord, nSize )        &
    &                                                             result( phi )
    ! ---------------------------------------------------------------------------
    !> vector size
    integer, intent(in) :: nSize
    !> source values of the square corners
    real(kind=rk), intent(in) :: srcVal(nSize, 7)
    !> interpolation location within the square
    real(kind=rk), intent(in) :: targetCoord(3)
    !> interpolated value
    real(kind=rk) :: phi( nSize )
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: phi_northFront, phi_southFront
    real(kind=rk) :: phi_northBack, phi_southBack
    real(kind=rk) :: phi_front, phi_back
    integer :: iSize
    ! ---------------------------------------------------------------------------
    do iSize = 1, nSize
      ! Linear interpolation on the cube front side (z = 0 )
      phi_northFront = (1._rk - targetCoord(1))* srcVal(iSize,3)               &
        &                     + targetCoord(1) * srcVal(iSize,4)
      phi_southFront = (1._rk - targetCoord(1))* srcVal(iSize,1)               &
        &                     + targetCoord(1) * srcVal(iSize,2)
      ! Linear interpolation on the cube back side (z = 1 )
      phi_northBack  = srcVal(iSize,7) !(1._rk - targetCoord(1))* srcVal(7) ! + targetCoord(1)*srcVal(8)
      phi_southBack  = (1._rk - targetCoord(1))* srcVal(iSize,5)               &
        &                     + targetCoord(1) * srcVal(iSize,6)
      ! Linear interpolation on the cube front side (z = 0 )
      phi_front = (1._rk - targetCoord(2))* phi_southFront                     &
        &                + targetCoord(2) * phi_northFront
      ! Linear interpolation on the cube back side (z = 1 )
      phi_back  = (1._rk - targetCoord(2))* phi_southBack                      &
        &                + targetCoord(2) * phi_northBack
      phi( iSize ) = (1._rk - targetCoord(3))* phi_front                       &
        &                   + targetCoord(3) * phi_back
    enddo

  end function tem_intp_trilinearReduced_vect
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function calculate the cross product of two 3D vector
  !!
  pure function cross_product3D(a, b) result( cross )
    ! ---------------------------------------------------------------------------
    !> resulting cross produkt
    real(kind=rk) :: cross(3)
    !> input vector a
    real(kind=rk), intent(in) :: a(3)
    !> input vector b
    real(kind=rk), intent(in) :: b(3)
    ! ---------------------------------------------------------------------------

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)

  end function cross_product3D
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function calculates the delta function used in the paper of Ota et al.
  !! [7] (bibliography of treelm) for a single value.
  !! \[
  !!     w(r) = \left\{ \begin{array}{ll} \frac{1}{8}\left ( 3-2|r|+\sqrt{1+4|r|
  !!                                           -4r^{2}} \right )& , |r| \le 1 \\
  !!                    \frac{1}{8}\left ( 5-2|r|+\sqrt{-7+12|r|-4r^{2}} \right
  !!                                                     )& , 1 \le |r| \le 2 \\
  !!                    0& , else \end{array} \right.
  !! \]
  !!
  function inamuroDelta1D(r) result( res )
    ! ---------------------------------------------------------------------------
    !> input point coordinate
    real(kind=rk), intent(in) :: r
    !> resulting value of the 1D delta function
    real(kind=rk) :: res
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: r_abs
    ! ---------------------------------------------------------------------------
    r_abs = abs(r)
    if( floor(r_abs) .eq. 0) then
      res = (3._rk-2._rk*r_abs+sqrt(1._rk+4._rk*r_abs-4._rk*r_abs**2))/8._rk
    else if ( floor(r_abs) .eq. 1 ) then
      res = (5._rk-2._rk*r_abs-sqrt(-7._rk+12._rk*r_abs-4._rk*r_abs**2))/8._rk
    else
      res = 0._rk
    end if

  end function inamuroDelta1D
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function calculates the delta function used in the paper of Ota et al.
  !! [7] (bibliography of treelm)
  !! for a vector by multiplying the results of the 1D version.
  !! \[
  !!     W(r) = w\left(\frac{r_1}{\Delta x}\right) \cdot
  !!            w\left(\frac{r_2}{\Delta x}\right) \cdot
  !!            w\left(\frac{r_3}{\Delta x}\right) \cdot
  !!            \frac{1}{\Delta x^3}
  !! \]
  !!
  function inamuroDelta3D(r, dx) result( res )
    ! ---------------------------------------------------------------------------
    !> input point coordinates
    real(kind=rk), intent(in) :: r(3)
    !> spatial discretization
    real(kind=rk), intent(in) :: dx
    !> resulting value of the 3D delta function
    real(kind=rk) :: res
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: tmpVal_X
    real(kind=rk) :: tmpVal_Y
    real(kind=rk) :: tmpVal_Z
    ! ---------------------------------------------------------------------------
    tmpVal_X = inamuroDelta1D(r(1)/dx)
    tmpVal_Y = inamuroDelta1D(r(2)/dx)
    tmpVal_Z = inamuroDelta1D(r(3)/dx)

    res = tmpVal_X * tmpVal_Y * tmpVal_Z / dx**3

  end function inamuroDelta3D
! ****************************************************************************** !


end module tem_math_module
! ****************************************************************************** !
