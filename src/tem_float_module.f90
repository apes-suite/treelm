! Copyright (c) 2016-2017, 2019-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> Floating point numbers need a special treatment due to their storage format.
!!
!! Some people already made up their minds about those things, thus that we can
!! make use of their findings.
!! See further:
!!   - http://www.floating-point-gui.de/errors/comparison/
!!   - https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
!!   - http://www.ssec.wisc.edu/~paulv/Fortran90/Utility/Compare_Float_Numbers.f90.html
!!
!! Most of this modules code traces back to the work of Paul van Delst,
!! CIMSS/SSEC, paul.vandelst@ssec.wisc.edu
module tem_float_module

  use env_module,       only: rk

  implicit none

  private

  public :: operator(.feq.)
  public :: operator(.fne.)
  public :: operator(.fgt.)
  public :: operator(.fge.)
  public :: operator(.flt.)
  public :: operator(.fle.)

  interface operator (.feq.)
    module procedure equal
    module procedure equal_array
  end interface
  interface operator (.fne.)
    module procedure notEqual
    module procedure notEqual_array
  end interface
  interface operator (.fgt.)
    module procedure greaterThan
  end interface
  interface operator (.fge.)
    module procedure greaterThanOrEqual
  end interface
  interface operator (.flt.)
    module procedure lessThan
  end interface
  interface operator (.fle.)
    module procedure lessThanOrEqual
  end interface

contains

  !> Relational operator to test the equality of floating point numbers.
  elemental function equal(a,b) result(res)
    !> Floating point value to be compared.
    real(kind=rk), intent(in) :: a
    !> Floating point value to be compared.
    real(kind=rk), intent(in) :: b
    !> The result is a logical value indicating whether the operands are equal
    !! to within numerical precision.
    logical :: res
    res = abs( a - b ) < spacing( max( abs( a ), abs( b ) ) )
  end function

  !> Relational operator to test the equality of two arrays of floating point
  !! numbers.
  pure function equal_array(a,b) result(res)
    !> Floating point array to be compared.
    real(kind=rk), intent(in) :: a(:)
    !> Floating point array to be compared.
    real(kind=rk), intent(in) :: b(:)
    !> The result is a logical value indicating whether the operands are equal
    !! to within numerical precision.
    logical :: res
    integer :: i
    res = size(a) == size(b)
    if (res) then
      do i = 1, size(a)
        res = res .and. equal(a(i), b(i))
      end do
    end if
  end function

  !> Relational operator to test the not-equality of floating point numbers.
  elemental function notEqual(a,b) result(res)
    !> Floating point value to be compared.
    real(kind=rk), intent(in) :: a
    !> Floating point value to be compared.
    real(kind=rk), intent(in) :: b
    !> The result is a logical value indicating whether the operands are not
    !! equal to within numerical precision.
    logical :: res
    res = .not. equal(a,b)
  end function

  !> Relational operator to test the not-equality of two floating point arrays.
  pure function notEqual_array(a,b) result(res)
    !> Floating point array to be compared.
    real(kind=rk), intent(in) :: a(:)
    !> Floating point array to be compared.
    real(kind=rk), intent(in) :: b(:)
    !> The result is a logical value indicating whether the operands are not
    !! equal to within numerical precision.
    logical :: res
    res = .not. equal_array(a,b)
  end function

  !> Relational operator to test if one operand is greater than another.
  elemental function greaterThan(a,b) result(res)
    !> Floating point value to be compared.
    real(kind=rk), intent(in) :: a
    !> Floating point value to be compared.
    real(kind=rk), intent(in) :: b
    !> The result is a logical value indicating whether the operand a is greater
    !! than b by more than the spacing between representable floating point
    !! numbers.
    logical :: res
    res = ( a - b ) >= spacing( max( abs( a ), abs( b ) ) )
  end function

  !> Relational operator to test if one operand is greater than or equal to
  !! another.
  elemental function greaterThanOrEqual(a,b) result(res)
    !> Floating point value to be compared.
    real(kind=rk), intent(in) :: a
    !> Floating point value to be compared.
    real(kind=rk), intent(in) :: b
    !> The result is a logical value indicating whether the operand a is greater
    !! than or equal to b by more than the spacing between representable
    !! floating point numbers.
    logical :: res
    res = greaterThan(a,b) .or. equal(a,b)
  end function

  !> Relational operator to test if one operand is less than another.
  elemental function lessThan(a,b) result(res)
    !> Floating point value to be compared.
    real(kind=rk), intent(in) :: a
    !> Floating point value to be compared.
    real(kind=rk), intent(in) :: b
    !> The result is a logical value indicating whether the operand a is less
    !! than b by more than the spacing between representable floating point
    !! numbers.
    logical :: res
    res = ( b - a ) >= spacing( max( abs( a ), abs( b ) ) )
  end function

  !> Relational operator to test if one operand is less than or equal to
  !! another.
  elemental function lessThanOrEqual(a,b) result(res)
    !> Floating point value to be compared.
    real(kind=rk), intent(in) :: a
    !> Floating point value to be compared.
    real(kind=rk), intent(in) :: b
    !> The result is a logical value indicating whether the operand a is less
    !! than or equal to b by more than the spacing between representable
    !! floating point numbers.
    logical :: res
    res = lessThan(a,b) .or. equal(a,b)
  end function

end module
