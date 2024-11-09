! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014 Rishabh Chandola <rishabh.chandola@student.uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
!> This module gathers information for the scalar cylindrical wave.
!!
module tem_cylindricalWave_module

  ! include treelm modules
  use env_module,         only: rk
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit

  ! include aotus modules
  use aotus_module,     only: aoterr_NonExistent, flu_state
  use aot_table_module, only: aot_table_open, aot_table_close,                 &
   &                          aot_table_length, aot_get_val

  implicit none

  public :: tem_load_cylindricalWave, tem_cylindricalWave_type, &
          & tem_eval_cylindricalWave

  !> This type contains datas to define the scalar cylindrical wave.
  type tem_cylindricalWave_type
    !> The wave order
    integer :: order
    !> The radial constant
    real(kind=rk) :: radialConstant
    !> The cut radius of the cylindrical wave
    real(kind=rk) :: radius
  end type tem_cylindricalWave_type


contains

  !> Load definition of the scalar cylindrical wave.
  subroutine tem_load_cylindricalWave(conf, thandle, me)
    ! ---------------------------------------------------------------------------
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: thandle
    type(tem_cylindricalWave_type), intent(out) :: me
    ! ---------------------------------------------------------------------------
    integer :: iError
    ! ---------------------------------------------------------------------------

    write(logUnit(0),*)'ERROR: the cylindrical wave function is not available' &
      &                // ' in this'
    write(logUnit(0),*)'       executable, due to the compiler lacking ' &
      &                // 'support of'
    write(logUnit(0),*)'       Bessel functions!'
    call tem_abort()

  end subroutine tem_load_cylindricalWave


  !> Calculate the function values for the cylindrical wave.
  function tem_eval_cylindricalWave(me, coord, time, n)  result(res)
    !> Spacetime function to evaluate
    type(tem_cylindricalWave_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> The current physical time
    real(kind=rk), intent( in )  :: time
    !> return value
    real(kind=rk) :: res(n)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: polar(n,2)
    integer :: i
    ! ---------------------------------------------------------------------------

    res = 0.0_rk

  end function tem_eval_cylindricalWave


end module tem_cylindricalWave_module
