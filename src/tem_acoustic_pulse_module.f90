! Copyright (c) 2016, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
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
module tem_acoustic_pulse_module
  use env_module,         only: rk
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit

  use aotus_module,     only: flu_state, aot_get_val

  implicit none

  private

  !> Analytical solution for an acoustic wave emitted by a Gaussian pulse as
  !! described in Tam: Computational Acoustics, a wave number approach.
  !! Appendix G.3.
  !!
  !! The pulse is described geometrical by an origin and a halfwidth, where
  !! center is the location of the maximum in the spherical Gaussian pulse, and
  !! halwidth is the radius at which the pulse has half the value of the
  !! maximum.
  !! The height of the pulse is given by amplitude.
  !! Note, this definition matches the gausspulse spatial function, which can be
  !! used as an initial condition to this problem.
  !! Additionally the speed of sound and the background pressure are required
  !! to compute the solution.
  type tem_acoustic_pulse_type
    !> Pulse height in the initial condition at the center over the background
    !! value.
    real(kind=rk) :: amplitude

    !> Location of the pulse center in 3D.
    real(kind=rk) :: center(3)

    !> Radius of the initial pulse, where half the amplitude is reached.
    real(kind=rk) :: halfwidth

    !> A background value to use (result is given by background+pulse).
    real(kind=rk) :: background

    !> Speed of sound is the velocity by which the acoustic wave is to be
    !! transported.
    real(kind=rk) :: speed_of_sound
  end type tem_acoustic_pulse_type

  public :: tem_acoustic_pulse_type
  public :: tem_load_acoustic_pulse
  public :: tem_eval_acoustic_pulse


contains


  ! ------------------------------------------------------------------------ !
  !> Load the definition of an acoustic pulse from a configuration Lua script.
  subroutine tem_load_acoustic_pulse(conf, thandle, me)
    type(flu_State) :: conf
    integer, intent(in) :: thandle
    type(tem_acoustic_pulse_type), intent(out) :: me
    ! -------------------------------------------------------------------- !
    integer :: iError, vError(3)
    ! -------------------------------------------------------------------- !

    write(logunit(1),*) 'Loading predefined function for acoustic pulse:'

    call aot_get_val( L       = conf,         &
      &               thandle = thandle,      &
      &               key     = 'amplitude',  &
      &               val     = me%amplitude, &
      &               ErrCode = iError        )

    if (iError /= 0) then
      write(logUnit(1),*) 'ERROR in tem_load_acoustic_pulse: not able '
      write(logUnit(1),*) 'to read amplitude from config file.'
      call tem_abort()
    end if

    write(logUnit(1),*) ' * amplitude =', me%amplitude

    call aot_get_val( L       = conf,         &
      &               thandle = thandle,      &
      &               key     = 'halfwidth',  &
      &               val     = me%halfwidth, &
      &               ErrCode = iError        )

    if (iError /= 0) then
      write(logUnit(1),*) 'ERROR in tem_load_acoustic_pulse: not able '
      write(logUnit(1),*) 'to read halfwidth from config file.'
      call tem_abort()
    end if

    write(logUnit(1),*) ' * halfwidth =', me%halfwidth

    call aot_get_val( L       = conf,              &
      &               thandle = thandle,           &
      &               key     = 'speed_of_sound',  &
      &               val     = me%speed_of_sound, &
      &               ErrCode = iError             )

    if (iError /= 0) then
      write(logUnit(1),*) 'ERROR in tem_load_acoustic_pulse: not able '
      write(logUnit(1),*) 'to read speed_of_sound from config file.'
      call tem_abort()
    end if

    write(logUnit(1),*) ' * speed_of_sound =', me%speed_of_sound

    call aot_get_val( L       = conf,                     &
      &               thandle = thandle,                  &
      &               key     = 'center',                 &
      &               val     = me%center,                &
                      default = [0.0_rk, 0.0_rk, 0.0_rk], &
      &               ErrCode = vError                    )

    if ( any(vError /= 0) ) then
      write(logUnit(1),*) 'ERROR in tem_load_acoustic_pulse: not able '
      write(logUnit(1),*) 'to read center from config file.'
      call tem_abort()
    end if

    write(logUnit(1),*) ' * center =', me%center

    call aot_get_val( L       = conf,          &
      &               thandle = thandle,       &
      &               key     = 'background',  &
      &               val     = me%background, &
                      default = 0.0_rk,        &
      &               ErrCode = iError         )

    if (iError /= 0) then
      write(logUnit(1),*) 'ERROR in tem_load_acoustic_pulse: not able '
      write(logUnit(1),*) 'to read background from config file.'
      call tem_abort()
    end if

    write(logUnit(1),*) ' * background =', me%background

  end subroutine tem_load_acoustic_pulse
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Evaluate the acoustic pulse at given points in space for one point in time.
  !!
  !! Exact solution for an acoustic wave from a Gaussian pulse in pressure.
  !! See Tam: Computational Acoustics, a wave number approach. Appendix G.3.
  !! Any point may be probed, the solution at the center is properly defined
  !! with a finite value.
  function tem_eval_acoustic_pulse(me, coord, time, n) result(res)
    !> Definition of the acoustic pulse to evaluate
    type(tem_acoustic_pulse_type), intent(in) :: me

    !> Number of different points to evaluate the acoustic pulse at.
    integer, intent(in) :: n

    !> 3D Coordinates of all points.
    real(kind=rk), intent(in) :: coord(n,3)

    !> Point in time to evaluate the points at.
    real(kind=rk), intent(in) :: time

    !> Analytical solution in all n points.
    real(kind=rk) :: res(n)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: radius(n)
    real(kind=rk), parameter :: zero_rad = 16.0_rk * tiny(time)
    real(kind=rk) :: wavepos
    real(kind=rk) :: ampfact
    real(kind=rk) :: expfact
    ! -------------------------------------------------------------------- !

    radius = sqrt( (coord(:,1)-me%center(1))**2  &
      &           + (coord(:,2)-me%center(2))**2 &
      &           + (coord(:,3)-me%center(3))**2 )

    wavepos = me%speed_of_sound * time
    ampfact = 0.5_rk * me%amplitude
    expfact = -log(2.0_rk) / me%halfwidth**2

    where (radius > zero_rad)
      res = me%background &
        & + (ampfact/radius) * ( (radius-wavepos)                   &
        &                        * exp(expfact*(radius-wavepos)**2) &
        &                      + (radius+wavepos)                   &
        &                        * exp(expfact*(radius+wavepos)**2) )
    elsewhere
      res = me%background + me%amplitude * exp(expfact*wavepos**2)       &
        &                                * (1.0_rk + 2*expfact*wavepos**2)
    end where

  end function tem_eval_acoustic_pulse
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


end module tem_acoustic_pulse_module
