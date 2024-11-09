! Copyright (c) 2012-2013 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2013, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> This module gathers the various predefined initial conditions
!!
module tem_ic_predefs_module

  ! include treelm modules
  use env_module,         only: rk
  use tem_param_module,   only: PI
  use tem_logging_module, only: logUnit
  use tem_aux_module,     only: tem_abort

  ! include aotus modules
  use aotus_module,     only: aoterr_NonExistent, flu_state
  use aot_table_module, only: aot_table_open, aot_table_close,                 &
   &                          aot_table_length, aot_get_val

  implicit none

  private

  public :: ic_gausspulse_type
  public :: load_ic_gausspulse
  public :: ic_gausspulse_for

  public :: ic_2dcrvp_type
  public :: load_ic_2dcrvp
  public :: ic_2dcrvpX_for, ic_2dcrvpY_for, ic_2dcrvpPressure_for

  public :: ic_tgv_type
  public :: load_ic_tgv
  public :: ic_tgv_pressure_for
  public :: ic_tgv_ux_for, ic_tgv_uy_for
  public :: ic_tgv_Sxx_for, ic_tgv_Syy_for, ic_tgv_Sxz_for, ic_tgv_Syz_for

  !> cutoff radius definition
  type cutoff_type
    !> cutoff is active?
    logical :: active
    !> cutoff values
    real(kind=rk) :: length
    !> cutoff start
    real(kind=rk) :: r_min
    !> cutoff end
    real(kind=rk) :: r_max
    !> linear behavior
    logical :: linear
    !> quadratic behavior
    logical :: quadratic
  end type cutoff_type

  !> This type contains datas to define gauss pulse
  type ic_gausspulse_type
    !> Gauss pulse center
    real(kind=rk) :: center(3)
    !> half width of gauss pulse from center
    real(kind=rk) :: halfwidth
    !> height or magnitude of gauss pulse
    real(kind=rk) :: amplitude
    !> reference value. In case of density, it is reference density
    real(kind=rk) :: background
    !> spatial step size
    real(kind=rk) :: dx
    !> time step size
    real(kind=rk) :: dt
  end type ic_gausspulse_type

  !> This type contains datas to define 2d co-rotating vortex pair
  type ic_2dcrvp_type
    !> spinning center
    real(kind=rk) :: center(3)
    !> distance of vortex centers / 2
    real(kind=rk) :: radius_rot
    !> core radius = radius_rot/3
    real(kind=rk) :: radius_C
    !> circulation of vortices
    real(kind=rk) :: circulation
    !> reference pressure
    real(kind=rk) :: p0
    !> reference density
    real(kind=rk) :: rho0
    !> adiabatic exponent
    real(kind=rk) :: kappa
    !> speed of sound
    real(kind=rk) :: cs
    !> rotating Mach number
    real(kind=rk) :: Ma
    !> position in time
    real(kind=rk) :: t
    type(cutoff_type) :: cutoff
    !> Approximation of the pressure distribution inside the core
    !! radius with a gaussian pulse model
    logical :: pressGaussModel
    !> vortex core velocity model: rankine
    logical :: rankineModel
    !> to match the gauss model to the pressure distribution
    !! Set to 2.2
    real(kind=rk) :: matchFactor
  end type ic_2dcrvp_type

  type ic_tgv_type
    !> Origin points
    real(kind=rk) :: x0(3)

    !> Ref velocity (X and Y)
    real(kind=rk) :: u0(2)

    !> Ref pressure
    real(kind=rk) :: p0

    !> Rate of decay coefficient
    real(kind=rk) :: tD

    !> Reynolds number
    real(kind=rk) :: Re

  end type ic_tgv_type

contains

! ****************************************************************************** !
  !> load gauss pulse variables to set initial condition
  !!
  subroutine load_ic_gausspulse(conf, thandle, me)
    ! ---------------------------------------------------------------------------
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: thandle
    !> Global gauss pulse data type
    type(ic_gausspulse_type), intent(out) :: me
    ! ---------------------------------------------------------------------------
    integer :: cent_handle
    integer :: i
    integer :: iError
    ! ---------------------------------------------------------------------------

    !center
    call aot_table_open( L       = conf,                                       &
      &                  parent  = thandle,                                    &
      &                  thandle = cent_handle,                                &
      &                  key     = 'center' )
    if (aot_table_length(L=conf, thandle=cent_handle) == 3) then
      do i=1,3
        call aot_get_val( L       = conf,                                      &
          &               thandle = cent_handle,                               &
          &               pos     = i,                                         &
          &               val     = me%center(i),                              &
          &               ErrCode = iError )
      end do
    else
      write(*,*) 'ERROR while reading the center of a gauss pulse,'
      write(*,*) 'should have 3 entries for each coordinate!'
      STOP
    end if
    call aot_table_close(conf, cent_handle)
    write(logUnit(1),*) ' * center =', me%center

    !halfwidth
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'halfwidth',                                   &
      &               val     = me%halfwidth,                                  &
      &               ErrCode = iError )
    write(logUnit(1),*) ' * halfwidth =', me%halfwidth

    !amplitude
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'amplitude',                                   &
      &               val     = me%amplitude,                                  &
      &               ErrCode = iError )
    write(logUnit(1),*) ' * amplitude =', me%amplitude

    !backgroud
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'background',                                  &
      &               val     = me%background,                                 &
      &               ErrCode = iError )
    write(logUnit(1),*) ' * background=', me%background

  end subroutine load_ic_gausspulse
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function defines gauss pulse
  !!
  !> This function computes gauss pulse for given array
  !! co-ordinate points and defined gauss parameters in LUA file.
  !! Gauss function:
  !!\( f(x) = a e^{\frac{(x-b)^2}{2c^2}} \)
  !! where,
  !!```
  !! a - pulse height,
  !! x - pulse center,
  !! c - pulse widthn,
  !! b - position along x
  !!```
  function ic_gausspulse_for(me, coord, n) result(res)
    ! ---------------------------------------------------------------------------
    !> number of return values
    integer, intent(in) :: n
    !> global gauss pulse data
    type(ic_gausspulse_type), intent(in) :: me
    !> coordinate of an element
    real(kind=rk), intent(in) :: coord(n, 3)
    !> return value which is sent to state variable
    real(kind=rk) :: res(n)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: distsquare(n)
    real(kind=rk) :: fact
    ! ---------------------------------------------------------------------------

    fact = -log(2.0_rk)/(me%halfwidth**2)

    distsquare = (coord(:,1) - me%center(1))**2                                &
      &        + (coord(:,2) - me%center(2))**2                                &
      &        + (coord(:,3) - me%center(3))**2

    distsquare = fact*distsquare

    res = me%background + (me%amplitude * exp(distsquare))

  end function ic_gausspulse_for
! ****************************************************************************** !


! ****************************************************************************** !
  !> load crvp variables to set initial condition
  !!
  subroutine load_ic_2dcrvp( conf, thandle, me )
    ! ---------------------------------------------------------------------------
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: thandle
    !> Global gauss pulse data type
    type(ic_2dcrvp_type), intent(out) :: me
    ! ---------------------------------------------------------------------------
    integer :: cent_handle
    integer :: i
    integer :: iError
    ! ---------------------------------------------------------------------------

    !center
    call aot_table_open( L       = conf,                                       &
      &                  parent  = thandle,                                    &
      &                  thandle = cent_handle,                                &
      &                  key     = 'center' )
    if (aot_table_length(L=conf, thandle=cent_handle) == 3) then
      do i=1,3
        call aot_get_val( L       = conf,                                      &
          &               thandle = cent_handle,                               &
          &               pos     = i,                                         &
          &               val     = me%center(i),                              &
          &               ErrCode = iError)
      end do
    else
      write(*,*) 'ERROR while reading the center of a vortices,'
      write(*,*) 'should have 3 entries for each coordinate!'
      STOP
    end if
    call aot_table_close(conf, cent_handle)
    write(logUnit(1),*) ' * center =', me%center

    ! get physical time steps dx and dt
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 't',                                           &
      &               val     = me%t,                                          &
      &               ErrCode = iError,                                        &
      &               default = 0._rk )
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'kappa',                                       &
      &               val     = me%kappa,                                      &
      &               ErrCode = iError,                                        &
      &               default = 1._rk )
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'p0',                                          &
      &               val     = me%p0,                                         &
      &               ErrCode = iError,                                        &
      &               default = 0._rk )
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'rho0',                                        &
      &               val     = me%rho0,                                       &
      &               ErrCode = iError,                                        &
      &               default = 1._rk )
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'cs',                                          &
      &               val     = me%cs,                                         &
      &               ErrCode = iError )
    if(( iError .gt. 0 )) then
      if( me%p0 .le. 0._rk ) then
        write(logUnit(1),*)' No speed of sound is given and the reference  '
        write(logUnit(1),*)' pressure for CRVP is 0. Choose > 0! '
        write(logUnit(1),*)' Either give the speed of sound cs = ... or use '
        write(logUnit(1),*)' a reference pressure > 0'
        call tem_abort
      end if
      me%cs = sqrt( me%kappa*me%p0/me%rho0 )
    end if
    if( me%cs .le. 0._rk ) then
      write(logUnit(1),*)' Error: Speed of sound <= 0  '
      call tem_abort
    end if

    ! Gaussian vortex core model for pressure
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'gpmodel',                                     &
      &               val     = me%pressGaussModel,                            &
      &               ErrCode = iError,                                        &
      &               default = .true. )
    me%matchFactor = 2.2_rk
    if( me%pressGaussModel ) then
      write(logUnit(1),*) ' * gaussian model for pressure in rC'
    else
      write(logUnit(1),*) ' * no model for pressure in rC'
    end if

    ! Rankine vortex core model for velocity
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'rankine',                                     &
      &               val     = me%rankineModel,                               &
      &               ErrCode = iError,                                        &
      &               default = .true. )
    if( me%rankineModel ) then
      write(logUnit(1),*) ' * Rankine vortex core model'
    end if

    !halfwidth
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'radius_rot',                                  &
      &               val     = me%radius_rot,                                 &
      &               ErrCode = iError )
    write(logUnit(1),*) ' * radius_rot =', me%radius_rot
    me%radius_C = me%radius_rot / 3._rk

    !amplitude
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'circulation',                                 &
      &               val     = me%circulation,                                &
      &               ErrCode = iError )
    write(logUnit(1),*) ' * circulation =', me%circulation

    me%Ma = me%circulation/(4._rk*pi*me%radius_rot*me%cs)
    ! some more info to the crvp initialization
    write(logUnit(1),*) ' * Ma_rot      =', me%Ma


    ! Set the cutoff definition
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'cutoff_length',                               &
      &               val     = me%cutoff%length,                              &
      &               ErrCode = iError )
    if ( btest( iError, aotErr_NonExistent ) ) then
      me%cutoff%active = .false.
    else
      me%cutoff%active = .true.
    endif

    !me%cutoff%r_min = 0.6_rk
    !me%cutoff%r_max = 0.95_rk
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'cutoff_rmin',                                 &
      &               val     = me%cutoff%r_min,                               &
      &               ErrCode = iError,                                        &
      &               default = huge(me%cutoff%r_min))
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'cutoff_rmax',                                 &
      &               val     = me%cutoff%r_max,                               &
      &               ErrCode = iError,                                        &
      &               default = huge(me%cutoff%r_max))
    if( me%cutoff%active ) then
      me%cutoff%linear    = .false.
      me%cutoff%quadratic = .true.
      write(logUnit(1),*)' * Cutoff length = ', me%cutoff%length
      if( me%cutoff%linear ) then
        write(logUnit(1),*)' * Cutoff progress linear'
      end if
      if( me%cutoff%quadratic ) then
        write(logUnit(1),*)' * Cutoff progress quadratic'
      end if
    endif

  end subroutine load_ic_2dcrvp
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function defines the y-velocity component of the
  !!   spinning (= co-rotating) vortex pair
  !! Source: complex velocity potential of both vortices
  !! complex coordinates:
  !! `z = x+i*y`
  !! Gamma ... circulation
  !!```
  !! b = r0*exp(i*omega*t)
  !! w(z,t) = Gamma/(2Pi*i)*ln(z^2-b^2)
  !! dw/dz = Gamma/(2Pi*i)*z/(z^2-b^2)
  !! u =  Re( dw/dz( z, t=0 )
  !! v = -Im( dw/dz( z, t=0 )
  !!```
  !! Unit of the result is in m/s, as the coordinates are given in physical
  !! coordinates and hence all other parameters also have to be physical ones
  !! As the potential induces a singularity inside the vortex,
  !! a vortex core model is employed. Here we use the Rankine vortex, where the
  !! velocity field inside the core radius `rC = 1/3` `r0` is approximated with a
  !! linear profile.
  !!
  function ic_2dcrvpX_for(me, coord, n) result(res)
    ! ---------------------------------------------------------------------------
    !> number of return values
    integer, intent(in) :: n
    !> global gauss pulse data
    type(ic_2dcrvp_type), intent(in) :: me
    !> coordinate of an element
    real(kind=rk), intent(in) :: coord(n, 3)
    !> return value which is sent to state variable
    real(kind=rk) :: res(n)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: r         ! radius from center of rotation
    real(kind=rk) :: rC, x_max ! core radius
    real(kind=rk) :: omega     ! angular speed
    real(kind=rk) :: x, y
    complex(kind=rk) :: z_coord, z_compl, b_coord, vortex1, vortex2, umax, z_l, z_r
    complex(kind=rk),parameter :: i_complex = (0.,1.)
    integer :: iElem
    ! ---------------------------------------------------------------------------

    omega   = me%circulation / (Pi*me%radius_rot*me%radius_rot*4._rk)
    b_coord = me%radius_rot*( cos(omega*me%t) + i_complex*sin(omega*me%t))
    rC = me%radius_C

    do iElem = 1, n
      x = coord(iElem,1) - me%center(1)
      y = coord(iElem,2) - me%center(2)
      z_coord = x + i_complex * y
      ! complex coordinate with respect to left  vortex
      z_l = z_coord + b_coord
      ! complex coordinate with respect to right vortex
      z_r = z_coord - b_coord

      ! left vortex. compare potential_single_w_z
      vortex1 = me%circulation/(2._rk*PI*i_complex*z_l)
      ! right vortex
      vortex2 = me%circulation/(2._rk*PI*i_complex*z_r)
      ! Get the velocity on the core radius for the core model

      if( me%rankineModel ) then
        x_max = -real(b_coord + rC, kind=rk )
        umax  =  me%circulation/(2._rk*PI*i_complex*(x_max + b_coord))
        z_compl = x - i_complex * y
        ! Rankine vortex model, if inside the core radius
        if(     real(z_r*(z_compl-b_coord), kind=rk) .lt. rC*rC ) then
          ! right vortex
          vortex2 = - (+umax*z_compl/rC - umax*b_coord/rC)
        elseif( real(z_l*(z_compl+b_coord), kind=rk) .lt. rC**2 ) then
          ! left vortex
          vortex1 = - (+umax*z_compl/rC + umax*b_coord/rC)
        end if
      end if

      r = abs(z_coord)
      res( iElem ) = cutoff_factor(me%cutoff, r)         &
        &            * (real( vortex1 + vortex2, kind=rk))

    end do

  end function ic_2dcrvpX_for
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function defines the y-velocity component of the
  !!   spinning (= co-rotating) vortex pair
  !! Source: complex velocity potential of both vortices
  !! complex coordinates:
  !! `z = x+i*y`
  !! Gamma ... circulation
  !!```
  !! b = r0*exp(i*omega*t)
  !! w(z,t) = Gamma/(2Pi*i)*ln(z^2-b^2)
  !! dw/dz = Gamma/(2Pi*i)*z/(z^2-b^2)
  !! u =  Re( dw/dz( z, t=0 )
  !! v = -Im( dw/dz( z, t=0 )
  !!```
  !! see ic_2dcrvpX_for as a documentation reference.
  !! This routine is exactly the same, except for that in the end, instead of
  !! evaluating the Re of the potential function, we have to evalute -Im
  !!
  function ic_2dcrvpY_for(me, coord, n) result(res)
    ! ---------------------------------------------------------------------------
    !> number of return values
    integer, intent(in) :: n
    !> global gauss pulse data
    type(ic_2dcrvp_type), intent(in) :: me
    !> coordinate of an element
    real(kind=rk), intent(in) :: coord(n, 3)
    !> return value which is sent to state variable
    real(kind=rk) :: res(n)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: r      ! radius from center of rotation
    real(kind=rk) :: rC, x_max ! core radius
    real(kind=rk) :: omega     ! angular speed
    complex(kind=rk) :: z_coord, b_coord, z_compl, vortex1, vortex2, umax
    complex(kind=rk),parameter :: i_complex = (0._rk,1._rk)
    integer :: iElem
    ! ---------------------------------------------------------------------------
    do iElem = 1, n
      z_coord    = ( coord(iElem,1) - me%center(1))                            &
        &        + i_complex*(coord(iElem,2) - me%center(2))
      z_compl    = ( coord(iElem,1) - me%center(1))                            &
        &        - i_complex*(coord(iElem,2) - me%center(2))

      r = abs(z_coord)
      rC = me%radius_C

      omega = me%circulation/(Pi*me%radius_rot*me%radius_rot*4._rk)
      b_coord    = me%radius_rot*( cos( omega*me%t )                           &
        &        + i_complex*sin( omega*me%t))

      ! left vortex. compare potential_single_w_z
      vortex1 =  me%circulation/(2._rk*PI*i_complex*(z_coord+b_coord))
      ! right vortex
      vortex2 =  me%circulation/(2._rk*PI*i_complex*(z_coord-b_coord))
      ! Get the velocity on the core radius for the core model
      x_max   = -real(b_coord + rC, kind=rk )
      umax    =  me%circulation/(2._rk*PI*i_complex*(x_max + b_coord))
      if( me%rankineModel ) then
      ! Rankine vortex model, if inside the core radius
      ! right vortex
        if( real(((z_coord) - b_coord)*(z_compl-b_coord), kind=rk)             &
          &                                                   .lt. rC*rC ) then
          vortex2 = -(umax*z_compl/rC - umax*b_coord/rC)
        elseif( real((z_coord + b_coord)*(z_compl + b_coord), kind=rk)         &
          &                                                   .lt. rC**2 ) then
          vortex1 = -(umax*z_compl/rC + umax*b_coord/rC)
        endif
      endif
      res( iElem ) = cutoff_factor(me%cutoff, r)                               &
        &          * (-aimag( vortex1 + vortex2))
    end do

  end function ic_2dcrvpY_for
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function defines the density of the
  !!   spinning (= co-rotating) vortex pair
  !! See the matlab file where the pressure is plot
  !! in the ase-testcases/ repo in
  !! musubi/crvp/matlab/crvp_velPress_plot.m
  !!
  !! As a reference, see
  !! [1] R. Fortenbach, 'Mehrskalenmodellierung von aeroakustischen Quellen in
  !! schwach kompressiblen Stroemungen,' pp. 1-151, Jul. 2006.
  !!
  !! Source: complex velocity potential of both vortices
  !! complex coordinates:
  !! `z = x+i*y`
  !! Gamma ... circulation
  !!```
  !! b = r0*exp(i*omega*t) ... rotation orbit
  !! w(z,t) = Gamma/(2Pi*i)*ln(z^2-b^2) ... potential function
  !! dw/dz = Gamma/(2Pi*i)*z/(z^2-b^2)  ... derivative of potential
  !! u =  Re( dw/dz( z, t=0 )           ... x -velocity
  !! v = -Im( dw/dz( z, t=0 )           ... y -velocity
  !! u0 = Gamma/(4Pi*r0)                ... rotation velocity at center of each
  !!                                        vortice
  !! omega = 2Pi/T0 = u0/r0 = Gamma/(4Pi*ro^2) ... rotation angular frequency
  !! Ma = u0/cs        ... rotation Mach number
  !! rho = rho0 - Ma^2*rho0/cs^2*( Re{ -omega*Gamma/Pi * b^2/(z^2-b^2)}
  !!                                                   + (u^2+v^2)/2 )
  !!```
  !! Unit of the result is in kg/m^3, as the coordinates are given in physical
  !! coordinates and hence all other parameters also have to be physical ones
  !!
  function ic_2dcrvpPressure_for( me, coord, n ) result(pressure)
    ! ---------------------------------------------------------------------------
    !> number of return values
    integer, intent(in) :: n
    !> global gauss pulse data
    type(ic_2dcrvp_type), intent(in) :: me
    !> coordinate of an element
    real(kind=rk), intent(in) :: coord(n, 3)
    !> return value which is sent to state variable
    real(kind=rk) :: pressure(n)
    ! ---------------------------------------------------------------------------
    complex(kind=rk) :: z_coord, b_coord, z_compl, b_coordChg, z_Core
    real(kind=rk) :: omega
    real(kind=rk) :: u(1),v(1)
    complex(kind=rk),parameter :: i_complex = (0._rk,1._rk)
    real(kind=rk) :: r, rC
    real(kind=rk) :: dist2Left, dist2Right, pMin, phiTC, uxC, uyC
    integer :: iElem
    ! ---------------------------------------------------------------------------

    do iElem = 1, n
      ! complex coordinate
      z_coord    = ( coord(iElem,1) - me%center(1))                            &
        &        + i_complex*(coord(iElem,2) - me%center(2))
      z_compl    = ( coord(iElem,1) - me%center(1))                            &
        &        - i_complex*(coord(iElem,2) - me%center(2))

      omega = me%circulation/(Pi*me%radius_rot*me%radius_rot*4._rk)
      ! Current radial coordinate
      r = abs( z_coord )
      rC = me%radius_C
      b_coord    = me%radius_rot*( cos( omega*me%t )                           &
        &        + i_complex*sin( omega*me%t))
      b_coordChg = b_coord ! *3._rk
      z_Core = b_coord + rC

      ! Calculate velocity
      u(:)   = ic_2dcrvpX_for( me, coord( iElem, :), 1 )
      v(:)   = ic_2dcrvpY_for( me, coord( iElem, :), 1 )

      ! Gaussian pulse approimxation inside core radius
      ! Compute the distance^2
      dist2Left  = real((z_coord + b_coord)*(z_compl + b_coordChg), kind=rk)
      dist2Right = real(((z_coord) - b_coord)*(z_compl-b_coordChg), kind=rk)
      if( dist2Right .lt. rC*rC .or. dist2Left .lt. rC*rC ) then
        uxC = real( me%circulation/(2._rk*PI*i_complex*(z_Core-b_coord)),      &
          &                                                            kind=rk)
        uyC = real(-aimag( me%circulation/(2._rk*PI*i_complex*                 &
          &                                        (z_Core-b_coord))), kind=rk)
        phiTC = real( -omega*me%circulation/PI*b_coord**2/                     &
          &                                   (z_Core**2-b_coord**2), kind=rk )
        pMin = me%rho0*( phiTC + (uxC**2+uyC**2)*0.5_rk);
        pMin = pMin * me%matchFactor
      end if

      if( dist2Right .lt. rC*rC ) then
        dist2Right = -0.5_rk/(rC*rC)*dist2Right
        if( me%pressGaussModel ) then
          ! velocity at the core radius: real(b_coord) + rC
          ! right vortex
          pressure( iElem ) = me%p0 - cutoff_factor( me%cutoff, r )            &
            &               * pMin *exp( dist2Right )
        else
          pressure( iElem ) = me%p0
        end if

      elseif( dist2Left .lt. rC*rC ) then
        dist2Left  = -0.5_rk/(rC*rC)*dist2Left
        if( me%pressGaussModel ) then
          ! velocity at the core radius: real(b_coord) + rC
          pressure( iElem ) = me%p0 - cutoff_factor( me%cutoff, r )            &
            &               * pMin*exp( dist2Left )
        else
          pressure( iElem ) = me%p0
        end if
      else
        pressure( iElem ) = ( - cutoff_factor(me%cutoff, r)                    &
          &               * ( real( -omega*me%circulation/PI*b_coord*b_coord/  &
          &                                 (z_coord*z_coord-b_coord*b_coord)) &
          &               + (u(1)*u(1) + v(1)*v(1))*0.5_rk ))*me%rho0 + me%p0
      end if
    end do

  end function ic_2dcrvpPressure_for
! ****************************************************************************** !


! ****************************************************************************** !
  !> return the cutoff multiplication factor
  !! This routine returns the cutoff factor for a circle of size r_min.
  !! Outside r_min, the quantity is
  !!
  !! - for me%linear == .true. :
  !!   linearly reduced to 0 until r_max.
  !! - for me%quadratic == .true. :
  !!   quadratically reduced to 0 until r_max.
  !! outside the radius r_max, the cutoff factor is set to zero
  !!
  function cutoff_factor(me, radius) result(cutoff_fac)
    ! ---------------------------------------------------------------------------
    !> global gauss pulse data
    type(cutoff_type), intent(in) :: me
    !> coordinate of an element
    real(kind=rk), intent(in) :: radius
    !> return value which is sent to state variable
    real(kind=rk) :: cutoff_fac
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: r_min, r_max ! minimum and maximum absolute radius
    real(kind=rk) :: a0, a1, a2   ! polynomial coefficients
    ! ---------------------------------------------------------------------------

    ! Define no cutoff as the default
    cutoff_fac = 1._rk

    if( me%active ) then
      ! If the cutoff is active ...
      ! first compute the absolute radius from the domain center:
      ! min for where to start cutting of
      r_min = me%length*me%r_min
      ! max for where to end cutting of
      r_max = me%length*me%r_max
      if( radius .le. r_min ) then
        cutoff_fac = 1._rk
      elseif( radius .gt. r_max ) then
        cutoff_fac = 0._rk
      else
        if( me%linear ) then
          ! Linear progress from r_min towards r_max
          cutoff_fac = 1._rk - (radius - r_min) / (r_max-r_min)
        elseif( me%quadratic ) then
          ! Quadratic progress from r_min towards r_max,
          ! where the derivative at r_min is zero for a smooth progression from
          ! the domain inside
          a0 = 1._rk / ( -r_min*r_min -r_max*r_max + 2._rk*r_min*r_max)
          a1 = -2._rk*a0*r_min
          a2 = 1._rk - a0*r_min*r_min - a1*r_min
          cutoff_fac = a0*radius*radius + a1*radius + a2
        end if
      endif
    endif

  end function cutoff_factor
! ****************************************************************************** !


! ****************************************************************************** !
  !> load gauss pulse variables to set initial condition
  !!
  subroutine load_ic_tgv(conf, thandle, me)
    ! ---------------------------------------------------------------------------
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: thandle
    !> TGV data type
    type(ic_tgv_type), intent(out) :: me
    ! ---------------------------------------------------------------------------
    integer :: cent_handle
    integer :: i
    integer :: iError
    ! ---------------------------------------------------------------------------

    ! Load x0
    call aot_table_open( L       = conf,                                       &
      &                  parent  = thandle,                                    &
      &                  thandle = cent_handle,                                &
      &                  key     = 'x0' )
    if (aot_table_length(L=conf, thandle=cent_handle) == 3) then
      do i=1,3
        call aot_get_val( L       = conf,         &
          &               thandle = cent_handle,  &
          &               pos     = i,            &
          &               val     = me%x0(i),     &
          &               ErrCode = iError )
      end do
    else
      me%x0 = [0.0_rk, 0.0_rk, 0.0_rk]
    end if
    call aot_table_close(conf, cent_handle)

    ! Load x0
    call aot_table_open( L       = conf,                                       &
      &                  parent  = thandle,                                    &
      &                  thandle = cent_handle,                                &
      &                  key     = 'u0' )
    if (aot_table_length(L=conf, thandle=cent_handle) == 2) then
      do i=1,2
        call aot_get_val( L       = conf,         &
          &               thandle = cent_handle,  &
          &               pos     = i,            &
          &               val     = me%u0(i),     &
          &               ErrCode = iError )
      end do
    else
      me%u0 = [1.0_rk, 1.0_rk]
    end if
    call aot_table_close(conf, cent_handle)

    ! Load p0
    call aot_get_val( L       = conf,    &
      &               thandle = thandle, &
      &               key     = 'p0',    &
      &               val     = me%p0,   &
      &               ErrCode = iError,  &
      &               default = 0._rk )

    ! Load Re
    call aot_get_val( L       = conf,    &
      &               thandle = thandle, &
      &               key     = 'Re',    &
      &               val     = me%Re,   &
      &               ErrCode = iError,  &
      &               default = 25._rk )

    me%tD = me%Re * 0.5_rk

    write(logUnit(5),"(A,3F8.4)") ' x0 = ', me%x0
    write(logUnit(5),"(A,2F8.4)") ' u0 = ', me%u0
    write(logUnit(5),"(A,1F8.4)") ' p0 = ', me%p0
    write(logUnit(5),"(A,1F8.4)") ' Re = ', me%Re
    write(logUnit(5),"(A,1F8.4)") ' tD = ', me%tD

  end subroutine load_ic_tgv
! ****************************************************************************** !

! ****************************************************************************** !
  pure function ic_tgv_pressure_for( me, coord, n ) result( pressure )
    ! ---------------------------------------------------------------------------
    !> global gauss pulse data
    type(ic_tgv_type), intent(in) :: me
    !> number of return values
    integer, intent(in) :: n
    !> coordinate of an element
    real(kind=rk), intent(in) :: coord(n, 3)
    !> return value which is sent to state variable
    real(kind=rk) :: pressure(n)
    ! ---------------------------------------------------------------------------
    integer :: iElem
    real(kind=rk) :: xc, yc, p1, p2 !, zc
    ! ---------------------------------------------------------------------------

    do iElem = 1, n
      xc = 2._rk * ( coord(iElem,1) - me%x0(1) )
      yc = 2._rk * ( coord(iElem,2) - me%x0(2) )
      !zc = 2._rk * ( coord(iElem,3) - me%x0(3) )
      ! p1 = cos(xc) * cos(zc) + 2._rk * cos(yc)
      ! p2 = cos(yc) * cos(zc) + 2._rk * cos(xc)
      p1 = cos(xc)
      p2 = cos(yc)
      pressure(iElem) = me%p0 - (p1+p2) * 0.25_rk
    end do

  end function ic_tgv_pressure_for
! ****************************************************************************** !

! ****************************************************************************** !
  pure function ic_tgv_ux_for( me, coord, n ) result( ux )
    ! ---------------------------------------------------------------------------
    !> global gauss pulse data
    type(ic_tgv_type), intent(in) :: me
    !> number of return values
    integer, intent(in) :: n
    !> coordinate of an element
    real(kind=rk), intent(in) :: coord(n, 3)
    !> return value which is sent to state variable
    real(kind=rk) :: ux(n)
    ! ---------------------------------------------------------------------------
    integer :: iElem
    ! ---------------------------------------------------------------------------

    do iElem = 1, n
      ux(iElem) =-me%u0(1) * ( cos(coord(iElem,1)-me%x0(1)) &
        &                    * sin(coord(iElem,2)-me%x0(2)) )
    end do
  end function ic_tgv_ux_for
! ****************************************************************************** !

! ****************************************************************************** !
  pure function ic_tgv_uy_for( me, coord, n ) result( uy )
    ! ---------------------------------------------------------------------------
    !> global gauss pulse data
    type(ic_tgv_type), intent(in) :: me
    !> number of return values
    integer, intent(in) :: n
    !> coordinate of an element
    real(kind=rk), intent(in) :: coord(n, 3)
    !> return value which is sent to state variable
    real(kind=rk) :: uy(n)
    ! ---------------------------------------------------------------------------
    integer :: iElem
    ! ---------------------------------------------------------------------------

    do iElem = 1, n
      uy(iElem) =  me%u0(2) * ( sin(coord(iElem,1)-me%x0(1)) &
        &                     * cos(coord(iElem,2)-me%x0(2)) )
    end do
  end function ic_tgv_uy_for
! ****************************************************************************** !

! ****************************************************************************** !
  pure function ic_tgv_Sxx_for( me, coord, n ) result( s )
    ! ---------------------------------------------------------------------------
    !> global gauss pulse data
    type(ic_tgv_type), intent(in) :: me
    !> number of return values
    integer, intent(in) :: n
    !> coordinate of an element
    real(kind=rk), intent(in) :: coord(n, 3)
    !> return value which is sent to state variable
    real(kind=rk) :: s(n)
    ! ---------------------------------------------------------------------------
    integer :: iElem
    ! ---------------------------------------------------------------------------

    do iElem = 1, n
      s(iElem) = me%u0(1) * ( sin(coord(iElem,1)-me%x0(1)) &
        &                   * sin(coord(iElem,2)-me%x0(2)) )
    end do
  end function ic_tgv_sxx_for
! ****************************************************************************** !

! ****************************************************************************** !
  pure function ic_tgv_Syy_for( me, coord, n ) result( s )
    ! ---------------------------------------------------------------------------
    !> global gauss pulse data
    type(ic_tgv_type), intent(in) :: me
    !> number of return values
    integer, intent(in) :: n
    !> coordinate of an element
    real(kind=rk), intent(in) :: coord(n, 3)
    !> return value which is sent to state variable
    real(kind=rk) :: s(n)
    ! ---------------------------------------------------------------------------
    integer :: iElem
    ! ---------------------------------------------------------------------------

    do iElem = 1, n
      s(iElem) = -me%u0(2) * ( sin(coord(iElem,1)-me%x0(1)) &
        &                    * sin(coord(iElem,2)-me%x0(2)) )
    end do
  end function ic_tgv_syy_for
! ****************************************************************************** !

! ****************************************************************************** !
  pure function ic_tgv_Sxz_for( me, coord, n ) result( s )
    ! ---------------------------------------------------------------------------
    !> global gauss pulse data
    type(ic_tgv_type), intent(in) :: me
    !> number of return values
    integer, intent(in) :: n
    !> coordinate of an element
    real(kind=rk), intent(in) :: coord(n, 3)
    !> return value which is sent to state variable
    real(kind=rk) :: s(n)
    ! ---------------------------------------------------------------------------
    integer :: iElem
    ! ---------------------------------------------------------------------------

    do iElem = 1, n
      s(iElem) = -me%u0(1) * ( sin(coord(iElem,1)-me%x0(1)) &
        &                    * cos(coord(iElem,2)-me%x0(2)) &
        &                    * sin(coord(iElem,3)-me%x0(3)) ) * 0.5_rk
    end do
  end function ic_tgv_Sxz_for
! ****************************************************************************** !

! ****************************************************************************** !
  pure function ic_tgv_Syz_for( me, coord, n ) result( s )
    ! ---------------------------------------------------------------------------
    !> global gauss pulse data
    type(ic_tgv_type), intent(in) :: me
    !> number of return values
    integer, intent(in) :: n
    !> coordinate of an element
    real(kind=rk), intent(in) :: coord(n, 3)
    !> return value which is sent to state variable
    real(kind=rk) :: s(n)
    ! ---------------------------------------------------------------------------
    integer :: iElem
    ! ---------------------------------------------------------------------------

    do iElem = 1, n
      s(iElem) =  me%u0(2) * ( cos(coord(iElem,1)-me%x0(1)) &
        &                    * sin(coord(iElem,2)-me%x0(2)) &
        &                    * sin(coord(iElem,3)-me%x0(3)) ) * 0.5_rk
    end do
  end function ic_tgv_Syz_for
! ****************************************************************************** !

end module tem_ic_predefs_module
! ****************************************************************************** !
