! Copyright (c) 2013 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
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
module tem_heaviside_gibbs_fun_module
  use env_module,         only: rk
  use tem_aux_module,     only: tem_abort
  use tem_param_module,   only: PI
  use tem_logging_module, only: logUnit

  use aotus_module,     only: flu_State,  aot_get_val

  implicit none
  private

  !> Defines a Heaviside function, including Gibbs oscillations.
  type tem_heaviside_gibbs_type
    !> The location of the jump
    real(kind=rk) :: center
    !> Approximation order
    integer :: order
    !> Asymptotic function value left of the jump
    real(kind=rk) :: left
    !> Asymptotic function value right of the jump
    real(kind=rk) :: right
  end type tem_heaviside_gibbs_type

  public :: tem_heaviside_gibbs_type, tem_load_heaviside_gibbs, tem_eval_heaviside_gibbs

  contains

  ! ****************************************************************************** !
  !> This subroutine loads the definition of a spatial Heaviside function
  !! including Gibbs oscillations occuring for a high order approximation.
  subroutine tem_load_heaviside_gibbs( conf, thandle, me )
    ! ---------------------------------------------------------------------------
    !> Heaviside function data
    type(tem_heaviside_gibbs_type),intent(out) :: me
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: thandle
    ! ---------------------------------------------------------------------------
    integer :: iError
    ! ---------------------------------------------------------------------------

    ! Center
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'center',                                      &
      &               val     = me%center,                                     &
      &               ErrCode = iError )
    ! Order
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'order',                                       &
      &               val     = me%order,                                      &
      &               ErrCode = iError )
    ! Left value
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'left',                                        &
      &               val     = me%left,                                       &
      &               ErrCode = iError )
    ! Right value
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'right',                                       &
      &               val     = me%right,                                      &
      &               ErrCode = iError )

    write(logUnit(1),*) ' Data for Heaviside function (including Gibbs'//      &
      &                'oscillations):'
    write(logUnit(1),*) ' * center =', me%center
    write(logUnit(1),*) ' * order  =', me%order
    write(logUnit(1),*) ' * left   =', me%left
    write(logUnit(1),*) ' * right  =', me%right

  end subroutine tem_load_heaviside_gibbs
! ****************************************************************************** !

  function tem_eval_heaviside_gibbs( me, coord, n) result(res)
    ! ---------------------------------------------------------------------------
    !> Description of the Heaviside function
    type(tem_heaviside_gibbs_type) :: me
    !> number of return values
    integer, intent( in ) :: n
    !> Coordinates to evaluate the function for
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value of the function
    real( kind=rk ) :: res(n)
    ! ---------------------------------------------------------------------------
    integer :: iPoint
    real(kind=rk) :: dist, z
    ! ---------------------------------------------------------------------------

    ! Loop over all the points
    do iPoint = 1, n

      ! Distance from the center
      dist = me%center - coord(iPoint,1)

      z = 0.5_rk * dist * sqrt(1.0_rk-me%center**2) * (2.0_rk * real(me%order,rk) + 1.0_rk)

      ! Calculate the point value
      res(iPoint) = 0.5 * (me%left + me%right) + ( dsinint(z)/PI ) * ( me%left - me%right )

    end do

  end function tem_eval_heaviside_gibbs

  !> Calculate sine integral of xvalue.
  !!   AUTHOR: Allan MacLeod
  !!           Dept. of Mathematics and Statistics
  !!           University of Paisley
  !!           Scotland
  !!           (e-mail: macl_ms0@paisley.ac.uk)
  function dsinint(xvalue) result(fn_val)
    ! ---------------------------------------------------------------------------
    real(kind=rk), intent(in) :: xvalue
    real(kind=rk)             :: fn_val
    ! ---------------------------------------------------------------------------
    INTEGER              :: i, indsgn
    real(kind=rk)             :: cx, fival, gival, sumden, sumnum, sx, x, xhigh, xsq
    real(kind=rk), parameter :: zero = 0.0_rk, one = 1.0_rk, six = 6.0_rk, &
                           twelve = 12.0_rk
    real(kind=rk), parameter :: piby2 = 1.5707963267948966192_rk
    real(kind=rk), parameter :: xlow = 4.47E-8_rk, xhigh1 = 2.32472E8_rk
    real(kind=rk), parameter :: xhigh2 = 9.0072E15_rk, xhigh3 = 1.4148475E16_rk
    real(kind=rk), parameter :: asintn(0:7) = (/ 1.0_rk,  &
             -0.44663998931312457298E-1_rk, 0.11209146443112369449E-2_rk,  &
             -0.13276124407928422367E-4_rk, 0.85118014179823463879E-7_rk,  &
             -0.29989314303147656479E-9_rk, 0.55401971660186204711E-12_rk, &
             -0.42406353433133212926E-15_rk /)
    real(kind=rk), parameter :: asintd(0:7) = (/ 1.0_rk,  &
              0.10891556624243098264E-1_rk, 0.59334456769186835896E-4_rk,  &
              0.21231112954641805908E-6_rk, 0.54747121846510390750E-9_rk,  &
              0.10378561511331814674E-11_rk, 0.13754880327250272679E-14_rk,&
              0.10223981202236205703E-17_rk /)
    real(kind=rk), parameter :: afn1(0:7) = (/ 0.99999999962173909991_rk,   &
             0.36451060338631902917E3_rk, 0.44218548041288440874E5_rk, &
             0.22467569405961151887E7_rk, 0.49315316723035561922E8_rk, &
             0.43186795279670283193E9_rk, 0.11847992519956804350E10_rk,&
             0.45573267593795103181E9_rk /)
    real(kind=rk), parameter :: afd1(0:7) = (/ 1.0_rk, 0.36651060273229347594E3_rk,  &
                     0.44927569814970692777E5_rk, 0.23285354882204041700E7_rk,  &
                     0.53117852017228262911E8_rk, 0.50335310667241870372E9_rk,  &
                     0.16575285015623175410E10_rk, 0.11746532837038341076E10_rk /)
    real(kind=rk), parameter :: agn1(0:8) = (/ 0.99999999920484901956_rk,     &
             0.51385504875307321394E3_rk, 0.92293483452013810811E5_rk,   &
             0.74071341863359841727E7_rk, 0.28142356162841356551E9_rk,   &
             0.49280890357734623984E10_rk, 0.35524762685554302472E11_rk, &
             0.79194271662085049376E11_rk, 0.17942522624413898907E11_rk /)
    real(kind=rk), parameter :: agd1(0:8) = (/ 1.0_rk, 0.51985504708814870209E3_rk,  &
                     0.95292615508125947321E5_rk, 0.79215459679762667578E7_rk,  &
                     0.31977567790733781460E9_rk, 0.62273134702439012114E10_rk,  &
                     0.54570971054996441467E11_rk, 0.18241750166645704670E12_rk,  &
                     0.15407148148861454434E12_rk /)
    real(kind=rk), parameter :: afn2(0:7) = (/ 0.19999999999999978257E1_rk,   &
             0.22206119380434958727E4_rk, 0.84749007623988236808E6_rk,   &
             0.13959267954823943232E9_rk, 0.10197205463267975592E11_rk,  &
             0.30229865264524075951E12_rk, 0.27504053804288471142E13_rk, &
             0.21818989704686874983E13_rk /)
    real(kind=rk), parameter :: afd2(0:7) = (/ 1.0_rk, 0.11223059690217167788E4_rk,  &
                     0.43685270974851313242E6_rk, 0.74654702140658116258E8_rk,  &
                     0.58580034751805687471E10_rk, 0.20157980379272098841E12_rk,  &
                     0.26229141857684496445E13_rk, 0.87852907334918467516E13_rk /)
    real(kind=rk), parameter :: agn2(0:8) = (/ 0.59999999999999993089E1_rk,   &
             0.96527746044997139158E4_rk, 0.56077626996568834185E7_rk,   &
             0.15022667718927317198E10_rk, 0.19644271064733088465E12_rk, &
             0.12191368281163225043E14_rk, 0.31924389898645609533E15_rk, &
             0.25876053010027485934E16_rk, 0.12754978896268878403E16_rk /)
    real(kind=rk), parameter :: agd2(0:8) = (/ 1.0_rk, 0.16287957674166143196E4_rk,  &
                      0.96636303195787870963E6_rk, 0.26839734750950667021E9_rk,  &
                      0.37388510548029219241E11_rk, 0.26028585666152144496E13_rk,  &
                      0.85134283716950697226E14_rk, 0.11304079361627952930E16_rk,  &
                      0.42519841479489798424E16_rk /)
    ! ---------------------------------------------------------------------------

    !   START COMPUTATION

    x = xvalue
    indsgn = 1
    IF ( x < zero ) THEN
      x = -x
      indsgn = -1
    END IF

    !   CODE FOR 0 <= |X| <= 6

    IF ( x <= six ) THEN
      IF ( x < xlow ) THEN
        fn_val = x
      ELSE
        sumnum = zero
        sumden = zero
        xsq = x * x
        DO i = 7 , 0 , -1
          sumnum = sumnum * xsq + asintn(i)
          sumden = sumden * xsq + asintd(i)
        END DO
        fn_val = x * sumnum / sumden
      END IF

    !   CODE FOR 6 < |X| <= 12

    ELSE IF ( x > six .AND. x <= twelve ) THEN
      sumnum = zero
      sumden = zero
      xsq = one / ( x * x )
      DO i = 7 , 0 , -1
        sumnum = sumnum * xsq + afn1(i)
        sumden = sumden * xsq + afd1(i)
      END DO
      fival = sumnum / ( x * sumden )
      sumnum = zero
      sumden = zero
      DO i = 8 , 0 , -1
        sumnum = sumnum * xsq + agn1(i)
        sumden = sumden * xsq + agd1(i)
      END DO
      gival = xsq * sumnum / sumden
      fn_val = piby2 - fival * COS(x) - gival * SIN(x)

    !   CODE FOR |X| > 12

    ELSE
      xhigh = MIN(xhigh2, xhigh3)
      IF ( x > xhigh ) THEN
        fn_val = piby2
      ELSE
        cx = COS(x)
        sx = SIN(x)
        xsq = one / ( x * x )
        IF ( x > xhigh1 ) THEN
          fn_val = piby2 - cx / x - sx * xsq
        ELSE
          sumnum = zero
          sumden = zero
          DO i = 7 , 0 , -1
            sumnum = sumnum * xsq + afn2(i)
            sumden = sumden * xsq + afd2(i)
          END DO
          fival =  ( one - xsq * sumnum / sumden ) / x
          sumnum = zero
          sumden = zero
          DO i = 8 , 0 , -1
            sumnum = sumnum * xsq + agn2(i)
            sumden = sumden * xsq + agd2(i)
          END DO
          gival =  ( one - xsq * sumnum / sumden ) * xsq
          fn_val = piby2 - cx * fival - sx * gival
        END IF
      END IF
    END IF
    IF ( indsgn == -1 ) fn_val = -fn_val
    RETURN
  end function dsinint


end module tem_heaviside_gibbs_fun_module
