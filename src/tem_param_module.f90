! Copyright (c) 2011-2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2012, 2014-2016, 2018-2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2015, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2021 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
!> author: Kannan Masilamani, Harald Klimach, Manuel Hasert, Jems Zudrop, Simon Zimny, ...
!! Parameter lists for Treelm
!!
!! This module holds the numbering of elements and subelements acoording
!! to the treelm numbering (Morton-curve).
!!
module tem_param_module
  use env_module, only:rk

  implicit none

  integer, parameter, public :: c_x  = 1
  integer, parameter, public :: c_y  = 2
  integer, parameter, public :: c_z  = 3
  !> Neighboring convention
  integer,parameter :: qQQQ     = 26    !< number of directions

  integer,parameter :: q__W     = 1     !< west               x-
  integer,parameter :: q__S     = 2     !< south              y-
  integer,parameter :: q__B     = 3     !< bottom             z-
  integer,parameter :: q__E     = 4     !< east               x+
  integer,parameter :: q__N     = 5     !< north              y+
  integer,parameter :: q__T     = 6     !< top                z+
  integer,parameter :: q_BS     = 7     !< bottom south       z-,y-
  integer,parameter :: q_TS     = 8     !< top south          z+,y-
  integer,parameter :: q_BN     = 9     !< bottom north       z-,y+
  integer,parameter :: q_TN     = 10    !< top north          z+,y+
  integer,parameter :: q_BW     = 11    !< bottom west        x-,z-
  integer,parameter :: q_BE     = 12    !< bottom east        x+,z-
  integer,parameter :: q_TW     = 13    !< top west           x-,z+
  integer,parameter :: q_TE     = 14    !< top east           x+,z+
  integer,parameter :: q_SW     = 15    !< south west         y-,x-
  integer,parameter :: q_NW     = 16    !< north west         y+,x-
  integer,parameter :: q_SE     = 17    !< south east         y-,x+
  integer,parameter :: q_NE     = 18    !< north east         y+,x+
  integer,parameter :: qBSW     = 19    !< bottom south west  x-,y-,z-
  integer,parameter :: qBSE     = 20    !< bottom south east  x+,y-,z-
  integer,parameter :: qBNW     = 21    !< bottom north west  x-,y+,z-
  integer,parameter :: qBNE     = 22    !< bottom north east  x+,y+,z-
  integer,parameter :: qTSW     = 23    !< top south west     x-,y-,z+
  integer,parameter :: qTSE     = 24    !< top south east     x+,y-,z+
  integer,parameter :: qTNW     = 25    !< top north west     x-,y+,z+
  integer,parameter :: qTNE     = 26    !< top north east     x+,y+,z+

  ! General direction index
  ! Each direction is uniquely represented by xyz
  ! where x,y,z can take one of the values: { -1, 0, 1 }
  ! here we use N to indicate -1, thus all direction can be writen as following
  !                                   X,  Y,  Z
  integer,parameter :: qN00 = 1   !< -1,  0,  0
  integer,parameter :: q0N0 = 2   !<  0, -1,  0
  integer,parameter :: q00N = 3   !<  0,  0, -1
  integer,parameter :: q100 = 4   !<  1,  0,  0
  integer,parameter :: q010 = 5   !<  0,  1,  0
  integer,parameter :: q001 = 6   !<  0,  0,  1

  integer,parameter :: q0NN = 7   !<  0, -1, -1
  integer,parameter :: q0N1 = 8   !<  0, -1,  1
  integer,parameter :: q01N = 9   !<  0,  1, -1
  integer,parameter :: q011 = 10  !<  0,  1,  1
  integer,parameter :: qN0N = 11  !< -1,  0, -1
  integer,parameter :: q10N = 12  !<  1,  0, -1
  integer,parameter :: qN01 = 13  !< -1,  0,  1
  integer,parameter :: q101 = 14  !<  1,  0,  1
  integer,parameter :: qNN0 = 15  !< -1, -1,  0
  integer,parameter :: qN10 = 16  !< -1,  1,  0
  integer,parameter :: q1N0 = 17  !<  1, -1,  0
  integer,parameter :: q110 = 18  !<  1,  1,  0

  integer,parameter :: qNNN = 19  !< -1, -1, -1
  integer,parameter :: qNN1 = 20  !< -1, -1,  1
  integer,parameter :: qN1N = 21  !< -1,  1, -1
  integer,parameter :: qN11 = 22  !< -1,  1,  1
  integer,parameter :: q1NN = 23  !<  1, -1, -1
  integer,parameter :: q1N1 = 24  !<  1, -1,  1
  integer,parameter :: q11N = 25  !<  1,  1, -1
  integer,parameter :: q111 = 26  !<  1,  1,  1

  integer,parameter :: q000 = 27  !< rest density is last

  integer, dimension(qQQQ),parameter  :: qInvDir = (/                          &
    &     q__E, q__N, q__T, q__W, q__S, q__B,                                  &
    &     q_TN, q_BN, q_TS, q_BS,                                              &
    &     q_TE, q_TW, q_BE, q_BW,                                              &
    &     q_NE, q_SE, q_NW, q_SW,                                              &
    &     qTNE, qTNW, qTSE, qTSW,                                              &
    &     qBNE, qBNW, qBSE, qBSW                    /)

  ! The reversed inverted direction map :)
  integer, dimension(qQQQ),parameter  :: qDir = (/                             &
    &     q__W, q__S, q__B, q__E, q__N, q__T,                                  &
    &     q_BS, q_TS, q_BN, q_TN,                                              &
    &     q_BW, q_BE, q_TW, q_TE,                                              &
    &     q_SW, q_NW, q_SE, q_NE,                                              &
    &     qBSW, qBSE, qBNW, qBNE,                                              &
    &     qTSW, qTSE, qTNW, qTNE                 /)

  !> mapping the first six neihbor directions to x, y, z (1,2,3)
  integer, parameter :: qAxis(6) = [1,2,3,1,2,3]

  !> Spatial offset according to the directions above for q__W, ...
  integer, dimension(qQQQ,3),parameter  :: qOffset =                           &
  reshape((/-1, 0, 0, 1, 0, 0, 0, 0, 0, 0,-1,+1,-1,+1,-1,-1,+1,+1,-1,+1,-1,+1,-1,+1,-1,+1,      &
    &        0,-1, 0, 0, 1, 0,-1,-1,+1,+1, 0, 0, 0, 0,-1,+1,-1,+1,-1,-1,+1,+1,-1,-1,+1,+1,      &
    &        0, 0,-1, 0, 0, 1,-1,+1,-1,+1,-1,-1,+1,+1, 0, 0, 0, 0,-1,-1,-1,-1,+1,+1,+1,+1/),(/qQQQ,3/))

  !> direction names of q__W, ... very useful for debuggin
  character(len=3), parameter           :: qDirName(qQQQ) = [                  &
     &        '  W', '  S', '  B', '  E', '  N', '  T',                        &
     &        ' BS', ' TS', ' BN', ' TN', ' BW', ' BE',                        &
     &        ' TW', ' TE', ' SW', ' NW', ' SE', ' NE',                        &
     &        'BSW', 'BSE', 'BNW', 'BNE',                                      &
     &        'TSW', 'TSE', 'TNW', 'TNE'                                       &
     &        ]

  !> Offset character bit to encode qOffset
  !! offset_bit = achar((qOffset(1)+1) + (qOffset(2)+1)*4 + (qOffset(3)+1)*16)
  !! Bit can be converted back to qOffset using
  !! qOffset(1) = mod(ichar(offset_bit),4) - 1
  !! qOffset(2) = mod(ichar(offset_bit),16)/4 - 1
  !! qOffset(3) = ichar(offset_bit)/16 - 1
  character, dimension(q000), parameter :: qOffset_inChar &
    & = (/ achar(20), achar(17), achar(5 ), &
    &      achar(22), achar(25), achar(37), &
    &      achar(1 ), achar(33), achar(9 ), &
    &      achar(41), achar(4 ), achar(6 ), &
    &      achar(36), achar(38), achar(16), &
    &      achar(24), achar(18), achar(26), &
    &      achar(0 ), achar(2 ), achar(8 ), &
    &      achar(10), achar(32), achar(34), &
    &      achar(40), achar(42), achar(21) /)

! -----------------------------------------------------------------------!
  ! Children conventions
  ! Numbering of children, starting in the lower, left, bottom corner.
  !> Spatial position of the child relative to the parent barycenter
  !! in integer coordinates
  integer, dimension(8,3),parameter  :: childPosition =                  &
                                 reshape((/-1,+1,-1,+1,-1,+1,-1,+1,      &
                                   &       -1,-1,+1,+1,-1,-1,+1,+1,      &
                                   &       -1,-1,-1,-1, 1, 1, 1, 1/),(/8,3/))


  !> Coordinates of children starting in the lower, left, bottom corner.
  !! This ordering follows the Z-curve
  integer, dimension(8,3), parameter :: childCoord = &
    &            reshape((/0, 1, 0, 1, 0, 1, 0, 1,   &
    &                      0, 0, 1, 1, 0, 0, 1, 1,   &
    &                      0, 0 ,0, 0, 1, 1, 1, 1/),(/8,3/))


  ! Transformation matrices for mapping the children
  integer, dimension(3,3,8),parameter  :: tem_rotationMatrix =             &
     reshape((/ 0,-1, 0,      &  ! 0 0 0
       &       -1, 0, 0,      &
       &        0, 0,-1,      &
       !        2nd entry in 3rd dimension
       &         1, 0, 0,     &  ! 1 0 0
       &         0,-1, 0,     &
       &         0, 0,-1,     &
       !        3rd entry in 3rd dimension
       &        -1, 0, 0,     &  ! 0 1 0
       &         0, 1, 0,     &
       &         0, 0,-1,     &
       !        4th entry in 3r!        4th entry in 3rd dimension
       &         0, 1, 0,     &!&         0, 0,-1,     &  ! 1 1 0
       &         1, 0, 0,     &!&         0, 1, 0,     &
       &         0, 0,-1,     &!&         1, 0, 0,     &
       !        5th entry in 3rd dimension
       &        -1, 0, 0,     &  ! 0 0 1
       &         0,-1, 0,     &
       &         0, 0, 1,     &
       !        6th entry in 3rd dimension
       &         0,-1, 0,     &!&         1, 0, 0,     &  ! 1 0 1
       &         1, 0, 0,     &!&         0, 0, 1,     &
       &         0, 0, 1,     &!&         0,-1, 0,     &
       !        7th entry in 3rd dimension
       &         0, 1, 0,     &  ! 0 1 1
       &        -1, 0, 0,     &
       &         0, 0, 1,     &
       !        8th entry  1 1 1
       &         1, 0, 0,     &  ! 1 1 1
       &         0, 1, 0,     &
       &         0, 0, 1      /),(/3,3,8/))

  ! MH: Why is this numbering scheme not according to the treelm general
  ! numbering?
  integer, dimension(7,3),parameter  :: parentNeighbor =                 &
                                 reshape((/ 1, 0, 0, 1, 0, 1, 1,         &
                                   &        0, 1, 0, 1, 1, 0, 1,         &
                                   &        0, 0, 1, 1, 1, 1, 0/),(/7,3/))
  ! map the apes numbering of the vertices to the one used by vtk
  ! element type 12 (VTK_HEXAHEDRON), 8 vertices
  ! our numbering:                  vtk numbering:
  !   3       4                       4       3
  !   ---------                       ---------
  !   |       |                       |       |
  !   |       |                       |       |
  !   ---------                       ---------
  !   1       2                       1       2
  !
  !                                                     0 1 2 3 4 5 6 7
  integer, dimension(8), parameter :: vtk8_vrtxMap = (/ 1,2,4,3,5,6,8,7 /)

  ! map the apes numbering of the vertices to the one used by vtk
  ! element type 25 (VTK_QUADRATIC_HEXAHEDRON), 20 vertices
  integer, dimension(20), parameter :: vtk20_vrtxMap =                         &
  !                        0 1 2 3 4 5 6 7 8  9 10 11 12 13 14 15 16 17 18 19
    &                   (/ 1,2,6,5,3,4,8,7,9,15,17,13,12,16,20,14,10,11,19,18 /)
!
! ------------------------------------------------------------------------------!
!
  ! State variable parameters
  integer, parameter :: prp_state     = 0         !< the mesoscopic state (pdfs)
  integer, parameter :: prp_density   = 1         !< fluid density
  integer, parameter :: prp_pressure  = 2         !< fluid pressure
  integer, parameter :: prp_wss       = 3         !< Wall Shear Stress
  integer, parameter :: prp_velMag    = 10        !< velocity magnitude
  integer, parameter :: prp_velX      = 11        !< velocity_X
  integer, parameter :: prp_velY      = 12        !< velocity_Y
  integer, parameter :: prp_velZ      = 13        !< velocity_Z


  ! Mathematical constants

  !> constant PI value
  real(kind=rk), parameter :: PI = 3.14159265358979323846_rk

  ! some other constants

  ! musubi related
  !> the speed of sound cs in various ways
  real(kind=rk), parameter :: cs       = 0.57735026918962576451_rk
  real(kind=rk), parameter :: csInv    = 1.73205080756887729352_rk
  real(kind=rk), parameter :: cs2inv   = 3._rk
  real(kind=rk), parameter :: cs2      = 1._rk/3._rk
  real(kind=rk), parameter :: cs4      = 1._rk/9._rk
  real(kind=rk), parameter :: t2cs4inv = 4.5_rk
  real(kind=rk), parameter :: t2cs2inv = 1.5_rk
  real(kind=rk), parameter :: cs4inv   = 9._rk
  real(kind=rk), parameter :: cs6inv   = 27._rk
  real(kind=rk), parameter :: cs8inv   = 81._rk
  real(kind=rk), parameter :: cs10inv  = 243._rk
  real(kind=rk), parameter :: cs12inv  = 729._rk
  real(kind=rk), parameter :: rho0     = 1._rk
  real(kind=rk), parameter :: rho0Inv  = 1._rk/rho0

  ! mathematical constants
  real(kind=rk), parameter :: sqrt2      = 1.4142135623731_rk
  real(kind=rk), parameter :: divSqrt2_2 = 0.707106781186548_rk
  real(kind=rk), parameter :: sqrt3      = 1.73205080756887729352_rk
  real(kind=rk), parameter :: two_sqrt2  = 2._rk * sqrt2

  real(kind=rk), parameter :: div1_2   = 1._rk /   2._rk
  real(kind=rk), parameter :: div1_3   = 1._rk /   3._rk
  real(kind=rk), parameter :: div1_4   = 1._rk /   4._rk
  real(kind=rk), parameter :: div1_6   = 1._rk /   6._rk
  real(kind=rk), parameter :: div1_7   = 1._rk /   7._rk
  real(kind=rk), parameter :: div1_8   = 1._rk /   8._rk
  real(kind=rk), parameter :: div1_9   = 1._rk /   9._rk
  real(kind=rk), parameter :: div1_12  = 1._rk /  12._rk
  real(kind=rk), parameter :: div1_16  = 1._rk /  16._rk
  real(kind=rk), parameter :: div1_18  = 1._rk /  18._rk
  real(kind=rk), parameter :: div1_21  = 1._rk /  21._rk
  real(kind=rk), parameter :: div1_24  = 1._rk /  24._rk
  real(kind=rk), parameter :: div1_27  = 1._rk /  27._rk
  real(kind=rk), parameter :: div1_36  = 1._rk /  36._rk
  real(kind=rk), parameter :: div1_42  = 1._rk /  42._rk
  real(kind=rk), parameter :: div1_48  = 1._rk /  48._rk
  real(kind=rk), parameter :: div1_54  = 1._rk /  54._rk
  real(kind=rk), parameter :: div1_72  = 1._rk /  72._rk
  real(kind=rk), parameter :: div1_108 = 1._rk / 108._rk
  real(kind=rk), parameter :: div1_216 = 1._rk / 216._rk

  real(kind=rk), parameter :: div2_3   = 2._rk /  3._rk
  real(kind=rk), parameter :: div2_8   = 2._rk /  8._rk
  real(kind=rk), parameter :: div2_9   = 2._rk /  9._rk
  real(kind=rk), parameter :: div2_27  = 2._rk / 27._rk

  real(kind=rk), parameter :: div4_3   = 4._rk /  3._rk
  real(kind=rk), parameter :: div4_9   = 4._rk /  9._rk
  real(kind=rk), parameter :: div4_21  = 4._rk / 21._rk
  real(kind=rk), parameter :: div4_27  = 4._rk / 27._rk
  
  real(kind=rk), parameter :: div3_2   = 3._rk /  2._rk
  real(kind=rk), parameter :: div3_4   = 3._rk /  4._rk
  real(kind=rk), parameter :: div3_4h  = 3._rk /  4.5_rk
  real(kind=rk), parameter :: div3_7   = 3._rk /  7._rk
  real(kind=rk), parameter :: div3_8   = 3._rk /  8._rk
  real(kind=rk), parameter :: div3_16  = 3._rk / 16._rk

  real(kind=rk), parameter :: div5_9   = 5._rk /  9._rk
  real(kind=rk), parameter :: div5_21  = 5._rk / 21._rk
  real(kind=rk), parameter :: div5_42  = 5._rk / 42._rk

  real(kind=rk), parameter :: div8_3   = 8._rk /   3._rk 
  real(kind=rk), parameter :: div8_7   = 8._rk /   7._rk
  real(kind=rk), parameter :: div8_27  = 8._rk /  27._rk
  real(kind=rk), parameter :: div9_16  = 9._rk /  16._rk

end module tem_param_module
! ****************************************************************************** !
