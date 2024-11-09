! Copyright (c) 2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
!> This module provides rotations of a child cell and
!! its relations to direct neighbors
!!
module tem_rotation_module

  ! include treelm modules
  use tem_param_module,    only: childPosition, tem_rotationMatrix

  implicit none

  private
  public :: tem_assign_rotatedNeighbors
  public :: tem_assign_rotatedNeighbors_ref2d

  contains

! ****************************************************************************** !
  !> Sets the source coordinate for the rotated unit cube
  !!
  function tem_assign_rotatedNeighbors( iCase ) result( rotCoord )
    ! ---------------------------------------------------------------------------
    !>
    integer, intent(in)    :: iCase
    !> rotated coordinates
    integer   :: rotCoord(3,8)
    ! ---------------------------------------------------------------------------
    integer   :: curCoord(3)
    integer   :: iNeigh
    ! ---------------------------------------------------------------------------
    do iNeigh = 1,8
      curCoord = childPosition( iNeigh, : )
      rotCoord( :, iNeigh ) = matmul(transpose(tem_rotationMatrix(:,:,iCase)), &
        &                            curCoord )
      rotCoord = (rotCoord+1)/2
    enddo
   end function tem_assign_rotatedNeighbors
! ****************************************************************************** !


! ****************************************************************************** !
  !> Sets the source coordinate for the rotated unit cube
  !!
  function tem_assign_rotatedNeighbors_ref2d( iCase ) result( sourceCoord )
    ! ---------------------------------------------------------------------------
    !>
    integer, intent(in)    :: iCase
    !> source coordinates
    integer  :: sourceCoord(3,4)
    ! ---------------------------------------------------------------------------

    sourceCoord = 0
    select case( iCase )
    case( 1)
      ! Child num 1: 0.75 0.75
      sourceCoord(:, 1) = (/ 1, 1, 0 /)
      sourceCoord(:, 2) = (/ 0, 1, 0 /)
      sourceCoord(:, 3) = (/ 1, 0, 0 /)
      sourceCoord(:, 4) = (/ 0, 0, 0 /)
    case( 2)
      ! Child num 2: 0.25 0.75
      sourceCoord(:, 1) = (/ 0, 1, 0 /)
      sourceCoord(:, 2) = (/ 0, 0, 0 /)
      sourceCoord(:, 3) = (/ 1, 1, 0 /)
      sourceCoord(:, 4) = (/ 1, 0, 0 /)
    case( 3)
      ! Child num 3: 0.75 0.25
      sourceCoord(:, 1) = (/ 1, 0, 0 /)
      sourceCoord(:, 2) = (/ 1, 1, 0 /)
      sourceCoord(:, 3) = (/ 0, 0, 0 /)
      sourceCoord(:, 4) = (/ 0, 1, 0 /)
    case( 4)
      ! Child num 4: 0.25 0.25
      sourceCoord(:, 1) = (/ 0, 0, 0 /)
      sourceCoord(:, 2) = (/ 1, 0, 0 /)
      sourceCoord(:, 3) = (/ 0, 1, 0 /)
      sourceCoord(:, 4) = (/ 1, 1, 0 /)
    case( 5)
      ! Child num 1: 0.75 0.75
      sourceCoord(:, 1) = (/ 1, 1, 1 /)
      sourceCoord(:, 2) = (/ 0, 1, 1 /)
      sourceCoord(:, 3) = (/ 1, 0, 1 /)
      sourceCoord(:, 4) = (/ 0, 0, 1 /)
    case( 6)
      ! Child num 2: 0.25 0.75
      sourceCoord(:, 1) = (/ 0, 1, 1 /)
      sourceCoord(:, 2) = (/ 0, 0, 1 /)
      sourceCoord(:, 3) = (/ 1, 1, 1 /)
      sourceCoord(:, 4) = (/ 1, 0, 1 /)
    case( 7)
      ! Child num 3: 0.75 0.25
      sourceCoord(:, 1) = (/ 1, 0, 1 /)
      sourceCoord(:, 2) = (/ 1, 1, 1 /)
      sourceCoord(:, 3) = (/ 0, 0, 1 /)
      sourceCoord(:, 4) = (/ 0, 1, 1 /)
    case( 8)
      ! Child num 4: 0.25 0.25
      sourceCoord(:, 1) = (/ 0, 0, 1 /)
      sourceCoord(:, 2) = (/ 1, 0, 1 /)
      sourceCoord(:, 3) = (/ 0, 1, 1 /)
      sourceCoord(:, 4) = (/ 1, 1, 1 /)
    end select

   end function tem_assign_rotatedNeighbors_ref2d
! ****************************************************************************** !


end module tem_rotation_module
! ****************************************************************************** !
