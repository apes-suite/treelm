! Copyright (c) 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
! ******************************************************************************
!> author: Kannan Masilamani
!! This module contains plane boundary definition and routines
! ******************************************************************************
module tem_plane_module
  use env_module,          only: rk
  use tem_math_module,     only: cross_product3D
  use tem_cube_module,     only: tem_cube_type
  use tem_triangle_module, only: tem_triangle_type,      &
    &                            tem_triangleCubeOverlap

  implicit none

  private

  public :: tem_plane_type
  public :: tem_planeCubeOverlap
  public :: tem_createPlane

  !> type contains intrinsic plane information
  type tem_plane_type
    !> origin of the plane in x,y,z coordinate system
    real(kind=rk) :: origin(3)
    !> two  vectors defining the plane length
    !! 1st index x,y,z coordinate, 2nd index vector number
    real(kind=rk) :: vec(3,2)
    !> Convert plane into two triangles
    type(tem_triangle_type) :: triangle(2)
    !> unit normal direction of a plane
    real(kind=rk) :: unitNormal(3)
  end type tem_plane_type

contains

  ! ***************************************************************************!
  !> This routine creates plane definition from given origin and two vectors
  !! \verbatim
  !! vecB____________
  !! /\              |
  !! |               |
  !! |               |
  !! |               |
  !! |               |
  !! |-------------->|
  !! origin      vecA
  !! \endverbatim
  subroutine tem_createPlane(me, origin, vecA, vecB)
    !--------------------------------------------------------------------------!
    type(tem_plane_type), intent(out) :: me !< plane to return
    real(kind=rk), intent(in) :: origin(3) !< plane origin
    real(kind=rk), intent(in) :: vecA(3) !< 1st vector
    real(kind=rk), intent(in) :: vecB(3) !< 2nd vector
    !--------------------------------------------------------------------------!
    integer :: iTri
    !--------------------------------------------------------------------------!
    me%origin = origin
    me%vec(:,1) = vecA
    me%vec(:,2) = vecB
    ! convert plane into triangle
    do iTri = 1,2
      me%triangle(iTri)%nodes(1:3,1) = me%origin
      me%triangle(iTri)%nodes(1:3,2) = me%origin + me%vec(:,iTri)
      me%triangle(iTri)%nodes(1:3,3) = me%origin + me%vec(:,1) + me%vec(:,2)
    end do

    ! normal direction = crossproduct(vecA, vecB)
    me%unitNormal = cross_product3D(vecA, vecB)
    me%unitNormal = me%unitNormal                       &
      & / sqrt(dot_product(me%unitNormal, me%unitNormal))

  end subroutine tem_createPlane
  ! ***************************************************************************!

  ! ***************************************************************************!
  !> This function checks for intersection of plane and cube by
  !! checking two triangles of plane with a cube
  function tem_planeCubeOverlap( plane, cube ) result(overlap)
    !--------------------------------------------------------------------------!
    type(tem_plane_type), intent(in) :: plane
    type(tem_cube_type), intent(in) :: cube
    logical :: overlap
    !--------------------------------------------------------------------------!
    integer :: iTri

    overlap = .false.
    do iTri = 1, 2
      overlap = overlap .or.                                        &
        &       tem_triangleCubeOverlap( plane%triangle(iTri), cube )
    end do

  end function tem_planeCubeOverlap
  ! ***************************************************************************!

end module tem_plane_module

!> \page plane Plane
!! Planes are defined by an origin and two vectors.
!! Internally in the code, planes are converted into two triangles.
!! \verbatim
!! vec2____________
!! /\              |
!! |               |
!! |               |
!! |               |
!! |               |
!! |-------------->|
!! origin      vec1
!! \endverbatim
!!
!!Defintion:
!!\verbatim
!!spaial_object={
!!     attribute = {
!!       kind = 'boundary',
!!       label = 'plane1',
!!     },
!!     geometry = {
!!       kind = 'canoND',
!!       object = {
!!         origin={2.0,-1.0,-1.0},
!!         vec={{0.0,2.0,0.0},
!!              {0.0,0.0,2.0}},
!!       }
!!     }
!!}
!!\endverbatim
!! Seeder file to generate the mesh with plane is given below:
!! \include testsuite/plane/seeder.lua
!! \n\n
!! The mesh generated with plane inside mesh:
!! \image html plane.png
!! \n
!! \image html plane_withedges.png
!! \n
!! Cut view of plane
!! \image html plane_cut.png
!! \n\n
!! Example lua file is available at \link testsuite/plane/seeder.lua
!! \example testsuite/plane/seeder.lua
