! Copyright (c) 2011 Konstantin Kleinheinz <k.kleinheinz@grs-sim.de>
! Copyright (c) 2011 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011-2013, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
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
!> @todo: Add description
!!
module tem_intersection_module

  ! include treelm modules
  use env_module,          only: rk, long_k
  use tem_float_module,    only: operator(.feq.)
  use tem_param_module,    only: qOffset, q__W, q_BS, qBSW
  use treelmesh_module,    only: treelmesh_type
  use tem_topology_module, only: tem_CoordOfId
  use tem_geometry_module, only: tem_BaryOfId

  implicit none

  private

  public :: tem_line
  public :: tem_plane
  public :: tem_intersec
  public :: tem_intersec_elem
  public :: exit_element
  public :: tem_intersec_line_plane
  public :: tem_intersec_line_line
  public :: tem_intersec_ray_point

  type tem_line
    !> real world coordinates of the start point of the ray
    real(kind=rk) :: coordStart(3)
    !> direction vector of the line
    real(kind=rk) :: direction(3)
  end type tem_line

  type tem_plane
    !> real world coordinates of one point of the plane
    real(kind=rk) :: coord(3)
    !> normal vector of the plane
    real(kind=rk) :: normal(3)
  end type tem_plane

  type tem_intersec
    !> value of lamda (parameter of the line) of the intersection
    real(kind=rk) :: lambda
    !> real world coordinates of the intersection
    real(kind=rk) :: coord(3)
    !> distance between the center of the face (element) and the intersection
    real(kind=rk) :: distance(3)
  end type tem_intersec

  type tem_intersec_elem
    !>
    real(kind=rk) :: center(3)
    !>
    real(kind=rk) :: length
  end type tem_intersec_elem

  contains
! ****************************************************************************** !
  !> This subroutine checks at which face, edge or corner the line leaves the
  !> element and calculates the next element.
  subroutine exit_element( TreeID, line, tree )
    ! ---------------------------------------------------------------------------
    integer(kind=long_k),intent(in) :: TreeID
    type(tem_line) :: line
    type(treelmesh_type), intent(in) :: tree
    ! ---------------------------------------------------------------------------
    type(tem_plane) :: face
    type(tem_line) :: edge
    real(kind=rk) :: corner(3)
    type(tem_intersec) :: intersection
    logical :: intersects
    type(tem_intersec_elem) :: elem
    integer       :: coord(4)        ! spatial index triple for a given ID
    integer :: iDir
    ! ---------------------------------------------------------------------------


    ! @todo: calculate the size of the element
    ! Kann ich von hier auf diese funktion zugreifen? wie ist das bei einem
    ! gridrefinement?
    ! auch in zeile 108 in tem_geometry.f90
    coord = tem_CoordOfId(TreeID)
    elem%length = tree%global%BoundingCubeLength / real(2**coord(4), kind=rk)
    ! calculate the center of the element
    elem%center = tem_BaryOfId(tree, TreeID)

    ! loop over the faces of the element.
    ! 1 - positive x-direction, 2 - negative x-direction
    ! 3 - positive y-direction, 4 - negative y-direction
    ! 5 - positive z-direction, 6 - negative z-direction
    dirLoop: do iDir = 1, 26

      ! FACE
      if( iDir .eq. q__W )then!.or. q__E .or. q__N .or. q__S .or. q__T .or. q__B ) then
        ! check whether or not it is possible that this face is the exit face
        if( line%direction(1) .lt. 0.0_rk ) then
          ! we need the normal of the face
          face%normal = qOffset( iDir, : )

          !< calculate the center of the face
          face%coord = elem%center + 0.5_rk * elem%length * face%normal

          call tem_intersec_line_plane( face, line, intersects, intersection )
          !> found intersection, so exit loop
          if( intersects ) exit dirLoop
        end if

      ! EDGE
      elseif( iDir .eq. q_BS ) then
        ! check whether or not it is possible that this edge is the exit edge
        if ( line%direction(2) .lt. 0.0_rk .and.                               &
          &  line%direction(3) .lt. 0.0_rk ) then
          ! we need the direction and one point of the edge
          edge%direction = qOffset( iDir, : )
          edge%coordStart = elem%center + 0.5_rk * elem%length                 &
            &             * qOffset( iDir, : )

          call tem_intersec_line_line( edge, line, intersects, intersection )

          !< found intersection, so exit loop
          if( intersects ) exit dirLoop
        end if

      elseif( iDir .eq. qBSW ) then
        ! check whether or not it is possible that this corner is the
        ! exit corner
        if( line%direction(1) .lt. 0.0_rk .and.                                &
          & line%direction(2) .lt. 0.0_rk .and.                                &
          & line%direction(3) .lt. 0.0_rk ) then
          !< calculate the coordinates of the corner
          corner = elem%center + 0.5_rk * elem%length * qOffset( iDir, : )
          !< check whether or not the corner intersects with the ray
          call tem_intersec_ray_point( corner, line, intersects, intersection )

          if( intersects ) exit dirLoop

        end if
      end if
    enddo dirLoop

  end subroutine exit_element
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine calculates the intersection between a plane and a line.
  !! It gives back the coordinates of the intersection, the multiple of the
  !! direction vector of the intersection and the distance of the intersection
  !! to the center point of the plan.
  !!
  subroutine tem_intersec_line_plane( plane, line, intersects, intersection )
    ! ---------------------------------------------------------------------------
    type(tem_plane), intent(in) :: plane
    type(tem_line), intent(in) :: line
    type(tem_intersec), intent(out) :: intersection
    logical, intent(out) :: intersects
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: alignment, dist
    ! ---------------------------------------------------------------------------

    alignment = dot_product(plane%normal, line%direction)
    dist = dot_product( plane%normal, ( plane%coord - line%coordStart ))
    intersects = (alignment > epsilon(alignment))

    if (intersects) then
      intersection%lambda = dist / alignment
      intersection%coord = line%coordStart + intersection%lambda               &
        &                * line%direction
    else
      if (dist < 16*tiny(dist)) then
        ! Line is parallel to the plane, but on the plane
        intersects = .true.
        intersection%lambda = 0.0_rk
      else
        ! Line is parallel to the plane, no intersection!
        intersection%lambda = huge(intersection%lambda)
      end if
      intersection%coord = line%coordStart
    end if
    intersection%distance = plane%coord - intersection%coord

  end subroutine tem_intersec_line_plane
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine calculates the intersection between a line and a line.
  !! It gives back the coordinates of the intersection, the multiple of the
  !! direction vector of the intersection and the distance of the intersection
  !! to the center point of the line.
  !!
  subroutine tem_intersec_line_line( edge, line, intersects, intersection )
    ! ---------------------------------------------------------------------------
    type(tem_line), intent(in) :: edge
    type(tem_line), intent(in) :: line
    type(tem_intersec), intent(out) :: intersection
    logical, intent(out) :: intersects
    ! ---------------------------------------------------------------------------
    real(kind=rk), dimension(3) :: diff_vector, normal
    real(kind=rk), dimension(3) :: enormal
    real(kind=rk) :: alignment
    real(kind=rk) :: dist_line, dist_edge
    ! ---------------------------------------------------------------------------

    ! check whether the two lines intersect
    ! They have to be in a common plane for an
    ! intersection, compute the normal of this
    ! plane
    normal(1) = edge%direction(2)*line%direction(3)                            &
      &       - edge%direction(3)*line%direction(2)
    normal(2) = edge%direction(3)*line%direction(1)                            &
      &       - edge%direction(1)*line%direction(3)
    normal(3) = edge%direction(1)*line%direction(2)                            &
      &       - edge%direction(2)*line%direction(1)

    alignment = normal(1)*normal(1)                                            &
      &       + normal(2)*normal(2)                                            &
      &       + normal(3)*normal(3)

    if ((alignment > 16*tiny(alignment))) then
      ! The lines are not colinear, they might intersect, compute
      ! the distance of parallel planes through the lines, if
      ! this is 0, they actually intersect.
      dist_line = dot_product(normal, line%coordStart)
      dist_edge = dot_product(normal, edge%coordStart)
      intersects = (abs(dist_line - dist_edge) < epsilon(dist_line))

      if (intersects) then
        ! They intersect, get the point of intersection
        diff_vector = edge%coordStart - line%coordStart
        enormal(1) = edge%direction(2)*normal(3)                               &
          &        - edge%direction(3)*normal(2)
        enormal(2) = edge%direction(3)*normal(1)                               &
          &        - edge%direction(1)*normal(3)
        enormal(3) = edge%direction(1)*normal(2)                               &
          &        - edge%direction(2)*normal(1)
        intersection%lambda = dot_product(diff_vector, enormal)                &
          &                 / dot_product(line%direction, enormal)
        intersection%coord = line%coordStart + intersection%lambda             &
          &                * line%direction
        intersection%distance = edge%coordStart - intersection%coord
      end if
    else
      ! Lines are colinear
      intersects = .false.
    end if

  end subroutine tem_intersec_line_line
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine checks whether a line intersects with a point
  !!
  subroutine tem_intersec_ray_point( point, line, intersects, intersection )
    ! ---------------------------------------------------------------------------
    type(tem_line) :: line
    real(kind=rk),dimension(3) :: point
    type(tem_intersec) :: intersection
    logical :: intersects
    ! ---------------------------------------------------------------------------
    real(kind=rk),dimension(3) :: test_lambda
    ! ---------------------------------------------------------------------------

    test_lambda = ( point - line%coordStart ) / line%direction
    if ( ( test_lambda(1) .feq. test_lambda(2) ) .and.                         &
      & ( test_lambda(1) .feq. test_lambda(3) ) ) then
      intersection%lambda = test_lambda(1)
      intersects = .true.
      intersection%coord = point
    end if

  end subroutine tem_intersec_ray_point
! ****************************************************************************** !


end module tem_intersection_module
! ****************************************************************************** !
