! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
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
!> author: Jens Zudrop
!! Module with routines to decide whether a faces is a compute face or not.
module tem_compteFaceRules_module
  ! include treelm modules
  use tem_aux_module,          only: tem_abort
  use tem_logging_module,      only: logUnit
  use tem_faceData_module,     only: tem_face_descriptor_type, &
    &                                tem_notExist_prp, tem_fluidFace_prp, &
    &                                tem_fromFinerFace_prp, &
    &                                tem_fromCoarserFace_prp, &
    &                                tem_remoteFace_prp, tem_bndFace_prp

  implicit none
  private

  public :: tem_isComputeFace

contains

! ****************************************************************************** !
  !> summary: Function to decide if a certain face is computed on the given rank
  !! and level or not.
  !!
  function tem_isComputeFace(facePos, faces, nEligibleChildren) result( isCompute )
    ! ---------------------------------------------------------------------------
    !> The position of the face you want to check.
    integer, intent(in) :: facePos
    !> The description of the faces.
    type(tem_face_descriptor_type), intent(in) :: faces
    !> Logical determine if given face is compute face or not.
    logical :: isCompute
    !> The number of eligible children for the vertical face dependency
    integer, intent(in) :: nEligibleChildren
    ! ---------------------------------------------------------------------------
    integer :: leftPrp, rightPrp
    ! ---------------------------------------------------------------------------

    ! In general we apply the following rule to determine which rank is
    ! computing a certain face:
    ! The owner (i.e. the rank) of the left element of a face computes the face.
    ! There are some exceptions for this rule, namely:
    ! If the left/right element is from finer, the finer level will compute
    ! the face.
    ! Furthermore, if the from finer is from remote, the remote process has to
    ! compute this face.

    leftPrp = faces%faceList%leftPrp%val(facePos)
    rightPrp = faces%faceList%rightPrp%val(facePos)

    ! So, start making our decision by considering all the cases explicitly.

    ! Fluid-fluid face (purely local)
    if( leftPrp.eq.tem_fluidFace_prp .and. rightPrp.eq.tem_fluidFace_prp ) then
      isCompute = .true.
    ! Fluid-ghostFromFiner face (purley local, no remote property set)
    elseif( leftPrp.eq.tem_fluidFace_prp .and.                                 &
      &     rightPrp.eq.tem_fromFinerFace_prp) then
      isCompute = .false.
      if(nEligibleChildren.eq.8) then
        isCompute = .true.
      end if
    ! GhostFromFiner-fluid face (purley local, no remote property set)
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_fluidFace_prp) then
      isCompute = .false.
      if(nEligibleChildren.eq.8) then
        isCompute = .true.
      end if
    ! GhostFromFiner-GhostFromFiner face (purley local, no remote property set)
    ! --> two ghost cells meet each other somehow. However, we do not have to
    !     compute these faces.
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_fromFinerFace_prp) then
      isCompute = .false.
    ! Fluid-ghostFromCoarser (purely local, no remote property set)
    elseif( leftPrp.eq.tem_fluidFace_prp .and.                                 &
      &     rightPrp.eq.tem_fromCoarserFace_prp) then
      isCompute = .true.
    ! GhostFromCoarser-fluid face (purely local, no remote property set)
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_fluidFace_prp) then
      isCompute = .true.
    ! GhostFromCoarser-ghostFromCoarser face (purely local, no remote property
    ! set)
    ! --> Two ghost from coarser elements meet each other somehow. However, we
    !     do not compute them.
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_fromCoarserFace_prp) then
      isCompute = .false.
    ! Fluid-remote face (the halo is not refined)
    elseif( leftPrp.eq.tem_fluidFace_prp .and.                                 &
      &     rightPrp.eq.tem_remoteFace_prp) then
      isCompute = .true.
    ! Remote-fluid face (the halo is not refined)
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_fluidFace_prp) then
      isCompute = .false.
    ! Remote-notExisting face (the remote is not refined)
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_notExist_prp) then
      isCompute = .false.
    ! notExisting-remote face (the remote is not refined)
    elseif( leftPrp.eq.tem_notExist_prp .and.                                  &
      &     rightPrp.eq.tem_remoteFace_prp) then
      isCompute = .false.
    ! notExisting-GhostFromFiner face (purely local, no remote property set)
    elseif( leftPrp.eq.tem_notExist_prp .and.                                  &
      &     rightPrp.eq.tem_fromFinerFace_prp) then
      isCompute = .false.
    ! GhostFromFiner-notExisting face (purely local, no remote property set)
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_notExist_prp) then
      isCompute = .false.
    ! Remote-remote face (sometimes, these faces can meet each other)
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_remoteFace_prp) then
      isCompute = .false.
    ! GhostFromCoarser-remote face (is computed on the remote rank, as
    ! my rank is living on the coarser level)
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_remoteFace_prp) then
      isCompute = .false.
    ! Remote-GhostFromCoarser face (is computed on the remote rank, as
    ! my rank is living on the coarser level)
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_fromCoarserFace_prp) then
      isCompute = .false.
    ! fluid-remote+fromCoaser face (the remote elements are actually coming
    ! from a higher level, so everything will be computed on my rank).
    elseif( leftPrp.eq.tem_fluidFace_prp .and.                                 &
      &     rightPrp.eq.ior(tem_remoteFace_prp, tem_fromCoarserFace_prp)) then
      isCompute = .true.
    ! remote+fromCoarser-fluid face (the remote elements are actually coming
    ! from a higher level, so everything will be computed on my rank).
    elseif( leftPrp.eq.ior(tem_remoteFace_prp,tem_fromCoarserFace_prp) .and.   &
      &     rightPrp.eq.tem_fluidFace_prp) then
      isCompute = .true.
    ! notExist-remote+fromCoaser face (nothing has to be done for this face,
    ! since no element of this face is located on my rank)
    elseif( leftPrp.eq.tem_notExist_prp .and.                                  &
      &     rightPrp.eq.ior(tem_remoteFace_prp, tem_fromCoarserFace_prp)) then
      isCompute = .false.
    ! remote+fromCoarser-remote face (nothing has to be done for this face,
    ! since no element of this face is located on my rank)
    elseif( leftPrp.eq.ior(tem_remoteFace_prp, tem_fromCoarserFace_prp) .and.  &
      &     rightPrp.eq.tem_notExist_prp) then
      isCompute = .false.
    ! remote+fromCoarser-remote+fromCoaser face (nothing has to be done for
    ! this face, since no element of this face is located on my rank)
    elseif( leftPrp.eq.ior(tem_remoteFace_prp, tem_fromCoarserFace_prp) .and.  &
      &     rightPrp.eq.ior(tem_remoteFace_prp, tem_fromCoarserFace_prp)) then
      isCompute = .false.
    ! fromFiner - notExist face
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_notExist_prp) then
      isCompute = .false.
    ! notExist - fromFiner face
    elseif( leftPrp.eq.tem_notExist_prp .and.                                  &
      &     rightPrp.eq.tem_fromFinerFace_prp) then
      isCompute = .false.
    ! fromCoarser - notExist face
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_notExist_prp) then
      isCompute = .false.
    ! notExist - fromCoarser face
    elseif( leftPrp.eq.tem_notExist_prp .and.                                  &
      &     rightPrp.eq.tem_fromCoarserFace_prp) then
      isCompute = .false.

    ! fluid-remote+fromFiner face (nothing is computed here, it will we be
    ! computed on the finer level remotely)
    elseif( leftPrp.eq.tem_fluidFace_prp .and.                                 &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp, tem_remoteFace_prp)) then
      isCompute = .false.
    ! remote+fromFiner-fluid face (nothing is computed here, it will we be
    ! computed on the finer level remotely)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp, tem_remoteFace_prp) .and.    &
      &     rightPrp.eq.tem_fluidFace_prp) then
      isCompute = .false.
    ! remote+fromFiner-remote face (everything is remote, so nothing to compute
    ! here)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp, tem_remoteFace_prp) .and.    &
      &     rightPrp.eq.tem_remoteFace_prp) then
      isCompute = .false.
    ! remote-remote+fromFiner face (everything is remote, so nothing to compute
    ! here)
    elseif( leftPrp.eq. tem_remoteFace_prp .and.                               &
      &     rightPrp.eq.ior(tem_remoteFace_prp, tem_fromFinerFace_prp)) then
      isCompute = .false.
    ! notExist-remote+fromFiner face (nothing is on my rank)
    elseif( leftPrp.eq. tem_notExist_prp .and.                                 &
      &     rightPrp.eq.ior(tem_remoteFace_prp, tem_fromFinerFace_prp)) then
      isCompute = .false.
    ! remote+fromFiner-notExist face (nothing is on my rank)
    elseif( leftPrp.eq.ior(tem_remoteFace_prp, tem_fromFinerFace_prp) .and.    &
      &     rightPrp.eq.tem_notExist_prp) then
      isCompute = .false.
    ! fromFiner-fromCoarser face (everything is local, this is an intermediate
    ! face level (i.e. leveljumps of difference larger than 1 occure here) )
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_fromCoarserFace_prp) then
      isCompute = .false.

    ! fromCoarser-fromFiner face (everything is local, this is an intermediate
    ! face level (i.e. leveljumps of difference larger than 1 occure here) )
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_fromFinerFace_prp) then
      isCompute = .false.

    ! fromFiner+remote - fromCoarser face (we do not have to recv info here, but
    ! we send it)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp) .and.     &
      &     rightPrp.eq.tem_fromCoarserFace_prp) then
      isCompute = .false.

    ! fromCoarser - remote+fromFiner face (we do not recv info here, but we send
    ! it)
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp)) then
      isCompute = .false.

    ! remote+fromFiner - remote+fromFiner face (everything is remote, nothing to
    ! be done)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp) .and.     &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp)) then
      isCompute = .false.

    ! remote+fromFiner - fromFiner face (everything is from finer, nothing to be
    ! done)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp) .and.     &
      &     rightPrp.eq.tem_fromFinerFace_prp) then
      isCompute = .false.

    ! fromFiner - remote+fromFiner face (everything is from finer, nothing to be
    ! done)
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp)) then
      isCompute = .false.

    ! fromFiner - remote face (nothing to be done)
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_remoteFace_prp) then
      isCompute = .false.

    ! remote - fromFiner face (nothing to be done)
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_fromFinerFace_prp) then
      isCompute = .false.

    ! remote+fromCoarser - fromCoarser face (elements are on coarser
    ! level, here and remote so nothing to be done)
    elseif( leftPrp.eq.ior(tem_fromCoarserFace_prp,tem_remoteFace_prp) .and.   &
      &     rightPrp.eq.tem_fromCoarserFace_prp) then
      isCompute = .false.

    ! fromCoarser - remote+fromCoarser face (elements are on coarser
    ! level, here and remote so nothing to be done)
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.ior(tem_fromCoarserFace_prp,tem_remoteFace_prp)) then
      isCompute = .false.

    ! fluid - boundary face
    elseif(leftPrp.eq.tem_fluidFace_prp .and. rightPrp.eq.tem_bndFace_prp) then
      isCompute = .true.

    ! boundary - fluid face
    elseif(leftPrp.eq.tem_bndFace_prp .and. rightPrp.eq.tem_fluidFace_prp) then
      isCompute = .true.

    ! remote - boundary face (nothing to be done)
    elseif(leftPrp.eq.tem_remoteFace_prp .and. rightPrp.eq.tem_bndFace_prp) then
      isCompute = .false.

    ! boundary - remote face (nothing to be done)
    elseif(leftPrp.eq.tem_bndFace_prp .and. rightPrp.eq.tem_remoteFace_prp) then
      isCompute = .false.

    ! fromFiner - boundary face (nothing to be done)
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_bndFace_prp) then
      isCompute = .false.

    ! boundary - fromFiner face (nothing to be done)
    elseif( leftPrp.eq.tem_bndFace_prp .and.                                   &
      &     rightPrp.eq.tem_fromFinerFace_prp) then
      isCompute = .false.

    ! fromFiner+remote - boundary face (nothing to be done)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp) .and.     &
      &     rightPrp.eq.tem_bndFace_prp) then
      isCompute = .false.

    ! boundary - fromFiner+remote face (nothing to be done)
    elseif( leftPrp.eq.tem_bndFace_prp .and.                                   &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp)) then
      isCompute = .false.

    ! fromCoarser - boundary face (nothing to be done)
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_bndFace_prp) then
      isCompute = .false.

    ! boundary - fromCoarser face (nothing to be done)
    elseif( leftPrp.eq.tem_bndFace_prp .and.                                   &
      &     rightPrp.eq.tem_fromCoarserFace_prp) then
      isCompute = .false.

    else
      write(logUnit(1),*) 'ERROR in tem_isComputeFace: not able to decide whether '//   &
        &        'face is computed locally or not, stopping'
      write(logUnit(1),*) 'left elem pos: ', faces%faceList%leftElemPos%val(facePos)
      write(logUnit(1),*) 'right elem pos: ', faces%faceList%rightElemPos%val(facePos)
      write(logUnit(1),*) 'left face prp: ', leftPrp
      write(logUnit(1),*) 'right face prp: ', rightPrp
      call tem_abort()
    end if

  end function tem_isComputeFace
! ****************************************************************************** !

end module tem_compteFaceRules_module
