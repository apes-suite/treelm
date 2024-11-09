! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2016, 2018-2019 Harald Klimach <harald.klimach@uni-siegen.de>
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
!! Datatypes and routines for face descriptions.
!!
!! Faces complement elements in the mesh.
!! Each element is enclosed by two faces in each direction. (For cubical
!! elements we get 6 faces enclosing the element).
!! We organize the faces by looking at each direction individually.
!! Thus, each face can be identified by its normal direction and an ID.
!! With this approach it is possible to use the same identification numbers
!! as for the elements and also to utilize the level descriptor that is
!! available for the elements, by just using it for each direction separately.
!!
!! For each face we have two sides (left in negative axis direction and right
!! in positive axis direction): left < | > right.
!! Accordingly we have two adjacent elements.
!! See [[tem_face_module]] for some more details about the internals and
!! the construction of the face information.
!!
!! Faces can have properties attached to them, similar to the properties
!! we have for elements.
!! It is possible to combine the properties by adding the corresponding
!! properties, i.e.
!! [[tem_faceData_module:tem_fromFinerFace_prp]]
!! + [[tem_faceData_module:tem_remoteFace_prp]]
!! represents a face which is communicated, but at the same time a face from
!! a finer level (i.e. in this situation is is a halo which is refined on
!! the remote partition).
module tem_faceData_module
  use env_module,              only: long_k
  use tem_dyn_array_module,    only: dyn_longArray_type, dyn_intArray_type
  use tem_grow_array_module,   only: grw_intArray_type, grw_longArray_type
  use tem_comm_module,         only: tem_communication_type
  use tem_construction_module, only: tem_levelDesc_type

  implicit none
  private

  !> Parameter to index left (1)
  integer, parameter :: tem_left = 1

  !> Parameter to index right (2)
  integer, parameter :: tem_right = 2

  !> Inverse face mapping.
  !!
  !! This simply maps 1 to 2 and 2 to 1.
  !! So it returns left when right is given, and right when left is given.
  integer, dimension(2), parameter :: tem_invFace_map = [ tem_right, &
    &                                                     tem_left   ]

  !> Mapping from spatial direction to left and right face.
  !!
  !! Each cubic element has 6 faces. In many solvers the different
  !! spatial directions are treated similarly, but just rotated. Therefore,
  !! in this module, we have left and right faces of an element.
  !! The following array provides a mapping from spatial direction
  !! to left or right face.
  !! In 2D the situation is as follows:
  !!```
  !!                                      right face in y
  !!                           ------------------------------
  !!                          /             / q__N         /
  !!       y dir             /             /              /
  !!        /               /  q__W       /      q__E    /
  !!       /     left face / <---------elem---------->  /   right face in x
  !!      /               /             /              /
  !!     /               /             /              /
  !!    /               /             / q__S         /
  !!   /               ------------------------------
  !!                            left face in y
  !!  ------------> x dir
  !!```
  !! In a similar way, we can generalize the upper situation to 3D. To convert
  !! direction a face direction, e.g.
  !! [[tem_param_module:q__e]] q__E, to a face direction access the
  !! following array by tem_dirToFace_map(q__E).
  integer, dimension(6), parameter :: tem_dirToFace_map =       &
    & [ tem_left, tem_left, tem_left,   &
    &   tem_right, tem_right, tem_right ]


  !> Left/right property of a face if the left/right element does not exist
  integer, parameter :: tem_notExist_prp = 1

  !> Left/right property of a face if the left/right element of the face is a
  !! fluid element.
  integer, parameter :: tem_fluidFace_prp = 2

  !> Left/right property of a face if the left/right element of the face is a
  !! from finer element.
  integer, parameter :: tem_fromFinerFace_prp = 4

  !> Left/right property of a face if the left/right element of the face is a
  !! from coarser element.
  integer, parameter :: tem_fromCoarserFace_prp = 8

  !> Left/right property of a face if the left/right element of the face is a
  !! halo element.
  integer, parameter :: tem_remoteFace_prp = 16

  !> Left/right property of a face if the left/right element of the face is a
  !! bnd element.
  integer, parameter :: tem_bndFace_prp = 32

  !> A unique list (i.e. guaranteees that no duplicates occur) collecting a
  !! set of faces in one direction.
  !!
  !! Each face is uniquely identified by its faceID and has two sides:
  !! left and right with corresponding neighbors and there may be different
  !! properties attached to each side.
  !!
  !! The number of faces in this descriptor is found in `faceID%nVals`.
  type tem_faceList_type
    !> The unique identifiers of the faces, i.e. the treeID of the elements
    !! on the left side of the face.
    !!
    !! treeID(leftelement) < face > right
    !! faceID(face) = treeID(leftelement)
    type(dyn_longArray_type) :: faceId

    !> The element id of the element on the right side of the face.
    type(grw_longArray_type) :: rightElemId
    !> Index of the element left of the face.
    type(grw_intArray_type) :: leftElemPos
    !> Index of the element right of the face.
    type(grw_intArray_type) :: rightElemPos
    !> Properties on the left side of the face.
    type(grw_intArray_type) :: leftPrp
    !> Properties on the right side of the face.
    type(grw_intArray_type) :: rightPrp
  end type tem_faceList_type

  !> Datatype to represent dependencies (in up- and downward direction)
  !! of the faces in a TREELM mesh. We only store face dependencies
  !! between the two adjacent refinement levels (i.e. vertical dependencies
  !! consider only a refinement level jump of 1 in both refinement directions).
  type tem_faceDep_type
    !> @todo JZ: currently we store vertical dependencies for all faces in the
    !! face,even if a face does not require a vertical dependency. Maybe we can
    !! replace this by a dynamic strucutre to avoid exessive memory consumption.

    !> The upwards (current level -> coarser level) facial dependencies. The
    !! size of this array is the number of faces. The array holds the face ids
    !! of the parent elements.
    !! If any element has no parent face (i.e. no upward dependency is
    !! necessary) this array holds a -1 at the element's position.
    integer(kind=long_k), allocatable :: parentFaceId(:)
    !> The parent elements position in the face description of the next coarser
    !! level. If any element does not have a parent face this array will hold a
    !! -1 entry at its position.
    integer, allocatable :: parentFacePos(:)
    !> The downwards (current level -> finer level) facial dependencies. The
    !! size of this array is 4 (as we have 4 child faces) and the number of
    !! faces on the current level. If a face does not have/require child faces
    !! its entries are set to -1.
    integer(kind=long_k), allocatable :: childFaceId(:,:)
    integer(kind=long_k), allocatable :: childFaceIdOp(:,:)
    !> The child element positions in the face desciption of the next finer
    !! level. If any element does not have child faces, its entries will be set
    !! to -1.
    integer, allocatable :: childFacePos(:,:)
    integer, allocatable :: childFacePosOp(:,:)
  end type tem_faceDep_type

  !> Iterator for a certain type of faces
  type tem_faceIterator_type
    !> Index of the element on the left side of the face.
    integer, allocatable :: leftPos(:)
    !> Index of the element on the right side of the face.
    integer, allocatable :: rightPos(:)
    !> Index of the face in the (overall) face description
    !! ([[tem_faceList_type]]).
    integer, allocatable :: facePos(:)
  end type tem_faceIterator_type

  !> Container to store the face that require interpolations.
  type tem_faceInterpolation_type
    !> The element position. Since interpolations are done
    !! facewise (left or right) we only need one element position here.
    integer, allocatable :: elemPos(:)
    !> The element positions of the finer child elements. As a face has
    !! exactly 4 children the first dimension is 4. The second dimension
    !! is the number of faces in this container.
    integer, allocatable :: childPos(:,:)
    !> The element position opposite to the refined face. E.g. if the
    !! left element is refined, the right face of the left element position
    !! is stored here.
    integer, allocatable :: elemPosOp(:)
    !> The child element positions opposite to the refined face.
    integer, allocatable :: childPosOp(:,:)
    !> Position of the face in the face description.
    integer, allocatable :: facePos(:)
  end type tem_faceInterpolation_type

  !> Datatype to provide face information for all faces
  !! in one direction on the same level.
  type tem_face_descriptor_type
    !> A dynamic list of faces in this descriptor.
    type(tem_faceList_type) :: faceList

    !> For each element in
    !! [[tem_faceData_module:tem_facelist_type]]
    !! we store the up- and downward facial dependencies.
    type(tem_faceDep_type) :: faceDep
    !> Datatype which holds iterators for the computable faces.
    type(tem_faceIterator_type) :: computeFace

    !> Datatype which holds iterators for the from finer faces.
    !! We separate left and right faces. To determine which direction
    !! corresponds to left or right face, have a look at this
    !! [[tem_faceData_module:tem_dirtoface_map]].
    type(tem_faceInterpolation_type) :: fromFinerFace(2)
    !> Buffer for state data which will be received before each compute step
    !! One buffer for the left faces of each cell and another one
    !! for the right faces of the cells. To determine which direction
    !! corresponds to left or right face, have a look at this
    !! [[tem_faceData_module:tem_dirtoface_map]].
    type(tem_communication_type) :: recvBuffer_state(2)
    !> Buffer for state data which will be send before each compute step.
    !! One buffer for the left faces of each cell and another one
    !! for the right faces of the cells. To determine which direction
    !! corresponds to left or right face, have a look at this
    !! [[tem_faceData_module:tem_dirtoface_map]].
    type(tem_communication_type) :: sendBuffer_state(2)
    !> Buffer for flux data which will be received before each compute step
    !! One buffer for the left faces of each cell and another one
    !! for the right faces of the cells. To determine which direction
    !! corresponds to left or right face, have a look at this
    !! [[tem_faceData_module:tem_dirtoface_map]].
    type(tem_communication_type) :: recvBuffer_flux(2)
    !> Buffer for flux data which will be send before each compute step.
    !! One buffer for the left faces of each cell and another one
    !! for the right faces of the cells. To determine which direction
    !! corresponds to left or right face, have a look at this
    !! [[tem_faceData_module:tem_dirtoface_map]].
    type(tem_communication_type) :: sendBuffer_flux(2)
  end type tem_face_descriptor_type

  !> Datatype for all faces in the mesh per level.
  type tem_face_type
    !> Face information: one descriptor for each direction (x,y,z).
    type(tem_face_descriptor_type) :: faces(3)

    !> Dimension-by-dimension level descriptors (these descriptors
    !! are necessary to build up the face descriptions).
    !! The first of this descriptor is build with (-1,0,0) and (+1,0,0)
    !! stencil for the x-direction.
    !! The second is build with (0,-1,0) and (0,+1,0)
    !! stencil for the y-direction.
    !! The third is build with (0,0,-1) and (0,0,+1)
    !! stencil for the z-direction.
    type(tem_levelDesc_type) :: dimByDimDesc(3)
  end type tem_face_type

  public :: tem_face_type, tem_faceIterator_type,                       &
    &       tem_face_descriptor_type, tem_faceList_type,                &
    &       tem_dirToFace_map,                                          &
    &       tem_left, tem_right, tem_invFace_map,       &
    &       tem_notExist_prp, tem_fluidFace_prp, tem_fromFinerFace_prp, &
    &       tem_fromCoarserFace_prp, tem_remoteFace_prp, tem_bndFace_prp

end module tem_faceData_module
