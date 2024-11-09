! Copyright (c) 2012, 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2012 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2014, 2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2013 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014 Rishabh Chandola <rishabh.chandola@student.uni-siegen.de>
! Copyright (c) 2014, 2016, 2019-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016-2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2018 Robin Weihe <robin.weihe@student.uni-siegen.de>
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
! **************************************************************************** !
!> summary: Module providing datatypes and routines regarding face information.
!! author: Jens Zudrop
!!
!! This module provides information about the existence and properties of
!! faces in an Octree based mesh. It is able to deliver info about the
!! horizontal (i.e. on the same refinement level) and vertical (i.e. between
!! the different levels) dependencies.
!!
!! If you want to make use of face descriptions in your code, you just have to
!! call [[tem_face_module:tem_build_face_info]] to buid up the face
!! descriptions in your solver. You have to make use of the created iterators
!! and face descriptions.
!! The `faces` structure created by [[tem_face_module:tem_build_face_info]]
!! is an array of [[tem_face_type]] with the face information for each level.
!! On each level there are three [[tem_face_descriptor_type]] for the three
!! spatial dimensions `faces(iLevel)%faces(1:3)`. Each of these descriptors
!! holds the complete list of faces on that level for this direction in
!! the `faceList`.
!!
!! In general, the construction of the face description follows the following
!! procedure:
!!
!! 1. Generate level descriptor in a dimension-by-dimension approach.
!! 2. Collect all the faces and attach properties to the faces.
!! 3. Extend the collected properties of the remote elements by from finer or
!!    from coarser properties (to indicate if a halo element is from coarser or
!!    from finer).
!! 4. Create vertical dependencies between the faces.
!! 5. Create list of faces which are relevant for computation and interpolation.
!! 6. Create communication buffers for send and receive of face information.
!!
!! Each face has two neighboring elements. To identify each face uniquely, we
!! use the tree ID of its left element.
!!```
!!       ---------------------------------------
!!      /                  /                  /
!!     /                  /                  /
!!    /   left elem.     /   right elem.    /
!!   /                  /                  /
!!  /                  /                  /
!!  --------------------------------------
!!```
!!
!! This module provides a basic description of the faces per spatial direction.
!! Building this description is achieved in several steps:
!!
!! 1. Collect all the faces (each level independently).
!! 2. Determine the properties of the faces.
!!    a. Check if halo elements are refined on the remote partition.
!! 3. Build the vertical dependencies between the faces of the different levels.
!! 4. Build buffers and lists of faces from the previous steps.
!!
!! To define a correct list of faces for the solvers, we attach left and right
!! properties to each face (in part 1 of the upper algorithm). We differentiate
!! between a left and right property of each face. The left one is determined
!! by it's left element and the right one by it's right element respectively.
!!```
!!              tem_left / / / tem_right
!!       -----------------------------------------------------
!!      /                      / / /                        /
!!     /                      / / /                        /
!!    /   left elem. --->    / / /    <--- right elem.    /
!!   /                      / / /                        /
!!  /                      / / /                        /
!!  ----------------------------------------------------
!!                       / / /
!!```
!!
!! We use one or more of the following properties to the left and right
!! properties of the faces:
!!
!! a. [[tem_faceData_module:tem_fluidface_prp]]
!! b. [[tem_faceData_module:tem_fromfinerface_prp]]
!! c. [[tem_faceData_module:tem_fromcoarserface_prp]]
!! d. [[tem_faceData_module:tem_remoteface_prp]]
!!
!! These properties are the basis for building the correct list of faces in step
!! 4 of the above algorithm.
!!
!! See [[Face Details]] for examples on face configurations.

module tem_face_module

  ! include treelm modules
  use env_module,              only: long_k
  use tem_dyn_array_module,    only: dyn_longArray_type, dyn_intArray_type, &
    &                                init, append, PositionOfVal
  use tem_grow_array_module,   only: grw_intArray_type, grw_longArray_type, &
    &                                init, append
  use tem_construction_module, only: tem_find_allElements,             &
    &                                tem_init_elemLevels,              &
    &                                tem_build_VerticalDependencies,   &
    &                                tem_build_HorizontalDependencies, &
    &                                tem_cleanupDependencyArrays,      &
    &                                tem_levelDesc_type, tem_treeIdInTotal
  use tem_param_module,        only: q__E, q__N, q__T, q__W, q__S, q__B, &
    &                                qOffset, qInvDir
  use tem_topology_module,     only: tem_directChildren, tem_ParentOf, &
    &                                tem_coordOfID, tem_IDofCoord
  use tem_logging_module,      only: logUnit
  use tem_aux_module,          only: tem_abort
  use tem_geometry_module,     only: tem_eligibleChildren
  use tem_stencil_module,      only: tem_stencilHeader_type, &
    &                                tem_stencil_dump
  use treelmesh_module,        only: treelmesh_type
  use tem_bc_prop_module,      only: tem_bc_prop_type
  use tem_comm_module,         only: tem_commpattern_type,tem_communication_type
  use tem_comm_env_module,     only: tem_comm_env_type
  use tem_faceData_module,     only: tem_face_type, tem_face_descriptor_type, &
    &                                tem_left, tem_right, &
    &                                tem_faceList_type,                   &
    &                                tem_notExist_prp, tem_fluidFace_prp, &
    &                                tem_fromFinerFace_prp,               &
    &                                tem_fromCoarserFace_prp,             &
    &                                tem_remoteFace_prp, tem_bndFace_prp
  use tem_compteFaceRules_module, only: tem_isComputeFace
  use tem_element_module,    only: eT_fluid, eT_halo,  &
    &                              eT_ghostFromCoarser, eT_ghostFromFiner
  use tem_spacetime_fun_module,  only: tem_st_fun_listElem_type

  implicit none

  private

  public :: tem_build_face_info


contains


  ! ************************************************************************** !
  !> summary: Routine to get face information for your whole mesh.
  !!
  !! This routine builds a description of the faces on each level. The
  !! resulting array of [[tem_face_type]] ranges from the minimal level
  !! of the mesh to the maximal level of the mesh.
  subroutine tem_build_face_info( tree, boundary, commPattern, proc, faces, &
                                & nEligibleChildren )
    ! --------------------------------------------------------------------------
    !> Tree representation of the mesh to buid the faces for.
    type(treelmesh_type), intent(inout) :: tree

    !> The boundaries of your simulation domain.
    type(tem_bc_prop_type), intent(in) :: boundary

    !> The communication pattern you use for the buffer.
    type(tem_commpattern_type), intent(in) :: commPattern

    !> Process description to use.
    type(tem_comm_env_type), intent(in) :: proc

    !> The number of eligible children for the vertical face dependency
    integer, intent(in) :: nEligibleChildren

    !> The created face descriptors (one for each level).
    type(tem_face_type), allocatable, intent(out) :: faces(:)
    ! --------------------------------------------------------------------------
    ! Level descriptor for each spatial direction and each level of your mesh.
    ! The level descriptor have to be constructed with the dimension by
    ! dimension stencils (+1, 0, -1) for each spatial direction.
    type(tem_levelDesc_type) :: levelDesc( 1:3, tree%global%minLevel  &
      &                                         :tree%global%maxLevel )
    type(tem_levelDesc_type), allocatable :: levelDescX(:)
    type(tem_levelDesc_type), allocatable :: levelDescY(:)
    type(tem_levelDesc_type), allocatable :: levelDescZ(:)
    integer :: iLevel
    ! --------------------------------------------------------------------------


    allocate(faces(tree%global%minLevel:tree%global%maxLevel))

    ! create the level desriptor in a dim-by-dim manner.
    write(logUnit(3),*)'Building dim-by-dim level descriptors for faces ...'
    call tem_dimByDim_construction( tree, boundary, commPattern, proc, &
      &                             levelDescX, levelDescY, levelDescZ )
    do iLevel = tree%global%minLevel,tree%global%maxLevel
      levelDesc(1,iLevel) = levelDescX(iLevel)
      levelDesc(2,iLevel) = levelDescY(iLevel)
      levelDesc(3,iLevel) = levelDescZ(iLevel)
    end do
    deallocate(levelDescX, levelDescY, levelDescZ)

    ! First we collect all the faces. A face is always identified by its left
    ! element tree id.
    ! We collect all the faces, even if only one element next to the face
    ! exists.
    ! Therefore, we run over fluids, ghostFromCoarser, ghostFromFiner and halo
    ! elements.
    ! We attach two-sided properties to the faces, too (depending on the element
    ! left and right to the face).
    call tem_collect_faces( levelDesc, tree%global%minLevel, &
      &                     tree%global%maxLevel, faces      )

    ! Now, we have to check all the faces which have the communication property
    ! again.
    ! We have to check if the halo element next to it is refined/from coarser on
    ! the remote rank. So, we extend the faces' properties.
    call tem_extend_remotePrp( levelDesc, tree%global%minLevel, &
      &                        tree%global%maxLevel, faces      )

    ! Now, we build the vertical dependency of the faces.
    call tem_faceDep_vertical( tree%global%minLevel, tree%global%maxLevel, &
      &                        faces, nEligibleChildren                    )

    !> \todo Remove all the unnecessary faces from the face descriptor

    ! Build the buffers for sent/received faces.
    call tem_build_faceBuffers( tree%global%minLevel, tree%global%maxLevel, &
      &                         faces, levelDesc                            )

    ! Build the iterators for the faces. These ones will be used by the solvers.
    call tem_build_faceLists(tree%global%minLevel, tree%global%maxLevel, &
      &                       nEligibleChildren, faces                   )

    !> \todo JZ: this is very memory consuming since we store them in different
    !! arrays several times.
    ! Buffer the dimension by dimension descriptor
    do iLevel = tree%global%minLevel,tree%global%maxLevel
      faces(iLevel)%dimByDimDesc(1) = levelDesc(1,iLevel)
      faces(iLevel)%dimByDimDesc(2) = levelDesc(2,iLevel)
      faces(iLevel)%dimByDimDesc(3) = levelDesc(3,iLevel)
    end do

  end subroutine tem_build_face_info
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> summary: Subroutine to build communication buffers for the faces.
  !!
  subroutine tem_build_faceBuffers( minLevel, maxLevel, faces, levelDesc )
    ! --------------------------------------------------------------------------
    !> The minimum refinement level of your mesh.
    integer, intent(in) :: minLevel
    !> The maximum refinement level of your mesh.
    integer, intent(in) :: maxLevel
    !> The communication pattern you want use for the buffer.
    ! type(tem_commpattern_type), intent(in) :: commPattern
    !> The created face descriptor.
    type(tem_face_type),intent(inout)  :: faces(minLevel:maxLevel)
    !> Dimension-by-dimension level descriptors
    type(tem_levelDesc_type), intent(in) :: levelDesc(1:3,minLevel:maxLevel)
    ! --------------------------------------------------------------------------
    integer :: iLevel, iDir
    ! --------------------------------------------------------------------------

    ! Create the face buffers for which I will receive information before
    ! my rank can do the timestep. Attention: this is consistent with the
    ! definition of a compute face.
    do iLevel = minLevel, maxLevel
      do iDir = 1, 3
        ! buffer for the state
        call tem_build_faceRecvBuffers(                                 &
          &    faces       = faces(iLevel)%faces(iDir),                 &
          &    levelDesc   = levelDesc(iDir, iLevel),                   &
          &    buf         = faces(iLevel)%faces(iDir)%recvBuffer_state )

        ! buffer for the flux
        call tem_build_faceRecvBuffers(                                &
          &    faces       = faces(iLevel)%faces(iDir),                &
          &    levelDesc   = levelDesc(iDir, iLevel),                  &
          &    buf         = faces(iLevel)%faces(iDir)%recvBuffer_flux )
      end do
    end do

    ! Create the face buffers for which I will send information before
    ! the other ranks can the timestep.
    do iLevel = minLevel, maxLevel
      do iDir = 1, 3
      ! buffer for the state
      call tem_build_faceSendBuffers(                                 &
        &    faces       = faces(iLevel)%faces(iDir),                 &
        &    levelDesc   = levelDesc(iDir, iLevel),                   &
        &    buf         = faces(iLevel)%faces(iDir)%sendBuffer_state )
      ! buffer for the flux
      call tem_build_faceSendBuffers(                                &
        &    faces       = faces(iLevel)%faces(iDir),                &
        &    levelDesc   = levelDesc(iDir, iLevel),                  &
        &    buf         = faces(iLevel)%faces(iDir)%sendBuffer_flux )
      end do
    end do

  end subroutine tem_build_faceBuffers
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> summary: Subroutine to build communication buffers for the faces the
  !! process will receive information for (before it can make the compute step).
  subroutine tem_build_faceRecvBuffers( levelDesc, faces, buf )
    ! -------------------------------------------------------------------------
    !> Dimension-by-dimension level descriptor for the current level and
    !! direction.
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> The communication pattern you want use for the buffer.
    ! type(tem_commpattern_type), intent(in) :: commPattern
    !> The created face descriptor.
    type(tem_face_descriptor_type),intent(inout)  :: faces
    !> The created receive buffer. Size is two, due to left and right limes of
    !! the face values. To access it use tem_left and tem_right.
    type(tem_communication_type), intent(out) :: buf(2)
    ! --------------------------------------------------------------------------
    integer :: iFace, iProc, iSide, elemPos, faceSide, tIdPos, rank, rankPos,  &
      &        elemAddPos, faceIndex
    logical :: wasAdded
    integer(kind=long_k) :: elemId
    ! The elements ids I will receive data for. We have to differentiate the
    ! left and right face of a face here.
    type(dyn_longArray_type) :: elemIds(2)
    ! The elements positions I will receive data for. We have to differentiate
    ! the left and right face of a face here.
    type(grw_intArray_type) :: elemPositions(2)
    ! The source ranks.
    type(dyn_intArray_type) :: srcRank(2)
    ! The position of the source ranks.
    type(grw_intArray_type) :: srcRankPos(2)
    ! --------------------------------------------------------------------------

    ! Init our growing arrays.
    call init( me = elemIds(1) , length = 16 )
    call init( me = elemIds(2) , length = 16 )
    call init( me = elemPositions(1) , length = 16 )
    call init( me = elemPositions(2) , length = 16 )
    call init( me = srcRankPos(1) , length = 16)
    call init( me = srcRankPos(2) , length = 16)
    call init( me = srcRank(1) , length = 16 )
    call init( me = srcRank(2) , length = 16 )

    ! Iterate over all the faces and count the communicate faces.
    do iFace = 1, faces%faceList%faceId%nVals

      ! Check whether we will receive info of one of the adjacent elements
      ! of this face.
      faceSide = tem_isRecvFace(iFace, faces)
      if (faceSide /= 0) then

        ! Check for which adjacent element (left or right) the current rank
        ! will receive information about and count it.
        if (faceSide == tem_left) then
          elemPos = faces%faceList%leftElemPos%val(iFace)
          elemId = faces%faceList%faceId%val(iFace)
          ! If we receive information about the left element, we receive
          ! information about it's right face!
          faceSide = tem_right
        else
          elemPos = faces%faceList%rightElemPos%val(iFace)
          elemId = faces%faceList%rightElemId%val(iFace)
          ! If we receive information about the right element, we receive
          ! information about it's left face!
          faceSide = tem_left
        end if

        ! Get the source rank. Therefore we have a look at the level
        ! descriptor of the current spatial direction and level.
        tIdPos = PositionOfVal( me = levelDesc%elem%tId, val = elemId )
        if (tIdPos <= 0) then
          write(*,*) 'ERROR in tem_build_faceRecvBuffers: not able to '//  &
            &        'identify source rank for a certain face, stopping ...'
          call tem_abort()
        else
          rank = levelDesc%elem%sourceProc%val(tIdPos)
        end if

        ! Store element position and source proc in a dynamic array.
        call append( me       = elemIds(faceSide), &
          &          val      = elemId,            &
          &          pos      = elemAddPos,        &
          &          wasAdded = wasAdded           )
        if (wasAdded) then
          call append( me = elemPositions(faceSide), val = elemPos )
          call append( me = srcRank(faceSide), val = rank, pos = rankPos )
          call append( me = srcRankPos(faceSide), val = rankPos )
        end if

      end if
    end do

    ! Now, we build the buffer for the faces (left and right values).
    do iSide = 1, 2
      ! Init the datastructures in the buffer
      buf(iSide)%nProcs = srcRank(iSide)%nVals
      allocate(buf(iSide)%proc( srcRank(iSide)%nVals ))
      allocate(buf(iSide)%nElemsProc( srcRank(iSide)%nVals ))
      allocate(buf(iSide)%elemPos( srcRank(iSide)%nVals ))
      do iProc = 1, buf(iSide)%nProcs
        call init( me = buf(iSide)%elemPos(iProc), length = 16 )
      end do

      ! Set the positions of the elements I will receive.
      do iFace = 1, elemPositions(iSide)%nVals
        faceIndex = elemIds(iSide)%sorted(iFace)
        elemPos = elemPositions(iSide)%val(faceIndex)
        rankPos = srcRankPos(iSide)%val(faceIndex)
        call append( me = buf(iSide)%elemPos(rankPos), val = elemPos )
      end do

      ! Count the number of elements we recv from each proc.
      do iProc = 1, buf(iSide)%nProcs
        buf(iSide)%proc(iProc) = srcRank(iSide)%val(iProc) - 1
        buf(iSide)%nElemsProc(iProc) = buf(iSide)%elemPos(iProc)%nVals
      end do
    end do

  end subroutine tem_build_faceRecvBuffers
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> summary: Subroutine to check if the current rank will receive information
  !! for the a given face before the compute step. Since face data
  !!
  pure function tem_isRecvFace(facePos, faces ) result( isRecvFace )
    ! --------------------------------------------------------------------------
    !> The position of the face in faces to be checked.
    integer , intent(in) :: facePos
    !> The description of the faces.
    type(tem_face_descriptor_type),intent(in)  :: faces
    !> Returns 0 if no information has to be received from remote ranks
    !! before a compute step can be done. In case that information has to be
    !! received it will return a non-zero positive value. Furhtermore,
    !! it determines whether information for the left face value
    !! (left adjacent element)
    !! or the right face value (right adjacent element) has to be received.
    !! In the first case it will return
    !! [[tem_faceData_module:tem_left]] and in the second case
    !! it will return
    !! [[tem_faceData_module:tem_right]].
    !! If something went wrong, e.g. an unknown combination occurs,
    !! it will return -1.
    integer :: isRecvFace
    ! --------------------------------------------------------------------------
    integer :: leftPrp, rightPrp
    ! --------------------------------------------------------------------------

    ! get the left and right property of the faces.
    leftPrp = faces%faceList%leftPrp%val(facePos)
    rightPrp = faces%faceList%rightPrp%val(facePos)

    isrecvface = -1 ! default value (error)

    ! So, check for all conditions under which face info has to be received.
    ! Fluid-fluid face (purely local)
    ! -> nothing to be sent
    if ( (leftPrp == tem_fluidFace_prp)        &
      &  .and. (rightPrp == tem_fluidFace_prp) ) then
      isRecvFace = 0

    ! Fluid-remote face (the halo is not refined)
    ! -> left element is on my rank, so I will compute this face. Therefore,
    !    I will also receive the remote information. In this case I will
    !    receive info from the right adjacent element.
    elseif ( (leftPrp == tem_fluidFace_prp)        &
      &     .and. (rightPrp == tem_remoteFace_prp) ) then
      isRecvFace = tem_right

    ! Remote-fluid face (the halo is not refined)
    ! -> left element is NOT on my rank, so another rank will compute this face.
    !    So, I have to send info for the left adjacent element, but I do not
    !    not receive information.
    elseif( (leftPrp == tem_remoteFace_prp)       &
      &     .and. (rightPrp == tem_fluidFace_prp) ) then
      isRecvFace = 0

    ! Remote-notExisting face (the remote is not refined)
    ! -> remote is not located on my rank and the other element does not
    !    exist. This face does not have to be computed on any rank, so my rank
    !    is not recv any info about any element adjacent to this face.
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_notExist_prp ) then
      isRecvFace = 0

    ! notExisting-remote face (the remote is not refined)
    ! -> remote is not located on my rank and the other element does not
    !    exist. This face does not have to be computed on any rank, so my rank
    !    is not recv any info about any element adjacent to this face.
    elseif( leftPrp.eq.tem_notExist_prp .and.                                  &
      &     rightPrp.eq.tem_remoteFace_prp ) then
      isRecvFace = 0

    ! Fluid-fromFiner face (purely local)
    ! -> everything is local, so there is no need to send or recv information
    !    from another rank.
    elseif( leftPrp.eq.tem_fluidFace_prp .and.                                 &
      &     rightPrp.eq.tem_fromFinerFace_prp ) then
      isRecvFace = 0

    ! FromFiner-fluid face (purely local)
    ! -> everything is local, so there is no need to send or recv information
    !    from another rank.
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_fluidFace_prp ) then
      isRecvFace = 0

    ! FromFiner-FromFiner face (purely local)
    ! -> everything is local. Somehow two ghost elements meet each other, but
    !    we do not have to recv or send any data for these faces.
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_fromFinerFace_prp ) then
      isRecvFace = 0

    ! FromCoarser-fluid face (purely local)
    ! -> everything is local, so no data has to be send or received.
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_fluidFace_prp ) then
      isRecvFace = 0

    ! Fluid-fromCoarser face (purely local)
    ! -> everything is local, so no data has to be send or received.
    elseif( leftPrp.eq.tem_fluidFace_prp .and.                                 &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isRecvFace = 0

    ! FromCoarser-fromCoarser face (purely local)
    ! -> everything is local, so no data has to be send or received.
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isRecvFace = 0

    ! Remote-fromCoarser face
    ! -> My rank holds the right elements. Therefore, my element is on a coarser
    !    level.
    !    However, the remote elements are refined, but not located on a single
    !    rank (otherwise it would have the notExist property instead of the
    !    remoteFace property). Therefore, my rank will send data to the remote
    !    ranks (as they are on a finer level), but we do not receive any info.
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isRecvFace = 0

    ! FormCoarser-Remote face
    ! -> My rank holds the left element. Therefore, my element is on a coarser
    !    level. However, the remote elements are refined, but not located on a
    !    single rank (otherwise it would have the notExist property instead of
    !    the remoteFace property). Therefore, my rank will send data to the
    !    remote ranks (as they are on a finer level), but we do not receive any
    !    info.
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_remoteFace_prp ) then
      isRecvFace = 0

    ! Remote-Remote face
    ! -> My rank holds none of the elements located next to this face. So we do
    !    not receive any information for this face.
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_remoteFace_prp ) then
      isRecvFace = 0

    ! notExist - remote face
    elseif( leftPrp.eq.tem_notExist_prp .and.                                  &
      &     rightPrp.eq.tem_remoteFace_prp) then
      isRecvFace = 0

    ! remote - notExist face
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_notExist_prp ) then
      isRecvFace = 0

    ! fluid-remote+fromCoaser face
    elseif( leftPrp.eq.tem_fluidFace_prp .and.                                 &
      &     rightPrp.eq.ior( tem_remoteFace_prp, tem_fromCoarserFace_prp )) then
      isRecvFace = tem_right

    ! remote+fromCoarser-fluid face
    elseif( leftPrp.eq.ior(tem_remoteFace_prp,tem_fromCoarserFace_prp) .and.   &
      &     rightPrp.eq.tem_fluidFace_prp ) then
      isRecvFace = tem_left

    ! notExist-remote+fromCoaser face (nothing to be received, this rank does
    ! not have any element shared by this face.)
    elseif( leftPrp.eq.tem_notExist_prp .and.                                  &
      &     rightPrp.eq.ior(tem_remoteFace_prp, tem_fromCoarserFace_prp)) then
      isRecvFace = 0

    ! remote+fromCoarser-remote face (nothing to be received, everything is
    ! remote)
    elseif( leftPrp.eq.ior(tem_remoteFace_prp, tem_fromCoarserFace_prp) .and.  &
      &     rightPrp.eq.tem_notExist_prp ) then
      isRecvFace = 0

    ! remote+fromCoarser-remote+fromCoaser face (nothing to be received,
    ! everything is remote)
    elseif( leftPrp.eq.ior(tem_remoteFace_prp, tem_fromCoarserFace_prp) .and.  &
      &     rightPrp.eq.ior(tem_remoteFace_prp, tem_fromCoarserFace_prp)) then
      isRecvFace = 0

    ! fromFiner - not exist face (nothing to be done)
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_notExist_prp ) then
      isRecvFace = 0

    ! notExist - fromFiner face (nothing to be done)
    elseif( leftPrp.eq.tem_notExist_prp .and.                                  &
      &     rightPrp.eq.tem_fromFinerFace_prp ) then
      isRecvFace = 0

    ! fromFiner+remote - fluid face (we do not have to recv info here, but we
    ! send it)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp) .and.     &
      &     rightPrp.eq.tem_fluidFace_prp ) then
      isRecvFace = 0

    ! fluid - remote+fromFiner face (we do not recv info here, but we send it)
    elseif( leftPrp.eq.tem_fluidFace_prp .and.                                 &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp)) then
      isRecvFace = 0

    ! fromCoarser - notExist face (nothing to receive)
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_notExist_prp ) then
      isRecvFace = 0

    ! notExist - fromCoarser face (nothing to receive)
    elseif( leftPrp.eq.tem_notExist_prp .and.                                  &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isRecvFace = 0

    ! remote+fromFiner-remote face (everything is remote)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp, tem_remoteFace_prp) .and.    &
      &     rightPrp.eq.tem_remoteFace_prp) then
      isRecvFace = 0

    ! remote-remote+fromFiner face (everything is remote)
    elseif( leftPrp.eq. tem_remoteFace_prp .and.                               &
      &     rightPrp.eq.ior(tem_remoteFace_prp, tem_fromFinerFace_prp)) then
      isRecvFace = 0

    ! notExist-remote+fromFiner face (nothing is on my rank)
    elseif( leftPrp.eq. tem_notExist_prp .and.                                 &
      &     rightPrp.eq.ior(tem_remoteFace_prp, tem_fromFinerFace_prp)) then
      isRecvFace = 0

    ! remote+fromFiner-notExist face (nothing is on my rank)
    elseif( leftPrp.eq.ior(tem_remoteFace_prp, tem_fromFinerFace_prp) .and.    &
      &     rightPrp.eq.tem_notExist_prp ) then
      isRecvFace = 0

    ! fromFiner-fromCoarser face (everything is local, this is an intermediate
    ! face level (i.e. leveljumps of difference larger than 1 occure here) )
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isRecvFace = 0

    ! fromCoarser-fromFiner face (everything is local, this is an intermediate
    ! face level (i.e. leveljumps of difference larger than 1 occure here) )
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_fromFinerFace_prp ) then
      isRecvFace = 0

    ! fromFiner+remote - fromCoarser face (we do not have to recv info here,
    ! but we send it)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp) .and.     &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isRecvFace = 0

    ! fromCoarser - remote+fromFiner face (we do not recv info here, but we
    ! send it)
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp)) then
      isRecvFace = 0

    ! remote+fromFiner - remote+fromFiner face (everything is remote, nothing
    ! to be done)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp) .and.     &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp)) then
      isRecvFace = 0

    ! remote+fromFiner - fromFiner face (everything is from finer, nothing to
    ! be done)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp) .and.     &
      &     rightPrp.eq.tem_fromFinerFace_prp ) then
      isRecvFace = 0

    ! fromFiner - remote+fromFiner face (everything is from finer, nothing to
    ! be done)
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp)) then
      isRecvFace = 0

    ! fromFiner - remote face (nothing to be done)
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_remoteFace_prp ) then
      isRecvFace = 0

    ! remote - fromFiner face (nothing to be done)
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_fromFinerFace_prp ) then
      isRecvFace = 0

    ! remote+fromCoarser - fromCoarser face (elements are on coarser
    ! level, here and remote so nothing to be done)
    elseif( leftPrp.eq.ior(tem_fromCoarserFace_prp,tem_remoteFace_prp) .and.   &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isRecvFace = 0

    ! fromCoarser - remote+fromCoarser face (elements are on coarser
    ! level, here and remote so nothing to be done)
    elseif( leftPrp.eq.tem_fromCoarserFace_prp  .and.                          &
      &     rightPrp.eq.ior(tem_fromCoarserFace_prp,tem_remoteFace_prp)) then
      isRecvFace = 0

    ! fluid - boundary face (nothing to be done)
    elseif(leftPrp.eq.tem_fluidFace_prp .and. rightPrp.eq.tem_bndFace_prp) then
      isRecvFace = 0

    ! boundary - fluid face (nothing to be done)
    elseif(leftPrp.eq.tem_bndFace_prp .and. rightPrp.eq.tem_fluidFace_prp) then
      isRecvFace = 0

    ! remote - boundary face (nothing to be done)
    elseif(leftPrp.eq.tem_remoteFace_prp .and. rightPrp.eq.tem_bndFace_prp) then
      isRecvFace = 0

    ! boundary - remote face (nothing to be done)
    elseif(leftPrp.eq.tem_bndFace_prp .and. rightPrp.eq.tem_remoteFace_prp) then
      isRecvFace = 0

    ! fromFiner - boundary face (nothing to be done)
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_bndFace_prp ) then
      isRecvFace = 0

    ! boundary - fromFiner face (nothing to be done)
    elseif( leftPrp.eq.tem_bndFace_prp .and.                                   &
      &     rightPrp.eq.tem_fromFinerFace_prp ) then
      isRecvFace = 0

    ! fromFiner+remote - boundary face (nothing to be done)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp) .and.     &
      &     rightPrp.eq.tem_bndFace_prp ) then
      isRecvFace = 0

    ! boundary - fromFiner+remote face (nothing to be done)
    elseif( leftPrp.eq.tem_bndFace_prp .and.                                   &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp)) then
      isRecvFace = 0

    ! fromCoarser - boundary face (nothing to be done)
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_bndFace_prp ) then
      isRecvFace = 0

    ! boundary - fromCoarser face (nothing to be done)
    elseif( leftPrp.eq.tem_bndFace_prp .and.                                   &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isRecvFace = 0

    end if

  end function tem_isRecvFace
! ****************************************************************************** !


  ! ************************************************************************** !
  !> summary: Subroutine to build communication buffers for the faces the
  !! current rank has to send information for (before the other ranks can make
  !! the compute step).
  subroutine tem_build_faceSendBuffers( levelDesc, faces, buf )
    ! --------------------------------------------------------------------------
    !> Dimension-by-dimension level descriptor for the current level and
    !! direction.
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> The communication pattern you want use for the buffer.
    ! type(tem_commpattern_type), intent(in) :: commPattern
    !> The created face descriptor.
    type(tem_face_descriptor_type),intent(inout)  :: faces
    !> The created send buffer. Size is two, due to left and right limes of
    !! the face values. To access it use tem_left and tem_right.
    type(tem_communication_type), intent(out) :: buf(2)
    ! --------------------------------------------------------------------------
    integer :: iFace, iProc, iSide, elemPos, faceSide, tIdPos, rank, rankPos,  &
      &        faceIndex, elemAddPos
    logical :: wasAdded
    integer(kind=long_k) :: elemId, trgElemId
    ! The element treeids I will send data for. We have to differentiate the
    ! left and right face of a face here.
    type(dyn_longArray_type) :: elemIds(2)
    ! The elements positions I will send data for. We have to differentiate the
    ! left and right face of a face here.
    type(grw_intArray_type) :: elemPositions(2)
    ! The source ranks.
    type(dyn_intArray_type) :: trgRank(2)
    ! The position of the source ranks.
    type(grw_intArray_type) :: trgRankPos(2)
    ! --------------------------------------------------------------------------

    ! Init our growing arrays.
    call init( me = elemIds(1) , length = 16 )
    call init( me = elemIds(2) , length = 16 )
    call init( me = elemPositions(1) , length = 16)
    call init( me = elemPositions(2) , length = 16)
    call init( me = trgRankPos(1) , length = 16)
    call init( me = trgRankPos(2) , length = 16)
    call init( me = trgRank(1) , length = 16 )
    call init( me = trgRank(2) , length = 16 )

    ! Iterate over all the faces and count the communicate faces.
    do iFace = 1, faces%faceList%faceId%nVals

      ! Check whether we will send info of one of the adjacent elements
      ! of this face.
      faceSide = tem_isSendFace(iFace, faces)
      if ( faceSide /= 0 ) then

        ! Check for which adjacent element (left or right) the current rank
        ! will send information about and count it.
        if (faceSide == tem_left) then
          ! The position and id of the elment I will send data for
          elemPos = faces%faceList%leftElemPos%val(iFace)
          elemId = faces%faceList%faceId%val(iFace)
          ! I will send information about my left element, so I will send
          ! the right face info of this element.
          faceSide = tem_right
          ! Get the element id in the opposite direction of the face. This
          ! will be used to determine the target rank.
          trgElemId = faces%faceList%rightElemId%val(iFace)
        else
          ! The position and id of the elment I will send data for
          elemPos = faces%faceList%rightElemPos%val(iFace)
          elemId = faces%faceList%rightElemId%val(iFace)
          ! I will send info about the right element, so I send info about
          ! it's left face.
          faceSide = tem_left
          ! Get the element id in the opposite direction of the face. This
          ! will be used to determine the target rank.
          trgElemId = faces%faceList%faceId%val(iFace)
        end if

        ! Get the source rank of the target element id. This will give us
        ! the target rank of the face.
        ! Therefore we have a look at the level
        ! descriptor of the current spatial direction and level.
        tIdPos = PositionOfVal( me = levelDesc%elem%tId, val = trgElemId )
        if (tIdPos <= 0) then
          write(*,*) 'ERROR in tem_build_faceSendBuffers: not able to '//  &
            &        'identify target rank for a certain face, stopping ...'
          call tem_abort()
        end if
        rank = levelDesc%elem%sourceProc%val(tIdPos)

        ! Store element position and soruce proc in a dynamic array.
        call append( me       = elemIds(faceSide), &
          &          val      = elemId,            &
          &          pos      = elemAddPos,        &
          &          wasAdded = wasAdded           )
        if(wasAdded) then
          call append( me = elemPositions(faceSide), val = elemPos )
          call append( me = trgRank(faceSide), val = rank, pos = rankPos )
          call append( me = trgRankPos(faceSide), val = rankPos )
        end if

      end if
    end do

    ! Now, we build the buffer for the faces (left and right values).
    do iSide = 1, 2
      ! Init the datastructures in the buffer
      buf(iSide)%nProcs = trgRank(iSide)%nVals
      allocate(buf(iSide)%proc( trgRank(iSide)%nVals ))
      allocate(buf(iSide)%nElemsProc( trgRank(iSide)%nVals ))
      allocate(buf(iSide)%elemPos( trgRank(iSide)%nVals ))
      do iProc = 1, buf(iSide)%nProcs
        call init( me = buf(iSide)%elemPos(iProc), length = 16 )
      end do

      ! Set the positions of the elements I will receive.
      do iFace = 1, elemPositions(iSide)%nVals
        faceIndex = elemIds(iSide)%sorted(iFace)
        elemPos = elemPositions(iSide)%val(faceIndex)
        rankPos = trgRankPos(iSide)%val(faceIndex)
        call append( me = buf(iSide)%elemPos(rankPos), val = elemPos )
      end do

      ! Count the number of elements we recv from each proc.
      do iProc = 1, buf(iSide)%nProcs
        buf(iSide)%proc(iProc) = trgRank(iSide)%val(iProc) - 1
        buf(iSide)%nElemsProc(iProc) = buf(iSide)%elemPos(iProc)%nVals
      end do
    end do

  end subroutine tem_build_faceSendBuffers
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> summary: Subroutine to check if the current rank will send information
  !! for the a given face before the compute step.
  !!
  pure function tem_isSendFace(facePos, faces) result( isSendFace )
    ! --------------------------------------------------------------------------
    !> The position of the face in faces to be checked.
    integer , intent(in) :: facePos
    !> The description of the faces.
    type(tem_face_descriptor_type),intent(in)  :: faces
    !> Returns 0 if no information has to be received from remote ranks
    !! before a compute step can be done. In case that information has to be
    !! received it will return a non-zero positive value. Furhtermore,
    !! it determines whether information for the left face value
    !! (left adjacent element)
    !! or the right face value (right adjacent element) has to be received.
    !! In the first case it will return
    !! [[tem_faceData_module:tem_left]] and in the second case
    !! it will return
    !! [[tem_faceData_module:tem_right]].
    !! If something went wrong, e.g. an unknown combination occurs,
    !! it will return -1.
    integer :: isSendFace
    ! --------------------------------------------------------------------------
    integer :: leftPrp, rightPrp
    ! --------------------------------------------------------------------------

    leftPrp = faces%faceList%leftPrp%val(facePos)
    rightPrp = faces%faceList%rightPrp%val(facePos)

    issendface = -1 ! default value (error)

    ! So, check for all conditions under which face info has to be send.

    ! Fluid-fluid face (purely local)
    ! -> nothing to send
    if( leftPrp.eq.tem_fluidFace_prp .and. rightPrp.eq.tem_fluidFace_prp ) then
      isSendFace = 0

    ! Fluid-remote face (the halo is not refined)
    ! -> left element is on my rank, so I will compute it. Therefore, I will
    !    receive information about the right element, but I do not have to send
    !    anything.
    elseif( leftPrp.eq.tem_fluidFace_prp .and.                                 &
      &     rightPrp.eq.tem_remoteFace_prp ) then
      isSendFace = 0

    ! Remote-fluid face (the halo is not refined)
    ! -> left element is on a remote rank, so I will not compute it. Therfore
    !    I will send the element to the remote rank.
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_fluidFace_prp ) then
      isSendFace = tem_right

    ! Remote-notExisting face (the remote is not refined)
    ! -> remote is not located on my rank and the other element does not
    !    exist. This face does not have to be computed on any rank, so my rank
    !    is not sending any info about any element adjacent to this face.
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_notExist_prp ) then
      isSendFace = 0

    ! notExisting-remote face (the remote is not refined)
    ! -> remote is not located on my rank and the other element does not
    !    exist. This face does not have to be computed on any rank, so my rank
    !    is not send any info about any element adjacent to this face.
    elseif( leftPrp.eq.tem_notExist_prp .and.                                  &
      &     rightPrp.eq.tem_remoteFace_prp ) then
      isSendFace = 0

    ! Fluid-fromFiner face (purely local)
    ! -> everything is local, so there is no need to send or recv information
    !    from another rank.
    elseif( leftPrp.eq.tem_fluidFace_prp .and.                                 &
      &     rightPrp.eq.tem_fromFinerFace_prp ) then
      isSendFace = 0

    ! FromFiner-fluid face (purely local)
    ! -> everything is local, so there is no need to send or recv information
    !    from another rank.
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_fluidFace_prp ) then
      isSendFace = 0

    ! FromFiner-FromFiner face (purely local)
    ! -> everything is local. Somehow two ghost elements meet each other, but
    !    we do not have to recv or send any data for these faces.
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_fromFinerFace_prp ) then
      isSendFace = 0

    ! FromCoarser-fluid face (purely local)
    ! -> everything is local, so no data has to be send or received.
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_fluidFace_prp ) then
      isSendFace = 0

    ! Fluid-fromCoarser face (purely local)
    ! -> everything is local, so no data has to be send or received.
    elseif( leftPrp.eq.tem_fluidFace_prp .and.                                 &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isSendFace = 0

    ! FromCoarser-fromCoarser face (purely local)
    ! -> everything is local, so no data has to be send or received.
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isSendFace = 0

    ! Remote-fromCoarser face
    ! -> The remote fluid elements are located on a finer level
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isSendFace = tem_right

    ! FromCoarser-remote face
    ! -> The remote fluid elements are located on a finer level
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_remoteFace_prp ) then
      isSendFace = tem_left

    ! Remote-remote face
    ! -> The halo elements meet each other. This can happen, but we do not have
    !    to send or receive information for this face, because we do not
    !    compute anything here.
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_remoteFace_prp ) then
      isSendFace = 0

    ! notExist - remote face (nothing to be send, we do not have info about
    ! adjacent elements)
    elseif( leftPrp.eq.tem_notExist_prp .and.                                  &
      &     rightPrp.eq.tem_remoteFace_prp ) then
      isSendFace = 0

    ! remote - notExist face (nothing to be send, we do not have any info about
    ! adjacent elements)
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_notExist_prp ) then
      isSendFace = 0

    ! fluid-remote+fromCoaser face (here, we receive data in some situations,
    ! but we never send info in this case)
    elseif( leftPrp.eq.tem_fluidFace_prp .and.                                 &
      &     rightPrp.eq.ior(tem_remoteFace_prp, tem_fromCoarserFace_prp)) then
      isSendFace = 0

    ! remote+fromCoarser-fluid face (here, we receive data in some situations,
    ! but we never send info for this faces)
    elseif( leftPrp.eq.ior(tem_remoteFace_prp,tem_fromCoarserFace_prp) .and.   &
      &     rightPrp.eq.tem_fluidFace_prp ) then
      isSendFace = 0

    ! notExist-remote+fromCoaser face (nothing to be send here)
    elseif( leftPrp.eq.tem_notExist_prp .and.                                  &
      &     rightPrp.eq.ior(tem_remoteFace_prp, tem_fromCoarserFace_prp)) then
      isSendFace = 0

    ! remote+fromCoarser-remote face (everything is remote, so nothing to be
    ! sent)
    elseif( leftPrp.eq.ior(tem_remoteFace_prp, tem_fromCoarserFace_prp) .and.  &
      &     rightPrp.eq.tem_notExist_prp ) then
      isSendFace = 0

    ! remote+fromCoarser-remote+fromCoaser face (everything is remote, so
    ! nothing to be sent)
    elseif( leftPrp.eq.ior(tem_remoteFace_prp, tem_fromCoarserFace_prp) .and.  &
      &     rightPrp.eq.ior(tem_remoteFace_prp, tem_fromCoarserFace_prp)) then
      isSendFace = 0

    ! fromFiner - not exist face (nothing to be send here).
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_notExist_prp ) then
      isSendFace = 0

    ! notExist - fromFiner face (nothing to be send here)
    elseif( leftPrp.eq.tem_notExist_prp .and.                                  &
      &     rightPrp.eq.tem_fromFinerFace_prp ) then
      isSendFace = 0

    ! fromFiner+remote - fluid face (will be send on the finer level)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp) .and.     &
      &     rightPrp.eq.tem_fluidFace_prp ) then
      isSendFace = 0

    ! fluid - remote+fromFiner face (will be send on the finer level)
    elseif( leftPrp.eq.tem_fluidFace_prp .and.                                 &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp)) then
      isSendFace = 0

    ! fromCoarser - notExist face (nothing to send)
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_notExist_prp ) then
      isSendFace = 0

    ! notExist - fromCoarser face (nothing to send)
    elseif( leftPrp.eq.tem_notExist_prp .and.                                  &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isSendFace = 0

    ! remote+fromFiner-remote face (everything is remote)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp, tem_remoteFace_prp) .and.    &
      &     rightPrp.eq.tem_remoteFace_prp ) then
      isSendFace = 0

    ! remote-remote+fromFiner face (everything is remote)
    elseif( leftPrp.eq. tem_remoteFace_prp .and.                               &
      &     rightPrp.eq.ior(tem_remoteFace_prp, tem_fromFinerFace_prp)) then
      isSendFace = 0

    ! notExist-remote+fromFiner face (nothing is on my rank)
    elseif( leftPrp.eq. tem_notExist_prp .and.                                 &
      &     rightPrp.eq.ior(tem_remoteFace_prp, tem_fromFinerFace_prp)) then
      isSendFace = 0

    ! remote+fromFiner-notExist face (nothing is on my rank)
    elseif( leftPrp.eq.ior(tem_remoteFace_prp, tem_fromFinerFace_prp) .and.    &
      &     rightPrp.eq.tem_notExist_prp ) then
      isSendFace = 0

    ! fromFiner-fromCoarser face (everything is local, this is an intermediate
    ! face level (i.e. leveljumps of difference larger than 1 occure here) )
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isSendFace = 0

    ! fromCoarser-fromFiner face (everything is local, this is an intermediate
    ! face level (i.e. leveljumps of difference larger than 1 occure here) )
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_fromFinerFace_prp ) then
      isSendFace = 0

    ! fromFiner+remote - fromCoarser face (send on a finer level)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp) .and.     &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isSendFace = 0

    ! fromCoarser - remote+fromFiner face (send on a finer level)
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp)) then
      isSendFace = 0

    ! remote+fromFiner - remote+fromFiner face (everything is remote, nothing
    ! to be done)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp) .and.     &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp)) then
      isSendFace = 0

    ! remote+fromFiner - fromFiner face (everything is from finer, nothing to
    ! be done)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp) .and.     &
      &     rightPrp.eq.tem_fromFinerFace_prp ) then
      isSendFace = 0

    ! fromFiner - remote+fromFiner face (everything is from finer, nothing to
    ! be done)
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp)) then
      isSendFace = 0

    ! fromFiner - remote face (nothing to be done)
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_remoteFace_prp ) then
      isSendFace = 0

    ! remote - fromFiner face (nothing to be done)
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_fromFinerFace_prp ) then
      isSendFace = 0

    ! remote+fromCoarser - fromCoarser face (elements are on coarser
    ! level, here and remote so nothing to be done)
    elseif( leftPrp.eq.ior(tem_fromCoarserFace_prp,tem_remoteFace_prp) .and.   &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isSendFace = 0

    ! fromCoarser - remote+fromCoarser face (elements are on coarser
    ! level, here and remote so nothing to be done)
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.ior(tem_fromCoarserFace_prp,tem_remoteFace_prp)) then
      isSendFace = 0

    ! fluid - boundary face (nothing to be done)
    elseif(leftPrp.eq.tem_fluidFace_prp .and. rightPrp.eq.tem_bndFace_prp) then
      isSendFace = 0

    ! boundary - fluid face (nothing to be done)
    elseif(leftPrp.eq.tem_bndFace_prp .and. rightPrp.eq.tem_fluidFace_prp) then
      isSendFace = 0

    ! remote - boundary face (nothing to be done)
    elseif(leftPrp.eq.tem_remoteFace_prp .and. rightPrp.eq.tem_bndFace_prp) then
      isSendFace = 0

    ! boundary - remote face (nothing to be done)
    elseif(leftPrp.eq.tem_bndFace_prp .and. rightPrp.eq.tem_remoteFace_prp) then
      isSendFace = 0

    ! fromFiner - boundary face (nothing to be done)
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_bndFace_prp ) then
      isSendFace = 0

    ! boundary - fromFiner face (nothing to be done)
    elseif( leftPrp.eq.tem_bndFace_prp .and.                                   &
      &     rightPrp.eq.tem_fromFinerFace_prp ) then
      isSendFace = 0

    ! fromFiner+remote - boundary face (nothing to be done)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp) .and.     &
      &     rightPrp.eq.tem_bndFace_prp ) then
      isSendFace = 0

    ! boundary - fromFiner+remote face (nothing to be done)
    elseif( leftPrp.eq.tem_bndFace_prp .and.                                   &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp)) then
      isSendFace = 0

    ! fromCoarser - boundary face (nothing to be done)
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_bndFace_prp ) then
      isSendFace = 0

    ! boundary - fromCoarser face (nothing to be done)
    elseif( leftPrp.eq.tem_bndFace_prp .and.                                   &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isSendFace = 0

    end if

  end function tem_isSendFace
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Creates dimension by dimension level descriptors.
  !!
  subroutine tem_dimByDim_construction( tree, boundary, commPattern, proc, &
    &                                   levelDescX, levelDescY, levelDescZ )
    ! --------------------------------------------------------------------------
    !> Tree representation of the mesh.
    type(treelmesh_type), intent(inout) :: tree
    !> The boundaries of your simulation domain
    type(tem_bc_prop_type), intent(in) :: boundary
    !> The communication pattern you use for the buffer.
    type(tem_commpattern_type), intent(in) :: commPattern
    !> Level descriptor for each spatial direction and each level of your mesh.
    !! The level descriptor have to be constructed with the dimension by
    !! dimension stencils (+1, 0, -1) for each spatial direction.
    type(tem_levelDesc_type), allocatable, intent(out) :: levelDescX(:)
    !> Level descriptor for each spatial direction and each level of your mesh.
    !! The level descriptor have to be constructed with the dimension by
    !! dimension stencils (+1, 0, -1) for each spatial direction.
    type(tem_levelDesc_type), allocatable, intent(out) :: levelDescY(:)
    !> Level descriptor for each spatial direction and each level of your mesh.
    !! The level descriptor have to be constructed with the dimension by
    !! dimension stencils (+1, 0, -1) for each spatial direction.
    type(tem_levelDesc_type), allocatable, intent(out) :: levelDescZ(:)
    !> Process description to use.
    type(tem_comm_env_type), intent(in) :: proc
    ! --------------------------------------------------------------------------
    ! The dimension-by-dimension stencil list.
    type(tem_stencilHeader_type) :: stencilList(3)
    ! --------------------------------------------------------------------------

    ! Create the dimension by dimension stencils
    write(logUnit(3),*) 'Creating dimension by dimension stencils for face ' &
      &                 //  'descriptor...'
    call tem_define_dimStencil( stencilList(1), tree%nElems, 1 )
    call tem_define_dimStencil( stencilList(2), tree%nElems, 2 )
    call tem_define_dimStencil( stencilList(3), tree%nElems, 3 )

    ! Dump the stencils we used.
    call tem_stencil_dump( stencilList( 1 ))
    call tem_stencil_dump( stencilList( 2 ))
    call tem_stencil_dump( stencilList( 3 ))

    ! Now, build the level descriptors for each of the stencils.
    call tem_create_levelDesc( tree, stencilList(1), boundary, commPattern, &
      &                        levelDescX, proc                             )
    call tem_create_levelDesc( tree, stencilList(2), boundary, commPattern, &
      &                        levelDescY, proc                             )
    call tem_create_levelDesc( tree, stencilList(3), boundary, commPattern, &
      &                        levelDescZ, proc                             )

  end subroutine tem_dimByDim_construction
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Create a face level descriptor
  !!
  subroutine tem_create_levelDesc( tree, stencil, boundary, commPattern, &
    &                              levelDesc, proc                       )
    ! --------------------------------------------------------------------------
    !> Tree representation of the mesh.
    type(treelmesh_type), intent(inout) :: tree
    !> The stencil to create the level descriptor for.
    type(tem_stencilHeader_type), intent(inout) :: stencil(1)
    !> The boundaries of your simulation domain
    type(tem_bc_prop_type), intent(in) :: boundary
    !> The communication pattern you use for the buffer.
    type(tem_commpattern_type), intent(in) :: commPattern
    !> The created level descriptor.
    type(tem_levelDesc_type), allocatable, intent(out) :: levelDesc(:)
    !> Process description to use.
    type(tem_comm_env_type), intent(in) :: proc
    ! --------------------------------------------------------------------------
    integer, allocatable :: levelPointer(:)
    ! --------------------------------------------------------------------------

    call tem_init_elemLevels( me             = levelDesc, &
      &                       boundary       = boundary,  &
      &                       tree           = tree,      &
      &                       stencils       = stencil    )

    call tem_find_allElements( tree           = tree,         &
      &                        levelDesc      = levelDesc,    &
      &                        levelPointer   = levelPointer, &
      &                        computeStencil = stencil,      &
      &                        commpattern    = commPattern,  &
      &                        proc           = proc          )

    ! JQ: allocated in tem_init_elemLevels
    ! do iLevel=tree%global%minlevel,tree%global%maxlevel
    !   allocate(levelDesc(iLevel)%neigh(1))
    ! end do
    call tem_build_HorizontalDependencies( iStencil       = 1,         &
      &                                    levelDesc      = levelDesc, &
      &                                    tree           = tree,      &
      &                                    computeStencil = stencil(1) )

    call tem_build_VerticalDependencies( levelDesc = levelDesc,            &
      &                                  minLevel  = tree%global%minLevel, &
      &                                  maxLevel  = tree%global%maxLevel  )

    call tem_cleanupDependencyArrays( levelDesc = levelDesc )

  end subroutine tem_create_levelDesc
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> summary: Creates the dimension-by-dimension stencils for a given spatial
  !! direction.
  !!
  subroutine tem_define_dimStencil( stencil, nElems, direction)
    ! --------------------------------------------------------------------------
    !> The stencil layout to set.
    type(tem_stencilHeader_type), intent(inout) :: stencil
    !> The number of elements which share this stencil
    integer, intent(in) :: nElems
    !> The spatial direction:
    !! 1 -> x direction
    !! 2 -> y direction
    !! 3 -> z direction
    integer, intent(in) :: direction
    ! --------------------------------------------------------------------------
    integer :: iElem, offsetIndex
    ! --------------------------------------------------------------------------

    stencil%useAll = .true.
    stencil%nElems = nElems
    call init( me = stencil%elem, length = 0 )
    stencil%QQN = 2
    stencil%QQ = 2
    allocate( stencil%cxDir(3, 2) )
    allocate( stencil%cxDirInv(2) )

    do iElem = 1,2
      offsetIndex = (iElem-1)*3 + direction
      stencil%cxDir(1:3, iElem) = qOffset(offsetIndex, 1:3)
      stencil%cxDirInv(iElem) = qInvDir(offsetIndex)
    end do

  end subroutine tem_define_dimStencil
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> summary: Routine to generate the separated lists of faces for
  !! communicated, refined/coarsened and local faces.
  !!
  subroutine tem_build_faceLists( minLevel, maxLevel, nEligibleChildren, faces )
    ! --------------------------------------------------------------------------
    !> Minimum level of your mesh.
    integer, intent(in) :: minLevel
    !> Maximum level of your mesh.
    integer, intent(in) :: maxLevel
    !> The number of eligible children for the vertical face dependency
    integer, intent(in) :: nEligibleChildren
    !> The created face descriptor.
    type(tem_face_type),intent(inout) :: faces(minLevel:maxLevel)
    ! --------------------------------------------------------------------------
    integer :: iLevel
    ! --------------------------------------------------------------------------

    levelLoop: do iLevel = minLevel, maxLevel

      ! Build the list of compute faces for the current level
      call tem_build_computeList( faces(iLevel), nEligibleChildren )

    end do levelLoop

    ! Build the list of from finer faces for the current level
    call tem_build_fromFinerList( minLevel, maxLevel, nEligibleChildren, faces )

  end subroutine tem_build_faceLists
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Extracts all the faces which will be computed by current rank.
  !!
  subroutine tem_build_fromFinerList( minLevel, maxLevel, nEligibleChildren, &
                                    & faces )
    ! --------------------------------------------------------------------------
    !> The min refinement level of the mesh.
    integer, intent(in) :: minLevel
    !> The max refinement level of the mesh.
    integer, intent(in) :: maxLevel
    !> The number of eligible children for the vertical face dependency
    integer, intent(in) :: nEligibleChildren
    !> The created face descriptor.
    type(tem_face_type),intent(inout) :: faces(minLevel:maxLevel)
    ! --------------------------------------------------------------------------
    integer :: nFromFinerFaces(3,2)
    integer :: iDir, iFace, iFromFinerFace(2), iDep, iLevel
    integer :: childFace, childFaceOp, childPos, childPosOp, isFromFiner
    ! --------------------------------------------------------------------------

    ! Loop over all the levels
    levelLoop: do iLevel = minLevel, maxLevel

      ! First, we count all the element which will be covered by this rank.
      nFromFinerFaces = 0
      do iDir = 1,3
        do iFace = 1, faces(iLevel)%faces(iDir)%faceList%faceId%nVals

          ! decide if the element will be computed on this level and this rank.
          isFromFiner = tem_isFromFinerFace(iFace, faces(iLevel)%faces(iDir))
          if( isFromFiner .eq. tem_left) then
            ! If left element is refined, its right face has to be interpolated.
            nFromFinerFaces(iDir,tem_right) =                          &
              &                     nFromFinerFaces(iDir,tem_right) + 1
          elseif( isFromFiner .eq. tem_right) then
            ! If right element is refined, its left face has to be interpolated.
            nFromFinerFaces(iDir,tem_left) =                           &
              &                      nFromFinerFaces(iDir,tem_left) + 1
          end if

        end do
      end do

      ! Allocate memory for the from finer faces
      do iDir = 1,3
        do iFace = 1,2
          allocate( faces(iLevel)%faces(iDir)%fromFinerFace(iFace)%            &
            &                          elemPos( nFromFinerFaces(iDir, iFace) ))
          allocate( faces(iLevel)%faces(iDir)%fromFinerFace(iFace)%            &
            &      childPos( nEligibleChildren, nFromFinerFaces(iDir, iFace) ))
          allocate( faces(iLevel)%faces(iDir)%fromFinerFace(iFace)%            &
            &                        elemPosOp( nFromFinerFaces(iDir, iFace) ))
          allocate( faces(iLevel)%faces(iDir)%fromFinerFace(iFace)%            &
            &    childPosOp( nEligibleChildren, nFromFinerFaces(iDir, iFace) ))
          allocate( faces(iLevel)%faces(iDir)%fromFinerFace(iFace)%            &
            &                          facePos( nFromFinerFaces(iDir, iFace) ))
        end do
      end do

      ! Now, we add all the compute elements to the iterator
      do iDir = 1,3
        iFromFinerFace(:) = 1
        do iFace = 1, faces(iLevel)%faces(iDir)%faceList%faceId%nVals

          ! Decide if the face has to be interpolated.
          isFromFiner = tem_isFromFinerFace(iFace, faces(iLevel)%faces(iDir))
          if( isFromFiner .eq. tem_left) then
            ! The left element is refined -> add the left element to the right
            ! faces that have to be interpolated.
            faces(iLevel)%faces(iDir)%fromFinerFace(tem_right)%        &
              &       elemPos(iFromFinerFace(tem_right)) =             &
              &       faces(iLevel)%faces(iDir)%faceList%leftElemPos%val(iFace)
            faces(iLevel)%faces(iDir)%fromFinerFace(tem_right)%        &
              &      elemPosOp(iFromFinerFace(tem_right)) =            &
              &      faces(iLevel)%faces(iDir)%faceList%rightElemPos%val(iFace)
            faces(iLevel)%faces(iDir)%fromFinerFace(tem_right)%        &
              &              facePos(iFromFinerFace(tem_right)) = iFace
            ! Now, we set the dependencies, i.e. we identify the positions of
            ! the child elements and set them in the from finer description.
            do iDep = 1,nEligibleChildren
              childFace = faces(iLevel)%faces(iDir)%faceDep%                   &
                &                                      childFacePos(iDep,iFace)
              childPos = faces(iLevel+1)%faces(iDir)%faceList%                 &
                &                                    leftElemPos%val(childFace)
              faces(iLevel)%faces(iDir)%fromFinerFace(tem_right)%      &
                &   childPos(iDep,iFromFinerFace(tem_right)) = childPos
              childFaceOp = faces(iLevel)%faces(iDir)%faceDep%                 &
                &                                      childFacePosOp(iDep,iFace)
              childPosOp = faces(iLevel+1)%faces(iDir)%faceList%               &
                &                                   rightElemPos%val(childFaceOp)
              faces(iLevel)%faces(iDir)%fromFinerFace(tem_right)%      &
                &            childPosOp(iDep,iFromFinerFace(tem_right))&
                & = childPosOp
            end do

            ! Count the number of right faces
            iFromFinerFace(tem_right) =                                &
              &                           iFromFinerFace(tem_right) + 1

          elseif( isFromFiner .eq. tem_right) then
            ! The right element is refined -> add the right element to the left
            ! faces that have to be interpolated.
            faces(iLevel)%faces(iDir)%fromFinerFace(tem_left)%         &
              &      elemPos(iFromFinerFace(tem_left)) =               &
              &      faces(iLevel)%faces(iDir)%faceList%rightElemPos%val(iFace)
            faces(iLevel)%faces(iDir)%fromFinerFace(tem_left)%         &
              &       elemPosOp(iFromFinerFace(tem_left)) =            &
              &       faces(iLevel)%faces(iDir)%faceList%leftElemPos%val(iFace)
            faces(iLevel)%faces(iDir)%fromFinerFace(tem_left)%         &
              &               facePos(iFromFinerFace(tem_left)) = iFace
            ! Now, we set the dependencies.
            ! Now, we set the dependencies, i.e. we identify the positions of
            ! the child elements and set them in the from finer description.
            do iDep = 1,nEligibleChildren
              childFace = faces(iLevel)%faces(iDir)%faceDep%                   &
                &                                      childFacePos(iDep,iFace)
              childPos = faces(iLevel+1)%faces(iDir)%faceList%                 &
                &                                   rightElemPos%val(childFace)
              faces(iLevel)%faces(iDir)%fromFinerFace(tem_left)%       &
                &    childPos(iDep,iFromFinerFace(tem_left)) = childPos
              childFaceOp = faces(iLevel)%faces(iDir)%faceDep%                 &
                &                                      childFacePosOp(iDep,iFace)
              childPosOp = faces(iLevel+1)%faces(iDir)%faceList%               &
                &                                    leftElemPos%val(childFaceOp)
              faces(iLevel)%faces(iDir)%fromFinerFace(tem_left)%       &
                & childPosOp(iDep,iFromFinerFace(tem_left)) = childPosOp
            end do

            ! Count the number of left faces
            iFromFinerFace(tem_left) =                                 &
              &                            iFromFinerFace(tem_left) + 1

          end if
        end do
      end do
    end do levelLoop

  end subroutine tem_build_fromFinerList
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> summary: Function to decide if a certain face is constructed from finer
  !! faces on the current rank.
  !!
  function tem_isFromFinerFace(facePos, faces) result( isFromFiner )
    ! --------------------------------------------------------------------------
    integer, intent(in) :: facePos
    type(tem_face_descriptor_type), intent(in) :: faces
    integer :: isFromFiner
    ! --------------------------------------------------------------------------
    integer :: leftPrp, rightPrp
    ! --------------------------------------------------------------------------

    leftPrp = faces%faceList%leftPrp%val(facePos)
    rightPrp = faces%faceList%rightPrp%val(facePos)

    ! Fluid-fluid face (purely local)
    if( leftPrp.eq.tem_fluidFace_prp .and. rightPrp.eq.tem_fluidFace_prp ) then
      isFromFiner = 0

    ! Fluid-ghostFromFiner face (purley local, no remote property set)
    elseif( leftPrp.eq.tem_fluidFace_prp .and.                                 &
      &     rightPrp.eq.tem_fromFinerFace_prp ) then
      isFromFiner = tem_right

    ! GhostFromFiner-fluid face (purley local, no remote property set)
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_fluidFace_prp ) then
      isFromFiner = tem_left

    ! GhostFromFiner-GhostFromFiner face (purley local, no remote property set)
    ! --> two ghost elements from finer meet each other somehow. However, we
    !     do not have to interpolate for these face combination.
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_fromFinerFace_prp ) then
      isFromFiner = 0

    ! Fluid-ghostFromCoarser face (purely local, no remote property)
    ! --> no from finer face.
    elseif( leftPrp.eq.tem_fluidFace_prp .and.                                 &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isFromFiner = 0

    ! GhostFromCoarser-fluid face (purely local, no remote property)
    ! --> no from finer face.
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_fluidFace_prp ) then
      isFromFiner = 0

    ! GhostFromCoarser-GhostFromCoarser face (purely local, no remote property)
    ! --> no from finer face. Somehow, two ghost from coarser meet each other,
    !     but it is not from finer face.
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isFromFiner = 0

    ! Fluid-remote face (no refinement property set)
    ! --> no refinement, so it is not from finer.
    elseif( leftPrp.eq.tem_fluidFace_prp .and.                                 &
      &     rightPrp.eq.tem_remoteFace_prp ) then
      isFromFiner = 0

    ! Remote-fluid face (no refinement property set)
    ! --> no refinement, so it is not from finer.
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_fluidFace_prp ) then
      isFromFiner = 0

    ! GhostFromCoarser-remote face
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_remoteFace_prp ) then
      isFromFiner = 0

    ! Remote-GhostFromCoarser face
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isFromFiner = 0

    ! Remote-remote face (somehow, two halo elements are meeting)
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_remoteFace_prp ) then
      isFromFiner = 0

    ! notExist - remote face (nothing from a finer level)
    elseif( leftPrp.eq.tem_notExist_prp .and.                                  &
      &     rightPrp.eq.tem_remoteFace_prp ) then
      isFromFiner = 0

    ! remote - notExist face (nothing from a finer level)
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_notExist_prp ) then
      isFromFiner = 0

    ! fluid-remote+fromCoaser face (nothing is from a finer level)
    elseif( leftPrp.eq.tem_fluidFace_prp .and.                                 &
      &     rightPrp.eq.ior(tem_remoteFace_prp, tem_fromCoarserFace_prp)) then
      isFromFiner = 0

    ! remote+fromCoarser-fluid face (nothing is from a finer level)
    elseif( leftPrp.eq.ior(tem_remoteFace_prp,tem_fromCoarserFace_prp) .and.   &
      &     rightPrp.eq.tem_fluidFace_prp) then
      isFromFiner = 0

    ! notExist-remote+fromCoaser face (neither my ranks hold an element for this
    ! face, nor elements are from finer. So, this face is not from finer).
    elseif( leftPrp.eq.tem_notExist_prp .and.                                  &
      &     rightPrp.eq.ior(tem_remoteFace_prp, tem_fromCoarserFace_prp)) then
      isFromFiner = 0

    ! remote+fromCoarser-remote face (neither my ranks hold an element for this
    ! face, nor elements are from finer. So, this face is not from finer).
    elseif( leftPrp.eq.ior(tem_remoteFace_prp, tem_fromCoarserFace_prp) .and.  &
      &     rightPrp.eq.tem_notExist_prp ) then
      isFromFiner = 0

    ! remote+fromCoarser-remote+fromCoaser face (neither my ranks hold an
    ! element for this face, nor elements are from finer. So, this face is not
    ! from finer).
    elseif( leftPrp.eq.ior(tem_remoteFace_prp, tem_fromCoarserFace_prp) .and.  &
      &     rightPrp.eq.ior(tem_remoteFace_prp, tem_fromCoarserFace_prp)) then
      isFromFiner = 0

    ! fromFiner - not exist face (this situation happens if a coarser element is
    ! located next to finer elements which are completely located on my rank. So
    ! my rank has to interpolate them upwards before we send them to the coarse
    ! element on the remote rank)
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_notExist_prp ) then
      isFromFiner = 0

    ! notExist - fromFiner face (this situation happens if a coarser element is
    ! located next to finer elements which are completely located on my rank. So
    ! my rank has to interpolate them upwards before we send them to the coarse
    ! element on the remote rank)
    elseif( leftPrp.eq.tem_notExist_prp .and.                                  &
      &     rightPrp.eq.tem_fromFinerFace_prp ) then
      isFromFiner = 0

    ! fromCoarser - notExist face (nothing is from finer)
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_notExist_prp ) then
      isFromFiner = 0

    ! notExist - fromCoarser face (nothing is from finer)
    elseif( leftPrp.eq.tem_notExist_prp .and.                                  &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isFromFiner = 0

    ! fluid-remote+fromFiner face (data is send to the remote rank on this
    ! level. The remote rank will also do the interpolation for me, so my rank
    ! does not have to do anyhting here.)
    elseif( leftPrp.eq.tem_fluidFace_prp .and.                                 &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp, tem_remoteFace_prp)) then
      isFromFiner = tem_right

    ! remote+fromFiner-fluid face (data is send to the remote rank on this
    ! level. The remote rank will also do the interpolation for me, so my rank
    ! does not have to do anyhting here.)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp, tem_remoteFace_prp) .and.    &
      &     rightPrp.eq.tem_fluidFace_prp ) then
      isFromFiner = tem_left

    ! remote+fromFiner-remote face (everything is remote)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp, tem_remoteFace_prp) .and.    &
      &     rightPrp.eq.tem_remoteFace_prp ) then
      isFromFiner = 0

    ! remote-remote+fromFiner face (everything is remote)
    elseif( leftPrp.eq. tem_remoteFace_prp .and.                               &
      &     rightPrp.eq.ior(tem_remoteFace_prp, tem_fromFinerFace_prp)) then
      isFromFiner = 0

    ! notExist-remote+fromFiner face (nothing is on my rank)
    elseif( leftPrp.eq. tem_notExist_prp .and.                                 &
      &     rightPrp.eq.ior(tem_remoteFace_prp, tem_fromFinerFace_prp)) then
      isFromFiner = 0

    ! remote+fromFiner-notExist face (nothing is on my rank)
    elseif( leftPrp.eq.ior(tem_remoteFace_prp, tem_fromFinerFace_prp) .and.    &
      &     rightPrp.eq.tem_notExist_prp) then
      isFromFiner = 0

    ! fromFiner-fromCoarser face (everything is local, this is an intermediate
    ! face level (i.e. leveljumps of difference larger than 1 occure here) )
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_fromCoarserFace_prp ) then
      isFromFiner = tem_left

    ! fromCoarser-fromFiner face (everything is local, this is an intermediate
    ! face level (i.e. leveljumps of difference larger than 1 occure here) )
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_fromFinerFace_prp ) then
      isFromFiner = tem_right

    ! fromFiner+remote - fromCoarser face (intermediate level for domain
    ! boundary and refinement boundary)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp) .and.     &
      &     rightPrp.eq.tem_fromCoarserFace_prp) then
      isFromFiner = tem_left

    ! fromCoarser - remote+fromFiner face (intermediate level for domain
    ! boundary and refinement boundary)
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp)) then
      isFromFiner = tem_right

    ! remote+fromFiner - remote+fromFiner face (everything is remote, nothing
    ! to be done)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp)  .and.    &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp)) then
      isFromFiner = 0

    ! remote+fromFiner - fromFiner face (everything is from finer, nothing to
    ! be done)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp) .and.     &
      &     rightPrp.eq.tem_fromFinerFace_prp) then
      isFromFiner = 0

    ! fromFiner - remote+fromFiner face (everything is from finer, nothing to
    ! be done)
    elseif( leftPrp.eq.tem_fromFinerFace_prp  .and.                            &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp)) then
      isFromFiner = 0

    ! fromFiner - remote face (nothing to be done)
    elseif(leftPrp.eq.tem_fromFinerFace_prp .and. rightPrp.eq.tem_remoteFace_prp) then
      isFromFiner = 0

    ! remote - fromFiner face (nothing to be done)
    elseif( leftPrp.eq.tem_remoteFace_prp .and.                                &
      &     rightPrp.eq.tem_fromFinerFace_prp) then
      isFromFiner = 0

    ! remote+fromCoarser - fromCoarser face (elements are on coarser
    ! level, here and remote so nothing to be done)
    elseif( leftPrp.eq.ior(tem_fromCoarserFace_prp,tem_remoteFace_prp) .and.   &
      &     rightPrp.eq.tem_fromCoarserFace_prp) then
      isFromFiner = 0

    ! fromCoarser - remote+fromCoarser face (elements are on coarser
    ! level, here and remote so nothing to be done)
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.ior(tem_fromCoarserFace_prp,tem_remoteFace_prp)) then
      isFromFiner = 0

    ! fluid - boundary face (nothing to be done)
    elseif(leftPrp.eq.tem_fluidFace_prp .and. rightPrp.eq.tem_bndFace_prp) then
      isFromFiner = 0

    ! boundary - fluid face (nothing to be done)
    elseif(leftPrp.eq.tem_bndFace_prp .and. rightPrp.eq.tem_fluidFace_prp) then
      isFromFiner = 0

    ! remote - boundary face (nothing to be done)
    elseif(leftPrp.eq.tem_remoteFace_prp .and. rightPrp.eq.tem_bndFace_prp) then
      isFromFiner = 0

    ! boundary - remote face (nothing to be done)
    elseif(leftPrp.eq.tem_bndFace_prp .and. rightPrp.eq.tem_remoteFace_prp) then
      isFromFiner = 0

    ! fromFiner - boundary face (nothing to be done)
    elseif( leftPrp.eq.tem_fromFinerFace_prp .and.                             &
      &     rightPrp.eq.tem_bndFace_prp) then
      isFromFiner = 0

    ! boundary - fromFiner face (nothing to be done)
    elseif( leftPrp.eq.tem_bndFace_prp .and.                                   &
      &     rightPrp.eq.tem_fromFinerFace_prp) then
      isFromFiner = 0

    ! fromFiner+remote - boundary face (nothing to be done)
    elseif( leftPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp) .and.     &
      &     rightPrp.eq.tem_bndFace_prp) then
      isFromFiner = 0

    ! boundary - fromFiner+remote face (nothing to be done)
    elseif( leftPrp.eq.tem_bndFace_prp .and.                                   &
      &     rightPrp.eq.ior(tem_fromFinerFace_prp,tem_remoteFace_prp)) then
      isFromFiner = 0

    ! fromCoarser - boundary face (nothing to be done)
    elseif( leftPrp.eq.tem_fromCoarserFace_prp .and.                           &
      &     rightPrp.eq.tem_bndFace_prp) then
      isFromFiner = 0

    ! boundary - fromCoarser face (nothing to be done)
    elseif( leftPrp.eq.tem_bndFace_prp .and.                                   &
      &     rightPrp.eq.tem_fromCoarserFace_prp) then
      isFromFiner = 0

!> \todo JZ: Add more cases for the different types of faces. Especially
!! combinations of finer/coarser elements and remote properties are still
!! missing.
!! For some combinations it might be necessary to determine if vertical
!! dependencies exist or not.

    else
      write(logUnit(1),*) 'ERROR in tem_isFromFinerFace: not able to decide whether '// &
        &        'face is locally from finer or not, stopping'
      write(logUnit(1),*) 'left elem pos: ', faces%faceList%leftElemPos%val(facePos)
      write(logUnit(1),*) 'right elem pos: ', faces%faceList%rightElemPos%val(facePos)
      write(logUnit(1),*) 'left face prp: ', leftPrp
      write(logUnit(1),*) 'right face prp: ', rightPrp
      call tem_abort()
    end if

  end function tem_isFromFinerFace
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Extracts all the faces which will be computed by current rank.
  !!
  subroutine tem_build_computeList( faces, nEligibleChildren )
    ! --------------------------------------------------------------------------
    !> The created face descriptor.
    type(tem_face_type),intent(inout) :: faces
    !> The number of eligible children for the vertical face dependency
    integer, intent(in) :: nEligibleChildren
    ! --------------------------------------------------------------------------
    integer :: nComputeFaces(3)
    integer :: iDir, iFace, iComputeFace
    ! --------------------------------------------------------------------------


    ! First, we count all the element which will be covered by this rank.
    nComputeFaces = 0
    do iDir = 1,3

      do iFace = 1, faces%faces(iDir)%faceList%faceId%nVals

        ! decide if the element will be computed on this level and this rank.
        if ( tem_isComputeFace(iFace,             &
          &                    faces%faces(iDir), &
          &                    nEligibleChildren) ) then
          nComputeFaces(iDir) = nComputeFaces(iDir) + 1
        end if

      end do

    end do

    ! Allocate memory for the compute faces
    do iDir = 1,3
      allocate( faces%faces(iDir)%computeFace%leftPos( nComputeFaces(iDir) ) )
      allocate( faces%faces(iDir)%computeFace%rightPos( nComputeFaces(iDir) ) )
      allocate( faces%faces(iDir)%computeFace%facePos( nComputeFaces(iDir) ) )
    end do

    ! Now, we add all the compute elements to the iterator
    do iDir = 1,3
      iComputeFace = 1
      do iFace = 1, faces%faces(iDir)%faceList%faceId%nVals

        ! Decide if the element will be computed on this level and this rank.
        ! If yes, we add it to the list of compute faces.
        if ( tem_isComputeFace(iFace,             &
          &                    faces%faces(iDir), &
          &                    nEligibleChildren) ) then
          faces%faces(iDir)%computeFace%leftPos(iComputeFace) =                &
            &                 faces%faces(iDir)%faceList%leftElemPos%val(iFace)
          faces%faces(iDir)%computeFace%rightPos(iComputeFace) =               &
            &                faces%faces(iDir)%faceList%rightElemPos%val(iFace)
          faces%faces(iDir)%computeFace%facePos(iComputeFace) = iFace
          iComputeFace = iComputeFace + 1
        end if

      end do
    end do
  end subroutine tem_build_computeList
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> summary: Builds all the vertical dependencies in the face description.
  !!
  subroutine tem_faceDep_vertical(minLevel, maxLevel, faces, nEligibleChildren )
    ! --------------------------------------------------------------------------
    !> Minimum level of your mesh.
    integer, intent(in) :: minLevel
    !> Maximum level of your mesh.
    integer, intent(in) :: maxLevel
    !> Face descriptor where the faces will be appended to.
    type(tem_face_type),intent(inout)  :: faces(minLevel:maxLevel)
    !> The number of eligible children for the vertical face dependency
    integer, intent(in) :: nEligibleChildren
    ! --------------------------------------------------------------------------
    integer :: iLevel, iDir
    ! --------------------------------------------------------------------------

    ! Initialize the container for the vertical face dependencies.
    levelLoop: do iLevel = minLevel, maxLevel
      do iDir = 1, 3
        call tem_init_faceDep( faces(iLevel)%faces(iDir), nEligibleChildren )
      end do
    end do levelLoop

    ! build the dependencies downward (coarse face -> finer faces)
    levelLoopDown: do iLevel = minLevel, maxLevel-1
      do iDir = 1, 3
        call tem_faceDep_verticalDown( faces(iLevel)%faces(iDir),         &
          &                            faces(iLevel+1)%faces(iDir), iDir, &
          &                            nEligibleChildren )
      end do
    end do levelLoopDown

    ! build the dependencies upward (fine face -> coarse face)
    levelLoopUp: do iLevel = maxLevel, minLevel+1, -1
      do iDir = 1, 3
        call tem_faceDep_verticalUp( faces(iLevel)%faces(iDir),                &
          &                          faces(iLevel-1)%faces(iDir) )
      end do
    end do levelLoopUp

  end subroutine tem_faceDep_vertical
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> summary: Subroutine to initialize the container for the vertical
  !! dependencies of a face descriptor.
  !!
  subroutine tem_init_faceDep( faces, nEligibleChildren )
    ! --------------------------------------------------------------------------
    !> The face description in which you want to initialize the vertical
    !! face dependency container.
    type(tem_face_descriptor_type), intent(inout) :: faces
    !> The number of eligible children for the vertical face dependency
    integer, intent(in) :: nEligibleChildren
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    allocate( faces%faceDep%parentFaceId(faces%faceList%faceId%nVals) )
    allocate( faces%faceDep%parentFacePos(faces%faceList%faceId%nVals) )
    allocate( faces%faceDep%childFaceId( nEligibleChildren, &
                                       & faces%faceList%faceId%nVals) )
    allocate( faces%faceDep%childFaceIdOp( nEligibleChildren, &
                                       & faces%faceList%faceId%nVals) )
    allocate( faces%faceDep%childFacePos( nEligibleChildren, &
                                       & faces%faceList%faceId%nVals) )
    allocate( faces%faceDep%childFacePosOp( nEligibleChildren, &
                                       & faces%faceList%faceId%nVals) )

    faces%faceDep%parentFaceId = -1_long_k
    faces%faceDep%parentFacePos = -1
    faces%faceDep%childFaceId = -1_long_k
    faces%faceDep%childFacePos(:,:) = -1

  end subroutine tem_init_faceDep
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> summary: Builds downward face depndencies (coarse face -> finer faces).
  !!
  subroutine tem_faceDep_verticalDown( coarseFaces, fineFaces, dir, &
                                     & nEligibleChildren )
    ! --------------------------------------------------------------------------
    !> Face description on the coarse level. The dependencies to the finer
    !! level will be appended to this face descriptor.
    type(tem_face_descriptor_type), intent(inout) :: coarseFaces
    !> Face description on the finer level.
    type(tem_face_descriptor_type), intent(in) :: fineFaces
    !> The spatial direction of the face to add the downward dependencies for.
    !! 1 --> x direction
    !! 2 --> y direction
    !! 3 --> z direction
    integer, intent(in) :: dir
    !> The number of eligible children for the vertical face dependency
    integer, intent(in) :: nEligibleChildren
    ! --------------------------------------------------------------------------
    integer :: iFace
    integer :: leftPrp, rightPrp
    ! --------------------------------------------------------------------------

    ! Run over all the faces and add the downward depndency (depending on the
    ! properties of the face).
    faceLoop: do iFace = 1, coarseFaces%faceList%faceId%nVals

      ! get its properties from left and right
      call tem_getFace_prp( coarseFaces, iFace, leftPrp, rightPrp )

      ! check if this face requires downwards dependencies. If yes, we just add
      ! this dependency to the coarse face descriptor.
      if( tem_reqDep_down(leftPrp, rightPrp) ) then
        call tem_addDep_down( iFace, coarseFaces, fineFaces, &
                            & dir, nEligibleChildren )
      end if

    end do faceLoop

  end subroutine tem_faceDep_verticalDown
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> summary: Builds upward face depndencies (fine face -> coarse faces).
  !!
  subroutine tem_faceDep_verticalUp( fineFaces, coarseFaces )
    ! --------------------------------------------------------------------------
    !> Face description on the fine level. The dependencies to the coarser
    !! level will be appended to this face descriptor.
    type(tem_face_descriptor_type), intent(inout) :: fineFaces
    !> Face description on the coarse level.
    type(tem_face_descriptor_type), intent(in) :: coarseFaces
    ! --------------------------------------------------------------------------
    integer :: iFace
    integer :: leftPrp, rightPrp
    ! --------------------------------------------------------------------------

    ! Run over all the faces and add the downward depndency (depending on the
    ! properties of the face).
    faceLoop: do iFace = 1, fineFaces%faceList%faceId%nVals

      ! get its properties from left and right
      call tem_getFace_prp( fineFaces, iFace, leftPrp, rightPrp )

      ! check if this face requires upwards dependencies. If yes, we just add
      ! this dependency to the coarse face descriptor.
      if( tem_reqDep_up(leftPrp, rightPrp) ) then
        call tem_addDep_up( iFace, fineFaces, coarseFaces )
      end if

    end do faceLoop

  end subroutine tem_faceDep_verticalUp
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Subroutine returns left and right propery of a face.
  !!
  subroutine tem_getFace_prp( faces, facePos, leftPrp, rightPrp)
    ! --------------------------------------------------------------------------
    !> The face descriptor the face is located in.
    type(tem_face_descriptor_type), intent(in) :: faces
    !> The position of the face in the face descriptor.
    integer, intent(in) :: facePos
    !> The left property of the face.
    integer, intent(out) :: leftPrp
    !> The right property of the face.
    integer, intent(out) :: rightPrp
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    ! read the left and right property
    leftPrp = faces%faceList%leftPrp%val(facePos)
    rightPrp = faces%faceList%rightPrp%val(facePos)

  end subroutine tem_getFace_prp
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Function to check if a face with given left and right property requires
  !! a downward (coarse->fine) dependency.
  function tem_reqDep_down(leftPrp, rightPrp) result( reqDep )
    ! --------------------------------------------------------------------------
    !> Left property of the face.
    integer, intent(in) :: leftPrp
    !> Right property of the face.
    integer, intent(in) :: rightPrp
    !> Does the face requires a vertical dependency.
    logical :: reqDep
    ! --------------------------------------------------------------------------

    reqDep = .false.

    ! Check if one of the face properties (left or right) has the from finer
    ! property.
    if ( iand(tem_fromFinerFace_prp,leftPrp).ne.0 &
      & .or.                                      &
      & iand(tem_fromFinerFace_prp,rightPrp).ne.0 ) then
      reqDep = .true.
    end if

  end function tem_reqDep_down
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Function to check if a face with given left and right property requires
  !! a upward (fine->coarse) dependency.
  !!
  function tem_reqDep_up(leftPrp, rightPrp) result( reqDep )
    ! --------------------------------------------------------------------------
    !> Left property of the face.
    integer, intent(in) :: leftPrp
    !> Right property of the face.
    integer, intent(in) :: rightPrp
    !> Does the face requires a vertical dependency.
    logical :: reqDep
    ! --------------------------------------------------------------------------

    reqDep = .false.

    ! Check if one of the face properties (left or right) has the from coarser
    ! property.
    if ( iand(tem_fromCoarserFace_prp,leftPrp)/=0  &
      &  .or.                                      &
      &  iand(tem_fromCoarserFace_prp,rightPrp)/=0 ) then
      reqDep = .true.
    end if

  end function tem_reqDep_up
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> summary: Adds an downward face dependency (coarse->fine) to the coarse face
  !! descriptor.
  !!
  subroutine tem_addDep_down( coarseFacePos, coarseFaces, fineFaces, dir, &
                            & nEligibleChildren )
    ! --------------------------------------------------------------------------
    !> The position of the coarse face in coarseFaces you want to add
    !! the child dependencies for.
    integer, intent(in) :: coarseFacePos
    !> The description of the faces on the coarse level. The dependency will be
    !! added to this face descriptor.
    type(tem_face_descriptor_type), intent(inout) :: coarseFaces
    !> The descriptor of the faces on the fine level.
    type(tem_face_descriptor_type), intent(in) :: fineFaces
    !> The spatial direction of the face to add the downward dependencies for.
    !! 1 --> x direction
    !! 2 --> y direction
    !! 3 --> z direction
    integer, intent(in) :: dir
    !> The number of eligible children for the vertical face dependency
    integer, intent(in) :: nEligibleChildren
    ! --------------------------------------------------------------------------
    integer(kind=long_k) :: faceId, rightElemId, childFaceId, childFaceIdOp, &
                          & leftOf_childFaceIdOp
    integer :: leftOf_coord(4)
    integer(kind=long_k) :: children(8), childrenOp(8)
    integer, allocatable :: eligibleChildren(:), eligibleChildrenOp(:)
    integer :: eligChildPar, eligChildParOp, iChild, childFacePos, childFacePosOp
    ! --------------------------------------------------------------------------

    ! Left element id of the face.
    faceId = coarseFaces%faceList%faceId%val(coarseFacePos)
    rightElemId = coarseFaces%faceList%rightElemId%val(coarseFacePos)

    ! Identify the all the eligible child elements.
    select case(nEligibleChildren)
    case(4)
      select case(dir)
      case(1) ! x - direction
        eligChildPar = q__E
        eligChildParOp = q__W
      case(2) ! y - direction
        eligChildPar = q__N
        eligChildParOp = q__S
      case(3) ! z - direction
        eligChildPar = q__T
        eligChildParOp = q__B
      case default
        write(*,*) 'ERROR in tem_addDep_down: unknown spatial direction, '//     &
          &        'stopping ...'
        call tem_abort()
      end select
      call tem_eligibleChildren(eligibleChildren, eligChildPar)
      call tem_eligibleChildren(eligibleChildrenOp, eligChildParOp)
      children = tem_directChildren(faceId)
      childrenOp = tem_directChildren(rightElemId)
    case(8)
      allocate(eligibleChildren(nEligibleChildren))
      allocate(eligibleChildrenOp(nEligibleChildren))
      eligibleChildren =  (/ 1, 2, 3, 4, 5, 6, 7, 8 /)
      eligibleChildrenOp =  (/ 1, 2, 3, 4, 5, 6, 7, 8 /)
      children = tem_directChildren(faceId)
      childrenOp = tem_directChildren(rightElemId)
    case default
      write(logUnit(1),*) 'ERROR in tem_addDep_down: Unknown number of'
      write(logUnit(1),*) ' eligible children, stopping ...'
      call tem_abort()
    end select

    ! Loop over the relevant, eligible children.
    childLoop: do iChild = 1,nEligibleChildren

      ! The tree id of the current eligible child
      childFaceId = children(eligibleChildren(iChild))
      childFaceIdOp = childrenOp(eligibleChildrenOp(iChild))

      ! Get the TreeID which is left of childFaceIdOp for face identification
      leftOf_coord = tem_CoordOfId( childFaceIdOp  )
      leftOf_coord(dir) = leftOf_coord(dir)-1
      leftOf_childFaceIdOp = tem_IdOfCoord(leftOf_coord)

      ! Get its position in the fine face descriptor.
      childFacePos = PositionOfVal( me  = fineFaces%faceList%faceId,           &
        &                           val = childFaceId )

      ! Get the op position in the fine face descriptor.
      childFacePosOp = PositionOfVal( me  = fineFaces%faceList%faceId,           &
        &                           val = leftOf_childFaceIdOp )

      ! If the child element exists in the fine face descriptor, add a
      ! dependency from coarse element to its child element.
      if(childFacePos<=0) then
        coarseFaces%faceDep%childFaceId(iChild, coarseFacePos) = -1_long_k
        coarseFaces%faceDep%childFacePos(iChild, coarseFacePos) = -1
      else
        coarseFaces%faceDep%childFaceId(iChild, coarseFacePos) = childFaceId
        coarseFaces%faceDep%childFacePos(iChild, coarseFacePos) = childFacePos
      end if

      if(childFacePosOp<=0) then
        coarseFaces%faceDep%childFaceIdOp(iChild, coarseFacePos) = -1_long_k
        coarseFaces%faceDep%childFacePosOp(iChild, coarseFacePos) = -1
      else
        coarseFaces%faceDep%childFaceIdOp(iChild, coarseFacePos) = leftOf_childFaceIdOp
        coarseFaces%faceDep%childFacePosOp(iChild, coarseFacePos) = childFacePosOp
      end if

    end do childLoop

  end subroutine tem_addDep_down
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> summary: Adds an upward face dependency (fine->coarse) to the fine face
  !! descriptor.
  !!
  subroutine tem_addDep_up( fineFacePos, fineFaces, coarseFaces )
    ! --------------------------------------------------------------------------
    !> The position of the face in the fineFaces descriptor you want to add.
    integer, intent(in) :: fineFacePos
    !> The description of the faces on the fine level. The dependency will be
    !! added to this face descriptor.
    type(tem_face_descriptor_type), intent(inout) :: fineFaces
    !> The descriptor of the faces on the coarse level.
    type(tem_face_descriptor_type), intent(in) :: coarseFaces
    ! --------------------------------------------------------------------------
    integer(kind=long_k) :: faceId, parentId
    integer :: parentPos
    ! --------------------------------------------------------------------------

    ! Left element id of the face.
    faceId = fineFaces%faceList%faceId%val(fineFacePos)

    ! Identify the parent element.
    parentId = tem_parentOf(faceId)

    ! Get its position in the coarse face descriptor.
    parentPos = PositionOfVal( me  = coarseFaces%faceList%faceId, &
      &                        val = parentId                     )

    ! If the parent element exists in the coarse face descriptor, add a
    ! dependency from fine element to its parent element.
    if( parentPos<=0 ) then
      fineFaces%faceDep%parentFaceId(fineFacePos) = -1_long_k
      fineFaces%faceDep%parentFacePos(fineFacePos) = -1
    else
      fineFaces%faceDep%parentFaceId(fineFacePos) = parentId
      fineFaces%faceDep%parentFacePos(fineFacePos) = parentPos
    end if

  end subroutine tem_addDep_up
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> summary: Collects all the faces which have at least one neighboring
  !! element that exists in the level descriptor.
  !!
  subroutine tem_collect_faces(levelDesc, minLevel, maxLevel, faces )
    ! --------------------------------------------------------------------------
    !> Minimum level of your mesh.
    integer, intent(in) :: minLevel
    !> Maximum level of your mesh.
    integer, intent(in) :: maxLevel
    !> Level descriptor for each spatial direction and each level of your mesh.
    type(tem_levelDesc_type), intent(in) :: levelDesc(1:3,minLevel:maxLevel)
    !> Face descriptor where the faces will be appended to.
    type(tem_face_type),intent(inout)  :: faces(minLevel:maxLevel)
    ! --------------------------------------------------------------------------
    integer :: iLevel, iDir
    ! --------------------------------------------------------------------------

    levelLoop: do iLevel = minLevel, maxLevel
      directionLoop: do iDir = 1, 3

        ! For this level and this direction we collect all the faces.
        call tem_get_faces( levelDesc = levelDesc(iDir, iLevel),  &
          &                 direction = iDir,                     &
          &                 faces     = faces(iLevel)%faces(iDir) )

      end do directionLoop
    end do levelLoop

  end subroutine tem_collect_faces
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Subroutine to initialize the dyanmic content of a face list.
  !!
  subroutine tem_init_faceList( faceList )
    ! --------------------------------------------------------------------------
    !> The face list to initialize.
    type(tem_faceList_type), intent(inout) :: faceList
    ! --------------------------------------------------------------------------

    call init( me = faceList%faceId, length = 32 )
    call init( me = faceList%rightElemId, length = 32 )
    call init( me = faceList%leftElemPos, length = 32 )
    call init( me = faceList%rightElemPos, length = 32 )
    call init( me = faceList%leftPrp, length = 32 )
    call init( me = faceList%rightPrp, length = 32 )

  end subroutine tem_init_faceList
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> summary: Collects all the faces on a given level of the mesh in a given
  !! direction. It attaches also two-sided properties to the face.
  !!
  subroutine tem_get_faces(levelDesc, direction, faces)
    ! --------------------------------------------------------------------------
    !> Level descriptor of the level of the mesh you want to collect the
    !! faces for.
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> The spatial direction to collect the faces for:
    !! 1 -> x direction
    !! 2 -> y direction
    !! 3 -> z direction
    integer, intent(in) :: direction
    !> Description of the faces on the current level.
    type(tem_face_descriptor_type), intent(inout) :: faces
    ! --------------------------------------------------------------------------
    integer :: iElem
    integer :: neighIndex
    integer :: neighPos
    integer :: elemPrp, neighPrp
    integer(kind=long_k) :: elemId, neighId
    ! --------------------------------------------------------------------------

    ! Initialize the facelist (with a dynamic array for the faceIDs)
    call tem_init_faceList( faces%faceList )

    ! Loop over all the elements (includes fluid, ghost and halos) have a look
    ! to the left and to the right face.
    ! Try to add both faces of the element to the list of all faces.
    elemLoop: do iElem = 1, levelDesc%nElems

      ! Get the current element property
      elemPrp = tem_get_elemPrp(levelDesc, iElem)
      elemId = levelDesc%total(iElem)

      ! Look for the left face of the current element.
      ! We use the compute stencil to determine the position of the neighboring
      ! left element.
      neighIndex = tem_left
      call tem_get_faceNeigh( levelDesc, iElem, direction, neighIndex, &
        &                     neighId, neighPos                        )
      neighPrp = tem_get_elemPrp(levelDesc, neighPos)

      ! For the left face, the current element will be on its right side.
      call tem_addFace( leftElemId   = neighId,  &
        &               leftElemPos  = neighPos, &
        &               leftElemPrp  = neighPrp, &
        &               rightElemId  = elemId,   &
        &               rightElemPos = iElem,    &
        &               rightElemPrp = elemPrp,  &
        &               faces        = faces     )

      ! Look for the right face of the current element.
      ! Get the right neighbor of the element.
      neighIndex = tem_right
      call tem_get_faceNeigh( levelDesc, iElem, direction, neighIndex, &
        &                     neighId, neighPos                        )
      neighPrp = tem_get_elemPrp(levelDesc, neighPos)

      ! For the right face, the current element will be on its left side.
      call tem_addFace( leftElemId   = elemId,   &
        &               leftElemPos  = iElem,    &
        &               leftElemPrp  = elemPrp,  &
        &               rightElemId  = neighId,  &
        &               rightElemPos = neighPos, &
        &               rightElemPrp = neighPrp, &
        &               faces        = faces     )

    end do elemLoop

  end subroutine tem_get_faces
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Function to get the face neighbor of a certain element in the
  !! level descriptor. Even non-existing face neighbors can be handled
  !! by this routine.
  !!
  subroutine tem_get_faceNeigh( levelDesc, elemPos, dir, leftOrRight, neighId, &
    &                           neighPos )
    ! --------------------------------------------------------------------------
    !> The level descriptor the element is located in.
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> The element position in the level descriptor.
    integer, intent(in) :: elemPos
    !> The spatial direction of the face.
    !!
    !! 1. -> x direction
    !! 2. -> y direction
    !! 3. -> z direction
    integer, intent(in) :: dir
    !> Find left or right face neighbor of the element.
    !! Use
    !! \ref tem_face_module::tem_left
    !! or
    !! \ref tem_face_module::tem_right.
    integer, intent(in) :: leftOrRight
    !> The treeid of the face neighbor element.
    integer(kind=long_k), intent(out) :: neighId
    !> Position of the neighbor in the total list, or 0 if element does not
    !! exist.
    integer, intent(out) :: neighPos
    ! --------------------------------------------------------------------------
    integer(kind=long_k) :: elemId
    integer :: elemCoord(4)
    integer :: dirOffset
    ! --------------------------------------------------------------------------

    ! We check the neighbor position of the element.
    neighPos = levelDesc%neigh(1)%nghElems( leftOrRight, elemPos )
    elemId = levelDesc%total(elemPos)

    ! Can we find the neighbor in the level descriptor. If not, we have
    ! to build the tree id our own.
    if (neighPos>0 .and. neighPos<=levelDesc%nElems) then
      neighId = levelDesc%total(neighPos)
    else
      ! So, the neighbor element does not exist in level descriptor.
      ! Let's create the tree id (of the neighbor on this level) by our own.
      elemCoord = tem_coordOfId(elemId)
      ! Initialize dirOffset just to silence compiler warning
      dirOffset = 0
      select case(leftOrRight)
      case(tem_left)
        dirOffset = -1
      case(tem_right)
        dirOffset = 1
      case default
        write(*,*) 'ERROR in tem_get_faceNeigh: unknown face neighbor '//      &
          &        'direction stopping ... '
        write(*,*) 'Given direction was: ', dir
        call tem_abort()
      end select
      elemCoord(dir) = elemCoord(dir) + dirOffset
      neighId = tem_idOfCoord(elemCoord)

      ! Check if it is boundary element or not
      if (neighPos<0) then
        ! In case of a boundary element, we keep the negative boundary
        ! unmodified.
      else
        ! If no, we try to locate the neighbor in the total list of the level
        ! descriptor
        neighPos = tem_treeIDinTotal( neighId, levelDesc)
      end if
    endif

  end subroutine tem_get_faceNeigh
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Adds a new face to the face description.
  !!
  !! If the face already exists in the face description,
  !! only the properties of the already existing face will be overwritten.
  subroutine tem_addFace( leftElemId, leftElemPos, leftElemPrp,          &
    &                     rightElemId, rightElemPos, rightElemPrp, faces )
    ! --------------------------------------------------------------------------
    !> Element id of the left element
    integer(kind=long_k), intent(in) :: leftElemId
    !> Element id of the right element
    integer(kind=long_k), intent(in) :: rightElemId
    !> Position of the left element in the level descriptor's total list.
    integer, intent(in) :: leftElemPos
    !> Properties of the left element.
    integer, intent(in) :: leftElemPrp
    !> Position of the right element in the level desriptor's total list.
    integer, intent(in) :: rightElemPos
    !> Properties of the right element
    integer, intent(in) :: rightElemPrp
    !> The face description the new face will be added to. If the face already
    !! exists in this face description. The existing face's property will be
    !! overwritten by the new ones.
    type(tem_face_descriptor_type), intent(inout) :: faces
    ! --------------------------------------------------------------------------
    logical :: wasAdded
    integer :: pos
    ! --------------------------------------------------------------------------


    ! Try to add the left tree id to the face id list of the faces
    call append( me       = faces%faceList%faceId, &
      &          val      = leftElemId,            &
      &          wasAdded = wasAdded,              &
      &          pos      = pos                    )

    ! Check if the face already existed in the face descriptor.
    ! If yes, we overwrite the properties. If no, we overwrite
    ! the existing properties.
    if(wasAdded) then ! New face
      ! Append right element, element positions and new properties
      call append( me = faces%faceList%rightElemId, val = rightElemId )
      call append( me = faces%faceList%leftElemPos, val = leftElemPos )
      call append( me = faces%faceList%rightElemPos, val = rightElemPos )
      call append( me = faces%faceList%leftPrp, val = leftElemPrp )
      call append( me = faces%faceList%rightPrp, val = rightElemPrp )
    else ! Face already existed
      ! Overwrite existing properies
      faces%faceList%leftPrp%val(pos) = leftElemPrp
      faces%faceList%rightPrp%val(pos) = rightElemPrp
    end if

  end subroutine tem_addFace
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> summary: Extends the properties of the faces which are marked as
  !! communication faces (at least by one of the neighboring elements). It
  !! checks if the communicated elements are from finer/coarser level and
  !! attaches the fromFiner/fromCoarser property to this face if necessary.
  !!
  subroutine tem_extend_remotePrp(levelDesc, minLevel, maxLevel, faces )
    ! --------------------------------------------------------------------------
    !> Minimum level of your mesh.
    integer, intent(in) :: minLevel
    !> Maximum level of your mesh.
    integer, intent(in) :: maxLevel
    !> Level descriptor for each level of your mesh.
    type(tem_levelDesc_type), intent(in) :: levelDesc(1:3,minLevel:maxLevel)
    !> The face descriptor to be corrected.
    type(tem_face_type),intent(inout)  :: faces(minLevel:maxLevel)
    ! --------------------------------------------------------------------------
    integer :: iLevel, iDir
    ! --------------------------------------------------------------------------

    ! We correct the properties for each level and each direction independently
    levelLoop: do iLevel = minLevel, maxLevel
      directionLoop: do iDir = 1, 3

        ! Check for from finer property of the communicated faces.
        call tem_extend_commFromFinerPrp( levelDesc(iDir, iLevel), iDir, &
          &                               faces(iLevel)%faces(iDir)      )

        ! Check for from coarser property of the communicated faces.
        call tem_extend_commFromCoarserPrp( levelDesc(iDir, iLevel), iDir, &
          &                                 faces(iLevel)%faces(iDir)      )

      end do directionLoop
    end do levelLoop

  end subroutine tem_extend_remotePrp
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> summary: Extend communication property for faces by the finer property if
  !! the neighboring halo element is from finer on the remote partition.
  !!
  subroutine tem_extend_commFromFinerPrp( levelDesc, direction, faces )
    ! --------------------------------------------------------------------------
    !> Level descriptor of the level of the mesh you want to collect the
    !! faces for.
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> The face direction.
    integer, intent(in) :: direction
    !> Description of the faces on the current level.
    type(tem_face_descriptor_type), intent(inout) :: faces
    ! --------------------------------------------------------------------------
    integer :: iProc, iElem
    integer :: neighIndex
    integer(kind=long_k) :: elemId, neighId
    integer :: elemPos, neighPos
    ! --------------------------------------------------------------------------

    ! We run over all the elements in the from finer buffer and attach
    ! the from finer property to the face the element is aligned to.
    processLoop: do iProc = 1, levelDesc%recvbufferFromFiner%nProcs

      elementLoop: do iElem = 1, levelDesc%recvbufferFromFiner%nElemsProc(iProc)

        ! The element position in the total list of the level descriptor.
        elemPos = levelDesc%recvbufferFromFiner%elemPos(iProc)%val(iElem)
        elemId = levelDesc%total(elemPos)

        ! Attach elem property for its right face (checks also if face exists).
        neighIndex = tem_right
        call tem_get_faceNeigh( levelDesc, elemPos, direction, neighIndex, &
          &                     neighId, neighPos                          )
        call tem_appendFace_prp( elemId, tem_fromFinerFace_prp, faces, &
          &                      tem_left                      )

        ! Attach elem property for its left face (checks also if face exists).
        neighIndex = tem_left
        call tem_get_faceNeigh( levelDesc, elemPos, direction, neighIndex, &
          &                     neighId, neighPos                          )
        call tem_appendFace_prp( neighId, tem_fromFinerFace_prp, faces, &
          &                      tem_right                      )

    end do elementLoop
  end do processLoop

  end subroutine tem_extend_commFromFinerPrp
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> summary: Extend communication property for faces by the coarser property if
  !! the neighboring halo element is from coarser on the remote partition.
  !!
  subroutine tem_extend_commFromCoarserPrp( levelDesc, direction, faces )
    ! --------------------------------------------------------------------------
    !> Level descriptor of the level of the mesh you want to collect the
    !! faces for.
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> The face direction.
    integer, intent(in) :: direction
    !> Description of the faces on the current level.
    type(tem_face_descriptor_type), intent(inout) :: faces
    ! --------------------------------------------------------------------------
    integer :: iElem, iProc
    integer :: neighIndex
    integer(kind=long_k) :: elemId, neighId
    integer :: elemPos, neighPos
    ! --------------------------------------------------------------------------

    ! We run over all the elements in the from coarser buffer and attach
    ! the from coarser property to the face the element is aligned to.
    processLoop: do iProc = 1, levelDesc%recvbufferFromCoarser%nProcs

      elementLoop: do iElem = 1, levelDesc%recvbufferFromCoarser%              &
        &                                                     nElemsProc(iProc)

        ! The element position in the total list of the level descriptor.
        elemPos = levelDesc%recvbufferFromCoarser%elemPos(iProc)%val(iElem)
        elemId = levelDesc%total(elemPos)

        ! Attach elem property for its right face (checks also if face exists).
        neighIndex = tem_right
        call tem_get_faceNeigh( levelDesc, elemPos, direction, neighIndex, &
          &                     neighId, neighPos                          )
        call tem_appendFace_prp( elemId, tem_fromCoarserFace_prp, faces, &
          &                      tem_left                        )

        ! Attach elem property for its left face (checks also if face exists).
        neighIndex = tem_left
        call tem_get_faceNeigh( levelDesc, elemPos, direction, neighIndex, &
          &                     neighId, neighPos                          )
        call tem_appendFace_prp( neighId, tem_fromCoarserFace_prp, faces, &
          &                      tem_right                        )

      end do elementLoop
    end do processLoop

  end subroutine tem_extend_commFromCoarserPrp
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Attaches another property to a given face (from left or right). If the face
  !! does not exist this routine will not do anything.
  !!
  subroutine tem_appendFace_prp( faceId, prp, faces, prp_dir )
    ! --------------------------------------------------------------------------
    !> The face identifier to append the face property to.
    integer(kind=long_k), intent(in) :: faceId
    !> The property to attach
    integer, intent(in) :: prp
    !> The face desriptor to update.
    type(tem_face_descriptor_type), intent(inout) :: faces
    !> Attach the property to the left or right property of the face.
    integer, intent(in) :: prp_dir
    ! --------------------------------------------------------------------------
    integer :: facePos, currentPrp, newPrp
    ! --------------------------------------------------------------------------

    ! Check if the face element exists and attach the property.
    facePos = PositionOfVal( faces%faceList%faceId, faceId )

    ! If the element is found: Check if the element has already the
    ! property we want to attach to it. If not, add it.
    if(facePos.gt.0) then

      ! Get the current face propery.
      if( prp_dir .eq. tem_left ) then
        currentPrp = faces%faceList%leftPrp%val(facePos)
      else
        currentPrp = faces%faceList%rightPrp%val(facePos)
      end if

      ! Melt the property we want to attach with the current one.
      newPrp = tem_melt_facePrp(currentPrp, prp)

      ! Set the melted face propery.
      if( prp_dir .eq. tem_left ) then
        faces%faceList%leftPrp%val(facePos) = newPrp
      else
        faces%faceList%rightPrp%val(facePos) = newPrp
      end if

    end if

  end subroutine tem_appendFace_prp
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Function to melt two properties together. The resulting property
  !! holds the union of firstPrp and secondPrp.
  !!
  function tem_melt_facePrp(firstPrp, secondPrp) result(meltedPrp)
    ! --------------------------------------------------------------------------
    !> The first property to be melted.
    integer, intent(in) :: firstPrp
    !> The second property to be melted.
    integer, intent(in) :: secondPrp
    !> The resulting property (the union of the first and second property)
    integer :: meltedPrp
    ! --------------------------------------------------------------------------

    ! Apply a bitwise OR operation to melt the two properties.
    meltedPrp = ior(firstPrp, secondPrp)

  end function tem_melt_facePrp
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> summary: Returns the property of a given element in the level descriptor
  !! which can be assigned to one of its face properties.
  !!
  function tem_get_elemPrp(levelDesc, elemPos) result(elemPrp)
    ! --------------------------------------------------------------------------
    !> The level descriptor that contains the investigated element.
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> The position of the element in the total list of the level descriptor.
    integer, intent(in) :: elemPos
    !> The element property of the element.
    integer :: elemPrp
    ! --------------------------------------------------------------------------

    elemPrp = 0

    ! If the element does not exist, we can only assign the non-existing
    ! property.
    if(elemPos == 0) then
      elemPrp = tem_notExist_prp
      return
    elseif(elemPos < 0) then
      elemPrp = tem_bndFace_prp
      return
    end if

    ! Check for fluid
    if( levelDesc%offset(1,1)+1 <= elemPos .and.                               &
      & elemPos<=levelDesc%offset(1,1)+levelDesc%elem%nElems( eT_fluid )  ) then
      elemPrp = tem_fluidFace_prp
    ! Check for ghost from coarser
    else if (levelDesc%offset(1,2)+1<=elemPos .and.                          &
      & elemPos<=levelDesc%offset(1,2)+levelDesc%elem%nElems( eT_ghostFromCoarser) ) then
      elemPrp = tem_fromCoarserFace_prp
    ! Check for ghost from finer
    else if (levelDesc%offset(1,3)+1<=elemPos .and.                          &
      & elemPos<=levelDesc%offset(1,3)+levelDesc%elem%nElems( eT_ghostFromFiner )) then
      elemPrp = tem_fromFinerFace_prp
    ! Check for halo
    else if (levelDesc%offset(1,4)+1<=elemPos .and.                          &
      & elemPos<=levelDesc%offset(1,4)+levelDesc%elem%nElems( eT_halo ) ) then
      elemPrp = tem_remoteFace_prp
    ! Everything else should not occur
    else
      write(*,*) 'ERROR in tem_get_elemPrp: Unknown element property, '//      &
        &        'stopping ...'
      call tem_abort()
    end if

  end function tem_get_elemPrp
  ! ************************************************************************** !


end module tem_face_module
! **************************************************************************** !
