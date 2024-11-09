! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2014-2015 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014 Julia Moos <julia.moos@student.uni-siegen.de>
! Copyright (c) 2014 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
module tem_adaptation_module

  ! include treelm modules
  use mpi
  use env_module,               only: rk, long_k, rk_mpi
  use treelmesh_module,         only: treelmesh_type, dump_treelmesh
  use tem_grow_array_module,    only: grw_longArray_type, init, append
  use tem_topology_module,      only: tem_directChildren, tem_levelOf
  use tem_stencil_module,       only: tem_stencilHeader_type,                  &
    &                                 tem_stencilElement_type
  use tem_construction_module,  only: tem_levelDesc_type
  use tem_element_module,       only: init, append, PositionOfVal, changeType, &
    &                                 eT_fluid, eT_sacrificed, eT_ghostFromFiner
  use tem_comm_env_module,      only: tem_comm_env_type
  use tem_construction_module,  only: identify_stencilNeigh
  use tem_property_module,      only: prp_chgElems
  ! use tem_adaptation_config_module, only: tem_adapt_type


  implicit none

  private

  public :: tem_adapt_dump_newMesh

  contains

! *******************************************************************************
  !> This routine prepares the ground work for dumping the adapted mesh to
  !! disk. The new treeIDs which were created while adaptive refinement
  !! are sorted within the levelDescriptor elem type, and then passed to the
  !! dump_treelmesh routine for dumping.
  !!
  subroutine tem_adapt_dump_newMesh(levelDesc, tree, proc )
    ! ---------------------------------------------------------------------------
    type(tem_levelDesc_type), intent(inout) :: levelDesc(:)
    type(treelmesh_type)    , intent(inout) :: tree
    type(tem_comm_env_type) , intent(in)    :: proc
    ! ---------------------------------------------------------------------------
    ! A temporary list of TreeIDs
    type( grw_longArray_type )  :: loc_treeID
    type( treelmesh_type )      :: newTree
    integer(kind=long_k) :: children(8)
    integer :: iElem, iChild
    integer :: elemPos, level
    integer :: iErr
    ! ---------------------------------------------------------------------------
    call init( me     = loc_treeID, length = tree%nElems )

    ! We loop over all the elements of OLD/ORIGINAL mesh
    ! this nElems will be updated later and has to be stored
    ! for this routine
    ! If an element was sacrificed, its 8 children are appended
    ! at its original position in the tree to preserve the space filling curve
    do iElem = 1, tree%nElems
      level   = tem_levelOf(tree%treeID(iElem) )
      elemPos = PositionOfVal( me  = levelDesc(level)%elem%tID,                &
        &                      val = tree%treeID(iElem) )
      if( .not. levelDesc(level)%elem%property%val(elemPos) == prp_chgElems ) then
        call append( me  = loc_treeID,           &
                     val = tree%treeID(iElem) )
      else
        children = tem_directChildren( tree%treeID(iElem) )
        do iChild = 1, 8
          call append( me  = loc_treeID,          &
            &          val = children(iChild) )
        end do
      end if
    end do
! @todo: The children need to be added recursively, and a separate routine
!! will be written for that
! @todo: A separate routine for coarsening needs to be implemented

    ! Initialize the newTree and inherit some data from old
    allocate( newTree%treeID          (loc_treeID%nVals) )
    allocate( newTree%Part_First      (loc_treeID%nVals) )
    allocate( newTree%Part_Last       (loc_treeID%nVals) )
    allocate( newTree%ElemPropertyBits(loc_treeID%nVals) )
    allocate( newTree%pathList        (loc_treeID%nVals) )

    newTree%treeID = loc_treeID%val
    newTree%nElems = loc_treeID%nVals
    newTree%global = tree%global

    ! Now dump the mesh to disk
    !KJ: This needs to be re-thought as where and when to do it
    call MPI_ALLREDUCE( loc_treeID%nVals, newTree%global%nElems, 1, MPI_INTEGER, &
      &                 MPI_SUM, proc%comm, iErr )
    call dump_treelmesh( me = newTree )

  end subroutine tem_adapt_dump_newMesh


! *******************************************************************************



end module tem_adaptation_module
