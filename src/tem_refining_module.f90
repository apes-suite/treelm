! Copyright (c) 2015-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Daniel Fleischer <daniel.fleischer@student.uni-siegen.de>
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
!> Module that allows the refinement of a given tree to create a mesh with a
!! different resolution.
!!
!! IMPORTANT NOTE: This is a module, that does not employ mesh adaptivity
!!                 strategies, instead, complete meshes are created anew!
!!                 Smarter, more efficient mesh adaptivity procedures are being
!!                 implemented, and should replace the functionality of this
!!                 module eventually.
module tem_refining_module
  use mpi

  use env_module, only: long_k, long_k_mpi
  use treelmesh_module, only: treelmesh_type
  use tem_bc_prop_module, only: tem_BC_prop_type
  use tem_geometry_module, only: tem_eligibleChildren
  use tem_grow_array_module, only: grw_long2dArray_type, init, append, destroy
  use tem_property_module, only: prp_hasbnd
  use tem_subtree_type_module, only: tem_subtree_type

  implicit none


contains


  ! ------------------------------------------------------------------------ !
  !> Refine all elements defined in subtree by one level in the original mesh,
  !! and create a new mesh.
  !!
  !! Orig_mesh needs to be a properly defined treelmesh, and subtree an
  !! accompanying subtree in that mesh.
  !! orig_bcs needs to be the boundary conditions accompanying the original
  !! mesh.
  !! New_mesh will be the resulting mesh after refining the elements at elempos.
  !! New_bcs will be the resulting boundary conditions for that new mesh.
  !! Per default, the elements from the orig_mesh that are not covered by the
  !! subtree are still part of new_mesh, just not refined.
  !! To change this behavior and include only the refined elements in the new
  !! mesh, set restrict_to_sub to true.
  !!
  !! The boundary properties are the only ones, that will be inherited to the
  !! new mesh.
  !! All other properties will get lost!!!
  subroutine tem_refine_global_subtree( orig_mesh, orig_bcs, subtree, ndims, &
    &                                   new_mesh, new_bcs, restrict_to_sub   )
    !> The original mesh to be refined.
    type(treelmesh_type), intent(in) :: orig_mesh

    !> Boundary conditions for the original mesh.
    type(tem_BC_prop_type), intent(in) :: orig_bcs

    !> (Process local) positions of elements to refine.
    type(tem_subtree_type), intent(in) :: subtree

    !> Number of dimensions for the refinement
    !!
    !! The dimensionality can restrict the elements to create.
    integer, intent(in) :: ndims

    !> Newly created refined mesh.
    type(treelmesh_type), intent(out) :: new_mesh

    !> Boundary conditions for the new mesh.
    type(tem_BC_prop_type), intent(out) :: new_bcs

    !> Flag to indicate, wether only refined elements should be put into the
    !! new mesh.
    logical, optional, intent(in) :: restrict_to_sub
    ! -------------------------------------------------------------------- !
    logical :: restrict
    logical :: refine_all
    integer :: iParent
    integer :: nParents
    integer :: parentpos
    integer :: childpos
    integer :: newpos
    integer :: subpos
    integer :: iError
    integer :: bcprop
    integer :: iBC_elem
    integer :: iElem
    integer :: iChild
    integer :: iProp
    integer :: iSide
    integer(kind=long_k) :: nNewElems
    integer(kind=long_k) :: nBCElems
    integer(kind=long_k) :: childID_off
    integer(kind=long_k) :: child_bcid(orig_bcs%nSides,2**ndims)
    integer, allocatable :: bc_child(:)
    type(grw_long2dArray_type) :: newbcid
    logical :: has_boundary
    ! -------------------------------------------------------------------- !

    refine_all = subtree%useGlobalMesh

    if (present(restrict_to_sub)) then
      restrict = restrict_to_sub
    else
      restrict = .false.
    end if

    nParents = subtree%nElems

    ! (local) Number of new elements:
    if (restrict) then
      nNewElems = nParents*(2**ndims)
    else
      nNewElems = nParents*(2**ndims-1) + orig_mesh%nElems
    end if

    new_mesh%nElems = int(nNewElems)

    new_mesh%global = orig_mesh%global

    ! Need to create new properties for the new mesh
    nullify(new_mesh%global%Property)
    nullify(new_mesh%Property)

    ! Find the boundary property
    bcprop = 0
    do iProp=1,orig_mesh%global%nProperties
      if (orig_mesh%global%Property(iProp)%bitpos == prp_hasBnd) then
        bcprop = iProp
        EXIT
      end if
    end do


    ! Only carry on boundary conditions
    if (bcprop > 0) then
      new_mesh%global%nProperties = 1
      new_bcs%nSides = orig_bcs%nSides
      new_bcs%nBCtypes = orig_bcs%nBCtypes
      new_bcs%BC_label = orig_bcs%BC_label
      call init( me     = newbcid,                             &
        &        width  = orig_bcs%nSides,                     &
        &        length = (ndims*2-1)*orig_bcs%property%nElems )
    else
      new_mesh%global%nProperties = 0
    end if

    ! Adapt the level range in the new mesh.
    new_mesh%global%maxlevel = new_mesh%global%maxlevel + 1
    if (refine_all .or. restrict) then
      new_mesh%global%minlevel = new_mesh%global%minlevel + 1
    end if

    ! Overall number of elements in the new mesh and offsets.
    call MPI_Exscan(nNewelems, new_mesh%ElemOffset, 1, long_k_mpi, &
      &             MPI_SUM, new_mesh%global%comm, iError          )
    ! Initialize for rank 0
    if (new_mesh%global%myPart == 0) new_mesh%ElemOffset = 0_long_k
    new_mesh%global%nElems = new_mesh%ElemOffset+nNewElems
    call MPI_Bcast( new_mesh%global%nElems, 1, long_k_mpi, &
      &             new_mesh%global%nParts-1,              &
      &             new_mesh%global%comm, iError           )

    allocate(new_mesh%treeID(nNewElems))
    allocate(new_mesh%ElemPropertyBits(nNewElems))
    allocate(new_mesh%Part_First(orig_mesh%global%nParts))
    allocate(new_mesh%Part_Last(orig_mesh%global%nParts))


    newpos = 0
    subpos = 1
    iBC_elem = 0

    ! Need to iterate through all elements of the original mesh to properly
    ! account for the boundary elements.
    origsall: do iParent=1,orig_mesh%nElems
      has_boundary = btest( orig_mesh%ElemPropertyBits(iParent), &
        &                   prp_hasbnd                           )
      if (has_boundary) iBC_elem = iBC_elem + 1
      if (refine_all) then
        parentPos = iParent
      else
        if (nParents >= subpos) then
          parentPos = subtree%map2global(subpos)
        else
          ! Element not part of the subtree!
          parentPos = -1
        end if
      end if

      dorefine: if (parentPos == iParent) then

        ! Element in subtree to be refined.

        childID_off = orig_mesh%treeID(iParent)*8_long_k
        childpos = newpos ! store current position for boundaries
        do iChild=1,2**ndims
          newpos = newpos+1
          new_mesh%treeID(newpos) = childID_off + iChild
          new_mesh%ElemPropertyBits(newpos) = 0_long_k
        end do
        if (has_boundary) then
          child_bcid = 0_long_k
          do iSide=1,2*ndims ! Only iterate over direct neighbors
            if (orig_bcs%boundary_ID(iSide, iBC_elem) > 0) then
              call tem_eligibleChildren( eligible_child = bc_child, &
                &                        direction      = iSide     )
              do iChild=1,2**(ndims-1)
                new_mesh%ElemPropertyBits(childpos+bc_child(iChild)) &
                  &  = ibset(0_long_k, prp_hasbnd)
                child_bcid(iSide, bc_child(iChild)) = orig_bcs               &
                  &                                   %boundary_ID( iSide,   &
                  &                                                 iBC_elem )
              end do
              deallocate(bc_child)
            end if
          end do

          do iChild=1,2**ndims
            if ( btest(new_mesh%ElemPropertyBits(childpos+iChild), &
              &        prp_hasBnd) ) then
              call append( me  = newbcid,             &
                &          val = child_bcid(:,iChild) )
            end if
          end do
        end if
        subpos = min(subpos + 1, subtree%nelems)

      else dorefine

        ! Element NOT in subtree to be refined.

        if (.not. restrict) then

          ! Add also non-refined elements to the new mesh.
          newpos = newpos+1
          new_mesh%treeID(newpos) = orig_mesh%treeID(iParent)
          new_mesh%ElemPropertyBits(newpos) = 0_long_k
          if (has_boundary) then
            new_mesh%ElemPropertyBits(newpos) = ibset(0_long_k, prp_hasbnd)

            call append( me  = newbcid,                         &
              &          val = orig_bcs%boundary_ID(:,iBC_elem) )
          end if
        end if

      end if dorefine

    end do origsall

    call MPI_Allgather( new_mesh%treeID(1), 1, long_k_mpi,  &
      &                 new_mesh%Part_First, 1, long_k_mpi, &
      &                 new_mesh%global%comm, iError        )

    call MPI_Allgather( new_mesh%treeID(new_mesh%nElems), 1, long_k_mpi, &
      &                 new_mesh%Part_Last, 1, long_k_mpi,               &
      &                 new_mesh%global%comm, iError                     )

    ! Set boundary property, if needed.
    if (bcprop > 0) then
      allocate(new_mesh%Property(1))
      allocate(new_mesh%global%Property(1))
      new_mesh%global%property(1)%label &
        &  = orig_mesh%global%Property(bcprop)%label
      new_mesh%global%property(1)%bitpos &
        &  = orig_mesh%global%Property(bcprop)%bitpos
      new_mesh%Property(1)%nElems = newbcid%nVals
      nBCElems = int(newbcid%nVals, kind=long_k)
      allocate(new_mesh%Property(1)%ElemID(newbcid%nVals))
      iBC_elem = 0
      do iElem=1,int(nNewElems)
        if ( btest(new_mesh%ElemPropertyBits(iElem), prp_hasBnd) ) then
          iBC_elem = iBC_elem + 1
          new_mesh%Property(1)%ElemID(iBC_elem) = iElem
        end if
      end do
      allocate(new_bcs%boundary_ID(new_bcs%nSides, newbcid%nVals))
      new_bcs%property => new_mesh%Property(1)
      new_bcs%header => new_mesh%global%Property(1)
      new_bcs%boundary_ID = newbcid%val(:,:newbcid%nVals)
      call destroy(newbcid)
      call MPI_Exscan(nBCElems, new_bcs%property%Offset, 1, long_k_mpi, &
        &             MPI_SUM, new_mesh%global%comm, iError             )
      nBCElems = new_bcs%property%Offset+nBCElems
      call MPI_Bcast(nBCElems, 1, long_k_mpi, new_mesh%global%nParts-1, &
        &            new_mesh%global%comm, iError                       )
      new_bcs%header%nElems = nBCElems
    else
      if (associated(new_mesh%property)) &
        &    deallocate(new_mesh%property)
      if (associated(new_mesh%global%property)) &
        &    deallocate(new_mesh%global%property)
      allocate(new_mesh%Property(0))
      allocate(new_mesh%global%Property(0))
    end if

  end subroutine tem_refine_global_subtree
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !

end module tem_refining_module
