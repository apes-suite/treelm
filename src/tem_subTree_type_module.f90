! Copyright (c) 2013-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013, 2016, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2016, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
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
!> author: Simon Zimny
!! Module for the subTree datatype.
!!
!! The subTree datatype consists of an array of pointers (map2global) to the
!! list of treeIDs of an actual tree in treelmesh format, global information on
!! the subTree,
!!
!!
module tem_subTree_type_module

  ! include treelm modules
  use mpi
  use env_module,          only: long_k, rk
  use tem_global_module,   only: tem_global_type
  use tem_property_module, only: tem_property_type
  use treelmesh_module,    only: treelmesh_type, dump_treelmesh
  use tem_aux_module,      only: tem_abort
  use tem_logging_module,  only: logUnit

  implicit none

  private

  public :: tem_subTree_type
  public :: tem_destroy_subTree
  public :: tem_local_subTree_from
  public :: tem_dump_subTree
  public :: tem_treeIDfrom_subTree

  !> Datatype for sub trees incl. local and global information as well as an
  !! array of pointers to the global list of treeIDs
  type tem_subTree_type
    !> This entry provides some global informations on the mesh, which is needed
    !! across all partitions.
    type(tem_global_type) :: global

    !> Total number of elements on this partition
    integer :: nElems

    !> Total number of points on this partition
    integer :: nPoints

    !> Total number of points on this subTree
    integer(kind=long_k) :: glob_nPoints

    !> (n, 3) array to store coordinate points for tracking
    real(kind=rk), allocatable :: points(:,:)

    !> does the subTree correspond with the corresponding one (valid for
    !! shape='all')
    logical :: useGlobalMesh = .False.

    !> Position of the treeID in the corresponding tree on the partition only
    !! set for useGlobalMesh == False and useLocalShape == False
    integer, allocatable :: map2global(:)

    !> Does the subTree correspond local partition?
    logical :: useLocalMesh = .False.

    !> Flag to indicate, whether a new communicator was created for this
    !! subtree.
    logical :: created_new_comm = .false.

    !> Actual array of treeIDs in the subTree. This is only defined for
    !! useGlobalMesh = .False. and logical :: useLocalMesh = .True.
    !! set for useGlobalMesh == False and useLocalShape == False
    integer(kind=long_k), allocatable :: treeID(:)

    !> array of bc_ids of boundary labels defined in tracking%shape boundary
    integer, allocatable :: bc_ID(:)

    !> Offset of this partition
    integer(kind=long_k) :: elemOffset

    !> Bit field of element properties for each element (the array has a length
    !! of nElems), each bit in the 8 byte integer indicates the presence or
    !! absence of a given property, which are further described in the
    !! Property component.
    integer(kind=long_k), allocatable :: ElemPropertyBits(:)

    !> declaration of additional elemental properties, the array has a length of
    !! global%nProperties, each property provides a list of local element
    !! indices of elements with this property.
    type(tem_property_type), pointer :: Property(:) => null()

  end type

contains

! ****************************************************************************** !
  !> This subroutine creates a subTree based on a provided map or list of
  !! treeids (in case a local shape is used) to the corresponding tree.
  !! This routine shall be called from a single process only, it is not intended
  !! to communicate and broadcast its information.
  !! requires the setting of me%global before hand.
  !!
  subroutine tem_local_subTree_from( me, map2global, treeID )
    ! ---------------------------------------------------------------------------
    !> subTree to be created from list of elements (map2global)
    type(tem_subTree_type), intent(inout)      :: me
    !> position of the treeID in the global treeID list
    integer, optional, intent(in)              :: map2global(:)
    !> list of treeIDs only use this in case a local shape is set
    integer(kind=long_k), optional, intent(in) :: treeID(:)
    ! ---------------------------------------------------------------------------
    integer :: nElemsList
    ! ---------------------------------------------------------------------------

    nElemsList = -9999
    ! determine number of elements in the provided map2global
    if(present(map2global) .and. .not. present(treeID))then
      nElemsList  = size(map2global)
      allocate(me%map2global( nElemsList ))
      ! copy the map2global
      me%map2global = map2global
    else if(present(treeID) .and. .not. present(map2global))then
      nElemsList  = size(treeID)
      allocate(me%treeID( nElemsList ))
      ! copy the list of treeID
      me%treeID = treeID
      ! in case a treeID is set we have to set uselocalmesh to true
      me%useLocalMesh = .true.
    else
      write(logUnit(1),*)'error:either none or both of treeID or map2global '
      write(logUnit(1),*)'      are passed to tem_subTree_from. '//            &
        &                'only one of them is allowed'
      call tem_abort()
    end if

    ! allocate the different arrays in the subTree type
    allocate(me%Property(me%global%nproperties))
    allocate(me%elemPropertyBits(nElemsList))

    ! assign treeID list to the new mesh. sort first
    me%nElems = nElemsList
    me%global%nElems = nElemsList

    ! reset properties. not needed
    me%elempropertybits = 0

  end subroutine tem_local_subTree_from
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine creates a new mesh in treelmesh format from a global and a
  !! sub tree and dumps it to disc in treelmesh format (needed for tracking in
  !! harvester format where the mesh has to be dumped to disc)
  !!
  subroutine tem_dump_subTree( me, globalTree, root_only )
    ! ---------------------------------------------------------------------------
    !> subTree to dump to disk
    type(tem_subTree_type), intent(in) :: me
    !> the global tree
    type(treelmesh_type), intent(in) :: globalTree
    !> root dump global mesh when true and
    !! all process dump its own mesh when false
    logical, intent(in), optional :: root_only
    ! ---------------------------------------------------------------------------
    type(treelmesh_type) :: local_tree
    integer :: iElem
    ! ---------------------------------------------------------------------------
    local_tree%global     = me%global
    local_tree%nElems     = me%nElems
    local_tree%elemOffset = me%elemOffset
    allocate(local_tree%treeID( me%nElems ))
    if (me%useLocalMesh) then
      local_tree%treeID = me%treeID
    else if (me%useGlobalMesh) then
      local_tree%treeID = globalTree%treeID
    else
      do iElem = 1, me%nElems
        local_tree%treeID( iElem ) = globalTree%treeID( me%map2global( iElem ))
      end do
    endif

    allocate( local_tree%ElemPropertyBits( size( me%ElemPropertyBits )))
    local_tree%ElemPropertyBits = me%ElemPropertyBits
    if (associated(me%property)) then
      allocate( local_tree%Property( size( me%Property )))
      local_tree%Property = me%Property
    end if

    call dump_treelmesh( me = local_tree, root_only = root_only )

  end subroutine tem_dump_subTree
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine derives all treeIDs of a subTree from the corresponding
  !! global tree and stores them in treeID.
  !!
  subroutine tem_treeIDfrom_subTree( me, glob_tree, treeID, bound )
    ! ---------------------------------------------------------------------------
    !> subTree of glob_tree
    type(tem_subTree_type), intent(in) :: me
    !> corresponding global tree
    type(treelmesh_type), intent(in) :: glob_tree
    !> complete treeID list
    integer(kind=long_k), intent(out), allocatable :: treeID(:)
    !> lower(1) and upper(2) bound of the element array
    integer, intent(in) :: bound(2)
    ! ---------------------------------------------------------------------------

    ! allocate the array of treeIDs
    allocate( treeID( bound(2)-bound(1)+1 ))
    ! set the number of elements in the subTree
    if( me%useGlobalMesh)then
      treeID = glob_tree%treeID( bound(1) : bound(2) )
    else if( me%useLocalMesh)then
      treeID = me%treeID( bound(1) : bound(2) )
    else
      treeID = glob_tree%treeID( me%map2global( bound(1) : bound(2) ) )
    end if

  end subroutine tem_treeIDfrom_subTree
! ****************************************************************************** !


  ! ------------------------------------------------------------------------ !
  !> This subroutine frees the resources used for a subtree.
  !!
  subroutine tem_destroy_subTree( me )
    ! -------------------------------------------------------------------- !
    !> subTree of glob_tree
    type(tem_subTree_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !
    integer :: iError
    ! -------------------------------------------------------------------- !

    ! deallocate the lists
    if (allocated(me%treeID)) deallocate(me%treeID)
    if (allocated(me%map2global)) deallocate(me%map2global)
    if (allocated(me%ElemPropertyBits)) deallocate(me%ElemPropertyBits)
    if (me%created_new_comm) then
      call MPI_Comm_free(me%global%comm, iError)
    end if
    me%created_new_comm = .false.
    me%useLocalMesh = .False.
    !me%useGlobalMesh = .False.

  end subroutine tem_destroy_subTree
  ! ------------------------------------------------------------------------ !


end module tem_subTree_type_module
! ****************************************************************************** !
