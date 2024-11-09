! Copyright (c) 2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
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
! ****************************************************************************** !
!> author: Daniel Harlacher
!! author: Manuel Hasert
!! author: Kartik Jain
!! author: Simon Zimny
!! The implementation is based on the SPartA Algorithm.
!!
!! Some elements which require further
!! treatment expose higher computational effort. This includes among others the
!! enforcing of boundary conditions, source terms, or the interpolation between
!! levels.  Depending on the selected behavior of elements, this result in
!! considerable higher effort on certain elements.
!! When running in parallel, it is important to balance the load between the
!! participating processing units equally in order to achieve a minimal time to
!! solution. It is always desired to minimize the waiting time of processing
!! units.
!! The initial disitribution of the elements among the processes assigns the
!! roughly the same amount of elements to each partition, regardless of their
!! cost.
!! This simple partitioning yields very good results for uniform grid
!! simulations, as the elements with additional behavior is of the order of
!! \(\mathcal{O}(n^2)\) while the total number of elements are of order
!! \(\mathcal{O}(n^3)\).
!!
!! Especially when using locally refined grids, the computational cost of
!! elements varies considerably. Besides the increased update interval of the
!! smaller elements, the interpolation routines are very expensive compared to
!! the pure computation of an element.
!!
!! To account for this variation of the computational cost, a load balancing
!! algorithm (see: [2] at Treelm bibliography) is employed to homogenize
!! the total cost among  the partitions. It is based on a space-filling curve,
!! which fits nicely into the TreElm framework. Additionally, each element of
!! the tree is assigned a weight representing the computational cost. The total
!! cost is then possibly equally distributed  among the partitions which results
!! in a re-partitioning of the mesh along the space-filling curve.
!!
module tem_Sparta_module

  ! include treelm modules
  use mpi
  use env_module, only: rk, long_k, rk_mpi, long_k_mpi, double_k
  use tem_logging_module,  only: logUnit
  use tem_aux_module,      only: tem_abort
  use tem_float_module,    only: operator(.feq.)

  implicit none
  private

  type tem_sparta_type
    integer, allocatable :: send_index(:)
    integer, allocatable :: recv_index(:)
    integer, allocatable :: send_count(:)
    integer, allocatable :: recv_count(:)

    integer :: old_size = 0
    integer :: new_size = 0
  end type tem_sparta_type

  interface tem_exchange_sparta
    module procedure tem_exchange_long
    module procedure tem_exchange_double
    module procedure tem_exchange_long2
    module procedure tem_exchange_double2
  end interface tem_exchange_sparta

  public :: tem_sparta_type
  public :: tem_balance_sparta
  public :: tem_init_sparta
  public :: tem_destroy_sparta
  public :: tem_exchange_sparta
  public :: tem_output_sparta
  public :: tem_derive_sparta

contains

! ****************************************************************************** !
  !> Return splitting positions based on the weights provided
  !! by each rank.
  !!
  !! This is the SPARTA algorithm which uses simple splitting based on
  !! given weights for all elements in the mesh.
  !!
  subroutine tem_balance_sparta(weight, myPart, nParts, comm, myElems, offset, &
    &                          sparta )
    ! ---------------------------------------------------------------------------
    !> Sorted list of weights corresponding to treeID order
    real(kind=rk),intent(in)  :: weight(:)
    integer,intent(in) :: myPart !< Rank of the calling process
    !> Number of procs the distribution should span
    integer, intent(in) :: nParts
    !> MPI Communicator
    integer, intent(in) :: comm
    !> number of elements
    integer, intent(inout) :: myElems
    !> Array of offsets with the size nParts. Offset index starts at 0.
    !! This Array needs to be allocate and deallocated outside
    integer(kind=long_k), intent(out) :: offset
    ! Count variables that state which rank gets how many elements
    ! *_count(rank0, ...)
    ! Right now the size is nParts but that will be changed in
    ! the near future to avoid O(p) allocations
    type( tem_sparta_type ), intent(inout) :: sparta
    ! ---------------------------------------------------------------------------
    integer :: iErr  ! MPI error variable
    integer :: iElem, iProc
    integer(kind=long_k) :: myElems_long
    real(kind=rk) :: w_sum, w_opt
    real(kind=rk) :: send, recv ! Send and receive buffers for MPI calls
    ! boundary values of the elements in which we search for splitters
    real(kind=rk) :: lower_boundary, upper_boundary
    ! local prefix sum array of myElems
    real(kind=rk), allocatable :: presum(:)
    integer :: rmin, rmax, lb, ub, left_off, mid
    real(kind=rk) :: opt_split, wsplit
    integer :: send_count(0:nParts-1)
    ! ---------------------------------------------------------------------------

    write(logUnit(5),*) "Balance by SpartA algorithm."

    send_count = 0

    ! Allocate Array for Prefix sum
    allocate(presum(myElems))
    ! Prefix sum over local weights. later on we will look for the splitter in
    ! this prefix sum
    presum(1) = weight(1)
    do iElem = 2,myElems
      presum(iElem) = presum(iElem-1) + weight(iElem)
    end do
    send = presum(myElems)
    ! sum up global total weight
    call MPI_ALLREDUCE(send, recv, 1, rk_mpi, mpi_sum, comm, iErr)

    w_sum = recv
    ! Calculate global optimum
    w_opt = w_sum / dble(nParts)

    ! Global prefix sum for weights
    call MPI_EXSCAN(send, recv, 1, rk_mpi, mpi_sum, comm, iErr)

    ! initialize splitter search
    lower_boundary=recv
    if (myPart == 0) lower_boundary = 0
    upper_boundary = lower_boundary + presum(myElems)

    rmin = max(floor(lower_boundary/w_opt),0)
    rmax = min(ceiling(upper_boundary / w_opt),nParts-1)

    ! Do splitter search
    left_off = 1
    do iProc = rmin,rmax
      lb = left_off
      ub = myelems
      opt_split = (iProc+1)*w_opt
      if (iProc*w_opt < upper_boundary) then
      do
        mid = (lb+ub)/2
        wsplit = presum(mid) + lower_boundary
        if (wsplit .feq. opt_split) exit
        if (wsplit < opt_split) then
           lb = mid
        else
           ub = mid
        end if
        ! exit if a single element was found, need to do this
        if (lb >= ub-1) exit
        !                    here, to have mid and wsplit set.
      end do
      if (ABS(wsplit - opt_split) > ABS(wsplit - opt_split - weight(mid))) then
        mid = mid - 1 ! return 0 if the splitter is left of the lower boundary
      else
        if (mid+1 <= myElems) then
          if (ABS(wsplit - opt_split)                                          &
            & > ABS(wsplit - opt_split + weight(mid+1))) then
            mid = mid + 1 ! return myElems at most
          end if
        else
          if (opt_split > upper_boundary) mid = myElems
        end if
      end if
      send_count(iProc) = mid - left_off + 1
      left_off = mid + 1
      end if
    end do
    ! finished splitter search. Communciate results.
    ! Each process needs to know how many elements to receive from which process
    call tem_set_sparta( sparta, comm, nParts, send_count )

    ! Calculate myElems and offset -----------------------------------
    ! total number of my elements after exchanging elements
    myElems = sparta%new_size
    myElems_long = int(myElems, kind=long_k)
    call mpi_exscan(myElems_long, offset, 1, long_k_mpi, mpi_sum, comm, ierr)
    if (myPart == 0) offset = 0
    ! write(*,"(3(A,I0))") 'myPart ', myPart, ' nElems: ', myElems, ' offset: ', offset
    ! Calculate myElems and offset -----------------------------------

  end subroutine tem_balance_sparta
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine tem_init_sparta( me, nParts )
    type( tem_sparta_type ), intent(inout) :: me
    integer, intent(in) :: nParts

    allocate( me%send_count(0:nParts-1) )
    allocate( me%recv_count(0:nParts-1) )
    allocate( me%send_index(0:nParts-1) )
    allocate( me%recv_index(0:nParts-1) )

  end subroutine tem_init_sparta
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine tem_exchange_long( me, val, nComponents, comm )
    ! ---------------------------------------------------------------------------
    type( tem_sparta_type ), intent(in) :: me
    integer, intent(in) :: nComponents
    integer(kind=long_k), allocatable, intent(inout) :: val(:)
    integer, intent(in) :: comm
    ! ---------------------------------------------------------------------------
    integer(kind=long_k), allocatable :: old_val(:)
    integer :: iError
    ! ---------------------------------------------------------------------------

    ! Assumption check
    if ( nComponents <= 0 ) then
      write(logUnit(0),*) 'When call tem_exchange_data, nComponents <= 0!'
      write(logUnit(0),*) 'Stop!'
      call tem_abort()
    end if

    call move_alloc( val, old_val )
    allocate( val(me%new_size*nComponents) )

    call mpi_alltoallv( old_val, me%send_count*nComponents, me%send_index*nComponents, long_k_mpi,&
                            val, me%recv_count*nComponents, me%recv_index*nComponents, long_k_mpi,&
                        comm, ierror )

    deallocate( old_val )
  end subroutine tem_exchange_long
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine tem_exchange_double( me, val, nComponents, comm )
    ! ---------------------------------------------------------------------------
    type( tem_sparta_type ), intent(in) :: me
    integer, intent(in) :: nComponents
    real(kind=double_k), allocatable, intent(inout) :: val(:)
    integer, intent(in) :: comm
    ! ---------------------------------------------------------------------------
    real(kind=double_k), allocatable :: old_val(:)
    integer :: iError
    ! ---------------------------------------------------------------------------

    ! Assumption check
    if ( nComponents <= 0 ) then
      write(logUnit(0),*) 'When call tem_exchange_data, nComponents <= 0!'
      write(logUnit(0),*) 'Stop!'
      call tem_abort()
    end if

    call move_alloc( val, old_val )
    allocate( val(me%new_size*nComponents) )

    call mpi_alltoallv( old_val, me%send_count*nComponents, me%send_index*nComponents, mpi_double_precision,&
                            val, me%recv_count*nComponents, me%recv_index*nComponents, mpi_double_precision,&
                        comm, ierror )

    deallocate( old_val )
  end subroutine tem_exchange_double
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine tem_exchange_long2( me, val, nComponents, comm )
    ! ---------------------------------------------------------------------------
    type( tem_sparta_type ), intent(in) :: me
    integer, intent(in) :: nComponents
    integer(kind=long_k), allocatable, intent(inout) :: val(:,:)
    integer, intent(in) :: comm
    ! ---------------------------------------------------------------------------
    integer(kind=long_k), allocatable :: old_val(:,:)
    integer :: iError
    ! ---------------------------------------------------------------------------

    call move_alloc( val, old_val )
    allocate( val(nComponents, me%new_size) )

    call mpi_alltoallv( old_val, me%send_count*nComponents, me%send_index*nComponents, long_k_mpi,&
                            val, me%recv_count*nComponents, me%recv_index*nComponents, long_k_mpi,&
                        comm, ierror )

    deallocate( old_val )
  end subroutine tem_exchange_long2
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine tem_exchange_double2( me, val, nComponents, comm )
    ! ---------------------------------------------------------------------------
    type( tem_sparta_type ), intent(in) :: me
    integer, intent(in) :: nComponents
    real(kind=double_k), allocatable, intent(inout) :: val(:,:)
    integer, intent(in) :: comm
    ! ---------------------------------------------------------------------------
    real(kind=double_k), allocatable :: old_val(:,:)
    integer :: iError
    ! ---------------------------------------------------------------------------

    call move_alloc( val, old_val )
    allocate( val(nComponents, me%new_size) )

    call mpi_alltoallv( old_val, me%send_count*nComponents, me%send_index*nComponents, mpi_double_precision,&
                            val, me%recv_count*nComponents, me%recv_index*nComponents, mpi_double_precision,&
                        comm, ierror )

    deallocate( old_val )
  end subroutine tem_exchange_double2
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine tem_destroy_sparta( me )
    type( tem_sparta_type ), intent(inout) :: me

    deallocate( me%send_count )
    deallocate( me%recv_count )
    deallocate( me%send_index )
    deallocate( me%recv_index )
  end subroutine tem_destroy_sparta
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine tem_output_sparta( me, outUnit )
    type( tem_sparta_type ), intent(in) :: me
    integer, intent(in) :: outUnit

    if ( allocated( me%send_count ) ) then
      write(outUnit, "(A,I0)") "Old size: ", me%old_size
      write(outUnit, "(A,I0)") "New size: ", me%new_size
      write(outUnit, "(A,I0)") "Size of send_count: ", size(me%send_count)
      write(outUnit, *) "Send_index: ", me%send_index(:)
      write(outUnit, *) "Send_count: ", me%send_count(:)
      write(outUnit, *) "Recv_index: ", me%recv_index(:)
      write(outUnit, *) "Recv_count: ", me%recv_count(:)
    else
      write(outUnit, *) "Data is NOT allocated!"
    end if

  end subroutine tem_output_sparta
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine tem_derive_sparta( origin, derived, nElems, elemPropertyBits, prpBit, &
    &                           comm, nParts )
    type( tem_sparta_type ), intent(in)  :: origin
    type( tem_sparta_type ), intent(inout) :: derived
    integer, intent(in) :: nElems
    integer(kind=long_k), intent(in) :: elemPropertyBits(nElems)
    integer, intent(in) :: prpBit, comm, nParts

    integer :: send_count(0:nParts-1), iPart, offset, iElem

    do iPart = 0, nParts - 1
      offset = origin%send_index(iPart)
      send_count(iPart) = 0
      do iElem = 1, origin%send_count(iPart)
        if (btest( elemPropertyBits( iElem+offset ), prpBit ) ) then
          send_count(iPart) = send_count(iPart) + 1
        end if
      end do ! iElem
    end do

    call tem_set_sparta( derived, comm, nParts, send_count )

  end subroutine tem_derive_sparta
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine tem_set_sparta( me, comm, nParts, send_count )
    type( tem_sparta_type ), intent(inout) :: me
    integer, intent(in) :: comm, nParts
    integer, intent(in) :: send_count(0:nParts-1)

    integer :: iProc, iErr

    ! Each process needs to know how many elements to receive from which process
    call mpi_alltoall(send_count, 1, mpi_integer, &
      &               me%recv_count, 1, mpi_integer,  &
      &               comm, iErr)

    ! Create additional structures for alltoallv ------------------------------
    me%send_index(0) = 0
    me%recv_index(0) = 0
    me%send_count(0) = send_count(0)
    do iProc = 1,nParts-1
      me%send_count(iProc) = send_count(iProc)
      me%send_index(iProc) = me%send_index(iProc-1) + send_count(iProc-1)
      me%recv_index(iProc) = me%recv_index(iProc-1) + me%recv_count(iProc-1)
    end do
    ! Create additional structures for alltoallv ------------------------------

    me%old_size = sum(me%send_count)
    me%new_size = sum(me%recv_count)

  end subroutine tem_set_sparta
! ****************************************************************************** !

end module tem_sparta_module
! ****************************************************************************** !
