! Copyright (c) 2011-2014 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2014 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
!> This module provides basic information on the parallel environment
module tem_comm_env_module

  !$ use omp_lib

  use mpi

  implicit none

  private

  public :: tem_comm_env_type
  public :: tem_comm_env_init, tem_comm_env_fin
  public :: tem_comm_env_init_empty

  !> Information about parallel runs
  type tem_comm_env_type
    !> size of MPI communicator
    integer :: comm_size
    !> MPI rank
    integer :: rank
    !> MPI root rank
    integer :: root
    !> MPI communicator
    integer :: comm
    !> Maximal Number of OpenMP threads
    integer :: nThreads

    !> Whether this process is the root
    logical :: isRoot

  end type tem_comm_env_type

  contains

! ****************************************************************************** !
  !> Initialize the environment. This routine is called by tem_start
  !! which should be the very first action in a program.
  !!
  subroutine tem_comm_env_init( proc, comm )
    ! --------------------------------------------------------------------------=
    !> The process communicator type
    type( tem_comm_env_type ) :: proc
    !> mpi communicator if it is predefined as in apesmate
    integer, intent(in), optional :: comm
    ! ---------------------------------------------------------------------------
    !> Error flag
    integer :: iError
    ! --------------------------------------------------------------------------

    proc%nThreads = 1
    !$ proc%nThreads = omp_get_max_threads()

    ! Init MPI rank, size and root
    ! if communicator is predefiend and passed use that one
    ! else default to mpi_comm_world
    if(present(comm)) then
      proc%comm = comm
    else
      !KM: with this call proc%comm must be freed by mpi_comm_free
      !so directly setting proc%comm = MPI_COMM_WORLD
      !call mpi_comm_dup( mpi_comm_world, proc%comm, iError )
      proc%comm = mpi_comm_world
    endif

    call mpi_comm_rank( proc%comm, proc%rank, iError )
    call mpi_comm_size( proc%comm, proc%comm_size, iError )
    proc%root = 0

    if ( proc%rank == proc%root ) then
      proc%isRoot = .true.
    else
      proc%isRoot = .false.
    end if

  end subroutine tem_comm_env_init
! ****************************************************************************** !


! ****************************************************************************** !
  !> Initialize a debug environment for just one process without envoking MPI
  !!
  subroutine tem_comm_env_init_empty( proc )
    ! ---------------------------------------------------------------------------
    !> The process communicator type
    type( tem_comm_env_type ) :: proc
    ! ---------------------------------------------------------------------------

    ! Init MPI rank, size and root
    proc%comm = 0
    proc%rank = 0
    proc%comm_size = 0
    proc%root = 0
    proc%isRoot = .false.

  end subroutine tem_comm_env_init_empty
! ****************************************************************************** !


! ****************************************************************************** !
  !> Finalize the environment. This routine is called by tem_finalize
  !! which should be the very last action in a program.
  !!
  subroutine tem_comm_env_fin( proc )
    ! ---------------------------------------------------------------------------
    !> The process communicator type
    type( tem_comm_env_type ) :: proc
    ! ---------------------------------------------------------------------------
    !> Error flag
    integer :: iError
    ! ---------------------------------------------------------------------------

    ! Free the communicator again
    call mpi_comm_free( proc%comm, iError )

  end subroutine tem_comm_env_fin
! ****************************************************************************** !

end module tem_comm_env_module
! ****************************************************************************** !
