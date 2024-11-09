! Copyright (c) 2013 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2013-2014, 2022 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2013 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013-2014, 2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2022 Jana Gericke <jana.gericke@dlr.de>
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
!> author: Kartik Jain
!! This module handles the various status bits and performs the relevant
!! communications. For example if one process ends, this information is stored
!! in a bit and communicated to all other processors.
!!
module tem_status_module

  use mpi
  use tem_comm_env_module, only: tem_comm_env_type

  implicit none

  private

  public :: tem_status_type
  public :: tem_stat_nFlags
  public :: tem_stat_interval, tem_stat_run_terminate, tem_stat_steady_state
  public :: tem_stat_max_sim, tem_stat_max_iter, tem_stat_max_clock
  public :: tem_stat_stop_file
  public :: tem_stat_nan_Detected
  public :: tem_stat_nonPhysical
  public :: tem_stat_global_error

  public :: tem_status_dump
  public :: tem_status_clear
  public :: tem_status_run_end
  public :: tem_status_run_terminate
  public :: tem_status_communicate
  public :: tem_status_communicate_delayed

  !> Number of available status flags.
  integer, parameter :: tem_stat_nFlags = 10

  ! Definition of the various status flags
  ! Remember to increase the counter of nFlags above if you add a new one.

  !> Indicator for the overall steering interval to being reached.
  integer, parameter :: tem_stat_interval      =  1

  !> A steady state condition has been reached.
  integer, parameter :: tem_stat_steady_state  =  2

  !> The application is to terminate abnormally.
  integer, parameter :: tem_stat_run_terminate =  3

  !> Some global error occurred.
  integer, parameter :: tem_stat_global_error  =  4

  !> A NaN was detected during the computation.
  integer, parameter :: tem_stat_nan_detected  =  5

  !> A non-physical state was detected during the computation.
  integer, parameter :: tem_stat_nonPhysical   =  6

  !> The maximal simulation time has been reached.
  integer, parameter :: tem_stat_max_sim       =  7

  !> The maximal number of iterations has been reached.
  integer, parameter :: tem_stat_max_iter      =  8

  !> The maximall wall-clock time was reached.
  integer, parameter :: tem_stat_max_clock     =  9

  !> A stop file was encountered.
  integer, parameter :: tem_stat_stop_file     = 10


  !> Define an array to hold all status flags
  type tem_status_type
    integer :: check_request = MPI_REQUEST_NULL
    logical :: bits(tem_stat_nFlags) = .false.
    logical :: oldbits(tem_stat_nFlags) = .false.
  end type


contains


  ! **************************************************************************** !
  !> Clear (unset) all status bits.
  subroutine tem_status_clear(me)
    ! --------------------------------------------------------------------------!
    !> Status type to initialize.
    type(tem_status_type), intent(inout) :: me
    ! --------------------------------------------------------------------------!

    me%bits = .false.

  end subroutine tem_status_clear
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Write the status in me to outUnit.
  !!
  !! Status will only be written if any of the flags is true.
  subroutine tem_status_dump(me, outUnit)
    ! --------------------------------------------------------------------------!
    !> Status to write on outunit.
    type(tem_status_type), intent(in) :: me

    !> The file unit to write to.
    integer, intent(in) :: outUnit
    ! --------------------------------------------------------------------------!

    if (any(me%bits)) then
      write(outUnit,*) '+--------------------------------------------+'
      write(outUnit,*) 'STATUS:'
      if (me%bits(tem_stat_interval)) &
        &  write(outUnit,*) ' * Status Interval is reached'

      if (me%bits(tem_stat_run_terminate)) &
        &  write(outUnit,*) ' * Simulation is to abnormally terminate'

      if (me%bits(tem_stat_steady_state)) &
        &  write(outUnit,*) ' * Simulation reached a steady state'

      if (me%bits(tem_stat_global_error)) &
        &  write(outUnit,*) ' * Encountered a global error'

      if (me%bits(tem_stat_nan_detected)) &
        &  write(outUnit,*) ' * Detected a NaN during computation'

      if (me%bits(tem_stat_nonPhysical)) &
        &  write(outUnit,*) ' * Ran into a non-physical state'

      if (me%bits(tem_stat_max_sim)) &
        &  write(outUnit,*) ' * Reached maximal simulation time'

      if (me%bits(tem_stat_max_iter)) &
        &  write(outUnit,*) ' * Reached maximal number of iterations'

      if (me%bits(tem_stat_max_clock)) &
        &  write(outUnit,*) ' * Reached maximal wall clock running time'

      if (me%bits(tem_stat_stop_file)) &
        &  write(outUnit,*) ' * Found a stop file'

      write(outUnit,*) '+--------------------------------------------+'
      write(outUnit,*)
    end if

  end subroutine tem_status_dump
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Decide if the simulation run should end based on the status flags.
  function tem_status_run_end(me) result(run_end)
    ! --------------------------------------------------------------------------!
    !> Status to check if the run has to end.
    type(tem_status_type), intent(in) :: me

    !> Result indicating if the run should come to a regular end.
    logical :: run_end
    ! --------------------------------------------------------------------------!

    run_end = (     me%bits(tem_stat_max_sim)      &
      &        .or. me%bits(tem_stat_max_iter)     &
      &        .or. me%bits(tem_stat_max_clock)    &
      &        .or. me%bits(tem_stat_stop_file)    &
      &        .or. me%bits(tem_stat_steady_state) )

  end function tem_status_run_end
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Decide if the simulation run should abnormally terminate based on the
  !! status flags.
  function tem_status_run_terminate(me) result(run_terminate)
    ! --------------------------------------------------------------------------!
    !> Status to check for an extraordinary termination of the run.
    type(tem_status_type), intent(in) :: me

    !> Result indicating if the run should com to an irregular end.
    logical :: run_terminate
    ! --------------------------------------------------------------------------!

    run_terminate = (     me%bits(tem_stat_run_terminate) &
      &              .or. me%bits(tem_stat_global_error)  &
      &              .or. me%bits(tem_stat_nan_detected)  &
      &              .or. me%bits(tem_stat_nonPhysical)   )

  end function tem_status_run_terminate
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Perform the communication of status bits
  subroutine tem_status_communicate(me, comm)
    ! --------------------------------------------------------------------------!
    !> Status to communicate.
    type(tem_status_type), intent(inout) :: me

    !> Communicator to use for the MPI reduction operation.
    integer, intent(in) :: comm
    ! --------------------------------------------------------------------------!
    integer :: iError
    logical :: local_bits(tem_stat_nFlags)
    ! --------------------------------------------------------------------------!

    local_bits = me%bits
    call mpi_allreduce( local_bits, me%bits, tem_stat_nFlags, MPI_LOGICAL, &
      &                 MPI_LOR, comm, iError                              )

  end subroutine tem_status_communicate
  ! **************************************************************************** !

  ! **************************************************************************** !
  !> Perform the communication of status bits with a nonblocking allreduce
  !! resulting in an delayed communication by one check_iter interval.
  subroutine tem_status_communicate_delayed(me, comm)
    ! --------------------------------------------------------------------------!
    !!> Status to communicate.
    type(tem_status_type), intent(inout) :: me

    !> Communicator to use for the MPI reduction operation.
    integer, intent(in) :: comm
    ! --------------------------------------------------------------------------!
    integer :: iError
    integer :: sync_status(MPI_STATUS_SIZE)
    logical :: local_bits(tem_stat_nFlags)
    ! --------------------------------------------------------------------------!

    local_bits = me%bits

    if (me%check_request /= MPI_REQUEST_NULL) then
      ! Wait on previous Iallreduce, to synchronize status bits.
      call MPI_Wait(me%check_request, sync_status, iError)
      me%bits = me%oldbits
    end if

    me%oldbits = local_bits

    ! Synchronize status from the current iteration (may only become known on
    ! other processes in next iteration, after the wait above)
    call MPI_Iallreduce( MPI_IN_PLACE, me%oldbits, tem_stat_nFlags,           &
      &                  MPI_LOGICAL, MPI_LOR, comm, me%check_request, iError )

  end subroutine tem_status_communicate_delayed
  ! **************************************************************************** !


end module tem_status_module
! ****************************************************************************** !
