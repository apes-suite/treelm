! Copyright (c) 2015,2021 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2015 Verena Krupp <verena.krupp@uni-siegen.de>
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
!> Module to track various relevant memory information during runtime.
!!
!! The values tracked are inspired by the suggestions on
!! the QT [documentation](https://doc.qt.io/archives/qtextended4.4/syscust-oom.html).
!!
!! Where it says:
!!
!! The purpose of the Main Memory Test is to decide how often to run the Memory
!! Monitor. The test involves reading certain values from /proc/meminfo and
!! summing them to estimate the amount of main memory available. The values
!! used are MemFree, Buffers, and Cached.
!! When the sum of these values drops below a certain percentage of the main
!! memory in the device (MemTotal), the timer that governs how often the Memory
!! Monitor should run is reset with a short timeout interval so the Memory
!! Monitor runs more often.
!! Otherwise, the timer is reset with a long timeout interval so the Memory
!! Monitor runs less often.
!!
!! Note that the sum of the values of MemFree, Buffers, and Cached is not the
!! amount of virtual memory that is actually available, i.e. it cannot be used
!! to determine the memory state. In practice, the amount of memory available is
!! much greater than this sum, because many of the memory pages consumed by a
!! process are from code space, and a code space page can become free as soon as
!! program execution leaves the page.
!!
!! The page fault test involves monitoring the number of major page faults that
!! occur in each timeout interval (long or short), normalizing them, and
!! averaging them over a specified number of timeout intervals (test cycle).
!! A major page fault is one that cannot be satisfied by a page already loaded
!! in main memory. It requires loading the requested page from secondary memory.
!! We expect that an abnormally high average number of major page faults per
!! test cycle implies that too many pages of main memory are occupied by data
!! pages of too many processes. These processes are competing for the few
!! remaining code space pages available, and page thrashing is the result.
!! Since Linux processes cannot be forced to relinquish data pages they no
!! longer need, an abnormally high average number of major page faults per test
!! cycle is a reliable sign that the system is running low on memory.
!!
!! The page fault test is run in three steps. Step one computes the number of
!! major page faults that occurred during the last timeout interval. The data
!! for this computation is read from the per process virtual memory status files
!! in /proc. The number of page faults is then stored in a circular list that
!! represents the last n timeout intervals, i.e. one test cycle. Step two
!! computes the average number of page faults over all the entries in the test
!! cycle.
module tem_trackmem_module
  use mpi

  use aotus_module, only: flu_state, aot_get_val

  use env_module, only: my_status_int, my_status_int_vec
  use tem_aux_module, only: tem_open
  use tem_logging_module, only: logunit

  implicit none


contains


  !> Write the current memory status into the memfile.
  !!
  !! The file will be opened and closed, so this might be slow.
  !! Only the root process writes to this file, but the data is gathered from
  !! all processes.
  !! The output will be prepended by the current date.
  !! We will track the VmHWM: from /proc/self/status,
  !! MemFree:, Buffers: and Cached: from /proc/meminfo
  !! and pgmajfault from /proc/vmstat
  !! The min, max and average across all processes will be recorded.
  subroutine tem_trackmem(memfile, iteration)
    character(len=*), intent(in) :: memfile
    integer, intent(in)          :: iteration

    integer :: myMem(5), minMem(5), maxMem(5)
    integer :: nProcs
    integer :: iError
    integer :: funit
    integer :: myRank
    logical :: fexists
    character(len=8) :: fstat
    character(len=8) :: today
    character(len=10) :: now
    real :: myScale(5)
    real :: sumMem(5)
    real :: avgMem(5)

    if (trim(memfile) /= '') then
      call MPI_Comm_Size(MPI_COMM_WORLD, nProcs, iError)

      myMem(1) = my_status_int('VmHWM:')
      myMem(2:4) = my_status_int_vec( key = ['MemFree:', 'Buffers:', &
        &                                    'Cached: '],            &
        &                             info = '/proc/meminfo'         )
      myMem(5) = my_status_int(key = 'pgmajfault', info = '/proc/vmstat')
      call MPI_Reduce(myMem, minMem, 5, MPI_INTEGER, MPI_MIN, 0, &
        &             MPI_COMM_WORLD, iError)
      call MPI_Reduce( myMem, maxMem, 5, MPI_INTEGER, MPI_MAX, 0, &
        &              MPI_COMM_WORLD, iError )

      myScale(1:4) = real(myMem(1:4))/1024.0
      myScale(5) = real(myMem(5))
      call MPI_Reduce( myScale, sumMem, 5, MPI_REAL, MPI_SUM, 0, &
        &              MPI_COMM_WORLD, iError                    )

      avgMem = sumMem / real(nProcs)

      inquire(file = trim(memfile), exist = fexists)
      if (fexists) then
        fstat = 'old'
      else
        fstat = 'new'
      end if

      call MPI_Comm_rank(MPI_COMM_WORLD, myRank, iError)
      if (myRank == 0) then
        call tem_open( file     = trim(memfile), &
          &            newunit  = funit,         &
          &            status   = trim(fstat),   &
          &            action   = 'write',       &
          &            position = 'append',      &
          &            form     = 'formatted'    )
        call date_and_time(date = today, time = now)
        write(funit,'(a)') today // ': ' // now
        write(funit, '(a,I5)') '  iteration ', iteration
        write(funit,'(a,3(1x,f11.3))') '   VmHWM   (min/max/avg):', &
          &                                  minMem(1)/1024.0, &
          &                                  maxMem(1)/1024.0, avgMem(1)
        write(funit,'(a,3(1x,f11.3))') '   MemFree (min/max/avg):', &
          &                                  minMem(2)/1024.0, &
          &                                  maxMem(2)/1024.0, avgMem(2)
        write(funit,'(a,3(1x,f11.3))') '   Buffers (min/max/avg):', &
          &                                  minMem(3)/1024.0, &
          &                                  maxMem(3)/1024.0, avgMem(3)
        write(funit,'(a,3(1x,f11.3))') '   Cached  (min/max/avg):', &
          &                                  minMem(4)/1024.0, &
          &                                  maxMem(4)/1024.0, avgMem(4)
        write(funit,'(a,i0,1x,i0,1x,f11.3)') '   pgmajflt(min/max/avg):', &
          &                                  minMem(5), maxMem(5), avgMem(5)
        write(funit,'(a)') ''
        close(funit)
      end if

    end if

  end subroutine tem_trackmem
  ! **************************************************************************** !
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Get a trackmem_file setting from the config script.
  subroutine tem_trackmem_load( me, conf )
    ! ---------------------------------------------------------------------------
    !> The trackmem file
    character(len=*), intent(out) :: me
    !> Handle to the Lua script containing the configuration.
    type(flu_state) :: conf
    ! ---------------------------------------------------------------------------
    integer :: iError
    ! ---------------------------------------------------------------------------

    ! load trackmem_file
    call aot_get_val( L       = conf,            &
      &               key     = 'trackmem_file', &
      &               val     = me,              &
      &               ErrCode = iError,          &
      &               default = ''               )

    if (trim(me) /= '') then
      write(logunit(1),*) ''
      write(logunit(1),*) 'Will track memory consumption into file:'
      write(logunit(1),*) trim(me)
      write(logunit(1),*) 'ATTENTION: This is a slow operation!!!'
      write(logunit(1),*) 'Deactivate the trackmem_file setting for production.'
      write(logunit(1),*) ''
    end if

  end subroutine tem_trackmem_load

end module tem_trackmem_module
