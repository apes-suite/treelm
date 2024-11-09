! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2014, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2013-2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2013 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
!! This module provides the datatypes for load balancing algorithms.
!!
module tem_balance_module

  ! include treelm modules
  use mpi
  use env_module,             only: rk, rk_mpi, labelLen
  use flu_binding,            only: flu_State
  use tem_time_module,        only: tem_time_never
  use tem_timeControl_module, only: tem_timeControl_type, tem_TimeControl_load
  use tem_tools_module,       only: tem_horizontalSpacer
  use tem_logging_module,     only: logUnit
  use tem_aux_module,         only: tem_abort

  ! include aotus modules
  use aotus_module,     only: aot_get_val
  use aot_table_module, only: aot_table_open, aot_table_close, aot_get_val

  implicit none
  private

  public :: tem_balance_type
  public :: tem_calc_imbalance
  public :: tem_balance_load

  !> Load balancing information and control
  type tem_balance_type

    !> load balancing type
    character(len=labelLen) :: kind = 'none'

    !> is dynamic load balancing activated?
    logical :: dynamic = .false.

    !> whether dump weight file
    logical :: weight = .false.

    !> control about when to do dynamic balancing
    type(tem_timeControl_type) :: timeControl

    !> last iteration the dynamic load balancing was called
    ! integer :: lastIter = 0

  end type tem_balance_type

contains

! ****************************************************************************** !
  !> Read all the configuration options for load balancing from the
  !! configuration.
  !!
  !! To activate the load balancing, use
  !! ~~~~~~~~~~~~~~~~~~~~~{.lua}
  !! balance = {
  !!   dynamic = true,          -- perform dynamic load balancing
  !!   kind = 'IBM',            -- kind = 'levelwise','boundary','IBM'
  !!   time_control = {         -- time control for slots when to perform load
  !!        min      = 0.01,    --   balancing
  !!        interval = 0.01,
  !!        max      = 1.0,
  !!   },
  !!   folder = 'balance/', -- directory where to store tmp data
  !! }
  !! ~~~~~~~~~~~~~~~~~~~~~
  !!
  subroutine tem_balance_load(me, conf)
    ! ---------------------------------------------------------------------------
    !> The load balancing type
    type(tem_balance_type), intent(inout) :: me
    !> lua configuration handle
    type(flu_state)                       :: conf
    ! ---------------------------------------------------------------------------
    integer :: iError, thandle
    character(len=labelLen) :: local_defaultString = 'simple'
    ! ---------------------------------------------------------------------------

    ! Open Balance table
    call aot_table_open(L       = conf,     &
      &                 thandle = thandle,  &
      &                 key     = "balance" )

    if( thandle /= 0 ) then
      ! table is open !
      call aot_get_val( L       = conf,       &
        &               thandle = thandle,    &
        &               key     = 'dynamic',  &
        &               val     = me%dynamic, &
        &               ErrCode = iError,     &
        &               default = .false. )

      call aot_get_val( L       = conf,               &
        &               thandle = thandle,            &
        &               key     = 'kind',             &
        &               val     = me%kind,            &
        &               ErrCode = iError,             &
        &               default = local_defaultString )

      call aot_get_val( L       = conf,       &
        &               thandle = thandle,    &
        &               key     = 'weight',   &
        &               val     = me%weight,  &
        &               ErrCode = iError,     &
        &               default = .false. )

      ! read the directory for dumping restart files of load balancing
      ! call aot_get_val( L       = conf,                                        &
      !   &               thandle = thandle,                                     &
      !   &               key     = 'folder',                                    &
      !   &               val     = me%restart%controller%writePrefix,           &
      !   &               ErrCode = iError )

      ! Load the time controls for dynamic load balancing
      call tem_timeControl_load( me     = me%timeControl,                      &
        &                        conf   = conf,                                &
        &                        parent = thandle )

    else ! NO balance table given

      me%kind = 'none'
      me%dynamic = .false.
      me%weight  = .false.
      me%timeControl%Interval = tem_time_never()
      me%timeControl%Min      = tem_time_never()
      me%timeControl%Max      = tem_time_never()

    end if

    call aot_table_close(L=conf, thandle=thandle)

  end subroutine tem_balance_load
! ****************************************************************************** !


! ****************************************************************************** !
  !> Evaluate the imbalance of all the processes by each rank.
  !!
  function tem_calc_imbalance( myCost, comm, nProcs, isRoot ) result ( imbalance )
    ! ---------------------------------------------------------------------------
    !> each process cpu cost. Basis to evaluate the imbalance
    real(kind=rk), intent(in) :: myCost
    !> MPI Communicator
    integer, intent(in) :: comm, nProcs
    !> in percentage
    real(kind=rk) :: imbalance
    !> Whether this rank is the root
    logical, intent(in) :: isRoot
    ! ---------------------------------------------------------------------------
    ! computation cost over all processes
    real(kind=rk) :: totalCost, maxCost
    integer :: iErr ! error handle
    ! ---------------------------------------------------------------------------

    call mpi_reduce( myCost,   maxCost, 1, rk_mpi, mpi_max, 0, comm, iErr)
    call mpi_reduce( myCost, totalCost, 1, rk_mpi, mpi_sum, 0, comm, iErr)

    imbalance = -9999._rk
    ! the value range of imbalance is [100%, nProcs]
    ! The smaller the value, the better the balance is
    if ( isRoot ) then
      totalCost = totalCost/dble(nProcs)
      imbalance = maxCost / totalCost * 100

      call tem_horizontalSpacer(fUnit=logUnit(5))
      write(logUnit(3),"(A)")         ' Load balance evaluation:'
      write(logUnit(3),"(A,F10.2,A)")   '  max CPU cost:', maxCost, ' (s)'
      write(logUnit(3),"(A,F10.2,A)")   '  ave CPU cost:', totalCost, ' (s)'
      write(logUnit(3),"(A,F10.2,A)") ' imbalance (i.e. max/ave):', imbalance, '%'
      call tem_horizontalSpacer(fUnit=logUnit(5))
    end if

  end function tem_calc_imbalance
! ****************************************************************************** !


end module tem_balance_module
! ****************************************************************************** !
