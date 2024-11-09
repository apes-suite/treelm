! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2013, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2015 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2014 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!
! Copyright (c) 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
!
! Parts of this file were written by Peter Vitt for University of Siegen.
!
! ****************************************************************************** !
!> Some auxilary routines, providing
!! frequently needed common tasks.
!!
module tem_solveHead_module

  ! include treelm modules
  use env_module,             only: LabelLen, PathLen
  ! This module has to be created on the fly by the waf environment
  use soi_revision_module,    only: soi_solver_revision

  ! include aotus module
  use aotus_module,    only: flu_State

  implicit none

  private

  public :: tem_solveHead_type
  public :: tem_solverTag
  public :: tem_init_solveHead

  !> Solver header information
  type tem_solveHead_type
    !> Name of current simulation
    character(len=PathLen) :: simName
    !> Mesh folder. set default to empty string
    character(len=PathLen)  :: meshFolder = ''
    !> main lua configfile full path.
    !! Default is set to empty string to do check on this
    !! before loading filename from command line
    character(len=PathLen)  :: configFile = ''
    character(len=LabelLen) :: solName    !< solver name
    character(len=LabelLen) :: version    !< version of the solver
    !> mercurial repository revision of the solver
    character(len=16) :: revision
    !> aotus lua conf state to load main configuration file
    !! and solver specific lua functions as chunk
    type(flu_state), allocatable :: conf(:)
    !> solver timer handle.
    !! use to lookup in timerData in timer_module to
    !! determine time taken for this solver
    integer :: timerHandle
    ! The filename to use for tracking the memory
    character(len=pathLen)       :: trackmem_file
  end type tem_solveHead_type

  contains

! ****************************************************************************** !
  !> Function to return a solver tag (combination of solver name and version)
  !!
  function tem_solverTag( solver ) result( tag )
    ! ---------------------------------------------------------------------------
    !> solver information
    type(tem_solveHead_type), intent(in) :: solver
    ! ---------------------------------------------------------------------------
    ! defining local variables
    character(len=LabelLen) :: tag  ! solver
    ! ---------------------------------------------------------------------------

    tag = trim(solver%solName(:labellen/2)) // '_' &
      & // trim(solver%version(:(labellen/2)-1))
  end function tem_solverTag
! ****************************************************************************** !

! ****************************************************************************** !
  !> Routine to initialize solver head with name, version and revision number
  subroutine tem_init_solveHead( me, solName, version )
    !> solver info
    type( tem_solveHead_type), intent(out) :: me
    !> name of the solver
    character(len=*), intent(in) :: solName
    !> version of the solver
    character(len=*), intent(in) :: version
    ! ---------------------------------------------------------------------------

    me%solName  = trim(solName)
    me%version = trim(version)
    me%revision = soi_solver_revision

  end subroutine tem_init_solveHead
! ****************************************************************************** !

end module tem_solveHead_module
! ****************************************************************************** !
