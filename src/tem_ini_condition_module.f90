! Copyright (c) 2012-2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012, 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2013, 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2013 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
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
!> This module provides the functionality to fill the initial values in the
!! computational domain.
module tem_ini_condition_module

  ! include treelm modules
  use env_module,            only: rk, labelLen, pathLen
  use tem_aux_module,        only: tem_abort
  use tem_tools_module,      only: tem_horizontalSpacer
  use tem_logging_module,    only: logUnit
  use tem_spatial_module,    only: tem_load_spatial, tem_spatial_type
  use tem_depend_module,     only: tem_depend_type, tem_load_depend

  ! include aotus modules
  use aot_table_module, only: aot_table_open, aot_table_close, aot_get_val
  use aotus_module,     only: flu_State

  implicit none

  private

  public :: tem_ini_condition_type
  public :: tem_load_ic

  !> Definition of the initial condition.
  !! The ini_state must exist for each variable
  type tem_ini_condition_type
    !> spatial variable names expected from config file
    character(len=labelLen), allocatable :: StateName(:)
    !> initial state of each variable. size is nVars
    type(tem_spatial_type), allocatable :: ini_state(:)
  end type tem_ini_condition_type

contains

! ****************************************************************************** !
  !> Load initial condition
  !!
  !! Check if restart file is defined. Then use restart file to initialize
  !! variables. If not then load state variables either as constant or function
  !! \author Kannan Masilamani
  !!
  subroutine tem_load_ic(me, conf, parent, StateName, key, errCode )
    ! ---------------------------------------------------------------------------
    !> Initial condition type
    type(tem_ini_condition_type), intent(out) :: me
    !> lua state type
    type(flu_State) :: conf
    !> array of defined initial state variables
    character(len=*), intent(in) :: StateName(:)
    !> errCode of variables
    !! Solver should take appropriate action according to errCode
    integer, intent(out) :: errCode(:)
    !>
    integer, intent(in), optional :: parent
    !> key for the initial cond.
    character(len=*), intent(in), optional :: key
    ! ---------------------------------------------------------------------------
    integer :: iState, nStates
    integer :: ic_table
    character(len=labelLen) :: lockey ! key for the initial cond.
    ! ---------------------------------------------------------------------------

    nStates = size(StateName)

    allocate(me%StateName(nStates))
    allocate(me%ini_state(nStates))
    me%StateName = StateName
    if ( present(key) ) then
      lockey = key
    else
      lockey = 'initial_condition'
    end if
    call aot_table_open(L=conf, parent=parent, thandle=ic_table, key=lockey)
    call tem_horizontalSpacer(fUnit=logUnit(1))
    write(logUnit(1),*) ' Loading initial conditions'
    if (ic_table == 0) then
      ! Interpret the initial_condition entry as name for the restart file,
      ! if it is not a table.
      call tem_abort('Initial condition table is not found!')
    else
      ! We are not doing a restart, need to find initial conditions for all
      ! requested states.
      do iState = 1, nStates
        call tem_load_spatial( me           = me%ini_state(iState),       &
          &                    conf         = conf,                       &
          &                    parent       = ic_table,                   &
          &                    key          = trim(me%stateName(iState)), &
          &                    errCode      = errCode( iState ),          &
          &                    defaultValue = 0._rk                       )
        ! check if the initial condition has kind = none
        ! if ( me%ini_state(iState)%kind == 'none' ) then
        !   write(*,*) 'Initial condition variable', trim(me%stateName(iState)), &
        !    & 'has kind none, stopping ...'
        !   call tem_abort()
        ! end if
      end do

    end if
    call aot_table_close(L=conf, thandle=ic_table)

  end subroutine tem_load_ic
! ****************************************************************************** !


end module tem_ini_condition_module
! ****************************************************************************** !
