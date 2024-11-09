! Copyright (c) 2012-2013 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2013, 2016, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Daniel Petró <daniel.petro@student.uni-siegen.de>
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
! This module contains routines to load the depend table and various components
! inside it. Depend table is usually defined in [[mus_scheme_module]] and
! [[mus_geomIncr_module]] module. Thus, the routines within this module are called
! in [[mus_scheme_module]] and [[mus_geomIncr_module]] module as required
!
module tem_depend_module
!! author: Simon Zimny, Kartik Jain

  ! include treelm modules
  use env_module,            only: labelLen
  use tem_varSys_module,     only: tem_varSys_type, tem_varSys_dump
  use tem_aux_module,        only: tem_abort
  use tem_condition_module,  only: tem_condition_type, tem_load_condition
  use tem_logging_module,    only: logUnit
  use tem_varMap_module,     only: tem_varMap_type, &
    &                              tem_create_varMap
  use tem_tools_module,      only: tem_horizontalSpacer

  ! include aotus modules
  use aotus_module,     only: flu_State, aot_get_val, aoterr_NonExistent, &
    &                         aoterr_Fatal
  use aot_table_module, only: aot_table_open, aot_table_close, &
    &                         aot_table_length, aot_get_val

  implicit none
  private

  public :: tem_depend_type
  public :: tem_load_depend
  public :: tem_init_depend

  !> Datatype containing information about dependency
  !! of geomIncrease header on other scheme(s)
  type tem_depend_type
    !> array of requested variable labels
    character(len=labelLen), allocatable :: varName(:)
    !> Contains name and position of variables to track in global varSys
    type(tem_varMap_type) :: varMap
    !> An instance of the condition type
    type(tem_condition_type), allocatable :: cond(:)
  end type tem_depend_type

  !> load depend table
  interface tem_load_depend
    module procedure tem_load_depend_vector
    module procedure tem_load_depend_single
  end interface tem_load_depend

  contains

! ****************************************************************************** !
  !> Load variables, parent scheme and conditions defined in lua file.
  !! This routine serves as a wrapper and calls the single routine which loads
  !! the various arguments
  !!
  subroutine tem_load_depend_vector(me, conf, parent, label, requireCond)
    !---------------------------------------------------------------------------
    !> list of depend types to be filled
    type(tem_depend_type), allocatable, intent(inout) :: me(:)
    !> lua state to read from
    type(flu_state), intent(in) :: conf
    !> parent table identifier
    integer, intent(in)        :: parent
    !> label to identify depend type
    character(len=*), intent(in) :: label
    !> if true? load condition table for each variable
    logical, optional, intent(in) :: requireCond
    !---------------------------------------------------------------------------
    integer :: iDep, nDeps
    integer :: dep_handle ! local depend table handle
    integer :: dep_sub_handle ! local depend table handle
    !---------------------------------------------------------------------------
    call tem_horizontalSpacer(fUnit = logUnit(1))
    write(logUnit(1),*) 'loading depend table of: '//trim(label)

    ! Open the depend table
    call aot_table_open( L       = conf,                                       &
      &                  parent  = parent,                                     &
      &                  thandle = dep_handle,                                 &
      &                  key     = 'depend' )

    if (dep_handle .ne. 0) then
      ! Check if there are multiple members inside depend
      call aot_table_open( L       = conf,                                     &
        &                  parent  = dep_handle,                               &
        &                  thandle = dep_sub_handle,                           &
        &                  pos     = 1 )
      ! If there is only one member, call the load routine once
      if (dep_sub_handle .eq. 0) then
        allocate(me(1))
        call aot_table_close( L = conf, thandle = dep_sub_handle )
        call tem_load_depend_single( me          = me(1),      &
          &                          conf        = conf,       &
          &                          parent      = dep_handle, &
          &                          label       = label,      &
          &                          requireCond = requireCond )
      else
        call aot_table_close( L = conf, thandle = dep_sub_handle )
        nDeps = aot_table_length( L = conf, thandle = dep_handle )
        allocate( me( nDeps ))

        ! Read the subtables individually
        do iDep = 1, nDeps
          call aot_table_open( L       = conf,                                 &
            &                  parent  = dep_handle,                           &
            &                  thandle = dep_sub_handle,                       &
            &                  pos     = iDep )

          call tem_load_depend_single( me          = me(iDep),       &
            &                          conf        = conf,           &
            &                          parent      = dep_sub_handle, &
            &                          label       = label,          &
            &                          requireCond = requireCond     )

          call aot_table_close( L = conf, thandle = dep_sub_handle )
        end do
      end if
    else
      write(logUnit(1),*) 'depend table is not defined'
    end if

    call aot_table_close( L = conf, thandle = dep_handle )
    call tem_horizontalSpacer(fUnit = logUnit(1))

  end subroutine tem_load_depend_vector
! ****************************************************************************** !


! ****************************************************************************** !
  !> Load single dependent variable of the scheme, in case of geomIncr: load the
  !! dependent variable as well as the conditions to be imposed for the geometry
  !! increase to take place.
  !!
  subroutine tem_load_depend_single( me, conf, parent, label, requireCond )
    !---------------------------------------------------------------------------
    !> depend type to be filled
    type(tem_depend_type), intent(inout) :: me
    !> lua state to read from
    type(flu_state), intent(in)          :: conf
    !> handle of parent table
    integer, intent(in)                  :: parent
    !> label to identify depend type
    character(len=*), intent(in) :: label
    !> if true? load condition table for each variable
    logical, optional, intent(in) :: requireCond
    !---------------------------------------------------------------------------
    integer, allocatable :: vError(:)
    integer :: cond_handle
    logical :: requireCond_loc
    !---------------------------------------------------------------------------

    ! if require condition is not present do not load condition table
    if (present(requireCond)) then
      requireCond_loc = requireCond
    else
      requireCond_loc = .false.
    end if

    ! Load the variable names
    call aot_get_val( val       = me%varName, &
      &               ErrCode   = vError,     &
      &               maxLength = 100,        &
      &               L         = conf,       &
      &               thandle   = parent,     &
      &               key       = 'variable'  )

    if ( any(btest(vError, aoterr_Fatal)) ) then
      write(logUnit(1),*) 'ERROR: could not load varnames for depend variable'
      call tem_abort()
    end if

    ! The below code snippet loads the conditions defined within depend table
    ! in the case of geometry increase, whereas in case this routine is called
    ! to load the depend table within schemes then conditions are not available
    ! and this load will not be effective in that case
    if (requireCond_loc) then
      call aot_table_open( L       = conf,                                     &
        &                  parent  = parent,                                  &
        &                  thandle = cond_handle,                              &
        &                  key     = 'condition')
      if (cond_handle /= 0) then
        ! condition table is defined
        ! close table and load with tem_load_condition
        call aot_table_close( L=conf, thandle = cond_handle )

        call tem_load_condition( me     = me%cond, &
          &                      conf   = conf,    &
          &                      parent = parent   )
      else
        write(logUnit(1),*) 'ERROR: Condition table is not defined in depend ' &
          &                 //'table of: '//trim(label)
        call tem_abort()
      end if ! Condition table defined

      ! check if there is condition for each variable
      if (size(me%cond) /= size(me%varname)) then
        write(logUnit(1),*) 'Error: Nr. of conditions \= Nr. of variables ' &
          &                 //'in depend table of: '//trim(label)
        write(logUnit(1),*) 'nCond: ', size(me%cond),  &
          &                 ' nVars: ', size(me%varname)
        call tem_abort()
      end if
    end if ! require condition

  end subroutine tem_load_depend_single
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine initializes the loaded depend table
  !!
  subroutine tem_init_depend( me, varSys )
    !---------------------------------------------------------------------------
    !> list of depend types to be filled
    type(tem_depend_type), intent(inout)  :: me(:)
    !> list of all global variable systems
    type(tem_varSys_type), intent(in)     :: varSys
    !---------------------------------------------------------------------------
    integer :: iDepend, nDepends
    !---------------------------------------------------------------------------
    nDepends = size(me)

    do iDepend = 1, nDepends
      ! map variables
      ! create depend variable position in the global varSys
      call tem_create_varMap( varname = me(iDepend)%varname, &
        &                     varSys  = varSys,              &
        &                     varMap  = me(iDepend)%varMap   )

    end do ! iDepend

  end subroutine tem_init_depend
! ****************************************************************************** !


end module tem_depend_module
! ****************************************************************************** !
