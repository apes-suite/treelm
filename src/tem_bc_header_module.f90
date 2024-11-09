! Copyright (c) 2011-2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2013 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2014 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
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
!> The tem_bc_header is providing an overview of the boundary conditions
!! defined in the configuration of the solver.
!!
!! It provides collects the
!! labels and the kinds of the configured BCs, and connects them to the
!! corresponding boundaries defined in the mesh, given in the BC_prop.
!!
module tem_bc_header_module

  ! include treelm modules
  use env_module,         only: LabelLen
  use tem_tools_module,   only: tem_horizontalSpacer
  use tem_aux_module,     only: tem_abort
  use tem_bc_prop_module, only: tem_bc_prop_type
  use tem_logging_module, only: logUnit

  ! include aotus modules
  use aot_table_module, only: aot_table_open, aot_table_close, aot_get_val,    &
    &                         aot_table_length
  use flu_binding,      only: flu_State

  implicit none

  private

  public :: tem_bc_header_type
  public :: tem_load_bc_header, tem_open_bc
  public :: tem_out_bc_header

  !> This type describes the general, not solver specific, header information
  !! given in the Lua configuration of the solvers.
  type tem_bc_header_type
    integer :: nBCs !< Number of Bounary conditions in the Lua config
    !> Label for each boundary condition
    character(len=LabelLen), allocatable :: label(:)
    !> The kind of each boundary condition
    character(len=LabelLen), allocatable :: BC_kind(:)
    !> ID of the boundary in the mesh property, as defined by the bc_prop
    !! If an entry is not positive, there is no corresponding boundary
    !! condition found.
    integer, allocatable :: BC_ID(:)
  end type


contains


! ****************************************************************************** !
  !> This subroutine reads in the boundary conditions specified in the
  !! configuration file, and connects them to the corresponding entries in the
  !! treelmesh.
  !! If there boundary conditions in the mesh, for which no configuration is
  !! found, the program is aborted!
  !!
  subroutine tem_load_bc_header( me, conf, parentHandle, BC_prop )
    ! ---------------------------------------------------------------------------
    !> The boundary conditions to fill
    type(tem_bc_header_type), intent(inout) :: me
    !> A handle to the Lua configuration script, to read the data from
    type(flu_State) :: conf
    !> The boundary properties of the treelmesh, to connect the header to
    type(tem_BC_prop_type), intent(in) :: BC_prop
    !> handle for schemes table
    integer, intent(in),optional       :: parentHandle
    ! ---------------------------------------------------------------------------
    integer :: iBC, nBCs
    integer :: bc_handle, myHandle
    integer :: iError
    ! ---------------------------------------------------------------------------

    call tem_horizontalSpacer(fUnit=logUnit(1))
    write(logUnit(1),*) ' Loading BC header'
    ! if fluid informations in scheme table parentHandle /= 0
    call aot_table_open( L       = conf,                                       &
      &                  parent  = parentHandle,                               &
      &                  thandle = bc_handle,                                  &
      &                  key     = 'boundary_condition' )
    nBCs = aot_table_length(L=conf, thandle=bc_handle)

    ! allocate the label, bc_kind and BC_ID arrays in me
    call tem_init_bc_header( me, nBCs )

    do iBC = 1, nBCs
      call aot_table_open( L       = conf,                                     &
        &                  parent  = bc_handle,                                &
        &                  thandle = myHandle,                                 &
        &                  pos     = iBC )
      if (myHandle /= 0) then
        call aot_get_val( L       = conf,                                      &
          &               thandle = myHandle,                                  &
          &               val     = me%label(iBC),                             &
          &               ErrCode = iError,                                    &
          &               key     = 'label',                                   &
          &               default = 'unnamed' )

        call aot_get_val( L       = conf,                                      &
          &               thandle = myHandle,                                  &
          &               val     = me%BC_kind(iBC),                           &
          &               ErrCode = iError,                                    &
          &               key     = 'kind',                                    &
          &               default = 'none' )
        write(logUnit(1),*)'  Found BC '//trim(me%Label(iBC))//' of kind '//   &
          &            trim(me%BC_kind(iBC))
      else
        write(logUnit(1),*)'WARNING: BC definition should be a table!'
        write(logUnit(1),*)'       : Assuming unnamed BC with kind none!'
        me%Label(iBC) = 'unnamed'
        me%BC_kind(iBC) = 'none'
      end if
      call aot_table_close(L=conf, thandle=myHandle)
    end do ! iBC

    call aot_table_close(L=conf, thandle=bc_handle)

    ! Now all the BCs defined in the config file are known,
    ! try to match them against the labels in the BC_prop.
    ! if a bc label can not be matched, treat it as a fatal error.
    call tem_match_bc_header( me = me, BC_prop = BC_prop )

    if ( nBCs .gt. 0 ) then
      write(logUnit(1),"(A, I0)") ' Number of Boundaries: ', me%nBCs
    else
      write(logUnit(1),*) ' No Boundary conditions!'
    end if

  end subroutine tem_load_bc_header
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine looks for a given label in the given boundary conditions
  !! table, and returns the according table handle.
  !! Note, that this should usually not be necessary, as the number of
  !! the header is given by the ordering in the bc_header_type, and you
  !! can use the desired position directly to look up a specific bc in the
  !! configuration script.
  !!
  subroutine tem_open_bc( label, bc_table, conf, thandle )
    ! ---------------------------------------------------------------------------
    !> The label to look for
    character(len=LabelLen), intent(in) :: label
    !> Handle to the boundary_condition table, to look in
    integer, intent(in) :: bc_table
    !> Handle of the Lua script to use
    type(flu_State) :: conf
    !> Returned handle of to the entry providing the requested label
    integer, intent(out) :: thandle
    ! ---------------------------------------------------------------------------
    character(len=LabelLen) :: bc_name
    integer :: myHandle
    integer :: iBC, iError
    integer :: nBCs
    ! ---------------------------------------------------------------------------

    nBCs = aot_table_length(L=conf, thandle=bc_table)
    ! Default thandle to nothing found.
    thandle = 0
    do iBC=1,nBCs
      call aot_table_open( L       = conf,                                     &
        &                  parent  = bc_table,                                 &
        &                  thandle = myHandle,                                 &
        &                  pos     = iBC )
      if (myHandle /= 0) then
        call aot_get_val( L       = conf,                                      &
          &               thandle = myHandle,                                  &
          &               val     = bc_name,                                   &
          &               ErrCode = iError,                                    &
          &               key     = 'label',                                   &
          &               default = 'unnamed' )
        if (adjustl(bc_name) == adjustl(label)) then
          ! Found the matching boundary table, returning it...
          thandle = myHandle
          exit
        end if
      end if
      call aot_table_close(L=conf, thandle=myHandle)
    end do

  end subroutine tem_open_bc
! ****************************************************************************** !


! ****************************************************************************** !
  !> This routine does the allocation job
  !!
  subroutine tem_init_bc_header( me, nBCs )
    ! ---------------------------------------------------------------------------
    !> The boundary header
    type(tem_bc_header_type), intent(out) :: me
    integer, intent(in) :: nBCs
    ! ---------------------------------------------------------------------------

    me%nBCs = nBCs
    allocate(me%Label(me%nBCs))
    allocate(me%BC_kind(me%nBCs))
    allocate(me%BC_ID(me%nBCs))

    ! Set initial values for BC_ID
    me%BC_ID = -1

  end subroutine tem_init_bc_header
! ****************************************************************************** !


! ****************************************************************************** !
  !> This routine match the labels in me(tem_bc_header_type) against the labels
  !! in the bcProp([[tem_BC_prop_type]]).
  !! If a bc label can not be matched, the code will STOP!
  !!
  subroutine tem_match_bc_header( me, BC_prop )
    ! ---------------------------------------------------------------------------
    !> The boundary header to fill
    type(tem_bc_header_type) :: me
    !> The boundary properties of the treelmesh, to connect the header to
    type(tem_BC_prop_type), intent(in) :: BC_prop
    ! ---------------------------------------------------------------------------
    integer :: iBC, iBCtype
    character(len=LabelLen) :: proplabel
    logical :: MatchedBC
    ! ---------------------------------------------------------------------------

    write(logUnit(3),"(A)") ' matching boundary defination in lua to mesh'
    ! Outside loop: loop over the bc labels in BC_prop
    do iBCtype = 1, BC_prop%nBCtypes

      proplabel = adjustl(BC_prop%BC_label(iBCtype))
      matchedBC = .false.
      ! Inside loop: loop over the bc labels in bc header
      do iBC = 1, me%nBCs
        if (me%BC_ID(iBC) == -1) then
          ! This ID has not been matched yet, check it now
          ! BC_ID is set to be -1 during initialization
          if (adjustl(me%label(iBC)) == proplabel) then
            ! Found a match, assign the configured BC, and leave the loop
            me%BC_ID(iBC) = iBCtype
            matchedBC = .true.
            exit
          end if
        end if
      end do ! iBC in me%label

      ! If this BC is not matched, give a error message and stop the code.
      if ( .not. matchedBC ) then
        write(logUnit(1),*)"ERROR:"
        write(logUnit(1),*)"No configuration found for boundary "//            &
          &             trim(proplabel)//" from the BC property in the mesh!"
        write(logUnit(1),*)"This could also be caused by double definition "   &
          &            //"of the same BC in the mesh."
        write(logUnit(1),*)"This is fatal, STOPPING execution now..."
        call tem_abort()
      end if
    end do ! iBCtype in bc_prop

    if (minval(me%BC_ID) == -1) then
      write(logUnit(1),*)"WARNING: There are more boundaries defined in the "  &
        &                //" config file, than in the mesh properties!"
    end if

  end subroutine tem_match_bc_header
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine tem_out_bc_header( me, outUnit )
    type( tem_bc_header_type ), intent(in) :: me
    integer, intent(in) :: outUnit
    integer :: iBC

    write( outUnit, "(A, I2)") ' in bc_header, number of BCs: ', me%nBCs
    write( outUnit, "(3A12)" ) 'label', 'BC kind', 'BC_ID'
    do iBC = 1, me%nBCs
      write( outUnit, "( A12, A12, I12)" ) &
        &               trim(me%label(iBC)), trim(me%BC_kind(iBC)), me%BC_ID(iBC)
    end do
  end subroutine tem_out_bc_header
! ****************************************************************************** !


end module tem_bc_header_module
! ****************************************************************************** !
