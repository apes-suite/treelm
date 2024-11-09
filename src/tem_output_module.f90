! Copyright (c) 2011-2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2012 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2011 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2011-2012 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2014, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2014 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> author: Manuel Hasert
!! This module provides the possibility of dumping the simulation status
!! (e.g. the state vector) to disc and reading the information later on to
!! restart the simulation at a certain point of time.
!!
!! Additionally informations like solver, version, number of elements,...
!! are stored in a seperate header file.
!!
module tem_output_module

  ! include treelm modules
  use tem_tools_module,       only: tem_horizontalSpacer
  use tem_timeControl_module, only: tem_timeControl_type, tem_timeControl_load
  use tem_logging_module,     only: logUnit

  ! include aotus module
  use aotus_module,     only: flu_State
  use aot_table_module, only: aot_table_open, aot_table_close, aot_get_val

  implicit none

  private

  public :: tem_load_output
  public :: tem_output_type

  !> Detailed information about output
  type tem_output_type
   !> is this output object active?
   logical :: active   = .false.
   !> is the output format vtk binary? MH: This is not very well defined here.
   logical :: vtk_bin  = .true.
   !> logical for cell based vtk out or not
   logical :: vertex   = .false.
   !> output file name
   character(len=40)  :: filename
   !> output folder name
   character(len=256) :: folder
   !> stores time control parameters
   type(tem_timeControl_type) :: timeControl
   !> fileformat
   character(len=40) :: fileformat
   !> additional quantity to dump
   logical :: shearstress = .false.
   !>
   logical :: wss = .false.
   ! Debug element output
   !> also write out the halo elements
   logical :: dumpHalos  = .false.
   !> also write out the ghost elements
   logical :: dumpGhosts = .false.
   !> write out in ascii format
   logical :: format_ascii = .false.
  end type tem_output_type

  contains

! ****************************************************************************** !
  !> Read in the output table to output VTK file from the Lua parameter file
  !!
  subroutine tem_load_output( me, conf, key, parent, default_active )
    ! ---------------------------------------------------------------------------
    !>
    type(tem_output_type), intent(inout) :: me
    !>
    type(flu_state) :: conf
    !>
    integer, optional, intent(in) :: parent
    !>
    character(len=*), optional, intent(in) :: key
    !>
    logical, optional, intent(in) :: default_active
    ! ---------------------------------------------------------------------------
    logical  :: def_active
    character(len=32) :: localKey
    integer :: out_table
    integer :: iError
    ! ---------------------------------------------------------------------------

    if( present( default_active )) then
      def_active = default_active
    else
      def_active = .false.
    endif

    me%active = def_active
    ! if the restart is part of another table, like e.g. with the debugStates
    if( present( key )) then
      localKey = key
    else
      localKey = 'output'
    endif

    call aot_table_open( L       = conf,                                       &
      &                  thandle = out_table,                                  &
      &                  parent  = parent,                                     &
      &                  key     = trim(localKey) )

    ! Need to look for:
    ! * active (mandatory)
    ! * boundary_type (default 1)
    ! * fileformat (default binary)
    if (out_table .ne. 0) then
      write(logUnit(1),*)"Loading DEPRECATED output information"
      call aot_get_val( L       = conf,                                        &
        &               thandle = out_table,                                   &
        &               key     = 'active',                                    &
        &               val     = me%active,                                   &
        &               default = def_active,                                  &
        &               ErrCode = iError )
      ! load output parameters only when active=true
      write(logUnit(1),*)'output: ', me%active
      if( me%active ) then
        write(logUnit(1),*)"Output parameters "
        call aot_get_val( L       = conf,                                      &
          &               thandle = out_table,                                 &
          &               key     = 'folder',                                  &
          &               val     = me%folder,                                 &
          &               default = 'output/',                                 &
          &               ErrCode = iError )
        write(logUnit(1),*)" where to store:       "//trim(me%folder)

        call aot_get_val( L       = conf,                                      &
          &               thandle = out_table,                                 &
          &               key     = 'ascii',                                   &
          &               val     = me%format_ascii,                           &
          &               default = .false.,                                   &
          &               ErrCode = iError )
        write(logUnit(1),*)" ascii:                ", me%format_ascii
        if(me%format_ascii)then
          me%vtk_bin=.false.
        end if

        call aot_get_val( L       = conf,                                      &
          &               thandle = out_table,                                 &
          &               key     = 'wss',                                     &
          &               val     = me%wss,                                    &
          &               default = .false.,                                   &
          &               ErrCode = iError )
        call aot_get_val( L       = conf,                                      &
          &               thandle = out_table,                                 &
          &               key     = 'shearstress',                             &
          &               val     = me%shearstress,                            &
          &               default = .false.,                                   &
          &               ErrCode = iError )
        call aot_get_val( L       = conf,                                      &
          &               thandle = out_table,                                 &
          &               key     = 'dumpGhosts',                              &
          &               val     = me%dumpGhosts,                             &
          &               default = .false.,                                   &
          &               ErrCode = iError )
        call aot_get_val( L       = conf,                                      &
          &               thandle = out_table,                                 &
          &               key     = 'dumpHalos',                               &
          &               val     = me%dumpHalos,                              &
          &               default = .false.,                                   &
          &               ErrCode = iError )

        if(me%dumpHalos) me%dumpGhosts = .true.

        call aot_get_val( L       = conf,                                      &
          &               thandle = out_table,                                 &
          &               key     = 'vertex',                                  &
          &               val     = me%vertex,                                 &
          &               default = .false.,                                   &
          &               ErrCode = iError )
        ! load time control to output tracking
        call tem_timeControl_load( conf           = conf,                      &
          &                        parent         = out_table,                 &
          &                        me             = me%timeControl )
      endif
    endif
    call aot_table_close(L=conf, thandle=out_table)
    call tem_horizontalSpacer(fUnit=logUnit(1))

  end subroutine tem_load_output
! ****************************************************************************** !


end module tem_output_module
! ****************************************************************************** !
