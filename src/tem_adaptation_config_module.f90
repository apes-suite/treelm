! Copyright (c) 2014 Julia Moos <julia.moos@student.uni-siegen.de>
! Copyright (c) 2015 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> This module contains routines to load the adaptation parameters
module tem_adaptation_config_module

  !include treelm modules
  use tem_timeControl_module,       only: tem_timeControl_type,                &
    &                                     tem_timeControl_load
  use tem_aux_module,               only: tem_abort
  use tem_logging_module,           only: logUnit

  !include aotus modules
  use aotus_module,           only: flu_State, aot_get_val,aoterr_NonExistent
  use aot_table_module,       only: aot_table_open, aot_table_close,           &
    &                               aot_table_length, aot_get_val

  implicit none
  private

  public :: tem_adapt_type
  public :: tem_load_adapt

  !> Datatype containig information about adaptation of the mesh
  type tem_adapt_type
    !>Contains variables which are requested in a lua file
    type(tem_timeControl_type) :: time

    !> is adaptation active?
    logical :: active = .false.
  end type tem_adapt_type

contains

! ****************************************************************************** !
  !> Load variables, parent scheme and conditions defined in lua file.
  subroutine tem_load_adapt( me, conf )
    ! ---------------------------------------------------------------------------
    !> adapt_type to be filled
    type( tem_adapt_type ), intent(inout) :: me
    !> lua state to read from
    type( flu_state ), intent(in)         :: conf
    ! ---------------------------------------------------------------------------
    integer :: adapt_handle
    ! ---------------------------------------------------------------------------

    ! Open the adapt table
    call aot_table_open( L       = conf,                                       &
      &                  thandle = adapt_handle,                               &
      &                  key     = 'adaptation' )

    if ( adapt_handle /= 0 ) then
      ! get timeControl conditions
      call tem_timeControl_load( me     = me%time,                             &
        &                        conf   = conf,                                &
        &                        parent = adapt_handle )

      me%active = .TRUE.
      write(logUnit(1),"(A)") 'Mesh Adaptation is active!'
    end if
    ! close the adaptation table
    call aot_table_close( L=conf, thandle=adapt_handle )

  end subroutine tem_load_adapt
! ****************************************************************************** !

end module tem_adaptation_config_module
! ****************************************************************************** !
