! Copyright (c) 2015-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2015-2016, 2018 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> Contains data types and routines needed for coupling via apesmate
module tem_coupling_module

  ! include treelm modules
  use env_module,           only: rk, labelLen, pathLen, globalMaxLevels

  use aotus_module,         only: flu_State, aot_get_val, aoterr_Fatal
  use aot_table_module,     only: aot_table_open,  &
    &                             aot_table_close, &
    &                             aot_table_length

  use tem_logging_module,    only: logUnit
  use tem_aux_module,        only: tem_abort
  use tem_grow_array_module, only: grw_intArray_type, grw_realArray_type
  use tem_comm_module,       only: tem_communication_type


  implicit none

  private

  public :: tem_aps_coupling_type
  public :: tem_aps_load_coupling

  ! type which include all exchange points information
  type cpl_value_type
    !> number of points per level
    integer :: nPnts = 0

    !> Global process ids to evaluate the points
    !! It is deallocated after recvBuffer is filled
    integer, allocatable :: pntRanks(:)

    !> Evaluated variable value on each point.
    !! If variable is time-independent then values are evaluated and stored
    !! at initialization stage, in this case point arrays are not stored.
    !! nComp = nScalars in the tem_coupling_type%varnames
    !! Access: (iVal-1)*nComp + iComp
    !type(grw_realArray_type) :: evalVal
    real(kind=rk), allocatable :: evalVal(:)

    !> Receive communication buffer to fill evalVal
    type(tem_communication_type) :: recvBuffer
  end type cpl_value_type

  !> Coupling description defined in config file from load space time function
  !! which is called from load boundary condition or load sources
  type tem_aps_coupling_type

    !> Remote domain label to get data from
    character(len=labelLen) :: rem_domLabel
    !> Domain ID of remote domain label
    integer :: rem_domID

    !> Number of variables to get from remote domain
    integer :: nVars
    !> List of variables to get from domain
    character(len=labelLen), allocatable :: varNames(:)

    !> nScalars of varNames
    !! Must be same as nComps in stFun
    integer :: nScalars

    !> Used to decided whether this spacetime functions are used
    !! for surface or volume i.e boundary or source.
    !! Boundary is treated as surface and source as volume
    !! coupling type can be rather surface or volume.
    !! For boundary. isSurface = 0
    !! For volume, isSurface = 1
    integer :: isSurface = -1
    !> store value on each level
    type(cpl_value_type) :: valOnLvl(globalMaxLevels)
  end type tem_aps_coupling_type


contains


  ! ****************************************************************************
  !> This routine loads coupling defintion from boundary condition table
  subroutine tem_aps_load_coupling(me, thandle, conf)
    ! -------------------------------------------------------------------------!
    !> Coupling description to be filled from config file
    type(tem_aps_coupling_type), intent(out) :: me
    !> Boundary condition sub table
    integer, intent(in) :: thandle
    !> Lua script to obtain the configuration data from.
    type(flu_State), intent(in) :: conf
    ! -------------------------------------------------------------------------!
    integer :: iError
    integer, allocatable :: vError(:)
    ! -------------------------------------------------------------------------!

    ! get the domain name with whom we want to couple
    call aot_get_val( L       = conf,            &
      &               thandle = thandle,         &
      &               key     = 'domain_from',   &
      &               val     = me%rem_domLabel, &
      &               ErrCode = iError           )
    if (iError .NE. 0) then
      write(logUnit(1),*) ' No coupling domain for the coupling boundary' &
        &                 //' is defined, coupling can not work, abort...'
      call tem_abort
    end if

    ! get the list of variables
    ! look for the request variable table
    call aot_get_val( val       = me%varNames,    &
      &               ErrCode   = vError,         &
      &               maxLength = 100,            &
      &               L         = conf,           &
      &               thandle   = thandle,        &
      &               key       = 'input_varname' )

    if ( any(btest(vError, aoterr_Fatal)) ) then
      write(logUnit(1),*) 'ERROR: could not input_varname for coupling variable'
      call tem_abort()
    end if

    me%nVars = size(me%varNames)

   end subroutine tem_aps_load_coupling
  ! ****************************************************************************


end module tem_coupling_module
! ****************************************************************************** !
