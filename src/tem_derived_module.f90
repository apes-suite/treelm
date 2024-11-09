! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2014 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012, 2014-2018 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2013, 2015-2016, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2013-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014-2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016-2017 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
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
! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015, 2018, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
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
!--------------------------------------------
!    A O S - Array of structures layout new
!-------------------------------------------
! Access to get_point value output
! Access to get_element value output
! ******************************************************************************
!> This module
!!
module tem_derived_module
  ! include treelm modules
  use tem_varSys_module,            only: tem_varSys_type,               &
    &                                     tem_varSys_solverData_evalElem_type
  use tem_variable_module,          only: tem_variable_type
  use tem_spacetime_fun_module,     only: tem_st_fun_linkedList_type
  use tem_spacetime_var_module,     only: tem_varSys_append_stFun
  use tem_logging_module,           only: logUnit
  use tem_aux_module,               only: tem_abort
  use tem_dyn_array_module,         only: PositionOfVal
  use tem_operation_var_module,     only: tem_varSys_append_operVar

  implicit none

  private

  public :: tem_varSys_append_luaVar

contains

  ! ****************************************************************************
  !> subroutine to add the variables from the input lua script to the varsys
  subroutine tem_varSys_append_luaVar( luaVar, varSys, st_funList, &
    &                                  solverData_evalElem         )
    !--------------------------------------------------------------------------
    !> variables defined in the lua file
    type(tem_variable_type), intent(in) :: luaVar(:)

    !> global variable system to which luaVar to be appended
    type(tem_varSys_type), intent(inout) :: varSys

    !> contains spacetime functions of all variables
    type(tem_st_fun_linkedList_type), intent(inout) :: st_funList

    !> A callback routine to allow the definition of solver specific
    !! element evaluation for space-time functions.
    !!
    !! This routine can be used to construct more than a single degree of
    !! freedom for a spacetime function in an element.
    type(tem_varSys_solverData_evalElem_type), &
      &  optional, intent(in) :: solverData_evalElem
    ! --------------------------------------------------------------------------
    integer :: iVar, varPos
    ! --------------------------------------------------------------------------
    if (size(luaVar) > 0) &
      & write(logUnit(5),*) 'Append variables defined in lua file to varSys'

    do iVar = 1, size(luaVar)
      write(logUnit(5),'(A,I2,A)') 'Appending variable ', iVar, ': ' &
        & // trim(luaVar(iVar)%label)
      varPos = PositionOfVal( me  = varSys%varName,          &
        &                     val = trim(luaVar(iVar)%label) )
      ! If variable already exist in varSys then do nothing
      if (varPos>0) then
        write(logUnit(5),*) 'Variable already exists!'
        cycle
      end if

      select case(trim(luaVar(iVar)%varType))
      case('st_fun')
        call tem_varSys_append_stFun(                    &
          &    stFunVar            = luaVar(iVar),       &
          &    varSys              = varSys,             &
          &    st_funList          = st_funList,         &
          &    solverData_evalElem = solverData_evalElem )

      case('operation')
        call tem_varSys_append_operVar(                       &
          &    operVar                  = luaVar(iVar),       &
          &    varSys                   = varSys,             &
          &    solverData_evalElem      = solverData_evalElem )

      case default
        if (associated(luaVar(iVar)%append_solverVar)) then
          call luaVar(iVar)%append_solverVar(                           &
            &                 varSys              = varSys,             &
            &                 solverData_evalElem = solverData_evalElem )
        else
          write(logUnit(1),*) 'WARNING: varType: '           &
            &                 // trim(luaVar(iVar)%varType)  &
            &                 // ' not supported. Variable ' &
            &                 // trim(luaVar(iVar)%label)    &
            &                 // ' is not appended.'
          cycle ! go to next variable
        end if
      end select

    end do !iVar

  end subroutine tem_varSys_append_luaVar
  ! ****************************************************************************


end module tem_derived_module
! ******************************************************************************
