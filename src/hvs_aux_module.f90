! Copyright (c) 2015 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2019 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> author: Kannan Masilamani
!! This module contains some aux routines for libharvesting
!!
module hvs_aux_module

  ! include treelm modules
  use env_module,           only: labelLen
  use tem_solveHead_module, only: tem_solveHead_type
  use tem_logging_module,   only: logUnit
  use tem_aux_module,       only: tem_print_execInfo, utc_date_string

  implicit none

  private

  public :: hvs_banner

contains

! ****************************************************************************** !
  !> Print harvester banner to screen
  !!
  subroutine hvs_banner( solveHead, supportSolName )
    ! ---------------------------------------------------------------------------
    type( tem_solveHead_type )  :: solveHead
    !> Name of the solver which invokes this routine
    character(len=*) :: supportSolName
    ! ---------------------------------------------------------------------------
    character(len=26) :: dat_string
    ! ---------------------------------------------------------------------------

    write(logunit(0),*)"                                                          "
    write(logunit(0),*)"        __                               __               "
    write(logunit(0),*)"       / /_  ____  ______   _____  _____/ /____  _____    "
    write(logunit(0),*)"      / __ \/ __ \/ ___/ | / / _ \/ ___/ __/ _ \/ ___/    "
    write(logunit(0),*)"     / / / / /_/ / /   | |/ /  __(__  ) /_/  __/ /        "
    write(logunit(0),*)"    /_/ /_/\____/_/    |___/\___/____/\__/\___/_"         &
      &                                                  //trim(solveHead%version)
    write(logunit(0),*)"                                                          "
    write(logunit(0),*)" (C) 2012 German Research School for Simulation Sciences  "
    write(logunit(0),*)"                                                          "
    call tem_print_execInfo()
    write(logunit(0),*)""
    write(logunit(0),*)"With "//trim(supportSolName)//" support."
    write(logunit(0),*)""
    dat_string = utc_date_string()
    write(logUnit(1),*)"Run at: "//trim(dat_string)
    write(logUnit(1),*)""
    write(logunit(0),*)"----------------------------------------------------------"

  end subroutine hvs_banner
! ****************************************************************************** !

end module hvs_aux_module
! ****************************************************************************** !
