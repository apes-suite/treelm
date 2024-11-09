! Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2017 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
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
module tem_operation_module

  use, intrinsic :: iso_c_binding,    only: c_ptr, c_f_pointer
  use env_module,                     only: globalMaxLevels
  use tem_grow_array_module,          only: grw_intArray_type
  use tem_reduction_transient_module, only: tem_reduction_transient_type

  implicit none

  private

  public :: tem_indexLvl_type
  public :: tem_varSys_op_data_type

  !> Type to store the index of points of the inputs
  !! levelwise, since points are levelwise
  type tem_indexLvl_type
    type( grw_intArray_type) :: indexLvl(globalMaxLevels)
  end type tem_indexLvl_type

  !> Type which is the method_data for derived variables, hence it
  !! consists of the point index for each input variable
  !! size: number of inputs variable
  type tem_varSys_op_data_type
    type(tem_indexLvl_type), allocatable :: input_pntIndex(:)

    !> A pointer to possibly additional solver data.
    !!
    !! This is for example used to keep a link to the projection data
    !! in Ateles to enable the construction of element data from the
    !! point data for the operation variables.
    type(c_ptr) :: solver_bundle

    !> time reduction data
    type(tem_reduction_transient_type) :: redTrans
  end type tem_varSys_op_data_type

contains


end module tem_operation_module
