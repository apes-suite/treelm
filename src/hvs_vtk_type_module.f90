! Copyright (c) 2015-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2018 Raphael Haupt <Raphael.Haupt@student.uni-siegen.de>
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
!> This module encapsulates the data type for VTK output format.
!
! Data type is defined in this file to avoid code dublication of data type in
! [[hvs_vtk_module]] and [[hvs_vtk_dummy]]
module hvs_vtk_type_module
  use env_module, only: PathLen, LabelLen

  implicit none

  !> Configuration of the VTK output.
  !!
  !! These are the settings, that can be configured by the
  !! user with respect to the VTK output.
  type hvs_vtk_config_type
    !> Format in which the data is present (either 'ascii' or 'binary')
    character(len=labelLen) :: dataform

    !> If a timestep output is done, the filename can either be built
    !! with the simulation time or the iteration.
    !!
    !! If this flag is set to true, iterations will be used, otherwise
    !! the simulation time (which is the default).
    logical :: iter_filename

    !> Flag to decided whether to write pvd file or not.
    !! Default is true.
    logical :: write_pvd

  end type hvs_vtk_config_type


  !> Description of the opened files for VTK output.
  type hvs_vtk_file_type
    !> File handle for the vtu file with the data.
    integer :: outunit

    !> Filehandle for the pvtu file for partitioned vtu data.
    integer :: punit

    !> Filehandle for the pvd file for time-series output
    integer :: pvdunit

    !> Basename of the VTK files to write
    character(len=pathLen) :: basename

    !> Name of the last opened file on this process
    !!
    !! If there is a pvtu written by the root process, the root process
    !! will store the name of the pvtu file here. Otherwise it contains
    !! the name of the vtu file.
    character(len=pathLen) :: last_opened_file

    !> Timestamp to construct the filename
    character(len=labelLen) :: timestamp

    !> Format in which the data is present (either 'ascii' or 'binary')
    !! todo: maybe remove this and use binary always?
    character(len=labelLen) :: dataform

    !> Flag to indicate, whether this process has to write the pvtu file.
    logical :: write_pvtu

    !> Indicator, wether celldata has been written to the VTK file.
    !!
    !! If true, the finalization will write a closing celldata tag.
    logical :: has_celldata = .false.

    !> Flag to decided whether to write pvd file or not.
    logical :: write_pvd

    !> number of vortices per Element
    !! depending on CellType
    !! check vtk-manual for more informations
    integer :: vtx_per_Elem = 8

    !> cell type used in VTK file
    !! 11 = voxel
    !!  8 = Pixel
    !!  4 = Poly_line
    !!  3 = Line
    integer :: CellType = 11

  end type hvs_vtk_file_type

end module hvs_vtk_type_module
