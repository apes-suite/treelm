! Copyright (c) 2015 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> This module encapsulates the dummy routines if executable is build with
!! --no_vtk.
!!
!! Actual routines are in hvs_vtk_module.f90
module hvs_vtk_module
  use, intrinsic :: iso_c_binding

  use aotus_module, only: flu_State

  use env_module, only: PathLen, LabelLen

  use treelmesh_module, only: treelmesh_type
  use tem_aux_module, only: tem_abort
  use tem_comm_env_module, only: tem_comm_env_type
  use tem_logging_module, only: logunit
  use tem_subtree_type_module, only: tem_subtree_type
  use tem_time_module, only: tem_time_type
  use tem_varsys_module, only: tem_varsys_type
  use tem_vrtx_module, only: tem_vrtx_type

  use hvs_vtk_type_module, only: hvs_vtk_config_type, hvs_vtk_file_type

  implicit none

  public :: hvs_vtk_config_load, hvs_vtk_init, hvs_vtk_open

contains

  ! ----------------------------------------------------------------------------!
  !> Read the VTK output configuration from a Lua script.
  subroutine hvs_vtk_config_load(me, conf, thandle)
    !> The VTK configuration settings to fill.
    type(hvs_vtk_config_type), intent(out) :: me

    !> Handle of the Lua script to load the configuration from.
    type(flu_state) :: conf

    !> Table handle to the table providing the VTK settings.
    integer, intent(in) :: thandle
    ! ----------------------------------------------------------------------!
    write(logUnit(0),*) 'Using dummy hvs_vtk_config_load because'
    write(logUnit(0),*) 'Executable is configured with --no_vtk'
    call tem_abort()

  end subroutine hvs_vtk_config_load
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  !> Initialize the type for VTK format
  subroutine hvs_vtk_init(vtk_file, vtk_config, basename, proc)
    ! --------------------------------------------------------------------------!
    !> The file description to open.
    type(hvs_vtk_file_type), intent(inout) :: vtk_file

    !> User specified settings for the output
    type(hvs_vtk_config_type), intent(in) :: vtk_config

    !> Basename for the output file, rank and suffix will be appended as
    !! needed.
    character(len=*), intent(in) :: basename

    !> Parallel environment to use for  the output.
    type(tem_comm_env_type), intent(in) :: proc
    ! --------------------------------------------------------------------------!

    write(logUnit(0),*) 'VTK output deactivated in this executable!'
    call tem_abort()

  end subroutine hvs_vtk_init
  ! ----------------------------------------------------------------------------!
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  !> Open the output files in VTK format.
  !!
  !! This will open VTU files and if multiple processes are used a PVTU file.
  !! We always write unstructured meshes, so we also write the header for the
  !! unstructured mesh here already.
  !! The actual mesh data is then to be written by hvs_vtk_write_meshdata.
  subroutine hvs_vtk_open(vtk_file, vtk_config, proc, time)
    !> The file description to open.
    type(hvs_vtk_file_type), intent(out) :: vtk_file

    !> User specified settings for the output
    type(hvs_vtk_config_type), intent(in) :: vtk_config

    !> Parallel environment to use for  the output.
    type(tem_comm_env_type), intent(in) :: proc

    !> Time information.
    !!
    !! If this is present, the filename will be built with a time stamp and
    !! the time point information is written into the vtu file.
    type(tem_time_type), intent(in), optional :: time
    ! ----------------------------------------------------------------------!
    character(len=PathLen) :: filename
    character(len=PathLen) :: headerline
    character :: linebreak
    ! ----------------------------------------------------------------------!

    write(logUnit(0),*) 'VTK output deactivated in this executable!'
    call tem_abort()

  end subroutine hvs_vtk_open
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  !> Write the mesh information into the VTK output files.
  !!
  !! Note: the unstructured meshes we write are consisting of voxels, this
  !!       assumption is built into the code and exploited.
  subroutine hvs_vtk_write_meshdata(vtk_file, vrtx, nElems)
    !> File handles to the files where the mesh data should be written to.
    type(hvs_vtk_file_type), intent(in) :: vtk_file

    !> Information on the vertices of the mesh
    type(tem_vrtx_type), intent(in) :: vrtx

    !> Number of elements in the mesh
    integer, intent(in) :: nElems
    ! ----------------------------------------------------------------------!
    write(logUnit(0),*) 'Using dummy hvs_vtk_write_meshdata because'
    write(logUnit(0),*) 'executable is configured with --no_vtk'
    call tem_abort()
  end subroutine hvs_vtk_write_meshdata
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  !> Convert the provided variable system data into celldata description in the
  !! given vtk files.
  !!
  subroutine hvs_vtk_write_varSys(vtk_file, varsys, varpos)
    !> Output info for vtu_output.
    type(hvs_vtk_file_type), intent(inout) :: vtk_file

    !> Description of the available variable system to get the given varnames
    !! from.
    type(tem_varSys_type), intent(in) :: varsys

    !> List of variable positions that should be written into the vtk output.
    !!
    !! If this is not provided, all variables from the varsys will be written
    !! to the vtk file.
    integer, optional, intent(in) :: varpos(:)
    ! ----------------------------------------------------------------------!
    write(logUnit(0),*) 'Using dummy hvs_vtk_write_varSys because'
    write(logUnit(0),*) 'executable is configured with --no_vtk'
    call tem_abort()

  end subroutine hvs_vtk_write_varSys
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  !> Dump the given data (input) with the given name in the given format (vtu)
  !! to the given unit.
  !!
  subroutine hvs_vtk_dump_data( vtk_file, varpos, varSys, mesh, time, subtree )
    !> VTK file to write data to.
    type(hvs_vtk_file_type), intent(in) :: vtk_file

    !> Position of the variable to write
    integer, intent(in) :: varpos

    !> Description of the available variable system to get the given varnames
    !! from.
    type(tem_varSys_type), intent(in) :: varsys

    !> Mesh to write the data on.
    type(treelmesh_type), intent(in) :: mesh

    !> Point in time to use for this data.
    !!
    !! Can be important for space-time function evaluations.
    type(tem_time_type), intent(in) :: time

    !> Optional restriction of the elements to output.
    type(tem_subtree_type), optional, intent(in) :: subtree
    ! ----------------------------------------------------------------------!
    write(logUnit(0),*) 'Using dummy hvs_vtk_dump_data because'
    write(logUnit(0),*) 'executable is configured with --no_vtk'
    call tem_abort()
  end subroutine hvs_vtk_dump_data
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  !> This routine finalizes the vtu file i.e closing cellData xml and
  !! creating pvtu file to combile all parallel vtu files
  subroutine hvs_vtk_close( vtk_file, proc )
    !> The file descriptor to close again.
    type(hvs_vtk_file_type), intent(in) :: vtk_file
    !> Communicator for the parallel environment.
    type(tem_comm_env_type), intent(in) :: proc
    ! ----------------------------------------------------------------------!
    write(logUnit(0),*) 'Using dummy hvs_vtk_close because'
    write(logUnit(0),*) 'executable is configured with --no_vtk'
    call tem_abort()
  end subroutine hvs_vtk_close

  ! ----------------------------------------------------------------------------!
  !> This routine finalizes the vtu file i.e closing cellData xml and
  !! creating pvtu file to combile all parallel vtu files
  subroutine hvs_vtk_closePVD(vtk_file, proc)
    !> The file descriptor to close again.
    type(hvs_vtk_file_type), intent(in) :: vtk_file
    !> Communicator for the parallel environment.
    type(tem_comm_env_type), intent(in) :: proc
    ! ----------------------------------------------------------------------!
    write(logUnit(0),*) 'Using dummy hvs_vtk_closePVD because'
    write(logUnit(0),*) 'executable is configured with --no_vtk'
    call tem_abort()
  end subroutine hvs_vtk_closePVD

end module hvs_vtk_module
