! Copyright (c) 2015-2016, 2018 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016, 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016, 2021 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2018 Raphael Haupt <Raphael.Haupt@student.uni-siegen.de>
! Copyright (c) 2021 Jana Gericke <jana.gericke@dlr.de>
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
!> This module encapsulates the output in VTK format.
!!
!! @todo HK: this stuff should probably be changed to use the VTK library stuff
!!           beneath, see [for example](http://www.vtk.org/doc/nightly/html/vtkXMLWriterF_8h_source.html)
!! Or maybe replace by https://github.com/szaghi/Lib_VTK_IO
!! Using the VTK API allows for greater flexibility and usage of more enhanced
!! features.
!!
!! The current workflow is as follows:
!! call [[hvs_vtk_config_load]] to get relevant user configurations.
!! Open the vtk files with [[hvs_vtk_open]].
!! Write the mesh data by calling [[hvs_vtk_write_meshdata]], this requires the
!! vertex information of the mesh, that needs to be obtained beforehand.
!! If there is data to be attached to the mesh data, describe it by calling
!! [[hvs_vtk_write_varSys]], which requires a variable system and a list of
!! variable names to write out of that variable system.
!! Then go on and call [[hvs_vtk_dump_data]] for each variable.
!! After all variables have been written, call [[hvs_vtk_close]] to close the
!! vtk files.
module hvs_vtk_module
  use, intrinsic :: iso_c_binding

  use aotus_module, only: flu_State, aot_get_val

  use env_module, only: PathLen, LabelLen, newUnit, double_k, rk, pathSep,&
    &                   isLittleEndian

  use hvs_base64_module, only: convert_to_base64

  use treelmesh_module,        only: treelmesh_type
  use tem_aux_module,          only: tem_abort, tem_open
  use tem_comm_env_module,     only: tem_comm_env_type
  use tem_logging_module,      only: logunit
  use tem_subtree_type_module, only: tem_subtree_type
  use tem_time_module,         only: tem_time_type,      &
    &                                tem_time_sim_stamp, &
    &                                tem_time_iter_stamp
  use tem_tools_module,        only: upper_to_lower
  use tem_varsys_module,       only: tem_varsys_type
  use tem_vrtx_module,         only: tem_vrtx_type

  use hvs_vtk_type_module, only: hvs_vtk_config_type, hvs_vtk_file_type

  implicit none

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
    character(len=labelLen) :: dformat
    integer :: iError
    ! ----------------------------------------------------------------------!

    call aot_get_val( L       = conf,       &
      &               thandle = thandle,    &
      &               key     = 'dataform', &
      &               val     = dformat,    &
      &               default = 'binary',   &
      &               ErrCode = iError      )

    dformat = adjustl(dformat)
    dformat = upper_to_lower(dformat)
    select case(trim(dformat))
    case('ascii','binary')
      me%dataform = trim(dformat)
    case default
      write(logunit(0),*) 'ERROR in hvs_vtk_config_load: ' &
        &                 // 'unknown data format '        &
        &                 // trim(dformat)
      write(logunit(0),*) 'dataform has to be one of:'
      write(logunit(0),*) '* ascii'
      write(logunit(0),*) '* binary'
      write(logunit(0),*)
      call tem_abort()
    end select

    call aot_get_val( L       = conf,             &
      &               thandle = thandle,          &
      &               key     = 'iter_filename',  &
      &               val     = me%iter_filename, &
      &               default = .false.,          &
      &               ErrCode = iError            )

    call aot_get_val( L       = conf,         &
      &               thandle = thandle,      &
      &               key     = 'write_pvd',  &
      &               val     = me%write_pvd, &
      &               default = .true.,       &
      &               ErrCode = iError        )

  end subroutine hvs_vtk_config_load
  ! ----------------------------------------------------------------------------!
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
    character(len=PathLen) :: headerline
    character :: linebreak
    character(len=labelLen) :: byte_order
    logical :: pvd_opened_already
    ! --------------------------------------------------------------------------!

    ! Copy the dataform for later usage without the vtk_config.
    vtk_file%dataform = vtk_config%dataform

    ! Assume no cell data, until actual celldata has been written.
    vtk_file%has_celldata = .false.

    ! Store the basename for later retrieval.
    vtk_file%basename = basename

    if ( isLittleEndian ) then
      byte_order = 'LittleEndian'
    else
      byte_order = 'BigEndian'
    end if

    ! Write pvtu file only when nProc > 1
    vtk_file%write_pvtu = ( (proc%comm_size > 1)   &
      &                     .and. (proc%rank == 0) )

    ! Write PVD file header from root and only if it is not
    ! deactivated in output table
    vtk_file%write_pvd = (proc%rank == 0) .and. vtk_config%write_pvd

    if ( vtk_file%write_pvd ) then
      inquire( file   = trim(basename)//'.pvd', &
        &      opened = pvd_opened_already      )

      if (.not. pvd_opened_already) then
        ! Only need to open pvd file if it was not already opened
        ! before
        call tem_open( newunit = vtk_file%pvdunit,       &
          &            file    = trim(basename)//'.pvd', &
          &            action  = 'write',                &
          &            status  = 'replace',              &
          &            form    = 'unformatted',          &
          &            access  = 'stream'                )

        linebreak = new_line('x')

        write(headerline,'(a)') '<?xml version="1.0"?>'
        write(vtk_file%pvdunit) trim(headerline)//linebreak

        write(headerline,'(a)') '<VTKFile type="Collection" ' &
          &                     // 'version="0.1" byte_order="'&
          &                     //trim(byte_order)//'">'
        write(vtk_file%pvdunit) trim(headerline)//linebreak

        write(headerline,'(a)') '<Collection>'
        write(vtk_file%pvdunit) trim(headerline)//linebreak
        flush(vtk_file%pvdunit)
      else
        ! PVD file already opened, get the unit it is connected to.
        inquire( file   = trim(basename)//'.pvd', &
          &      number = vtk_file%pvdunit        )
      end if
    end if

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
  subroutine hvs_vtk_open(vtk_file, use_iter, proc, time)
    !> The file description to open.
    type(hvs_vtk_file_type), intent(inout) :: vtk_file

    !> User specified settings for the output
    ! type(hvs_vtk_config_type), intent(in) :: vtk_config
    !> Whether to use iteration as part of filename
    logical, intent(in) :: use_iter

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
    integer :: pos
    character(len=labelLen) :: byte_order
    ! ----------------------------------------------------------------------!

    if ( isLittleEndian ) then
      byte_order = 'LittleEndian'
    else
      byte_order = 'BigEndian'
    end if

    if (proc%comm_size > 1) then
      write(filename,'(a,i6.6)') trim(vtk_file%basename) // '_p', proc%rank
    else
      write(filename,'(a)') trim(vtk_file%basename)
    end if

    vtk_file%timestamp = ''
    if (present(time)) then
      if ( use_iter ) then
        write(vtk_file%timestamp, '(a)') &
          & '_t' // trim(tem_time_iter_stamp(time))
      else
        write(vtk_file%timestamp, '(a)') &
          & '_t' // trim(tem_time_sim_stamp(time))
      end if
    end if

    write(filename,'(a)') trim(filename) // trim(vtk_file%timestamp) // '.vtu'

    if (.not.vtk_file%write_pvtu) then
      write(logunit(3),*) 'Opening VTU file: ' // trim(filename)
      vtk_file%last_opened_file = trim(filename)
    end if

    ! Open the file (always unformatted stream, for ascii output numbers will
    ! be converted to strings before writing).
    call tem_open( newunit = vtk_file%outunit, &
      &            file    = trim(filename),   &
      &            action  = 'write',          &
      &            status  = 'replace',        &
      &            form    = 'unformatted',    &
      &            access  = 'stream'          )

    linebreak = new_line('x')

    ! Writing the header to the VTU file:
    write(headerline,'(a)') '<?xml version="1.0"?>'
    write(vtk_file%outunit) trim(headerline)//linebreak

    write(headerline,'(a)') '<VTKFile type="UnstructuredGrid" version="0.1"' &
      &                     //  ' byte_order="'&
      &                     //trim(byte_order)//'">'
    write(vtk_file%outunit) trim(headerline)//linebreak

    write(headerline,'(a)') ' <UnstructuredGrid>'
    write(vtk_file%outunit) trim(headerline)//linebreak

    ! Add point in time information to the VTU file.
    if (present(time)) then
      write(headerline,'(a)') '<FieldData>'
      write(vtk_file%outunit) trim(headerline)//linebreak

      write(headerline,'(a)') '<DataArray type="Float64" Name="TIME" '// &
        &                     'NumberOfTuples="1" format="ascii">'
      write(vtk_file%outunit) trim(headerline)//linebreak

      write(headerline,*) time%sim
      write(vtk_file%outunit) trim(headerline)//linebreak

      write(headerline,'(a)') '</DataArray>'
      write(vtk_file%outunit) trim(headerline)//linebreak

      write(headerline,'(a)') '<DataArray type="Int32" Name="CYCLE" '// &
        &                     'NumberOfTuples="1" format="ascii">'
      write(vtk_file%outunit) trim(headerline)//linebreak

      write(headerline,*) time%iter
      write(vtk_file%outunit) trim(headerline)//linebreak

      write(headerline,'(a)') '</DataArray>'
      write(vtk_file%outunit) trim(headerline)//linebreak

      write(headerline,'(a)') '</FieldData>'
      write(vtk_file%outunit) trim(headerline)//linebreak
    end if

    ! Open the PVTU file if necessary
    if (vtk_file%write_pvtu) then
      write(filename,'(a)') trim(vtk_file%basename)         &
        &                 //trim(vtk_file%timestamp)//'.pvtu'
      write(logunit(3),*) 'Opening PVTU file: ' // trim(filename)
      vtk_file%last_opened_file = trim(filename)

      call tem_open( newunit = vtk_file%punit, &
        &            file    = trim(filename), &
        &            action  = 'write',        &
        &            status  = 'replace',      &
        &            form    = 'unformatted',  &
        &            access  = 'stream'        )

      write(headerline,'(a)') '<?xml version="1.0"?>'
      write(vtk_file%punit) trim(headerline)//linebreak

      write(headerline,'(a)') '<VTKFile type="PUnstructuredGrid" ' &
        &                     // 'version="0.1" byte_order="'&
        &                     //trim(byte_order)//'">'
      write(vtk_file%punit) trim(headerline)//linebreak

      write(headerline,'(a)') ' <PUnstructuredGrid GhostLevel="0">'
      write(vtk_file%punit) trim(headerline)//linebreak
      write(vtk_file%punit) linebreak
    end if

    ! Append current filename in pvd file.
    ! if pvtu is present, write pvtu filename else vtu filename
    if ( vtk_file%write_pvd ) then
      pos = INDEX(trim(filename), pathSep, .true.)
      write(headerline,'(a)') '  <DataSet timestep="' &
        &                     //trim(tem_time_sim_stamp(time))//'" file="' &
        &                     //trim(filename(pos+1:))//'"/>'
      write(vtk_file%pvdunit) trim(headerline)//linebreak
      flush(vtk_file%pvdunit)
    end if
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
    real(kind=c_double), target, allocatable :: indata(:)
    integer(kind=c_int), target, allocatable :: indata_int(:)
    integer(kind=c_int_least8_t), target, allocatable :: indata_int8(:)
    character(len=PathLen) :: headerline
    character :: linebreak
    logical :: use_ascii
    integer :: iVrtx, iElem
    ! ----------------------------------------------------------------------!

    linebreak = new_line('x')
    use_ascii = (trim(vtk_file%dataform) == 'ascii')

    write(headerline,'(a,i0,a,i0,a)') '  <Piece NumberOfPoints="',         &
      &                               vrtx%nVertices, '" NumberOfCells="', &
      &                               nElems, '">'
    write(vtk_file%outunit) trim(headerline)//linebreak
    write(vtk_file%outunit) linebreak

    ! Write the point data:
    write(headerline,'(a)') '   <Points>'
    write(vtk_file%outunit) trim(headerline)//linebreak
    write(headerline,'(a)') '    <DataArray Name="xyz" '                   &
      &                     // 'NumberOfComponents="3" type="Float64" '    &
      &                     // 'format="' // trim(vtk_file%dataform) // '">'
    write(vtk_file%outunit) trim(headerline)//linebreak
    if (use_ascii) then
      do ivrtx = 1, vrtx%nvertices
        write(headerline,*) real(vrtx%coord%val(:,ivrtx), kind = double_k)
        write(vtk_file%outunit) trim(headerline)//linebreak
      end do
    else
      allocate(indata(vrtx%nVertices*3))
      do ivrtx = 1, vrtx%nvertices
        indata(1+(ivrtx-1)*3:ivrtx*3) = vrtx%coord%val(:,ivrtx)
      end do

      call convert_to_base64( indata, vrtx%nVertices*3, vtk_file%outunit )
      write(vtk_file%outunit) linebreak

      deallocate(indata)
    end if
    write(headerline,'(a)') '    </DataArray>'
    write(vtk_file%outunit) trim(headerline)//linebreak
    write(headerline,'(a)') '   </Points>'
    write(vtk_file%outunit) trim(headerline)//linebreak
    write(vtk_file%outunit) linebreak

    ! Write the cell data:
    write(headerline,'(a)') '   <Cells>'
    write(vtk_file%outunit) trim(headerline)//linebreak
    write(headerline,'(a)') '    <DataArray type="Int32" Name="connectivity"' &
      &                     //' format="'//trim(vtk_file%dataform)//'">'
    write(vtk_file%outunit) trim(headerline)//linebreak

    if (use_ascii) then
      do iElem = 1, nElems
        write(headerline,*) vrtx%map2global(iElem,:) - 1
        write(vtk_file%outunit) trim(headerline)//linebreak
      end do
    else
      allocate(indata_int(vrtx%maxVertices))
      do iElem = 1, nElems
        indata_int(                                                        &
        &    (iElem-1)*vtk_file%vtx_per_Elem +1 : iElem*vtk_file%vtx_per_Elem) &
        &  = vrtx%map2global(iElem,:) - 1
      end do
      call convert_to_base64( indata_int, vrtx%maxVertices, vtk_file%outunit )

      deallocate(indata_int)
      write(vtk_file%outunit) linebreak
    end if

    write(headerline,'(a)') '    </DataArray>'
    write(vtk_file%outunit) trim(headerline)//linebreak

    ! Write offsets:
    write(headerline,'(a)') '    <DataArray type="Int32" Name="offsets"' &
      &                     //' format="'//trim(vtk_file%dataform)//'">'
    write(vtk_file%outunit) trim(headerline)//linebreak
    if (use_ascii) then
      do iElem = 1, nElems
        write(headerline,*) iElem*vtk_file%vtx_per_Elem
        write(vtk_file%outunit) trim(headerline)
      end do
      write(vtk_file%outunit) linebreak
    else
      allocate(indata_int(nElems))
      do iElem = 1, nElems
        indata_int(iElem) = iElem*vtk_file%vtx_per_Elem
      end do
      call convert_to_base64( indata_int, nElems, vtk_file%outunit )
      deallocate(indata_int)
      write(vtk_file%outunit) linebreak
    end if
    write(headerline,'(a)') '    </DataArray>'
    write(vtk_file%outunit) trim(headerline)//linebreak

    ! Write the cell type.
    write(headerline,'(a)') '    <DataArray type="UInt8" Name="types"' &
      &                     // ' format="'//trim(vtk_file%dataform)//'">'
    write(vtk_file%outunit) trim(headerline)//linebreak
    if (use_ascii) then
      do iElem = 1, nElems
        write(headerline,*) vtk_file%CellType
        write(vtk_file%outunit) trim(headerline)
      end do
      write(vtk_file%outunit) linebreak
    else
      allocate(indata_int8(nElems))
      indata_int8 = int(vtk_file%CellType, c_int_least8_t)
      call convert_to_base64( indata_int8, nElems, vtk_file%outunit )
      deallocate(indata_int8)
      write(vtk_file%outunit) linebreak
    end if
    write(headerline,'(a)') '    </DataArray>'
    write(vtk_file%outunit) trim(headerline)//linebreak

    write(headerline,'(a,i0,a)') '   </Cells>'
    write(vtk_file%outunit) trim(headerline)//linebreak
    write(vtk_file%outunit) linebreak

    ! Add mesh information to the pvtu file accordingly:
    if (vtk_file%write_pvtu) then
      write(headerline,'(a)') '  <PPoints>'
      write(vtk_file%punit) trim(headerline)//linebreak
      write(headerline,'(a)') '   <PDataArray Name="xyz" '                &
        &                     // 'NumberOfComponents="3" type="Float64" ' &
        &                     // 'format="' // trim(vtk_file%dataform) // '"/>'
      write(vtk_file%punit) trim(headerline)//linebreak
      write(headerline,'(a)') '  </PPoints>'
      write(vtk_file%punit) trim(headerline)//linebreak
      write(vtk_file%punit) linebreak

      write(headerline,'(a)') '  <PCells>'
      write(vtk_file%punit) trim(headerline)//linebreak
      write(headerline,'(a)') '   <PDataArray type="Int32" '    &
        &                     // 'Name="connectivity" format="' &
        &                     // trim(vtk_file%dataform) // '"/>'
      write(vtk_file%punit) trim(headerline)//linebreak
      write(headerline,'(a)') '   <PDataArray type="Int32" ' &
        &                     // 'Name="offsets" format="'   &
        &                     // trim(vtk_file%dataform) // '"/>'
      write(vtk_file%punit) trim(headerline)//linebreak
      write(headerline,'(a)') '   <PDataArray type="UInt8" ' &
        &                     // 'Name="types" format="'     &
        &                     // trim(vtk_file%dataform) // '"/>'
      write(vtk_file%punit) trim(headerline)//linebreak
      write(headerline,'(a)') '  </PCells>'
      write(vtk_file%punit) trim(headerline)//linebreak
      write(vtk_file%punit) linebreak
    end if

  end subroutine hvs_vtk_write_meshdata
  ! ----------------------------------------------------------------------------!
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  !> Write the cell data description into the VTK files.
  !!
  !! This will write the celldata string into the VTK files.
  subroutine hvs_vtk_write_cd_Header(vtk_file, celldata_str)
    !> Handles for the VTK files to write the celldata to.
    type(hvs_vtk_file_type), intent(in) :: vtk_file

    !> Actual string to put into the XML to describe the celldata.
    character(len=*), intent(in) :: celldata_str
    ! ----------------------------------------------------------------------!
    character :: linebreak
    character(len=PathLen) :: headerline
    ! ----------------------------------------------------------------------!

    linebreak = new_line('x')

    ! cell based vtk output
    write(headerline,'(a)') '   <CellData'//trim(celldata_str)//'>'
    write(vtk_file%outunit) trim(headerline)//linebreak
    if (vtk_file%write_pvtu) then
      write(headerline,'(a)') '  <PCellData'//trim(celldata_str)//'>'
      write(vtk_file%punit) trim(headerline)//linebreak
    end if

  end subroutine hvs_vtk_write_cd_Header
  ! ----------------------------------------------------------------------------!
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
    integer, intent(in) :: varpos(:)
    ! ----------------------------------------------------------------------!
    character :: linebreak
    character(len=PathLen) :: scalars
    character(len=PathLen) :: vectors
    integer :: scal_count, vect_count
    integer :: iVar, ipos
    ! ----------------------------------------------------------------------!

    linebreak = new_line('x')

    scalars = ''
    vectors = ''
    scal_count = 0
    vect_count = 0

    do iVar=1,size(varpos)
      ipos = varpos(iVar)
      if (varsys%method%val(ipos)%nComponents == 1) then
        write(scalars, '(a)') adjustl(trim(scalars)//' ' &
          &                   // trim(varSys%varname%val(iPos)))
        scal_count = scal_count + 1
      else
        write(vectors, '(a)') adjustl(trim(vectors)//' ' &
          &                   // trim(varSys%varname%val(iPos)))
        vect_count = vect_count + 1
      end if
    end do

    if (scal_count > 0) then
      write(scalars, '(a)') ' Scalars="'//trim(scalars)//'"'
    end if
    if (vect_count > 0) then
      write(vectors, '(a)') ' Vectors="'//trim(vectors)//'"'
    end if

    if (scal_count+vect_count > 0) then
      call hvs_vtk_write_cd_Header( vtk_file     = vtk_file,        &
        &                           celldata_str = trim(scalars)    &
        &                                          // trim(vectors) )

      vtk_file%has_celldata = .true.
    end if

  end subroutine hvs_vtk_write_varSys
  ! ----------------------------------------------------------------------------!
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  !> Dump the given data (input) with the given name in the given format (vtu)
  !! to the given unit.
  !!
  subroutine hvs_vtk_dump_data( vtk_file, varpos, varSys, mesh, time, subtree )
    ! --------------------------------------------------------------------------!
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
    ! --------------------------------------------------------------------------!
    integer :: nElems
    integer :: iElem
    integer :: nComp
    integer :: epos(1)
    character(len=PathLen) :: headerline
    character :: linebreak
    character(len=LabelLen) :: dataname
    real(kind=c_double), target, allocatable :: indata(:)
    integer :: spatialComp
    logical :: use_ascii
    ! --------------------------------------------------------------------------!

    nComp = varsys%method%val(varPos)%nComponents
    dataname = varsys%varname%val(varPos)
    if (present(subtree)) then
      nElems = subtree%nElems
    else
      nElems = mesh%nElems
    end if

    linebreak = new_line('x')
    use_ascii = (trim(vtk_file%dataform) == 'ascii')

    spatialComp = nComp

    ! Promote 2D data to 3D, as we always create 3D meshes.
    ! 3rd component will be set to 0.
    if (nComp == 2) then
      spatialComp = 3
    end if

    ! Open the data array.
    write(headerline,'(a,i0,a)') '    <DataArray type="Float64" Name="'        &
      &                          // trim(dataname)                             &
      &                          // '" NumberOfComponents="', spatialComp,     &
      &                          '" format="' // trim(vtk_file%dataform) // '">'
    write(vtk_file%outunit) trim(headerline)//linebreak

    if (vtk_file%write_pvtu) then
      write(headerline,'(a,i0,a)') '    <PDataArray type="Float64" Name="'   &
        &                          // trim(dataname)                         &
        &                          // '" NumberOfComponents="', spatialComp, &
        &                          '" format="' // trim(vtk_file%dataform)   &
        &                          //'"/>'
      write(vtk_file%punit) trim(headerline)//linebreak
    end if

    ! Dump the actual values into the data array.
    if (use_ascii) then

      allocate(indata(spatialComp))
      indata = 0.0_rk
      do iElem = 1, nElems
        if (present(subtree)) then
          epos = subtree%map2global(iElem)
        else
          epos = iElem
        end if

        call varsys%method%val(varpos)                    &
          &        %get_element( varsys = varsys,         &
          &                      elempos = epos,          &
          &                      time    = time,          &
          &                      tree    = mesh,          &
          &                      nElems  = 1,             &
          &                      nDofs   = 1,             &
          &                      res     = indata(:nComp) )
        write(headerline,*) indata
      end do

    else

      allocate(indata(nElems*spatialComp))
      if (nComp == 2) then
        do iElem = 1, nElems
          if (present(subtree)) then
            epos = subtree%map2global(iElem)
          else
            epos = iElem
          end if
          call varsys%method%val(varpos)                          &
            &        %get_element( varsys  = varsys,              &
            &                      elempos = epos,                &
            &                      time    = time,                &
            &                      tree    = mesh,                &
            &                      nElems  = 1,                   &
            &                      nDofs   = 1,                   &
            &                      res     = indata(1+(iElem-1)*3 &
            &                                       :iElem*3-1)   )
          indata(iElem*3) = 0.0_c_double
        end do
      else
        do iElem = 1, nElems
          if (present(subtree)) then
            epos = subtree%map2global(iElem)
          else
            epos = iElem
          end if
          call varsys%method%val(varpos)                              &
            &        %get_element( varsys  = varsys,                  &
            &                      elempos = epos,                    &
            &                      time    = time,                    &
            &                      tree    = mesh,                    &
            &                      nElems  = 1,                       &
            &                      nDofs   = 1,                       &
            &                      res     = indata(1+(iElem-1)*nComp &
            &                                       :iElem*nComp)     )
        end do
      end if
      call convert_to_base64( indata, nElems*spatialComp, vtk_file%outunit )
      write(vtk_file%outunit) linebreak
      deallocate(indata)

    end if

    ! Close the data array again.
    write(headerline,'(a)') '    </DataArray>'
    write(vtk_file%outunit) trim(headerline)//linebreak

  end subroutine hvs_vtk_dump_data
  ! ----------------------------------------------------------------------------!
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
    character(len=PathLen) :: headerline
    character :: linebreak
    character(len=PathLen) :: filename
    integer :: irank
    integer :: pos
    ! ----------------------------------------------------------------------!

    linebreak = new_line('x')

    if (vtk_file%has_celldata) then
      write(headerline,'(a)') '   </CellData>'
      write(vtk_file%outunit) trim(headerline)//linebreak
      write(vtk_file%outunit) linebreak
      if (vtk_file%write_pvtu) then
        write(headerline,'(a)') '  </PCellData>'
        write(vtk_file%punit) trim(headerline)//linebreak
        write(vtk_file%punit) linebreak
      end if
    end if

    write(headerline,'(a)') '  </Piece>'
    write(vtk_file%outunit) trim(headerline)//linebreak
    write(headerline,'(a)') ' </UnstructuredGrid>'
    write(vtk_file%outunit) trim(headerline)//linebreak
    write(headerline,'(a)') '</VTKFile>'
    write(vtk_file%outunit) trim(headerline)//linebreak

    close(vtk_file%outunit)

    ! write pvtu file from root process
    if (vtk_file%write_pvtu) then

      pos = INDEX(trim(vtk_file%basename), pathSep, .true.)
      do irank = 0, proc%comm_size-1
        write(filename,'(a,i6.6)') trim(vtk_file%basename(pos+1:)) &
          &                        // '_p', irank
        write(filename,'(a)') trim(filename) // trim(vtk_file%timestamp) &
          &                   // '.vtu'
        write(headerline,'(a)') '<Piece Source="'//trim(filename)//'"/>'
        write(vtk_file%punit) trim(headerline)//linebreak
      end do

      write(headerline,'(a)') '</PUnstructuredGrid>'
      write(vtk_file%punit) trim(headerline)//linebreak
      write(headerline,'(a)') '</VTKFile>'
      write(vtk_file%punit) trim(headerline)//linebreak

      !close the unit
      close(vtk_file%punit)

    end if

  end subroutine hvs_vtk_close
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  ! ----------------------------------------------------------------------------!
  !> This routine closes PVD file
  subroutine hvs_vtk_closePVD(vtk_file)
    ! --------------------------------------------------------------------------!
    !> The file description to open.
    type(hvs_vtk_file_type), intent(inout) :: vtk_file
    ! --------------------------------------------------------------------------!
    character(len=PathLen) :: headerline
    character :: linebreak
    ! --------------------------------------------------------------------------!
    if ( vtk_file%write_pvd ) then
      linebreak = new_line('x')

      write(headerline,'(a)') '</Collection>'
      write(vtk_file%pvdunit) trim(headerline)//linebreak

      write(headerline,'(a)') '</VTKFile>'
      write(vtk_file%pvdunit) trim(headerline)//linebreak
      close(vtk_file%pvdunit)
    end if
  end subroutine hvs_vtk_closePVD
  ! ----------------------------------------------------------------------------!


  ! *****************************************************************************
  !> Dumps the debug_data into a file.
  !!
  !! This routine takes a one-dimensional array (one value per element, ordered
  !! liek the treeIDs) and stores it into a vtk file. The name of the vtk-file
  !! is debug_dump, the values are called debug_value.
  !!
  !! The caller has to provide the vrtx, which can be created with the
  !! tem_calc_vrtx_coord-routine.
  subroutine hvs_dump_debug_array( proc, tree, time, vrtx, debug_data)
    !---------------------------------------------------------------------------
    type(tem_comm_env_type), intent(in) :: proc
    type(treelmesh_type), intent(in) :: tree
    type(tem_time_type), intent(in) :: time
    type(tem_vrtx_type) :: vrtx
    real(kind=rk) :: debug_data(tree%nElems)
    !---------------------------------------------------------------------------
    type(hvs_vtk_file_type) :: vtk_file
    type(hvs_vtk_config_type) :: vtk_config
    character :: linebreak
    character(len=PathLen) :: headerline
    !---------------------------------------------------------------------------

    vtk_config%dataform = 'binary'
    vtk_config%write_pvd = .false.

    call hvs_vtk_init( vtk_file   = vtk_file,     &
      &                vtk_config = vtk_config,   &
      &                basename   = 'debug_dump', &
      &                proc       = proc          )
    vtk_file%has_celldata = .true.


    call hvs_vtk_open( vtk_file = vtk_file,   &
      &                use_iter = .true.,     &
      &                proc     = proc,       &
      &                time     = time        )

    call hvs_vtk_write_meshdata( vtk_file = vtk_file,   &
      &                          vrtx     = vrtx,       &
      &                          nElems   = tree%nElems )


    linebreak = new_line('x')
    write(vtk_file%outunit) '<CellData Scalars="debug_value">'//linebreak

    write(headerline,'(a,i0,a)') '    <DataArray type="Float64" ' &
      & // ' Name="debug_value"'                                  &
      & // ' NumberOfComponents="', 1, '"'                        &
      & // ' format="' // trim(vtk_file%dataform) // '">'
    write(vtk_file%outunit) trim(headerline) // linebreak

    call convert_to_base64( debug_data, tree%nElems, vtk_file%outunit )
    write(vtk_file%outunit) linebreak

    ! Close the data array again.
    write(headerline,'(a)') '    </DataArray>'
    write(vtk_file%outunit) trim(headerline) // linebreak

    call hvs_vtk_close( vtk_file, proc)

  end subroutine hvs_dump_debug_array
  ! *****************************************************************************

end module hvs_vtk_module
