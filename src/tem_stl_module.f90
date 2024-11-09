! Copyright (c) 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
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
! ******************************************************************************
!> This module provides stl data type and routine to load all stl files and
!! store them in stl_node_type
!!
!! \author
!!
module tem_stl_module

  ! include treelm modules
  use env_module,                only: rk, PathLen
  use tem_aux_module,            only: tem_abort
  use tem_tools_module,          only: tem_horizontalSpacer
  use tem_logging_module,        only: logunit
  use tem_stlb_io_module,        only: tem_read_stlb, tem_size_stlb
  use stla_io,                   only: stla_read, stla_size, stla_check
  use tem_triangle_module,       only: tem_triangle_type,                      &
    &                                  grw_triangleArray_type, append,         &
    &                                  tem_triangleCubeOverlap
  use tem_transformation_module, only: tem_transformation_type
  use tem_cube_module,           only: tem_cube_type

  ! include aotus modules
  use flu_binding,               only: flu_State
  use aotus_module,              only: flu_State, aot_get_val,                 &
    &                                  aoterr_Fatal, aoterr_NonExistent,       &
    &                                  aoterr_WrongType
  use aot_table_module,          only: aot_table_open, aot_table_close,        &
    &                                  aot_table_length
  use aot_out_module,   only: aot_out_type, aot_out_val,                       &
    &                         aot_out_open_table, aot_out_close_table

  implicit none

  private

  public :: tem_load_stl, tem_stlCubeOverlap
  public :: tem_stlData_type
  public :: tem_stlHead_out

  integer, parameter:: stl_ascii = 1 !< stl ascii file identification
  integer, parameter:: stl_bin   = 2 !< stl bin file identification

  !> This type contains STL information read from the Lua configuration.
  type tem_stlHead_type
    character(len=PathLen) :: filename !< Filename of the STL File
    integer :: fileformat !< File format of the STL File
    integer :: nSolids !< number of solids
    integer :: nTris !< number of Tris
    integer :: nTexts !< number of texts
    integer :: nNodes !< number of nodes
    integer :: id !< stl header identification
  end type tem_stlHead_type


  !> Triangle information for all the STLs
  type tem_stlData_type
    !> Header information loaded from stl and config file
    type(tem_stlHead_type), allocatable :: head(:)
    !> List off all the nodes listed in the STL files
    !! size: 1st dimension = 3 (x,y,z) , 2nd dimension = nNodes
    real(kind=rk), allocatable :: nodes(:,:)
    !> Pointers to node array for each triangle
    !! tri_node(1:3, iTris) has three nodes position of each triangle
    !! size: 1st dimension = 3 (3 nodes), 2nd dimension = nTris
    integer, allocatable :: tri_node(:,:)
    !> Number of nodes
    integer :: nNodes
    !> Number of triangles
    integer :: nTris
  end type tem_stlData_type

  !> interface to write out stl filename in lua format to a file
  interface tem_stlHead_out
    module procedure tem_stlHead_out_scal
    module procedure tem_stlHead_out_vec
  end interface tem_stlHead_out

contains

  ! ****************************************************************************
  !> This routine loads stl_files table from configuration file
  !! Need to look for:
  !! * filename (mandatory)
  !! * fileformat (default binary)
  subroutine tem_load_stlHead( me, conf, thandle )
    ! --------------------------------------------------------------------------!
    !> contain stl files information from config file
    type(tem_stlHead_type), allocatable, intent(out) :: me(:)
    type(flu_State) :: conf !< lua state
    integer, intent(in) :: thandle !< stl object handle
    ! --------------------------------------------------------------------------!
    integer :: stl_handle, stl_subHandle
    integer :: nObjects, iObj
    ! --------------------------------------------------------------------------!


    call aot_table_open(L=conf, parent = thandle, thandle=stl_handle,          &
      & key='object')
    call aot_table_open(L=conf, parent = stl_handle, thandle =stl_subHandle,  &
      & pos = 1 )

    if( stl_subHandle .eq. 0 ) then
      !object is a single table
      call aot_table_close(L=conf, thandle=stl_subHandle)
      allocate(me (1) )
      call tem_load_stlHead_single( me = me(1), conf = conf,                   &
        & thandle = stl_handle )
    else
      !object is a multiple table
      call aot_table_close(L=conf, thandle=stl_subHandle)
      nObjects = aot_table_length(L=conf, thandle=stl_handle)
      allocate(me (nObjects) )
      do iObj = 1, nObjects
        call aot_table_open( L = conf, parent = stl_handle,                    &
          & thandle = stl_subhandle, pos = iObj )
        call tem_load_stlHead_single( me = me(iObj), conf = conf,              &
          & thandle = stl_subHandle )
        call aot_table_close(L=conf, thandle=stl_subHandle)
      end do
    end if

    call aot_table_close(L=conf, thandle=stl_Handle)

  end subroutine tem_load_stlHead
  ! ****************************************************************************

  ! ****************************************************************************
  !> This routine load stl_files table from configuration file
  !! Need to look for:
  !! * filename (mandatory)
  !! * fileformat (default binary)
  subroutine tem_load_stlHead_single( me, conf, thandle )
    ! --------------------------------------------------------------------------!
    !> contain stl files information from config file
    type(tem_stlHead_type), intent(inout) :: me
    type(flu_State) :: conf !< lua state
    integer, intent(in) :: thandle !< stl object handle
    ! --------------------------------------------------------------------------!
    integer :: iError
    character(len=6) :: formString
    ! --------------------------------------------------------------------------!

    me%fileformat = stl_bin

    call aot_get_val(L=conf, thandle=thandle, &
      &              val=me%filename, ErrCode=iError, &
      &              key='filename', pos=1)
    if (btest(iError, aoterr_Fatal)) then
      write(logunit(0),*) 'FATAL Error occured, while retrieving filename:'
      if (btest(iError, aoterr_NonExistent)) &
        &  write(logunit(0),*) 'Variable not existent!'
      if (btest(iError, aoterr_WrongType)) &
        &  write(logunit(0),*) 'Variable has wrong type!'
      call tem_abort()
    end if

    call aot_get_val(L=conf, thandle=thandle, &
      &              val=formString, ErrCode=iError, &
      &              key='format', pos=3, default='binary')
    if (btest(iError, aoterr_Fatal)) then
      write(logunit(0),*) 'FATAL Error occured, while retrieving fileformat:'
      if (btest(iError, aoterr_NonExistent)) &
        &  write(logunit(0),*) 'Variable not existent!'
      if (btest(iError, aoterr_WrongType)) &
        &  write(logunit(0),*) 'Variable has wrong type!'
      call tem_abort()
    else
      if (btest(iError, aoterr_NonExistent)) &
        & write(logunit(5),*) 'Variable fileformat not set in configuration, ' &
        &                     // 'using default value!'
    end if
    if (formString(1:1) == 'a') me%fileformat = stl_ascii


    select case(me%fileformat)
      case(stl_bin)
        formstring = 'binary'
      case(stl_ascii)
        formstring = 'ascii'
    end select

    write(logunit(1),"(A)") '       STL file: '//trim(me%filename)
    write(logunit(1),"(A)") '         format: '//trim(formstring)

  end subroutine tem_load_stlHead_single
  ! *****************************************************************************

  ! *****************************************************************************
  !> Read in all STL files from disk
  !! which are specified in the config file
  subroutine tem_read_stlFiles( stl_data )
    ! ---------------------------------------------------------------------------!
    !> stl data of current spatial object
    type( tem_stlData_type ), intent(inout) :: stl_data
    ! ---------------------------------------------------------------------------!
    !local variable
    integer :: nodeOffset !< Offset in the nodelist for multiple STLs
    integer :: triOffset !< Offset in the trilist for multiple STLs
    integer :: iFile, ierr
    integer :: nNodes
    integer :: nTris !< total number of triangles loaded
    ! ---------------------------------------------------------------------------!

    write(logunit(1),*) " Reading in STL Headers ..."
    ! Read headerfiles to determine number of nodes and tris
    do iFile = 1, size(stl_data%head)
      select case( stl_data%head(iFile)%fileformat )
      case( stl_ascii )
        call stla_size( trim(stl_data%head(iFile)%filename),  &
          &                   stl_data%head(iFile)%nSolids,   &
          &                   stl_data%head(iFile)%nNodes,    &
          &                   stl_data%head(iFile)%nTris,     &
          &                   stl_data%head(iFile)%nTexts)
      case( stl_bin )
        call tem_size_stlb( stl_data%head(iFile)%filename,                           &
          &                 stl_data%head(iFile)%nNodes,                             &
          &                 stl_data%head(iFile)%nTris)

      end select
    end do

    !Sum over all read in STL Files to generate one list of elements in the end
    nTris     = 0
    nNodes   = 0
    do iFile = 1, size(stl_data%head)
      nTris = nTris + stl_data%head(iFile)%nTris
      nNodes = nNodes + stl_data%head(iFile)%nNodes
    end do

    write(logunit(5),"(A,I0)") " Total number of triangles: ", nTris
    write(logunit(5),"(A,I0)") " Total number of nodes:     ", nNodes

    !allocate node and triangle arrays
    stl_data%nNodes = nNodes
    stl_data%nTris = nTris
    allocate(stl_data%nodes(1:3,nNodes))
    allocate(stl_data%tri_node(1:3,nTris))

    ! Read in node values from STL files to stl_data%nodes
    ! Store three nodes position of each triangle to stl_date%tri_node
    write(logunit(1),*) " Reading in STL Files ..."

    nodeOffset = 1
    triOffset = 1
    do iFile = 1, size(stl_data%head)
      select case(stl_data%head(iFile)%fileformat)
      case(stl_ascii)
        call stla_read(input_file_name = trim(stl_data%head(iFile)%filename), &
        &    node_num  = stl_data%head(iFile)%nNodes,  &
        &    face_num  = stl_data%head(iFile)%nTris, &
        &    node_xyz  = stl_data%nodes(1:3,nodeOffset:nodeOffset+stl_data%head(iFile)%nNodes-1), &
        &    face_node = stl_data%tri_node(1:3,triOffset:triOffset+stl_data%head(iFile)%nTris-1), &
        &    ierror    = ierr)

      case(stl_bin)
        call tem_read_stlb(filename = stl_data%head(iFile)%filename, &
        &    nNodesRead = stl_data%head(iFile)%nNodes,  &
        &    nTrisRead  = stl_data%head(iFile)%nTris, &
        &    nodes  = stl_data%nodes(1:3,nodeOffset:nodeOffset+stl_data%head(iFile)%nNodes-1), &
        &    tri_node = stl_data%tri_node(1:3,triOffset:triOffset+stl_data%head(iFile)%nTris-1), &
        &    ierror    = ierr)
      end select

      stl_data%tri_node(1:3,triOffset:triOffset+stl_data%head(iFile)%nTris-1) =  &
      & stl_data%tri_node(1:3,triOffset:triOffset+stl_data%head(iFile)%nTris-1) + (nodeOffset-1)

      nodeOffset = nodeOffset + stl_data%head(iFile)%nNodes
      triOffset = triOffset + stl_data%head(iFile)%nTris

    end do

    write(logunit(1),*) " Done."

  end subroutine tem_read_stlFiles
  ! *****************************************************************************

  ! *****************************************************************************
  !> This routine loads STL files from config and reads the triangles from the
  !! files into the dynamic array of triangles.
  subroutine tem_load_stl(stl_data, transform, conf, thandle)
    ! --------------------------------------------------------------------------!
    !> Array array of triangles in stlData
    type(tem_stlData_type), intent(out) :: stl_data
    !> transformation for spatial object
    type(tem_transformation_type), intent(in) :: transform
    !> Lua state
    type(flu_state) :: conf
    integer, intent(in) :: thandle !< handle for canonical objects
    ! --------------------------------------------------------------------------!
    integer :: iTri, iVer
    ! --------------------------------------------------------------------------!

    ! Load stl files from config file.
    call tem_load_stlhead(me = stl_data%head, conf = conf, thandle = thandle)

    ! Load triangles and nodes from stl files.
    call tem_read_stlFiles(stl_data = stl_data )

    ! if transformation is active apply transformation to all triangles
    ! first apply deformation and then translation
    if(transform%active) then
      if(transform%deform%active) then
        do iTri=1, stl_data%nTris
          do iVer=1,3
            stl_data%nodes(:, stl_data%tri_node(iVer, iTri)) = &
              &            matmul( transform%deform%matrix, &
              &            stl_data%nodes(:, stl_data%tri_node(iVer, iTri)) )
          enddo
        enddo
      endif
      if(transform%translate%active) then
        do iTri=1, stl_data%nTris
          do iVer=1,3
            stl_data%nodes(:, stl_data%tri_node(iVer, iTri)) = &
              &            transform%translate%vec + &
              &            stl_data%nodes(:, stl_data%tri_node(iVer, iTri))
          enddo
        enddo
      endif
    endif

  end subroutine tem_load_stl
  ! *****************************************************************************

  ! ****************************************************************************
  !> Compute, if the triangles in stl intersects the cube.
  function tem_stlCubeOverlap(stl_data, cube) result(overlaps)
    !--------------------------------------------------------------------------!
    type(tem_stlData_type), intent(in) :: stl_data
    type(tem_cube_type), intent(in) :: cube
    logical :: overlaps
    !--------------------------------------------------------------------------!
    integer :: iTri
    type(tem_triangle_type) :: triangle
    !--------------------------------------------------------------------------!

    overlaps = .false.
    do iTri=1, stl_data%nTris
      triangle%nodes(:,1) = stl_data%nodes( :, stl_data%tri_node(1,iTri) )
      triangle%nodes(:,2) = stl_data%nodes( :, stl_data%tri_node(2,iTri) )
      triangle%nodes(:,3) = stl_data%nodes( :, stl_data%tri_node(3,iTri) )
      overlaps = overlaps .or. tem_triangleCubeOverlap(triangle, cube)
      ! if intersection is found then terminate the loop and return overlaps
      if (overlaps) return
    end do

  end function tem_stlCubeOverlap
  ! ****************************************************************************

  ! ************************************************************************** !
  !> Write out an array of stlHeads in lua format
  !!
  subroutine tem_stlHead_out_vec( me, conf )
    ! --------------------------------------------------------------------------
    !> stlHead types to write out
    type( tem_stlHead_type ), intent(in) :: me(:)
    !> Aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    ! --------------------------------------------------------------------------
    ! counter
    integer :: i
    ! --------------------------------------------------------------------------

    ! create a table with name stlHead
    call aot_out_open_table( put_conf = conf, tname = 'object' )

    do i = 1, size(me)
      call tem_stlHead_out_scal( me(i), conf )
    end do

    call aot_out_close_table( put_conf = conf )

  end subroutine tem_stlHead_out_vec
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Write out a stlHead shape in lua format
  !!
  subroutine tem_stlHead_out_scal( me, conf )
    ! --------------------------------------------------------------------------
    !> stlHead types to write out
    type( tem_stlHead_type ), intent(in) :: me
    !> Aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    ! --------------------------------------------------------------------------
    ! ---------------------------------------------------------------------------

    ! create a table with name stlHead if not exist
    if( conf%level == 0 ) then
      call aot_out_open_table( put_conf = conf, tname = 'object' )
    else
      call aot_out_open_table( put_conf = conf )
    end if

    call aot_out_val( put_conf = conf, vname = 'filename', val = me%filename )

    call aot_out_close_table( put_conf = conf )

  end subroutine tem_stlHead_out_scal
  ! ************************************************************************** !


end module tem_stl_module

!> \page stl STL
!! stl files can be used as geometry kind. At first, seeder load the
!! triangles from stl files to the temporary stl_data type and
!! then each triangle in the stl_data is converted to tem_triangle_type.\n
!! stl geometry requires filename and stl file format. If file format
!! is not provided, default is set to binary.
!! Valid definition:
!! \li Single stl
!! \verbatim
!! geometry={
!!   kind='stl',
!!     object={
!!       filename='cube.stl',
!!       format = 'ascii' -- if not provided, default is binary
!!     }
!! }
!! \endverbatim
!!
!! \li Multiple stls
!! \verbatim
!! geometry={
!!   kind='stl',
!!     object={
!!       {
!!       filename = 'cube.stl'
!!       },
!!       {
!!       filename = 'cylinder.stl'
!!       }
!!     }
!! }
!! \endverbatim
!! \n\n
!! Seeder file to create mesh with single 'stl' geometry:
!! \include testsuite/stl/seeder.lua
!! \n\n
!! Mesh with 'stl' geometry created by seeder file:
!! \image html stl_cylinder.png
!! \endverbatim
!! Example lua file is available at \link testsuite/stl/seeder.lua
!! \example testsuite/stl/seeder.lua

