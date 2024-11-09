! Copyright (c) 2013-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013-2014, 2019, 2021 Harald Klimach <harald.klimach@dlr.de>
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
! **************************************************************************** !
!> author: Simon Zimny
!! This module provides functionality for reading and writing stl files.
!!
!!
module tem_stlb_io_module

  ! incude treelm modules
  use mpi
  use env_module,            only: single_k,rk, LabelLen, PathLen, rk_mpi, eps
  use tem_logging_module,    only: logUnit
  use tem_aux_module,        only: tem_open, tem_abort
  use tem_math_module,       only: cross_product3D
  use tem_time_module,       only: tem_time_type, tem_time_sim_stamp
  use tem_comm_env_module,   only: tem_comm_env_type
  use tem_grow_array_module, only: grw_int2darray_type, init, append, destroy

  implicit none

  private

  public :: tem_size_stlb
  public :: tem_read_stlb
  public :: tem_dump_stlb


contains


  ! ************************************************************************** !
  !> This subroutine reads the number of nodes and triangles from a given binary
  !! stl file.
  !!
  subroutine tem_size_stlb( filename, nNodes, nTris )
    ! --------------------------------------------------------------------------
    !> name of the binary stl file
    character(len=PathLen),intent(in) :: filename
    !> number of nodes in the file
    integer,intent(out) :: nNodes
    !> number of triangles in the file
    integer,intent(out) :: nTris
    ! --------------------------------------------------------------------------
    character(len=LabelLen) :: header
    integer :: stlUnit
    logical :: file_exists
    ! --------------------------------------------------------------------------

    inquire(file = trim(filename), exist = file_exists)
    if (.not. file_exists) then
      write(logUnit(0),*) trim(filename)//" file was not found. Dying..."
      call tem_abort()
    endif

    call tem_open( newunit = stlUnit,        &
      &            file    = trim(filename), &
      &            access  = 'stream',       &
      &            action  = 'read',         &
      &            status  = 'old'           )

    read(stlUnit) header
    read(stlUnit) nTris
    nNodes = nTris * 3

    close(stlUnit)

  end subroutine tem_size_stlb
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine actually reads the data (points, triangles, normals) from the
  !! binary file and stores them.
  !!
  subroutine tem_read_stlb( filename, nNodesRead, nTrisRead, nodes, tri_node,  &
    &                       ierror )
    ! --------------------------------------------------------------------------
    !> name of the binary stl file
    character(len=PathLen),intent(in) :: filename
    !> point coordinates read from the stl-file
    !! size: 3, nPoints_total
    real(kind=rk),intent(out) :: nodes(:,:)
    !> connectivity array for the triangles
    !! size: 3, nTriangles_total
    integer,intent(out) :: tri_node(:,:)
    !> Number of nodes read from the file header, to compare against the actual
    !! number of nodes read
    integer,intent(in)  :: nNodesRead
    !> Number of triangles read from the file header, to compare against the
    !! actual number of triangles read
    integer,intent(in)  :: nTrisRead
    !> error while openeing the file, or if the number of nodes/trias do not
    !! match to the ones read from the header (if error -> iError > 0)
    integer,intent(out) :: iError
    ! --------------------------------------------------------------------------
    ! buffer for reading the data from file
    real(kind=single_k)       :: temp(3) ! has to be single_k
    character(len=80) :: header ! has to be of length 80
    character(len=2) :: attribute
    integer :: nTriangles ! has to be 4 byte integer
    integer :: stlUnit
    integer :: i
    integer :: nNodes
    integer :: nTris
    ! --------------------------------------------------------------------------

    call tem_open( newunit = stlUnit,        &
      &            file    = trim(filename), &
      &            access  = 'stream',       &
      &            action  = 'read',         &
      &            status  = 'old'           )

    read(stlUnit) header
    read(stlUnit) nTriangles

    nNodes = 0
    nTris = 0
    do i=1,nTriangles
      ! assign node coordinates
      ! and nodes to tris
      nTris = nTris+1
      read(stlUnit) temp(:) ! normal vector

      read(stlUnit) temp(:) ! node 1 coords
      nNodes = nNodes+1
      nodes(1:3,nNodes) = real(temp(1:3),kind=rk)
      tri_node(1,nTris) = nNodes
      read(stlUnit) temp(:) ! node 2 coords
      nNodes = nNodes+1
      nodes(1:3,nNodes) = real(temp(1:3),kind=rk)
      tri_node(2,nTris) = nNodes
      read(stlUnit) temp(:) ! node 3 coords
      nNodes = nNodes+1
      nodes(1:3,nNodes) = real(temp(1:3),kind=rk)
      tri_node(3,nTris) = nNodes
      read(stlUnit) attribute
    end do

    close(stlUnit)
    if ( nTris .ne. nTrisRead) then
      write(logUnit(0),*) "Inconsistency found in number of triangles!"
      iError = 1
    endif
    if ( nNodes .ne. nNodesRead) then
      write(logUnit(0),*) "Inconsistency found in number of nodes!"
      iError = 1
    endif

  end subroutine tem_read_stlb
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine dumps a set of nodes and triangles to disc.
  !!
  !! The nodes and their connectivity are passed to the routine. The normals
  !! are passed optional or calculated internally. The outputfile name is
  !! composed of the $outprefix,$time,'.stl'.
  !!
  !! @note: This routine is a testing routine each process will dump all
  !!        of its nodes in a seperate file (this may include outdated data
  !!        as well)
  !!
  subroutine tem_dump_stlb( outprefix, nodes, triangles, proc, header,   &
    &                       normals, time )
    ! --------------------------------------------------------------------------
    !> output prefix for the filename
    character(len=*), intent(in) :: outprefix
    !> nodes to be dumped (size: 3*nNodes)
    real(kind=rk), intent(in) :: nodes(:)
    !> triangles to be dumped (size: 3, nTrias)
    integer, intent(in) :: triangles(:,:)
    !> process description to use
    type(tem_comm_env_type), intent(in) :: proc
    !> optional header to be dumped
    character(len=80), optional, intent(in) :: header
    !> optional array of normals, if not passed normals will be calculated
    !! internally
    real(kind=rk), optional, intent(in) :: normals(:,:)
    !> optional simulation time to be appended to the filename
    type(tem_time_type), optional, intent(in) :: time
    ! --------------------------------------------------------------------------
    real(kind=single_k) :: loc_normals( 3, size( triangles, 2 ))
    integer :: iTria
    ! temporary vectors to calculate the normals
    real(kind=rk) :: a(3), b(3)
    ! filename
    character(len=PathLen) :: filename
    ! output unit
    integer :: outUnit
    ! local header
    character(len=80) :: loc_header
    ! attribute to be dumped
    character(len=2) :: attribute
    ! error variable
    integer :: iError
    ! timestamp for the filename
    character(len=12) :: timeStamp
    ! temporary min and max position for X,Y,Z coordinates in the linearized
    ! array of nodes
    integer :: minPos1, maxPos1, minPos2, maxPos2, minPos3, maxPos3
    ! temporary array of node coordinates
    ! size = size(nodes)
    real(kind=rk), allocatable :: dump_nodes(:)
    ! number of entries in the nodes array
    integer :: nEntries
    real( kind=rk ) :: huge_real
    logical :: validTria
    ! temporary growing array of valid triangles to be dumped
    type(grw_int2darray_type) :: dump_trias
    integer :: iPoint
    ! --------------------------------------------------------------------------

    huge_real = huge(huge_real)

    nEntries = size( nodes )
    allocate( dump_nodes( nEntries ))

    ! first communicate all point coordinates to rank 0
    call mpi_reduce( nodes, dump_nodes, nEntries, rk_mpi, MPI_MIN, proc%root,  &
                     proc%comm, iError )

    ! now only root continues to calculate and dump information
    if( proc%rank .eq. proc%root )then

      if( present( normals ))then
        if( size( triangles, 2 ) .eq. size( normals, 2 ) )then
          loc_normals = real(normals, kind=single_k)
        else
          write(logUnit(0),*) " The number of triangles have to match the " &
            &              // "number of normals!!! This is not the case "  &
            &              // "(nTrias: ", size( triangles, 2 ),            &
            &                 " vs. nNormals: ", size( normals, 2 ),        &
            &                 " Stopping..."
          call tem_abort()
        end if
      else
        ! calculate the normals
        do iTria = 1, size( triangles, 2 )
          minPos1 = (triangles( 1, iTria )-1)*3+1
          maxPos1 = (triangles( 1, iTria )-1)*3+3
          minPos2 = (triangles( 2, iTria )-1)*3+1
          maxPos2 = (triangles( 2, iTria )-1)*3+3
          minPos3 = (triangles( 3, iTria )-1)*3+1
          maxPos3 = (triangles( 3, iTria )-1)*3+3
          a = dump_nodes( minPos2:maxPos2 ) - dump_nodes( minPos1:maxPos1 )
          b = dump_nodes( minPos3:maxPos3 ) - dump_nodes( minPos1:maxPos1 )
          loc_normals( :, iTria ) = real(cross_product3D( a, b ), kind=single_k)
        end do
      end if

      ! initialize the growing array of valid triangles to be dumped
      call init( me = dump_trias, width = 3 )

      ! now check the individual coordinates and remove all triangles that
      ! own a point coordinate which is not .lt. huge(rk)
      do iTria = 1, size( triangles, 2 )
        ! suppose the triangle is valid
        validTria = .true.
        point_loop: do iPoint = 1, 3
          minPos1 = (triangles( iPoint, iTria )-1)*3+1
          maxPos1 = (triangles( iPoint, iTria )-1)*3+3
          ! check if any coordinate for iPoint is not valid
          if( .not. all( dump_nodes( minPos1:maxPos1 ) .lt. huge_real ))then
            ! if this is the case set the validTria to false and
            ! exit the point loop
            validTria = .false.
            exit point_loop
          end if
        end do point_loop
        if( validTria )then
          ! append the valid triangle
          call append( me  = dump_trias,           &
            &          val = triangles( :, iTria ))
        end if
      end do

      if( present( header ))then
        loc_header = header
      else
        loc_header = ''
      end if

      ! initialize the attribute to be empty
      attribute = ''

      if( present( time ))then
        timestamp = adjustl(tem_time_sim_stamp(time))
        ! assemble the filename
        write(filename,'(a)')trim(outprefix)//'_t'//trim(timestamp)//".stl"
      else
        write(filename,'(a)')trim(outprefix)//".stl"
      end if

      ! open the output file and
      call tem_open( newunit = outUnit,        &
        &            file    = trim(filename), &
        &            access  = 'stream',       &
        &            action  = 'write',        &
        &            status  = 'replace'       )

      ! ... dump the header
      write(outUnit) loc_header
      ! ... dump the number of triangles
      write(outUnit) dump_trias%nVals

      ! .. for every triangle dump
      do iTria = 1, dump_trias%nVals
        ! ... calculate the correct positions in the linearized array
        minPos1 = (dump_trias%val( 1, iTria )-1)*3+1
        maxPos1 = (dump_trias%val( 1, iTria )-1)*3+3
        minPos2 = (dump_trias%val( 2, iTria )-1)*3+1
        maxPos2 = (dump_trias%val( 2, iTria )-1)*3+3
        minPos3 = (dump_trias%val( 3, iTria )-1)*3+1
        maxPos3 = (dump_trias%val( 3, iTria )-1)*3+3
        ! ... the normal vector
        write(outUnit)loc_normals( :, iTria )
        ! ... the vertices
        write(outUnit)real( dump_nodes( minPos1:maxPos1 ) ,kind=single_k )
        write(outUnit)real( dump_nodes( minPos2:maxPos2 ) ,kind=single_k )
        write(outUnit)real( dump_nodes( minPos3:maxPos3 ) ,kind=single_k )
        write(outUnit)attribute
      end do

      close(outunit)
      ! destroy the growing 2D array of valid triangles
      call destroy( me = dump_trias )

    end if

    deallocate( dump_nodes )

  end subroutine tem_dump_stlb
  ! ************************************************************************** !


end module tem_stlb_io_module
! **************************************************************************** !
