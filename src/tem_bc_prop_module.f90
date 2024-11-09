! Copyright (c) 2011-2014, 2016-2017, 2019, 2023 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2012, 2014-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011-2012 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2011-2012, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2015 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2017 Raphael Haupt <Raphael.Haupt@student.uni-siegen.de>
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
!> author: Harald Klimach
!! author: Kartik Jain
!! author: Jiaxing Qi
!! A module to provide informations on boundary condition [[tem_bc_module]]
!! property for elements.
!!
module tem_bc_prop_module

  ! include treelm modules
  use mpi
  use env_module,       only: long_k, rk, tem_create_EndianSuffix, &
    &                         labelLen, long_k_mpi, io_buffer_size, rk_mpi
  use tem_param_module, only: q__B, q__T, q__N, q__S, q__W, q__E, &
    &                         q_BS, q_BN, q_BW, q_BE,             &
    &                         q_TS, q_TN, q_TW, q_TE,             &
    &                         q_SW, q_NW, q_SE, q_NE,             &
    &                         qBSW, qBSE, qBNW, qBNE, qTSW, qTSE, qTNW, qTNE
  use treelmesh_module, only: treelmesh_type
  use tem_topology_module, only: tem_coordOfID, tem_IdOfCoord, &
    &                            tem_firstIDAtLevel
  use tem_prophead_module, only: tem_prophead_type
  use tem_property_module, only: tem_property_type, prp_hasBnd
  use tem_logging_module,  only: logUnit
  use tem_debug_module,    only: tem_print_array, dbgUnit
  use tem_aux_module,      only: tem_open, check_mpi_error

  ! include aotus modules
  use aot_out_module,   only: aot_out_type, aot_out_val, aot_out_open,         &
    &                         aot_out_close, aot_out_open_table,               &
    &                         aot_out_close_table
  use aot_table_module, only: aot_table_open, aot_table_close
  use aotus_module,     only: open_config_file, aot_get_val, close_config,     &
    &                         flu_State

  implicit none

  private

  public :: init_tem_bc_prop, load_tem_bc_prop
  public :: tem_bc_prop_sublist
  public :: tem_BC_prop_type
  public :: tem_empty_BC_prop, dump_tem_BC_propHeader
  public :: tem_unload_bc_prop
  public :: load_tem_BC_qVal
  public :: dump_tem_BC_prop
  public :: dump_tem_BC_qVal
  public :: dump_tem_BC_normal
  public :: load_tem_BC_normal
  public :: tem_debug_bc_prop
  public :: tem_bc_prop_pos

  type tem_BC_prop_type
    !> Pointer to treelmesh_type%global%property
    type(tem_prophead_type),  pointer :: header => null()

    !> Number of sides in the elements, each side might be assigned a
    !! boundary condition ID
    !! The first 6 are dedicated to "direct" neighbors, and shall always be
    !! present.
    !! The next 12 are dedicated to neighbors where two of the indices in the
    !! 3 dimensional address have to be changed (neighbors at the edges of
    !! the cube)
    !! After these 18 entries, additional 8 entries might be used for neighbors
    !! at the corners.
    !!
    !! By default we should store all 26 neighbors, to provide an as complete
    !! neighborhood as possible to arbitrary solvers.
    integer :: nSides

    !> Number of different Boundary conditions used in the complete domain
    !! e.g. wall, inflow, outflow
    integer :: nBCtypes

    !> Array of labels identifying each of the boundary conditions.
    !! This array has a length of nBCtypes
    character(len=LabelLen), allocatable :: BC_label(:)

    !> Logical array indicating whether each boundary is high order wall
    !! This array has a length of nBCtypes
    logical, allocatable :: hasQVal(:)

    !> Logical array indicating whether each boundary provides normal wall
    !! information.
    logical, allocatable :: hasNormal(:)

    !> Actual q-value array for high order wall boundary conditions
    !! The fist dimension has a length of nSides
    !! The second dimension has a length of the number of elements with the
    !! boundary condition property. tem_property_type%nElems
    real(kind=rk), allocatable :: qVal(:,:)

    !> Actual normals array for boundary conditions providing the wall normals
    !! The fist dimension has length 3
    !! The second dimension has a length of the number of elements with the
    !! hasNormal property. tem_property_type%nElems
    real(kind=rk), allocatable :: normal(:,:)

    !> Actual boundary condition identification for all sides of elements with
    !! this property.
    !! The first dimension has a length of nSides
    !! The second dimension has a length of the number of elements with the
    !! boundary condition property. tem_property_type%nElems
    !! An ID of 0 for a side indicates, that there is no boundary in that
    !! direction.
    !! A negative ID indicates a periodic boundary and that the given absolute
    !! of the given value is to be taken as neighboring ID in that direction.
    integer(kind=long_k), allocatable :: boundary_ID(:,:)

    !> Pointer to treelmesh_type%property
    type(tem_property_type),  pointer :: property => null()

  end type tem_BC_prop_type


contains


  ! ------------------------------------------------------------------------ !
  !> Get the boundary property position in the list of all properties in tree.
  !!
  !! There should be only one boundary property, and the first one found will
  !! be returned. If there is no boundary property at all a `-1` is returned.
  pure function tem_bc_prop_pos(tree) result(bcpos)
    ! -------------------------------------------------------------------- !
    !> Tree to find the boundary condition property for.
    type(treelmesh_type), intent(in) :: tree

    !> Position of the boundary condition property in the list of all
    !! properties.
    integer :: bcpos
    ! -------------------------------------------------------------------- !
    integer :: iProp
    ! -------------------------------------------------------------------- !

    bcpos = -1
    do iProp=1,tree%global%nProperties
      if (tree%global%Property(iprop)%bitpos == prp_hasBnd) then
        bcpos = iProp
        EXIT
      end if
    end do

  end function tem_bc_prop_pos
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


! ****************************************************************************** !
  !> Initialize boundary conditions of a given tree.
  !!
  subroutine init_tem_bc_prop( tree, mypart, comm, bc )
    ! ---------------------------------------------------------------------------
    type(treelmesh_type), intent(in)    :: tree
    integer, intent(in)                 :: mypart
    !> Communicator to use
    integer, intent(in) :: comm
    type(tem_bc_prop_type), intent(out) :: bc
    ! ---------------------------------------------------------------------------
    integer :: iprop
    logical :: found_bc
    ! ---------------------------------------------------------------------------

    iProp = tem_bc_prop_pos(tree)
    found_bc = iProp > 0

    if (found_bc) then
      bc%header => tree%global%Property(iprop)
      bc%property => tree%property(iprop)

      select case(trim(bc%header%label))
      case('internal 0D BC')
        call load_bc_intern_0D(tree = tree, me = bc)
      case('internal 1D BC')
        call load_bc_intern_1D(tree = tree, me = bc, xbounds=.false.)
      case('bounded internal 1D BC')
        call load_bc_intern_1D(tree = tree, me = bc, xbounds=.true.)
      case('internal 2D BC')
        call load_bc_intern_2D(tree = tree, me = bc)
      case default
        call load_tem_bc_prop(me       = bc,                               &
          &                   offset   = bc%property%Offset,               &
          &                   nElems   = bc%property%nElems,               &
          &                   basename = trim(tree%global%dirname)//'bnd', &
          &                   comm     = comm,                             &
          &                   mypart   = mypart                            )
      end select
    else
      call tem_empty_BC_prop(bc)
    end if

  end subroutine init_tem_bc_prop
! ****************************************************************************** !


! ****************************************************************************** !
  !> load bc property header from lua file, boundaryID from bnd.lsb
  !!
  subroutine load_tem_BC_prop( me, offset, nElems, basename, myPart, comm )
    ! ---------------------------------------------------------------------------
    !> Boundary condition construct to load the data into
    type(tem_BC_prop_type), intent(inout) :: me
    !> Offset of the local set of elements in the global list
    integer(kind=long_k), intent(in) :: offset
    !> Local number of elements with this property
    integer, intent(in) :: nElems
    !> Name of the file, the data is stored in, will be appended with
    !! ".ascii" for the header information and ".lsb" or ".msb" for the
    !! binary data.
    character(len=*), intent(in) :: basename
    !> Partition to load
    integer, intent(in) :: myPart
    !> Communicator to use
    integer, intent(in) :: comm
    ! ---------------------------------------------------------------------------
    type( flu_State ) :: conf ! lua flu state to read lua file
    integer :: i
    integer, parameter :: root = 0
    integer :: BCcomm
    integer :: color, iError
    logical :: participant !< If the local rank is a participant in BC
    character(len=4) :: EndianSuffix
    character(len=256) :: headerfile
    character(len=256) :: datafile
    integer :: thandle, typesize
    integer(kind=long_k), allocatable :: buffer(:)
    integer(kind=long_k), allocatable :: globbuffer(:)
    integer(kind=MPI_OFFSET_KIND)     :: displacement
    integer :: fh, etype, ftype, iostatus( MPI_STATUS_SIZE )
    integer :: file_rec_len
    integer :: myBCrank

    ! ---------------------------------------------------------------------------

    headerfile = trim(basename)//'.lua'
    EndianSuffix = tem_create_EndianSuffix()
    datafile = trim(basename)//trim(EndianSuffix)

    if (me%header%nElems > 0) then
      write(logUnit(1), *) 'Load boundary ID from file: '//trim(datafile)
    end if

    if (myPart == root) then
      ! Read the header only on the root process, broadcast to all others
      ! open mesh header file
      call open_config_file( L = conf, filename = headerfile )
      call aot_get_val( L       = conf,      &
        &               key     = 'nSides',  &
        &               val     = me%nSides, &
        &               ErrCode = iError     )
      call aot_get_val( L       = conf,        &
        &               key     = 'nBCtypes',  &
        &               val     = me%nBCtypes, &
        &               ErrCode = iError       )
    end if

    call MPI_Bcast(me%nSides, 1, MPI_INTEGER, root, comm, iError)
    call MPI_Bcast(me%nBCtypes, 1, MPI_INTEGER, root, comm, iError)

    allocate(me%BC_label(me%nBCtypes))
    allocate(me%hasQVal(me%nBCtypes))
    me%hasQVal(:) = .false.

    if (myPart == root) then
      ! Now read the list of boundary labels
      call aot_table_open( L = conf, thandle = thandle, key = 'bclabel' )
      do i=1,me%nBCtypes
        call aot_get_val( L       = conf,           &
          &               thandle = thandle,        &
          &               pos     = i,              &
          &               val     = me%BC_label(i), &
          &               ErrCode = iError          )
      end do
      call aot_table_close( L = conf, thandle = thandle )
      call close_config( conf )
    end if

    call MPI_Bcast( me%BC_label, LabelLen*me%nBCtypes, MPI_CHARACTER, &
      &             root, comm, iError                                )

    allocate(me%boundary_ID(me%nSides, nElems))

    participant = ( nElems > 0 )

    if (participant) then
      color = 1
    else
      color = MPI_UNDEFINED
    end if

    ! Split the communicator
    call MPI_COMM_SPLIT(comm, color, myPart, BCcomm, iError)

    if (nElems > 0) then

      allocate( buffer(me%nSides * nElems) )

      if (me%header%nElems*me%nSides > io_buffer_size) then
        write(logUnit(5), *) 'read with MPI'


        ! Create a contiguous type to describe the vector per element
        call MPI_TYPE_CONTIGUOUS( me%nSides, long_k_mpi, etype, iError )
        call check_mpi_error(iError,'contiguous etype in load_tem_BC_prop')
        call MPI_TYPE_COMMIT(     etype, iError )
        call check_mpi_error(iError,'commit etype in load_tem_BC_prop')
        call MPI_TYPE_SIZE(etype, typesize, iError )
        call check_mpi_error(iError,'typesize in load_tem_BC_prop')

        ! Calculate displacement for file view
        displacement = offset * typesize * 1_MPI_OFFSET_KIND

        ! Create a MPI CONTIGUOUS as ftype for file view
        call MPI_TYPE_CONTIGUOUS(nElems, etype, ftype, iError)
        call check_mpi_error(iError,'contiguous ftype in load_tem_BC_prop')
        call MPI_TYPE_COMMIT( ftype, iError )
        call check_mpi_error( iError, 'commit ftype in load_tem_BC_prop')

        ! Open the binary file for MPI I/O (read)
        call MPI_FILE_OPEN( BCcomm, trim(datafile), MPI_MODE_RDONLY,   &
          &                 MPI_INFO_NULL, fh, iError                  )
        call check_mpi_error( iError, 'Open File in load_tem_BC_prop')

        ! Set the view for each process on the file above
        call MPI_FILE_SET_VIEW( fh, displacement, etype, ftype, "native",  &
          &                     MPI_INFO_NULL, iError                      )
        call check_mpi_error( iError,'Set File view in load_tem_BC_prop')

        ! Read data from the file
        call MPI_FILE_READ_ALL( fh, buffer, nElems, etype, iostatus, iError )
        call check_mpi_error( iError,'Read All in load_tem_BC_prop')

        !Free the MPI_Datatypes which were created and close the file
        call MPI_TYPE_FREE (etype, iError)
        call check_mpi_error( iError,'free etype in load_tem_BC_prop')
        call MPI_TYPE_FREE (ftype, iError)
        call check_mpi_error( iError,'free ftype in load_tem_BC_prop')
        call MPI_FILE_CLOSE(fh,    iError)
        call check_mpi_error( iError,'close file in load_tem_BC_prop')
        ! END IO-part

      else

        ! File is so small, it probably is faster to read it on a single process
        ! and broadcast the data.
        call MPI_Comm_rank(BCcomm, myBCrank, iError)


        allocate(globbuffer(me%header%nElems*me%nSides))

        if (myBCrank == 0) then
          write(logUnit(5), *) 'read on a single process'
          inquire(iolength = file_rec_len) globbuffer
          call tem_open( newunit = fh,             &
            &            file    = trim(datafile), &
            &            recl    = file_rec_len,   &
            &            action  = 'read',         &
            &            access  = 'direct',       &
            &            form    = 'unformatted'   )

          read(fh, rec=1) globbuffer

          close(fh)
         end if

        call MPI_Bcast( globbuffer, int(me%nSides*me%header%nElems), &
          &             long_k_mpi, 0, BCcomm, iError                )

        buffer = globbuffer(offset*me%nSides+1:(offset+nElems)*me%nSides)
        deallocate(globbuffer)

      end if

      do i=1,nElems
        me%boundary_ID(:,i) = buffer( ((i-1)*me%nSides+1) : (i*me%nSides) )
      end do

      deallocate(buffer)

   end if

  end subroutine load_tem_BC_prop
! ****************************************************************************** !


! ****************************************************************************** !
  !> Load internal BC property for 1D Line.
  !!
  subroutine load_BC_intern_1D( tree, me, xbounds, nSides )
    ! ---------------------------------------------------------------------------
    type(treelmesh_type), intent(in) :: tree
    !> Boundary condition construct to load the data into
    type(tem_BC_prop_type), intent(inout) :: me
    !> Set boundaries east and west in X direction?
    logical, intent(in) :: xbounds
    !> Required sides to set, defaults to 26.
    integer, intent(in), optional :: nSides
    ! ---------------------------------------------------------------------------
    integer :: nelems
    integer :: iElem
    integer :: neigh_coord(4)
    integer :: my_coord(4)
    integer(kind=long_k) :: tID, west_ID, east_ID
    integer(kind=long_k) :: firstID, lastID
    ! ---------------------------------------------------------------------------

    if (present(nSides)) then
      me%nSides = nSides
    else
      me%nSides = 26
    end if

    if (xbounds) then
      me%nBCtypes = 2
    else
      me%nBCtypes = 0
    end if

    firstID = -1_long_k
    lastID = -1_long_k

    nElems = int(tree%global%nElems)
    if (nElems < 2**tree%global%minLevel .or. xbounds) then
      ! If the line does not span the complete cube axis, set the first and
      ! last ID here. This will be used for a proper periodicity at the
      ! truncated end. Without truncation, these two IDs are negative and no
      ! element will match the test below, thus resulting in usual full cube
      ! behavior.
      firstID = tem_firstIdAtLevel(tree%global%minlevel)
      lastID = tem_IdOfCoord(                               &
        & coord  = [ nelems-1, 0, 0, tree%global%minlevel], &
        & offset = firstID                                  )
    end if

    allocate(me%BC_label(me%nBCtypes))
    if (xbounds) then
      me%BC_label(1) = 'west'
      me%BC_label(2) = 'east'
    end if

    allocate(me%boundary_ID(me%nSides, me%Property%nElems))
    me%boundary_ID = 0_long_k
    do iElem=1,me%Property%nElems
      tID = tree%treeID(me%Property%ElemID(iElem))
      my_coord = tem_coordofID(tID)

      ! Neighbor in the west
      if (tID /= firstID) then
        neigh_coord = my_coord + [-1, 0, 0, 0]
        west_ID = -tem_IdOfCoord(neigh_coord)
      else
        if (xbounds) then
          west_ID = 1
        else
          west_ID = -lastID
        end if
        me%boundary_ID(q__W, iElem) = west_ID
      end if

      ! Neighbor in the east
      if (tID /= lastID) then
        neigh_coord = my_coord + [1, 0, 0, 0]
        east_ID = -tem_IdOfCoord(neigh_coord)
      else
        if (xbounds) then
          east_ID = 2
        else
          east_ID = -firstID
        end if
        me%boundary_ID(q__E, iElem) = east_ID
      end if

      me%boundary_ID(q__B, iElem) = -tID
      me%boundary_ID(q__T, iElem) = -tID
      me%boundary_ID(q__S, iElem) = -tID
      me%boundary_ID(q__N, iElem) = -tID
      if (me%nSides > 6) then
        me%boundary_ID(q_BS, iElem) = -tID
        me%boundary_ID(q_TS, iElem) = -tID
        me%boundary_ID(q_BN, iElem) = -tID
        me%boundary_ID(q_TN, iElem) = -tID
        ! Neighbors in the west
        me%boundary_ID(q_BW, iElem) = west_ID
        me%boundary_ID(q_TW, iElem) = west_ID
        me%boundary_ID(q_SW, iElem) = west_ID
        me%boundary_ID(q_NW, iElem) = west_ID
        ! Neighbors in the east
        me%boundary_ID(q_BE, iElem) = east_ID
        me%boundary_ID(q_TE, iElem) = east_ID
        me%boundary_ID(q_SE, iElem) = east_ID
        me%boundary_ID(q_NE, iElem) = east_ID
        if (me%nSides > 18) then
          ! Neighbors in the west
          me%boundary_ID(qBSW, iElem) = west_ID
          me%boundary_ID(qTSW, iElem) = west_ID
          me%boundary_ID(qBNW, iElem) = west_ID
          me%boundary_ID(qTNW, iElem) = west_ID
          ! Neighbors in the east
          me%boundary_ID(qBSE, iElem) = east_ID
          me%boundary_ID(qTSE, iElem) = east_ID
          me%boundary_ID(qBNE, iElem) = east_ID
          me%boundary_ID(qTNE, iElem) = east_ID
        end if
      end if
    end do

  end subroutine load_BC_intern_1D
! ****************************************************************************** !


! ****************************************************************************** !
  !> Load internal BC property for a single fully periodic element.
  !!
  subroutine load_BC_intern_0D( tree, me, nSides )
    ! ---------------------------------------------------------------------------
    type(treelmesh_type), intent(in) :: tree
    !> Boundary condition construct to load the data into
    type(tem_BC_prop_type), intent(inout) :: me
    !> Required sides to set, defaults to 26.
    integer, intent(in), optional :: nSides
    ! ---------------------------------------------------------------------------

    if (present(nSides)) then
      me%nSides = nSides
    else
      me%nSides = 26
    end if

    me%nBCtypes = 0

    allocate(me%BC_label(me%nBCtypes))
    allocate(me%boundary_ID(me%nSides, me%Property%nElems))
    ! All links point back to the element itself.
    me%boundary_ID = -tree%treeID(1)

  end subroutine load_BC_intern_0D
! ****************************************************************************** !


! ****************************************************************************** !
  !> Load internal BC property for 2D Slice.
  !!
  subroutine load_BC_intern_2D( tree, me, nSides )
    ! ---------------------------------------------------------------------------
    type(treelmesh_type), intent(in) :: tree
    !> Boundary condition construct to load the data into
    type(tem_BC_prop_type), intent(inout) :: me
    !> Required sides to set, defaults to 26.
    integer, intent(in), optional :: nSides
    ! ---------------------------------------------------------------------------
    integer :: iElem
    integer :: neigh_coord(4)
    integer :: my_coord(4)
    integer(kind=long_k) :: tID, neighID
    ! ---------------------------------------------------------------------------

    if (present(nSides)) then
      me%nSides = nSides
    else
      me%nSides = 26
    end if

    me%nBCtypes = 0

    allocate(me%BC_label(me%nBCtypes))
    allocate(me%boundary_ID(me%nSides, me%Property%nElems))
    me%boundary_ID = 0_long_k
    do iElem=1,me%Property%nElems
      tID = tree%treeID(me%Property%ElemID(iElem))
      me%boundary_ID(q__B, iElem) = -tID
      me%boundary_ID(q__T, iElem) = -tID
      if (me%nSides > 6) then
        my_coord = tem_coordofID(tID)
        ! Neighbor in the west
        neigh_coord = my_coord + [-1, 0, 0, 0]
        neighID = tem_IdOfCoord(neigh_coord)
        me%boundary_ID(q_BW, iElem) = -neighID
        me%boundary_ID(q_TW, iElem) = -neighID
        ! Neighbor in the east
        neigh_coord = my_coord + [+1, 0, 0, 0]
        neighID = tem_IdOfCoord(neigh_coord)
        me%boundary_ID(q_BE, iElem) = -neighID
        me%boundary_ID(q_TE, iElem) = -neighID
        ! Neighbor in the south
        neigh_coord = my_coord + [0, -1, 0, 0]
        neighID = tem_IdOfCoord(neigh_coord)
        me%boundary_ID(q_BS, iElem) = -neighID
        me%boundary_ID(q_TS, iElem) = -neighID
        ! Neighbor in the north
        neigh_coord = my_coord + [0, +1, 0, 0]
        neighID = tem_IdOfCoord(neigh_coord)
        me%boundary_ID(q_BN, iElem) = -neighID
        me%boundary_ID(q_TN, iElem) = -neighID
        if (me%nSides > 18) then
          ! Neighbor in the south-west
          neigh_coord = my_coord + [-1, -1, 0, 0]
          neighID = tem_IdOfCoord(neigh_coord)
          me%boundary_ID(qBSW, iElem) = -neighID
          me%boundary_ID(qTSW, iElem) = -neighID
          ! Neighbor in the south-east
          neigh_coord = my_coord + [+1, -1, 0, 0]
          neighID = tem_IdOfCoord(neigh_coord)
          me%boundary_ID(qBSE, iElem) = -neighID
          me%boundary_ID(qTSE, iElem) = -neighID
          ! Neighbor in the north-west
          neigh_coord = my_coord + [-1, +1, 0, 0]
          neighID = tem_IdOfCoord(neigh_coord)
          me%boundary_ID(qBNW, iElem) = -neighID
          me%boundary_ID(qTNW, iElem) = -neighID
          ! Neighbor in the north-east
          neigh_coord = my_coord + [+1, +1, 0, 0]
          neighID = tem_IdOfCoord(neigh_coord)
          me%boundary_ID(qBNE, iElem) = -neighID
          me%boundary_ID(qTNE, iElem) = -neighID
        end if
      end if
    end do

  end subroutine load_BC_intern_2D
! ****************************************************************************** !


! ****************************************************************************** !
  !> Initialize an empty boundary condition
  subroutine tem_empty_BC_prop( me )
    ! ---------------------------------------------------------------------------
    !> Boundary condition construct to load the data into
    type(tem_BC_prop_type), intent(out) :: me
    ! ---------------------------------------------------------------------------

    me%nSides = 0
    me%nBCtypes = 0
    allocate(me%BC_label(0))
    allocate(me%boundary_ID(0,0))
    allocate(me%hasQVal(0))
    allocate(me%qVal(0,0))

    write(logUnit(5), *) 'Cleaned boundary property type.'

  end subroutine tem_empty_BC_prop
! ****************************************************************************** !

  ! ------------------------------------------------------------------------ !
  !> Create the boundary property for a restricted set of elements given by
  !! sublist (position of elements in tree, usually from a subtree).
  !!
  !! After creating a new tree with [[tem_create_tree_from_sub]], this
  !! routine can be used to create the according boundary information on
  !! the restricted set of elements.
  subroutine tem_bc_prop_sublist(tree, bc, header, property, sublist, sub_bc)
    ! -------------------------------------------------------------------- !
    !> The original tree from which the subset is to be selected.
    type(treelmesh_type), intent(in) :: tree

    !> The boundary condition property in the original mesh (tree).
    type(tem_BC_prop_type), intent(in) :: bc

    !> Header description of the boundary condition property in the sublist.
    !!
    !! This information has to be gathered for the elements of the sublist
    !! beforehand.
    type(tem_prophead_type), target, intent(in) :: header

    !> Property description of the boundary condition property in the sublist.
    !!
    !! This information has to be gathered for the elements of the sublist
    !! beforehand.
    type(tem_property_type), target, intent(in) :: property

    !> List of elements to get the boundary information for.
    integer, intent(in) :: sublist(:)

    !> New boundary property description for just the elements provided in
    !! sublist.
    !!
    !! This may be used to correctly describe the boundary conditions in
    !! a subtree for example.
    type(tem_BC_prop_type), intent(out) :: sub_bc
    ! -------------------------------------------------------------------- !
    integer :: iElem
    integer :: iBCElem
    integer :: iSubBCElem
    integer :: iSubElem
    integer :: nSubElems
    ! -------------------------------------------------------------------- !

    if (associated(bc%property)) then

      sub_bc%nSides   = bc%nSides
      sub_bc%nBCtypes = bc%nBCtypes
      if (allocated(bc%BC_label)) then
        allocate(sub_bc%BC_label(size(bc%BC_label)))
        sub_bc%BC_label = bc%BC_label
      end if
      sub_bc%header   => header
      sub_bc%property => property

      allocate(sub_bc%boundary_ID(bc%nSides, property%nElems))

      nSubElems = size(sublist)
      iSubElem = 1
      iBCElem = 0
      iSubBCElem = 0
      do iElem=1,tree%nElems
        if (iSubElem > nSubElems) EXIT
        if ( btest(tree%elemPropertyBits(iElem), prp_hasBnd) ) then
          iBCElem = iBCElem + 1
          if (sublist(iSubElem) == iElem) then
            iSubBCelem = iSubBCelem + 1
            sub_bc%boundary_id(:,iSubBCelem) = bc%boundary_ID(:,iBCElem)
          end if
        end if
        if (sublist(iSubElem) == iElem) iSubElem = iSubElem+1
      end do

    else

      call tem_empty_BC_prop(sub_bc)

    end if

  end subroutine tem_bc_prop_sublist
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


! ****************************************************************************** !
  !> dump bc properties header information to lua file
  !!
  subroutine dump_tem_BC_propHeader( me, headerfile )
    ! ---------------------------------------------------------------------------
    !> Boundary condition construct to load the data into
    type(tem_BC_prop_type), intent(in) :: me
    !> name of the bc header lua file
    character(len=*), intent(in) :: headerfile
    ! ---------------------------------------------------------------------------
    integer :: i
    type(aot_out_type) :: conf ! aotus lua state to write output
    ! ---------------------------------------------------------------------------

    call aot_out_open( conf, headerfile )
    call aot_out_val( put_conf = conf, vname = 'nSides', val = me%nSides )
    call aot_out_val( put_conf = conf, vname = 'nBCtypes', val = me%nBCtypes )
    call aot_out_open_table( conf, 'bclabel' )
    do i = 1, me%nBCtypes
      call aot_out_val( conf, val = me%BC_label(i) )
    end do
    call aot_out_close_table(conf)
    ! close the mesh header file
    call aot_out_close(conf)

  end subroutine dump_tem_BC_propHeader
! ****************************************************************************** !


! ****************************************************************************** !
  !> dump bc properties
  subroutine dump_tem_BC_prop( me, offset, nElems, basename, myPart, comm )
    ! ---------------------------------------------------------------------------
    !> Boundary condition construct to load the data into
    type(tem_BC_prop_type), intent(in) :: me
    !> Offset of the local set of elements in the global list
    integer(kind=long_k), intent(in) :: offset
    !> Local number of elements with this property
    integer, intent(in) :: nElems
    !> Name of the file, the data is stored in, will be appended with
    !! ".lua" for the header information and ".lsb" or ".msb" for the
    !! binary data.
    character(len=*), intent(in) :: basename
    !> Partition to dump
    integer, intent(in) :: myPart
    !> Communicator to use
    integer, intent(in) :: comm
    ! ---------------------------------------------------------------------------
    integer :: root
    integer :: locomm
    character(len=256) :: headerfile
    character(len=256) :: datafile
    character(len=4) :: EndianSuffix
    ! ---------------------------------------------------------------------------
    integer(kind=MPI_OFFSET_KIND)     :: displacement
    integer :: fh, etype, ftype, ierror, iostatus( MPI_STATUS_SIZE ), typesize
    ! ---------------------------------------------------------------------------

    root = 0

    locomm = comm

    headerfile = trim(basename)//'.lua'
    EndianSuffix = tem_create_EndianSuffix()
    datafile = trim(basename)//trim(EndianSuffix)

    if (myPart == root) then
      ! Only root partition needs to write the header
      !open up the mesh header lua file to dump the stuff using aotus library
      call dump_tem_BC_propHeader( me, headerfile )
    end if

    if (nElems > 0) then
      !> Mpi IO
      write(logUnit(1),*) 'Write boundary ID to file: ' // trim(datafile)

      ! Open the binary file for MPI I/O (Write)
      call MPI_FILE_OPEN( comm, trim(datafile),                  &
        &                 ior(MPI_MODE_WRONLY,MPI_MODE_CREATE),  &
        &                 MPI_INFO_NULL, fh, iError              )
      call check_mpi_error( iError,'Open File in dump_tem_BC_prop')

      ! Create a contiguous type to describe the vector per element
      call MPI_TYPE_CONTIGUOUS( me%nSides, long_k_mpi, etype, iError )
      call check_mpi_error(iError,'contiguous etype in dump_tem_BC_prop')
      call MPI_TYPE_COMMIT( etype, iError )
      call check_mpi_error( iError,'commit etype in dump_tem_BC_prop')
      call MPI_TYPE_SIZE(etype, typesize, iError )
      call check_mpi_error(iError,'typesize in dump_tem_BC_prop')

      ! Calculate displacement for file view
      displacement = offset * typesize * 1_MPI_OFFSET_KIND

      ! Create a MPI CONTIGUOUS as ftype for file view
      call MPI_TYPE_CONTIGUOUS(nElems, etype, ftype, iError)
      call check_mpi_error(iError,'contiguous ftype in dump_tem_BC_prop')
      call MPI_TYPE_COMMIT( ftype, iError )
      call check_mpi_error( iError,'commit ftype in dump_tem_BC_prop')

      ! Set the view for each process on the file above
      call MPI_FILE_SET_VIEW( fh, displacement, etype, ftype, "native",  &
        &                     MPI_INFO_NULL, iError )
      call check_mpi_error( iError,'Set File view in dump_tem_BC_prop')

      ! Read data from the file
      call MPI_FILE_WRITE_ALL( fh, me%boundary_ID, nElems, etype, iostatus, iError )
      call check_mpi_error( iError,'File write all in dump_tem_BC_prop')

      !Free the MPI_Datatypes which were created and close the file
      call MPI_TYPE_FREE (etype, iError)
      call check_mpi_error(iError,'free etype in dump_tem_BC_prop')
      call MPI_TYPE_FREE (ftype, iError)
      call check_mpi_error(iError,'free ftype in dump_tem_BC_prop')
      call MPI_FILE_CLOSE(fh,    iError)
      call check_mpi_error(iError,'close file in dump_tem_BC_prop')
      ! END IO-part
    end if

  end subroutine dump_tem_BC_prop
! ****************************************************************************** !


! ****************************************************************************** !
  subroutine dump_tem_BC_logicalHeader( headerfile, propname, nBCtypes, &
    &                                   flag_per_BC                     )
    ! ---------------------------------------------------------------------------
    character(len=*), intent(in) :: headerfile
    character(len=*), intent(in) :: propname
    integer, intent(in) :: nBCtypes
    logical, intent(in) :: flag_per_BC(nBCtypes)
    ! ---------------------------------------------------------------------------
    integer :: i
    type(aot_out_type) :: conf ! aotus lua state to write output
    ! ---------------------------------------------------------------------------

    call aot_out_open( conf, headerfile )
    call aot_out_open_table( conf, propname )
    do i = 1, nBCtypes
      call aot_out_val( conf, val = flag_per_BC(i) )
    end do
    call aot_out_close_table(conf)
    ! close the mesh header file
    call aot_out_close(conf)
  end subroutine dump_tem_BC_logicalHeader
! ****************************************************************************** !


! ****************************************************************************** !
  !> dump bc properties
  subroutine dump_tem_BC_realArray( offset, arraylen, nElems, propdat, basename, &
    &                               comm                                         )
    ! ---------------------------------------------------------------------------
    !> Offset of the local set of elements in the global list
    integer(kind=long_k), intent(in) :: offset
    !> Length of real array for each element
    integer, intent(in) :: arraylen
    !> Local number of elements with this property
    integer, intent(in) :: nElems
    !> Real-valued property data for each element to write out
    real(kind=rk), intent(in) :: propdat(arraylen, nElems)
    !> Name of the file, the data is stored in, will be appended with
    !! ".lua" for the header information and ".lsb" or ".msb" for the
    !! binary data.
    character(len=*), intent(in) :: basename
    !> Communicator to use
    integer, intent(in) :: comm
    ! ---------------------------------------------------------------------------
    integer, parameter :: root = 0
    character(len=256) :: datafile
    character(len=4) :: EndianSuffix
    ! ---------------------------------------------------------------------------
    integer(kind=MPI_OFFSET_KIND)     :: displacement
    integer :: fh, etype, ftype, ierror, iostatus( MPI_STATUS_SIZE ), typesize
    ! ---------------------------------------------------------------------------

    EndianSuffix = tem_create_EndianSuffix()
    datafile = trim(basename)//trim(EndianSuffix)

    if (nElems > 0) then
      !> Mpi IO
      write(logUnit(1),*) 'Write qVal to file: ' // trim(datafile)

      ! Open the binary file for MPI I/O (Write)
      call MPI_FILE_OPEN( comm, trim(datafile),                  &
        &                 ior(MPI_MODE_WRONLY,MPI_MODE_CREATE),  &
        &                 MPI_INFO_NULL, fh, iError              )
      call check_mpi_error( iError,'File open in dump_tem_BC_realArray')

      ! Create a contiguous type to describe the vector per element
      call MPI_TYPE_CONTIGUOUS( arrayLen, rk_mpi, etype, iError )
      call check_mpi_error( iError,'contiguous etype in dump_tem_BC_realArray')
      call MPI_TYPE_COMMIT(     etype, iError )
      call check_mpi_error( iError,'commit etype in dump_tem_BC_realArray')
      call MPI_TYPE_SIZE(etype, typesize, iError )
      call check_mpi_error(iError,'typesize in dump_tem_BC_realArray')

      ! Calculate displacement for file view
      displacement = offset * typesize * 1_MPI_OFFSET_KIND

      ! Create a MPI CONTIGUOUS as ftype for file view
      call MPI_TYPE_CONTIGUOUS(nElems, etype, ftype, iError)
      call check_mpi_error( iError,'contiguous ftype in dump_tem_BC_realArray')
      call MPI_TYPE_COMMIT( ftype, iError )
      call check_mpi_error( iError,'commit ftype in dump_tem_BC_realArray')

      ! Set the view for each process on the file above
      call MPI_FILE_SET_VIEW( fh, displacement, etype, ftype, "native",  &
        &                     MPI_INFO_NULL, iError )
      call check_mpi_error( iError,'set File view in dump_tem_BC_realArray')

      ! Read data from the file
      call MPI_FILE_WRITE_ALL( fh, propdat, nElems, etype, iostatus, iError )
      call check_mpi_error( iError,'File write all in dump_tem_BC_realArray')

      !Free the MPI_Datatypes which were created and close the file
      call MPI_TYPE_FREE (etype, iError)
      call check_mpi_error( iError,'free etype in dump_tem_BC_realArray')
      call MPI_TYPE_FREE (ftype, iError)
      call check_mpi_error( iError,'free ftype in dump_tem_BC_realArray')
      call MPI_FILE_CLOSE(fh,    iError)
      call check_mpi_error( iError,'close file in dump_tem_BC_realArray')
      ! END IO-part
    end if

  end subroutine dump_tem_BC_realArray
! ****************************************************************************** !


! ****************************************************************************** !
  !> load bc realarray data from disk
  !!
  subroutine load_tem_BC_logicalHeader( nBCtypes, propname, basename, &
    &                                   flag_per_BC, myPart, comm     )
    ! ---------------------------------------------------------------------------
    !> Number of boundary condition types
    integer, intent(in)                   :: nBCtypes
    !> Name of the property to load
    character(len=*), intent(in)          :: propname
    !> Name of the file, the data is stored in, will be appended with
    !! ".lua" for the header information and ".lsb" or ".msb" for the
    !! binary data.
    character(len=*), intent(in)          :: basename
    !> The flags to set for each boundary condition type
    logical, intent(out)                  :: flag_per_BC(nBCtypes)
    !> Partition to load
    integer, intent(in)                   :: myPart
    !> Communicator to use
    integer, intent(in)                   :: comm
    ! ---------------------------------------------------------------------------
    type( flu_State ) :: conf ! lua flu state to read lua file
    integer :: i
    integer :: iError
    integer, parameter :: root = 0
    character(len=256) :: headerfile
    integer :: thandle
    ! ---------------------------------------------------------------------------

    ! set header file name
    headerfile = trim(basename)//'.lua'

    if (myPart == root) then
      ! Now read whether each boundary has q value
      call open_config_file(L = conf, filename = headerfile)
      call aot_table_open(L = conf, thandle = thandle, key = propName)
      do i=1,nBCtypes
        call aot_get_val( L       = conf,           &
          &               thandle = thandle,        &
          &               pos     = i,              &
          &               val     = flag_per_BC(i), &
          &               ErrCode = iError          )
      end do
      call aot_table_close( L = conf, thandle = thandle )
      call close_config( conf )
    end if

    call MPI_Bcast(flag_per_BC, nBCtypes, MPI_LOGICAL, root, comm, iError)
  end subroutine load_tem_BC_logicalHeader
! ****************************************************************************** !


! ****************************************************************************** !
  !> load bc realarray data from disk
  !!
  subroutine load_tem_BC_realArray( offset, propname, arraylen, nElems, &
    &                               propdat, basename, myPart, comm     )
    ! ---------------------------------------------------------------------------
    !> Offset of the local set of elements in the global list
    integer(kind=long_k), intent(in)      :: offset
    !> Name of the property to load
    character(len=*), intent(in)          :: propname
    !> Length of the real data array to read per element
    integer, intent(in)                   :: arraylen
    !> Local number of elements with this property
    integer, intent(in)                   :: nElems
    !> real array data to fill
    real(kind=rk), intent(out)            :: propdat(arraylen, nElems)
    !> Name of the file, the data is stored in, will be appended with
    !! ".lua" for the header information and ".lsb" or ".msb" for the
    !! binary data.
    character(len=*), intent(in)          :: basename
    !> Partition to load
    integer, intent(in)                   :: myPart
    !> Communicator to use
    integer, intent(in)                   :: comm
    ! ---------------------------------------------------------------------------
    integer :: i
    integer, parameter :: root = 0
    integer :: propcomm
    integer :: color, iError
    logical :: participant !< If the local rank is a participant in Qval
    character(len=4) :: EndianSuffix
    character(len=256) :: datafile
    real(kind=rk), allocatable :: buffer(:)
    integer(kind=MPI_OFFSET_KIND)     :: displacement
    integer :: fh, etype, ftype, iostatus( MPI_STATUS_SIZE ), typesize
    ! ---------------------------------------------------------------------------

    ! set binary file name
    EndianSuffix = tem_create_EndianSuffix()
    datafile = trim(basename)//trim(EndianSuffix)

    participant = ( nElems > 0 )

    If( participant ) then
      color = 1
    else
      color = MPI_UNDEFINED
    end if

    ! Split the communicator
    call MPI_COMM_SPLIT(comm, color, myPart, propcomm, iError)

    if (nElems > 0) then
      write(logUnit(1), *) 'Load '//propname//' from file: '//trim(datafile)

      allocate( buffer( arraylen * nElems ) )

      ! Open the binary file for MPI I/O (Write)
      call MPI_FILE_OPEN( propcomm, trim(datafile), MPI_MODE_RDONLY, &
        &                 MPI_INFO_NULL, fh, iError                  )
      call check_mpi_error( iError,'File open in load_tem_BC_realarray')

      ! Create a contiguous type to describe the vector per element
      call MPI_TYPE_CONTIGUOUS( arraylen, rk_mpi, etype, iError )
      call check_mpi_error( iError,'contiguous etype in load_tem_BC_realarray')
      call MPI_TYPE_COMMIT(     etype, iError )
      call check_mpi_error( iError,'commit etype in load_tem_BC_realarray')
      call MPI_TYPE_SIZE(etype, typesize, iError )
      call check_mpi_error(iError,'typesize in load_tem_BC_realarray')

      ! Calculate displacement for file view
      displacement = offset * typesize * 1_MPI_OFFSET_KIND


      ! Create a MPI CONTIGUOUS as ftype for file view
      call MPI_TYPE_CONTIGUOUS(nElems, etype, ftype, iError)
      call check_mpi_error( iError,'contiguous ftype in load_tem_BC_realarray')
      call MPI_TYPE_COMMIT( ftype, iError )
      call check_mpi_error( iError,'commit ftype in load_tem_BC_realarray')

      ! Set the view for each process on the file above
      call MPI_FILE_SET_VIEW( fh, displacement, etype, ftype, "native",  &
        &                     MPI_INFO_NULL, iError )
      call check_mpi_error( iError,'set File view in load_tem_BC_realarray' )

      ! Read data from the file
      call MPI_FILE_READ_ALL( fh, buffer, nElems, etype, iostatus, iError )
      call check_mpi_error( iError,'File read all in load_tem_BC_realarray')

      !Free the MPI_Datatypes which were created and close the file
      call MPI_TYPE_FREE (etype, iError)
      call check_mpi_error( iError,'free etype in load_tem_BC_realarray')
      call MPI_TYPE_FREE (ftype, iError)
      call check_mpi_error( iError,'free ftype in load_tem_BC_realarray')
      call MPI_FILE_CLOSE(fh,    iError)
      call check_mpi_error( iError,'close file in load_tem_BC_realarray')
      ! END IO-part

      do i=1,nElems
        propdat(:,i) = buffer( ((i-1)*arraylen+1) : (i*arraylen) )
      end do
      deallocate( buffer )

    end if

  end subroutine load_tem_BC_realarray
! ****************************************************************************** !

! ****************************************************************************** !
  !> dump qval header information to lua file
  !!
  subroutine dump_tem_BC_qValHeader( me, headerfile )
    ! ---------------------------------------------------------------------------
    !> Boundary condition construct to load the data into
    type(tem_BC_prop_type), intent(in) :: me
    !> name of the bc header lua file
    character(len=*), intent(in) :: headerfile
    ! ---------------------------------------------------------------------------

    call dump_tem_BC_logicalHeader( headerfile, 'hasQVal', me%nBCtypes, &
      &                             me%hasQVal                          )

  end subroutine dump_tem_BC_qValHeader
! ****************************************************************************** !


! ****************************************************************************** !
  !> dump bc properties
  subroutine dump_tem_BC_qVal( me, offset, nElems, basename, myPart, comm )
    ! ---------------------------------------------------------------------------
    !> Boundary condition construct to load the data into
    type(tem_BC_prop_type), intent(in) :: me
    !> Offset of the local set of elements in the global list
    integer(kind=long_k), intent(in) :: offset
    !> Local number of elements with this property
    integer, intent(in) :: nElems
    !> Name of the file, the data is stored in, will be appended with
    !! ".lua" for the header information and ".lsb" or ".msb" for the
    !! binary data.
    character(len=*), intent(in) :: basename
    !> Partition to dump
    integer, intent(in) :: myPart
    !> Communicator to use
    integer, intent(in) :: comm
    ! ---------------------------------------------------------------------------
    integer, parameter :: root = 0
    character(len=256) :: headerfile
    ! ---------------------------------------------------------------------------

    headerfile = trim(basename)//'.lua'

    if (myPart == root) then
      ! Only root partition needs to write the header
      !open up the mesh header lua file to dump the stuff using aotus library
      call dump_tem_BC_qValHeader( me, headerfile )
    end if

    call dump_tem_BC_realArray( offset   = offset,    &
      &                         arraylen = me%nSides, &
      &                         nElems   = nElems,    &
      &                         propdat  = me%qval,   &
      &                         basename = basename,  &
      &                         comm     = comm       )

  end subroutine dump_tem_BC_qVal
! ****************************************************************************** !


! ****************************************************************************** !
  !> load bc qVal header from lua file, qVal from qVal.lsb
  !!
  subroutine load_tem_BC_qVal( me, offset, nElems, basename, myPart, comm )
    ! ---------------------------------------------------------------------------
    !> Boundary condition construct to load the data into
    type(tem_BC_prop_type), intent(inout) :: me
    !> Offset of the local set of elements in the global list
    integer(kind=long_k), intent(in)      :: offset
    !> Local number of elements that have qVal
    integer, intent(in)                   :: nElems
    !> Name of the file, the data is stored in, will be appended with
    !! ".lua" for the header information and ".lsb" or ".msb" for the
    !! binary data.
    character(len=*), intent(in)          :: basename
    !> Partition to load
    integer, intent(in)                   :: myPart
    !> Communicator to use
    integer, intent(in)                   :: comm
    ! ---------------------------------------------------------------------------
    ! ---------------------------------------------------------------------------

    call load_tem_BC_logicalHeader( nBCtypes = me%nBCtypes,   &
      &                             propName = 'hasQVal',     &
      &                             basename = basename,      &
      &                             flag_per_BC = me%hasQVal, &
      &                             myPart   = myPart,        &
      &                             comm     = comm           )

    allocate(me%qVal(me%nSides, nElems))

    call load_tem_BC_realArray( offset   = offset,    &
      &                         propname = 'hasQVal', &
      &                         arraylen = me%nSides, &
      &                         nElems   = nElems,    &
      &                         propdat  = me%qVal,   &
      &                         basename = basename,  &
      &                         myPart   = myPart,    &
      &                         comm     = comm       )

  end subroutine load_tem_BC_qVal
! ****************************************************************************** !


! ****************************************************************************** !
  !> dump normal header information to lua file
  !!
  subroutine dump_tem_BC_NormalHeader( me, headerfile )
    ! ---------------------------------------------------------------------------
    !> Boundary condition construct to load the data into
    type(tem_BC_prop_type), intent(in) :: me
    !> name of the bc header lua file
    character(len=*), intent(in) :: headerfile
    ! ---------------------------------------------------------------------------
    ! ---------------------------------------------------------------------------

    call dump_tem_BC_logicalHeader( headerfile, 'hasNormal', me%nBCtypes, &
      &                             me%hasNormal                          )

  end subroutine dump_tem_BC_NormalHeader
! ****************************************************************************** !


! ****************************************************************************** !
  !> dump normal information
  subroutine dump_tem_BC_normal( me, offset, nElems, basename, myPart, comm )
    ! ---------------------------------------------------------------------------
    !> Boundary condition construct to load the data into
    type(tem_BC_prop_type), intent(in) :: me
    !> Offset of the local set of elements in the global list
    integer(kind=long_k), intent(in) :: offset
    !> Local number of elements with this property
    integer, intent(in) :: nElems
    !> Name of the file, the data is stored in, will be appended with
    !! ".lua" for the header information and ".lsb" or ".msb" for the
    !! binary data.
    character(len=*), intent(in) :: basename
    !> Partition to dump
    integer, intent(in) :: myPart
    !> Communicator to use
    integer, intent(in) :: comm
    ! ---------------------------------------------------------------------------
    integer, parameter :: root = 0
    character(len=256) :: headerfile
    ! ---------------------------------------------------------------------------

    headerfile = trim(basename)//'.lua'
    if (myPart == root) then
      ! Only root partition needs to write the header
      !open up the mesh header lua file to dump the stuff using aotus library
      call dump_tem_BC_NormalHeader( me, headerfile )
    end if

    call dump_tem_BC_realArray( offset   = offset,    &
      &                         arraylen = 3,         &
      &                         nElems   = nElems,    &
      &                         propdat  = me%normal, &
      &                         basename = basename,  &
      &                         comm     = comm       )

  end subroutine dump_tem_BC_normal
! ****************************************************************************** !


! ****************************************************************************** !
  !> load bc qVal header from lua file, qVal from qVal.lsb
  !!
  subroutine load_tem_BC_normal( me, offset, nElems, basename, myPart, comm )
    ! ---------------------------------------------------------------------------
    !> Boundary condition construct to load the data into
    type(tem_BC_prop_type), intent(inout) :: me
    !> Offset of the local set of elements in the global list
    integer(kind=long_k), intent(in)      :: offset
    !> Local number of elements that have qVal
    integer, intent(in)                   :: nElems
    !> Name of the file, the data is stored in, will be appended with
    !! ".lua" for the header information and ".lsb" or ".msb" for the
    !! binary data.
    character(len=*), intent(in)          :: basename
    !> Partition to load
    integer, intent(in)                   :: myPart
    !> Communicator to use
    integer, intent(in)                   :: comm
    ! ---------------------------------------------------------------------------
    ! ---------------------------------------------------------------------------

    call load_tem_BC_logicalHeader( nBCtypes = me%nBCtypes,     &
      &                             propName = 'hasNormal',     &
      &                             basename = basename,        &
      &                             flag_per_BC = me%hasNormal, &
      &                             myPart   = myPart,          &
      &                             comm     = comm             )

    allocate(me%normal(3, nElems))

    call load_tem_BC_realArray( offset   = offset,      &
      &                         propname = 'hasNormal', &
      &                         arraylen = 3,           &
      &                         nElems   = nElems,      &
      &                         propdat  = me%normal,   &
      &                         basename = basename,    &
      &                         myPart   = myPart,      &
      &                         comm     = comm         )

  end subroutine load_tem_BC_normal
! ****************************************************************************** !


! ****************************************************************************** !
  subroutine tem_debug_bc_prop( me )
    ! ---------------------------------------------------------------------------
    !> Boundary condition construct to load the data into
    type(tem_BC_prop_type), intent(in) :: me
    ! ---------------------------------------------------------------------------
    integer :: iBC, iElem, iVal
    character(len=256) :: writeBuf
    ! ---------------------------------------------------------------------------


    write(dbgUnit(1), "(A)") ''
    write(dbgUnit(1), "(A)") '-------------------------------------------------'
    write(dbgUnit(1), "(A)") '           BC prop information'
    ! debug output bcID and labels
    write(dbgUnit(1), "(A)")    ''
    write(dbgUnit(1), "(A,I0)") ' nSides: '  , me%nSides
    write(dbgUnit(1), "(A,I0)") ' nBCtypes: ', me%nBCtypes
    write(dbgUnit(1), "(A,I0)") ' nElems of bc_prop: ', me%property%nElems
    write(dbgUnit(1), "(A)")    ''

    do iBC = 1, me%nBCtypes
      write(dbgUnit(1), "(A,I2,A,L1)") ' bcID: ', iBC, &
        &                           ', label: '//trim(me%BC_label( iBC )), &
        &                           ', has qVal? ', me%hasQVal(iBC)
    end do
    write(dbgUnit(1), "(A)")    ''

    ! debug output boundary ID
    do iElem = 1, me%property%nElems
      write(writeBuf, "(A,I0,A)") ' boundary_ID(:,', iElem, ') = '
      do iVal = 1, me%nSides
        write(writeBuf, "(A, I3)") trim(writeBuf), me%boundary_ID( iVal, iElem )
      end do
      write(dbgUnit(1), "(A)") trim(writeBuf)
    end do

    write(dbgUnit(1), "(A)") '-------------------------------------------------'
    write(dbgUnit(1), "(A)") ''

  end subroutine tem_debug_bc_prop
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine tem_unload_BC_prop( me )
    ! ---------------------------------------------------------------------------
    !> Boundary property type
    type(tem_BC_prop_type), intent(inout) :: me
    ! ---------------------------------------------------------------------------

    deallocate( me%BC_label    )
    deallocate( me%boundary_ID )
    deallocate( me%hasQVal     )
    if (allocated(me%qVal)) deallocate( me%qVal        )

    write(logUnit(7), "(A)") 'Unload boundary property type.'

  end subroutine tem_unload_BC_prop
! ****************************************************************************** !

end module tem_bc_prop_module
! ****************************************************************************** !
