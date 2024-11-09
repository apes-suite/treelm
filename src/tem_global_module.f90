! Copyright (c) 2011-2016, 2019, 2022 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2011-2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2012 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011-2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2012, 2014, 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013 Kartik Jain <kartik.jain@uni-siegen.de>
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
! **************************************************************************** !
!> This module provides the global part of the TREe based ELemental Mesh, which
!! is needed on each partition.
!!
module tem_global_module

  ! include treelm modules
  use mpi
  use env_module,          only: rk, long_k, rk_mpi, long_k_mpi, labelLen, &
    &                            pathLen
  use tem_aux_module,      only: tem_abort
  use tem_logging_module,  only: logUnit
  use tem_prophead_module, only: tem_prophead_type, load_tem_prophead, &
    &                            dump_tem_prophead

  ! include aotus modules
  use aot_out_module,   only: aot_out_type, aot_out_val,            &
    &                         aot_out_open, aot_out_close,          &
    &                         aot_out_open_table, aot_out_close_table
  use aot_table_module, only: aot_table_open, aot_table_close
  use aotus_module,     only: open_config_file, aot_get_val, close_config, &
    &                         flu_State

  implicit none
  private

  !> A datatype describing the global mesh information, present in all
  !! partitions (on all processes).
  !!
  !! It includes the mesh name, the origin and extent of the universe, meaning
  !! the bounding cube, inside which the tree mesh was created.
  !! Further, the lowest and highest refinement levels are stored as well as the
  !! number of different properties, assigned to elements in the mesh.
  !! The properties are linked to in-detail descriptive data types themselves.
  type tem_global_type
    !> Length of the outermost bounding cube
    real(kind=rk) :: BoundingCubeLength
    !> Origin of the cube, extent of length in each direction from this point.
    real(kind=rk) :: Origin(3)
    !> Origin of the cube which is effectively covered by tree elements
    real(kind=rk) :: effOrigin(3) = 0._rk
    !> Length of the hexahedral domain which is effectively covered by
    !! tree elements
    real(kind=rk) :: effLength(3) = 0._rk
    !> name of the predefined mesh generated on the fly
    character(len=labelLen) :: predefined = ''
    !> total number of elements
    integer(kind=long_k) :: nElems

    ! @todo: these MPI data should NOT belong to global_type
    !        only tree related global information should be here
    !> Number of parts, the mesh is partitioned into !< process%comm_size
    integer :: nParts
    !> The part to be processed locally (0:nParts-1) !< process%rank
    integer :: myPart
    !> MPI communicator of the processes, the mesh is distributed on
    integer :: comm

    !> Minimal element level present in the mesh
    integer :: minLevel
    !> Maximal element level present in the mesh
    integer :: maxLevel
    !> Name of the mesh
    character(len=LabelLen) :: label
    !> Space for comments
    character(len=LabelLen) :: comment
    !> Total number of properties present in the mesh
    !! properties: fluid, solid, has boundaries, has elements, is deformed
    integer :: nProperties
    !> The information of real bounding cube which is read from the mesh header
    real(kind=rk) :: effboundingCube(3,2)

    !> declarations for each property (The bit position of the flag for the
    !! property, the total number of elements with it, and label for its
    !! declaration).
    !! For a detailed description see \ref tem_property_module
    type(tem_prophead_type), pointer :: Property(:) => null()

    !> Name of the directory, where the data for the mesh is read from.
    !! Has to have a trailing directory separator.
    character(len=pathlen) :: dirname

    !> did the mesh changed since last time dumping it??
    logical :: meshChange = .false.
  end type


  public :: tem_global_type
  public :: tem_global_mesh_free
  public :: tem_mesh_out
  public :: load_tem_global
  public :: dump_tem_global
  public :: tem_global_mesh_read


contains


  ! ************************************************************************ !
  !> write mesh information into lua file.
  !!
  !! This routine is specially need for write restart to write mesh info for
  !! predefined mesh such that read restart will generate predefined mesh
  !!
  subroutine tem_mesh_out( me, conf )
    ! -------------------------------------------------------------------- !
    !> Structure to store header in
    type(tem_global_type), intent(in) :: me
    !> aotus lua state to write output
    type(aot_out_type), intent(inout) :: conf
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: lineLength
    ! -------------------------------------------------------------------- !

    select case(trim(me%predefined))
      case('cube','slice','line')
        call aot_out_open_table(conf, 'mesh')

        call aot_out_val( put_conf = conf,                &
          &               val      = trim(me%predefined), &
          &               vname    = 'predefined'         )
        call aot_out_val( put_conf = conf,      &
          &               val      = me%Origin, &
          &               vname    = 'origin'   )
        if (me%nElems < 2**me%maxlevel) then
          lineLength = me%BoundingCubeLength*0.5_rk**me%maxLevel*me%nElems
          call aot_out_val( put_conf = conf,       &
            &               val      = lineLength, &
            &               vname    = 'length'    )
          call aot_out_val( put_conf = conf,           &
            &               val      = me%nElems,      &
            &               vname    = 'element_count' )
        else
          call aot_out_val( put_conf = conf,                  &
            &               val      = me%BoundingCubeLength, &
            &               vname    = 'length'               )
          call aot_out_val( put_conf = conf,             &
            &               val      = me%maxLevel,      &
            &               vname    = 'refinementLevel' )
        end if

        call aot_out_close_table(conf)

      case('line_bounded')
        lineLength = me%BoundingCubeLength*0.5_rk**me%maxLevel*me%nElems
        call aot_out_open_table(conf, 'mesh')

        call aot_out_val( put_conf = conf,                &
          &               val      = trim(me%predefined), &
          &               vname    = 'predefined'         )
        call aot_out_val( put_conf = conf,      &
          &               val      = me%Origin, &
          &               vname    = 'origin'   )
        call aot_out_val( put_conf = conf,       &
          &               val      = lineLength, &
          &               vname    = 'length'    )
        call aot_out_val( put_conf = conf,           &
          &               val      = me%nElems,      &
          &               vname    = 'element_count' )

        call aot_out_close_table(conf)

      case('single')
        call aot_out_open_table(conf, 'mesh')

        call aot_out_val( put_conf = conf,                &
          &               val      = trim(me%predefined), &
          &               vname    = 'predefined'         )
        call aot_out_val( put_conf = conf,      &
          &               val      = me%Origin, &
          &               vname    = 'origin'   )
        call aot_out_val( put_conf = conf,                  &
          &               val      = me%BoundingCubeLength, &
          &               vname    = 'length'               )
        call aot_out_close_table(conf)

      case default
        call aot_out_val( put_conf = conf,             &
          &               val      = trim(me%dirName), &
          &               vname    = 'mesh'            )

    end select

  end subroutine tem_mesh_out
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> @todo: Add description
  !!
  subroutine tem_global_mesh_read( me, conf, myPart, nParts, comm )
    ! -------------------------------------------------------------------- !
    !> Structure to store header in
    type(tem_global_type), intent(out) :: me
    !> lua flu state to read mesh from
    type( flu_State ) :: conf
    !> The process local part (= MPI Rank in comm)
    integer, intent(in) :: myPart
    !> Number of partitions, the mesh is partitioned into (= Number of MPI
    !! processes in comm).
    integer, intent(in) :: nParts
    !> MPI Communicator to use
    integer, intent(in) :: comm
    ! -------------------------------------------------------------------- !
    character(len=pathLen) :: dirname
    integer :: commLocal
    integer :: tem_handle
    integer :: iError
    ! -------------------------------------------------------------------- !

    ! Use the incoming communicator
    commLocal = comm

    write(logUnit(1),*) 'Obtaining HEADER of the configured mesh'
    call aot_table_open(L=conf, thandle=tem_handle, key='mesh')

    if (tem_handle /= 0) then
      ! The mesh is actually given as a table, parse it and
      ! generate a mesh internally.
      write(logUnit(1),*)'  Generating HEADER for an internally defined mesh'
      call tem_global_mesh_internal( me, conf, tem_handle, myPart, nParts,     &
        &                            commLocal )
    else
      ! The mesh is not a table, try to interpret it as a string.
      call aot_get_val( L       = conf,                                        &
        &               key     = 'mesh',                                      &
        &               val     = dirname,                                     &
        &               ErrCode = iError,                                      &
        &               default = 'mesh/' )
      call load_tem_global(me, trim(dirname), myPart, nParts, commLocal)
    end if
  end subroutine tem_global_mesh_read
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> A routine to load global informations from the header file in the given
  !! directory.
  !!
  subroutine load_tem_global( me, dirname, myPart, nParts, comm )
    ! -------------------------------------------------------------------- !
    !> Structure to store header in
    type(tem_global_type), intent(out) :: me
    !> Directory containing the mesh informations
    character(len=*), intent(in) :: dirname
    !> The process local part (= MPI Rank in comm)
    integer, intent(in) :: myPart
    !> Number of partitions, the mesh is partitioned into (= Number of MPI
    !! processes in comm).
    integer, intent(in) :: nParts
    !> MPI Communicator to use
    integer, intent(in) :: comm
    ! -------------------------------------------------------------------- !
    character(len=300) :: headname
    integer :: iError
    integer :: root
    integer :: i
    logical :: ex
    integer :: thandle, sub_handle
    type( flu_State ) :: conf ! lua flu state to read lua file
    ! -------------------------------------------------------------------- !

    root = 0

    me%comm = comm

    me%myPart = myPart
    me%nParts = nParts
    me%dirname = trim(adjustl(dirname))
    headname = trim(me%dirname)//'header.lua'
    write(logUnit(1), *) 'Load mesh header from file: '//trim(headname)

    if (myPart == root) then
      inquire(file=trim(headname), exist=ex)
      if (.not. ex) then
        write(*,*) 'File ',trim(headname),' not found. Aborting.'
        stop
      endif
      !! Read the header only on the root process, broadcast to all others
      ! open mesh header file
      call open_config_file(L = conf, filename = trim(headname))
      ! load label
      call aot_get_val( L       = conf,     &
        &               key     = 'label',  &
        &               val     = me%label, &
        &               ErrCode = iError    )
      call aot_get_val( L       = conf,       &
        &               key     = 'comment',  &
        &               val     = me%comment, &
        &               ErrCode = iError      )

      ! Open boundingbox table
      call aot_table_open( L = conf, thandle = thandle, key='boundingbox' )

      ! Read the origin
      call aot_table_open( L       = conf,       &
        &                  parent  = thandle,    &
        &                  thandle = sub_handle, &
        &                  key     = 'origin'    )
      do i = 1,3
        call aot_get_val( L       = conf,         &
          &               thandle = sub_handle,   &
          &               pos     = i,            &
          &               val     = me%origin(i), &
          &               ErrCode = iError        )
      end do
      call aot_table_close( L = conf, thandle = sub_handle )

      ! Read the bounding cube length
      call aot_get_val( L       = conf,                  &
        &               thandle = thandle,               &
        &               key     = 'length',              &
        &               val     = me%BoundingCubeLength, &
        &               ErrCode = iError                 )
      ! Close boundingbox table again
      call aot_table_close( L = conf, thandle = thandle )

      call aot_get_val( L       = conf,      &
        &               key     = 'nElems',  &
        &               val     = me%nElems, &
        &               ErrCode = iError     )

      call aot_get_val( L       = conf,        &
        &               key     = 'minLevel',  &
        &               val     = me%minLevel, &
        &               ErrCode = iError       )

      call aot_get_val( L       = conf,        &
        &               key     = 'maxLevel',  &
        &               val     = me%maxLevel, &
        &               ErrCode = iError       )

      call aot_get_val( L       = conf,           &
        &               key     = 'nProperties',  &
        &               val     = me%nProperties, &
        &               ErrCode = iError          )

      ! Read the effective bounding cube parameters
      ! Open the effective bounding box table
      call aot_table_open( L = conf, thandle = thandle, key='effBoundingbox' )
      ! Read the origin
      call aot_table_open( L       = conf,       &
        &                  parent  = thandle,    &
        &                  thandle = sub_handle, &
        &                  key     = 'origin'    )
      do i = 1,3
        call aot_get_val( L       = conf,                    &
          &               thandle = sub_handle,              &
          &               pos     = i,                       &
          &               val     = me%effboundingcube(i,1), &
          &               ErrCode = iError                   )
      end do
      call aot_table_close( L = conf, thandle = sub_handle )
      me%effOrigin = me%effboundingcube(:,1)
      ! Read the effective length (min and max)
      call aot_table_open( L       = conf,       &
        &                  parent  = thandle,    &
        &                  thandle = sub_handle, &
        &                  key     = 'effLength' )
      do i = 1,3
        call aot_get_val( L       = conf,            &
          &               thandle = sub_handle,      &
          &               pos     = i,               &
          &               val     = me%effLength(i), &
          &               ErrCode = iError           )
        me%effboundingcube(i,2) = me%effboundingcube(i,1) + me%effLength(i)
      end do
      call aot_table_close( L = conf, thandle = sub_handle )

    end if
    write(logUnit(1),*) 'The real bounding cube is...'
    write(logUnit(1),*) '  min: ',me%effBoundingCube(:,1)
    write(logUnit(1),*) '  max: ',me%effBoundingCube(:,2)

    !! Broadcast the header informations to all processes.
    call MPI_Bcast(me%nElems, 1, long_k_mpi, root, me%comm, iError)
    call MPI_Bcast(me%label, LabelLen, MPI_CHARACTER, root, me%comm, iError)
    call MPI_Bcast(me%comment, LabelLen, MPI_CHARACTER, root, me%comm, iError)
    call MPI_Bcast(me%BoundingCubeLength, 1, rk_mpi, root, me%comm, iError)
    call MPI_Bcast(me%Origin, 3, rk_mpi, root, me%comm, iError)
    call MPI_Bcast(me%minLevel, 1, MPI_INTEGER, root, me%comm, iError)
    call MPI_Bcast(me%maxLevel, 1, MPI_INTEGER, root, me%comm, iError)
    call MPI_Bcast(me%nProperties, 1, MPI_INTEGER, root, me%comm, iError)
    call MPI_Bcast(me%effBoundingCube, 6, rk_mpi, root, me%comm, iError)

    if (associated(me%Property)) deallocate(me%property)
    allocate(me%Property(me%nProperties))

    call load_tem_prophead( me     = me%Property, &
      &                     myPart = myPart,      &
      &                     comm   = me%comm,     &
      &                     conf   = conf,        &
      &                     root   = root         )

    if (myPart == root) then
      call close_config(conf)
    end if

  end subroutine load_tem_global
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> @todo: Add description
  !!
  subroutine tem_global_mesh_internal( me, conf, thandle, myPart, nParts, comm )
    ! -------------------------------------------------------------------- !
    !> Structure to load the mesh to
    type(tem_global_type), intent(out) :: me
    !> Directory containing the mesh informations
    type(flu_State) :: conf
    !> Handle for the table to read the description
    !! of the mesh from.
    integer, intent(in) :: thandle
    !> Partition to use on the calling process (= MPI Rank in comm)
    integer, intent(in) :: myPart
    !> Number of partitions, the mesh is partitioned into (= Number of MPI
    !! processes in comm).
    integer, intent(in) :: nParts
    !> MPI Communicator to use
    integer, intent(in) :: comm
    ! -------------------------------------------------------------------- !
    character(len=labelLen) :: meshtype
    integer :: iError
    ! -------------------------------------------------------------------- !

    call aot_get_val( L       = conf,         &
      &               thandle = thandle,      &
      &               val     = meshtype,     &
      &               ErrCode = iError,       &
      &               key     = 'predefined', &
      &               default = 'cube'        )

    select case(trim(meshtype))
    case('cube')
      ! Generate a single level full cube mesh
      call gen_treelm_cube_global( me, conf, thandle,   &
        &                          myPart, nParts, comm )
    case('slice')
      ! Generate a single level slice
      call gen_treelm_slice_global( me, conf, thandle,   &
        &                           myPart, nParts, comm )
    case('line', 'line_bounded')
      ! Generate a single level line
      call gen_treelm_line_global( me, conf, thandle,                   &
        &                          myPart, nParts, comm, trim(meshtype) )
    case('single')
      ! Generate a single level slice
      call gen_treelm_single_global( me, conf, thandle,   &
        &                            myPart, nParts, comm )
    case('default')
      write(logUnit(1),*)'Do not know how to generate mesh ' &
        &            //trim(meshtype)
      call tem_abort()
    end select
  end subroutine tem_global_mesh_internal
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Generate the header for the simple full cube mesh.
  !!
  subroutine gen_treelm_cube_global( me, conf, thandle, myPart, nParts, comm )
    ! -------------------------------------------------------------------- !
    !> Structure to load the mesh to
    type(tem_global_type), intent(out) :: me
    !> Directory containing the mesh informations
    type(flu_State) :: conf
    !> Handle for the table to read the description
    !! of the mesh from.
    integer, intent(in) :: thandle
    !> Partition to use on the calling process (= MPI Rank in comm)
    integer, intent(in) :: myPart
    !> Number of partitions, the mesh is partitioned into (= Number of MPI
    !! processes in comm).
    integer, intent(in) :: nParts
    !> MPI Communicator to use
    integer, intent(in) :: comm
    ! -------------------------------------------------------------------- !
    integer :: iError
    integer :: orig_err(3)
    ! -------------------------------------------------------------------- !

    write(logUnit(1),*) 'Creating HEADER for a full cubical mesh ' &
      &                 // 'without boundaries'

    me%nParts = nParts
    me%myPart = myPart
    me%comm = comm

    ! Get the origin of the cube:
    call aot_get_val( L       = conf,                    &
      &               thandle = thandle,                 &
      &               key     = 'origin',                &
      &               val     = me%origin,               &
      &               ErrCode = orig_err,                &
      &               default = [0.0_rk, 0.0_rk, 0.0_rk] )

    ! Get the length of the cube:
    call aot_get_val( L       = conf,                  &
      &               thandle = thandle,               &
      &               val     = me%BoundingCubeLength, &
      &               ErrCode = iError,                &
      &               key     = 'length',              &
      &               default = 1.0_rk                 )

    ! Get the refinement level:
    call aot_get_val( L       = conf,             &
      &               thandle = thandle,          &
      &               val     = me%minlevel,      &
      &               ErrCode = iError,           &
      &               key     = 'refinementLevel' )

    me%label = 'Generic_Cube'
    me%predefined = 'cube'

    if (me%minlevel == 0) then
      me%BoundingCubeLength = me%BoundingCubeLength*2
      me%minlevel = 1
      me%label = 'Generic_Single'
      me%predefined = 'single'
    end if

    me%maxLevel = me%minLevel

    write(me%comment,'(a15,i7,a16,i2,a1)') &
      &  'Generated with ', nParts, ' parts on Level ', me%minlevel, '.'
    me%dirname = './'

    ! No properties in this mesh
    me%nProperties = 0
    if (associated(me%Property)) deallocate(me%property)
    allocate(me%Property(me%nProperties))

  end subroutine gen_treelm_cube_global
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Generate the simple single level mesh of a slice in the full cube.
  !!
  subroutine gen_treelm_slice_global( me, conf, thandle, myPart, nParts, comm )
    ! -------------------------------------------------------------------- !
    !> Structure to load the mesh to
    type(tem_global_type), intent(out) :: me
    !> Directory containing the mesh informations
    type(flu_State) :: conf
    !> Handle for the table to read the description
    !! of the mesh from.
    integer, intent(in) :: thandle
    !> Partition to use on the calling process (= MPI Rank in comm)
    integer, intent(in) :: myPart
    !> Number of partitions, the mesh is partitioned into (= Number of MPI
    !! processes in comm).
    integer, intent(in) :: nParts
    !> MPI Communicator to use
    integer, intent(in) :: comm
    ! -------------------------------------------------------------------- !
    integer :: iError
    integer :: orig_err(3)
    ! -------------------------------------------------------------------- !

    write(logUnit(1),*)'Creating HEADER for a slice mesh without boundaries'

    me%nParts = nParts
    me%myPart = myPart
    me%comm = comm

    ! Get the origin of the cube:
    call aot_get_val( L       = conf,                    &
      &               thandle = thandle,                 &
      &               key     = 'origin',                &
      &               val     = me%origin,               &
      &               ErrCode = orig_err,                &
      &               default = [0.0_rk, 0.0_rk, 0.0_rk] )

    ! Get the length of the cube:
    call aot_get_val( L       = conf,                  &
      &               thandle = thandle,               &
      &               val     = me%BoundingCubeLength, &
      &               ErrCode = iError,                &
      &               key     = 'length',              &
      &               default = 1.0_rk                 )

    ! Get the refinement level:
    call aot_get_val( L       = conf,             &
      &               thandle = thandle,          &
      &               val     = me%minlevel,      &
      &               ErrCode = iError,           &
      &               key     = 'refinementLevel' )

    me%label = 'Generic_Slice'
    me%predefined = 'slice'

    if (me%minlevel == 0) then
      me%BoundingCubeLength = me%BoundingCubeLength*2
      me%minlevel = 1
      me%label = 'Generic_Single'
      me%predefined = 'single'
    end if

    me%maxLevel = me%minLevel

    write(me%comment,'(a15,i7,a16,i2,a1)') &
      &   'Generated with ', nParts, ' parts on Level ', me%minlevel, '.'
    me%dirname = './'

    ! Only boundary property in this mesh.
    me%nProperties = 1
    if (associated(me%Property)) deallocate(me%property)
    allocate(me%Property(me%nProperties))

  end subroutine gen_treelm_slice_global
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Generate the simple single level mesh of a line in the full cube.
  !!
  subroutine gen_treelm_line_global( me, conf, thandle, myPart, nParts, comm, &
    &                                predefined )
    ! -------------------------------------------------------------------- !
    !> Structure to load the mesh to
    type(tem_global_type), intent(out) :: me
    !> Directory containing the mesh informations
    type(flu_State) :: conf
    !> Handle for the table to read the description
    !! of the mesh from.
    integer, intent(in) :: thandle
    !> Partition to use on the calling process (= MPI Rank in comm)
    integer, intent(in) :: myPart
    !> Number of partitions, the mesh is partitioned into (= Number of MPI
    !! processes in comm).
    integer, intent(in) :: nParts
    !> MPI Communicator to use
    integer, intent(in) :: comm
    character(len=*), intent(in) :: predefined
    ! -------------------------------------------------------------------- !
    integer :: iError
    integer :: orig_err(3)
    integer :: level
    integer :: elementcount
    ! -------------------------------------------------------------------- !

    write(logUnit(1),*) 'Creating HEADER for a line mesh'

    me%nParts = nParts
    me%myPart = myPart
    me%comm = comm

    ! Get the origin of the cube:
    call aot_get_val( L       = conf,                    &
      &               thandle = thandle,                 &
      &               key     = 'origin',                &
      &               val     = me%origin,               &
      &               ErrCode = orig_err,                &
      &               default = [0.0_rk, 0.0_rk, 0.0_rk] )

    ! Get the length of the cube:
    call aot_get_val( L       = conf,                  &
      &               thandle = thandle,               &
      &               val     = me%BoundingCubeLength, &
      &               ErrCode = iError,                &
      &               key     = 'length',              &
      &               default = 1.0_rk                 )

    ! Get the refinement level:
    call aot_get_val( L       = conf,              &
      &               thandle = thandle,           &
      &               val     = level,             &
      &               ErrCode = iError,            &
      &               key     = 'refinementLevel', &
      &               default = -1                 )

    ! Get the element count:
    call aot_get_val( L       = conf,            &
      &               thandle = thandle,         &
      &               val     = elementcount,    &
      &               ErrCode = iError,          &
      &               key     = 'element_count', &
      &               default = -1               )

    me%label = 'Generic_Line'
    me%predefined = trim(predefined)

    if (predefined == 'line_bounded') then
      ! Need a padding of at least 1 element, to allow boundary definitions
      ! at both ends in Ateles.
      me%minlevel = max( ceiling(log(real(elementcount+1,kind=rk)) &
        &                        / log(2.0_rk)),                 1 )
      me%BoundingCubeLength = (me%BoundingCubeLength * 2**me%minlevel) &
        &                    / elementcount
    else if (elementcount > 1) then
      me%minlevel = ceiling(log(real(elementcount,kind=rk))/log(2.0_rk))
    else if (level > 0) then
      me%minlevel = level
    else if ((level == 0) .or. (elementcount == 1)) then
      me%BoundingCubeLength = me%BoundingCubeLength*2
      me%minlevel = 1
      me%label = 'Generic_Single'
      me%predefined = 'single'
      ! Reset elementcount, to avoid overwriting of the settings later on.
      elementcount = -1
    else
      write(logunit(1),*) 'For a line you need to state either refinementLevel'
      write(logunit(1),*) 'or element_count. None of them found!'
      write(logunit(1),*) 'STOPPING'
      call tem_abort()
    end if

    if (elementcount > 1) then
      me%BoundingCubeLength = (me%BoundingCubeLength * 2**me%minlevel) &
        &                    / elementcount
    end if

    me%maxLevel = me%minLevel

    write(me%comment,'(a15,i7,a16,i2,a1)') &
      &   'Generated with ', nParts, ' parts on Level ', me%minlevel, '.'
    me%dirname = './'

    ! Only boundary property in this mesh.
    me%nProperties = 1
    if (associated(me%Property)) deallocate(me%property)
    allocate(me%Property(me%nProperties))

  end subroutine gen_treelm_line_global
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Generate the single element mesh.
  !!
  subroutine gen_treelm_single_global( me, conf, thandle, myPart, nParts, comm )
    ! -------------------------------------------------------------------- !
    !> Structure to load the mesh to
    type(tem_global_type), intent(out) :: me
    !> Directory containing the mesh informations
    type(flu_State) :: conf
    !> Handle for the table to read the description
    !! of the mesh from.
    integer, intent(in) :: thandle
    !> Partition to use on the calling process (= MPI Rank in comm)
    integer, intent(in) :: myPart
    !> Number of partitions, the mesh is partitioned into (= Number of MPI
    !! processes in comm).
    integer, intent(in) :: nParts
    !> MPI Communicator to use
    integer, intent(in) :: comm
    ! -------------------------------------------------------------------- !
    integer :: iError
    integer :: orig_err(3)
    ! -------------------------------------------------------------------- !

    write(logUnit(1),*) 'Creating HEADER for a single element mesh without' &
      &                 //' boundaries'

    me%nParts = nParts
    me%myPart = myPart
    me%comm = comm

    ! Get the origin of the cube:
    call aot_get_val( L       = conf,                    &
      &               thandle = thandle,                 &
      &               key     = 'origin',                &
      &               val     = me%origin,               &
      &               ErrCode = orig_err,                &
      &               default = [0.0_rk, 0.0_rk, 0.0_rk] )

    ! Get the length of the cube:
    call aot_get_val( L       = conf,                  &
      &               thandle = thandle,               &
      &               val     = me%BoundingCubeLength, &
      &               ErrCode = iError,                &
      &               key     = 'length',              &
      &               default = 1.0_rk                 )


    me%minlevel = 1
    me%maxLevel = me%minLevel
    me%label = 'Generic_Single'
    me%predefined = 'single'
    write(me%comment,'(a15,i7,a16,i2,a1)') &
      &   'Generated with ', nParts, ' parts on Level ', me%minlevel, '.'
    me%dirname = './'

    ! Only boundary property in this mesh.
    me%nProperties = 1
    if (associated(me%Property)) deallocate(me%property)
    allocate(me%Property(me%nProperties))

  end subroutine gen_treelm_single_global
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Routine to free allocated memory for the header again.
  subroutine tem_global_mesh_free(me)
    type(tem_global_type), intent(inout) :: me

    if (associated(me%Property)) deallocate(me%Property)
    me%effOrigin = 0.0_rk
    me%effLength = 0.0_rk
    me%predefined = ''

  end subroutine tem_global_mesh_free
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> A routine to dump global informations into the mesh header in lua format
  !!
  subroutine dump_tem_global( me, root_only)
    ! -------------------------------------------------------------------- !
    !> global type to be dumped
    type(tem_global_type), intent(in) :: me
    !> root dump global mesh when true and
    !! all process dump its own mesh when false
    logical, intent(in), optional :: root_only
    ! -------------------------------------------------------------------- !
    character(len=300) :: headname
    integer :: version
    integer :: i
    type(aot_out_type) :: conf ! aotus lua state to write output
    integer :: locroot
    logical :: root_out
    ! -------------------------------------------------------------------- !
    if(present(root_only)) then
      root_out = root_only
    else
      root_out = .true.
    endif

    if (root_out) then
      locroot = 0
    else
      locroot = me%myPart
    end if

    version = 1 !! Provide a version string, for the file format
    headname = trim(me%dirname)//'header.lua'

    if (me%mypart == locroot) then
      ! Write the header only on the root process
      ! open up the mesh header lua file to dump the stuff using aotus library
      call aot_out_open( conf, headname )
      !write header version
      call aot_out_val( put_conf = conf,      &
        &               vname    = 'version', &
        &               val      = version    )
      !write mesh label/name
      call aot_out_val( put_conf = conf,     &
        &               vname    = 'label',  &
        &               val = trim(me%label) )
      !comments about the mesh
      call aot_out_val( put_conf = conf,            &
        &               vname    = 'comment',       &
        &               val      = trim(me%comment) )
      !write bounding box information in tha table
      call aot_out_open_table( conf, 'boundingbox' )
      !create another table for origin
      call aot_out_open_table( conf, 'origin' )
      !write three values of origin
      do i = 1,3
        call aot_out_val( conf, me%Origin(i) )
      end do
      call aot_out_close_table(conf)
      !write bounding box length
      call aot_out_val( put_conf = conf,                 &
        &               vname    = 'length',             &
        &               val      = me%BoundingCubeLength )
      call aot_out_close_table(conf)
      call aot_out_val( put_conf = conf,     &
        &               vname    = 'nElems', &
        &               val      = me%nElems )
      call aot_out_val( put_conf = conf,       &
        &               vname    = 'minLevel', &
        &               val      = me%minLevel )
      call aot_out_val( put_conf = conf,       &
        &               vname    = 'maxLevel', &
        &               val      = me%maxLevel )
      call aot_out_val( put_conf = conf,          &
        &               vname    = 'nProperties', &
        &               val      = me%nProperties )

      ! write effective boundary
      call aot_out_open_table( conf, 'effBoundingbox' )
      !create another table for origin
      call aot_out_open_table( conf, 'origin' )
      !write three values of origin
      do i = 1,3
        call aot_out_val( conf, me%effOrigin(i) )
      end do
      call aot_out_close_table(conf) ! close origin table
      !create another table for origin
      call aot_out_open_table( conf, 'effLength' )
      !write three values of origin
      do i = 1,3
        call aot_out_val( conf, me%effLength(i) )
      end do
      call aot_out_close_table(conf) ! close length table
      call aot_out_close_table(conf) ! close effBoundingBox table

      ! dump property head
      call dump_tem_prophead( me     = me%Property, &
        &                     conf   = conf         )

      ! close the mesh header file
      call aot_out_close(conf)
    end if

  end subroutine dump_tem_global
  ! ************************************************************************ !


end module tem_global_module
! **************************************************************************** !
