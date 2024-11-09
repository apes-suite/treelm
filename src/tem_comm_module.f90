! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2013, 2016, 2018-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2011-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012, 2015-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012 Sathish Krishnan P S <s.krishnan@grs-sim.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016-2017 Raphael Haupt <Raphael.Haupt@student.uni-siegen.de>
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
! **************************************************************************** !!
!> This module provides the data structure for the communication
!! during the simulation.
!!
!! Several exchange methods are implemented. CoCo is heavily used here to allow
!! for a concise definition of exchanges for various data types.
!! The basic idea is to initialize the buffers, use them in exchanges and
!! finalize them when not needed anymore.
!! In the definition of the buffers, an array of positions in the original
!! linearized data array is used to describe the origin or target positions for
!! communicated data.
!!
!! @note If you introduce a new type to exchange, you will need to introduce
!!      appropriate CoCo copy statements everywhere.
!!
module tem_comm_module

  use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer

  ! include treelm modules
  use mpi
  use env_module,            only: rk, rk_mpi, long_k, long_k_mpi
  use tem_aux_module,        only: tem_abort, check_mpi_error
  use tem_logging_module,    only: logUnit, tem_toStr
  use tem_grow_array_module, only: grw_intArray_type,     &
    &                              init, append, destroy
  use tem_dyn_array_module,  only: dyn_intArray_type, append, destroy, &
    &                              PositionOfVal
  use hvs_sizeof_module,     only: c_sizeof
  use mem_for_mpi_module,    only: alloc_mpif_mem, free_mpif_mem
  use tem_sparse_comm_module, only: use_sparse_alltoall, tem_sparse_alltoall_int

  ! include aotus modules
  use flu_binding,  only: flu_State
  use aotus_module, only: aot_get_val

  implicit none

  private

  public :: tem_communication_type
  public :: tem_commPattern_type
  public :: tem_comm_dumpType
  public :: tem_load_commPattern
  public :: tem_comm_init
  public :: tem_comm_count
  public :: tem_comm_createBuffer
  public :: tem_comm_alltoall_int
  public :: tem_comm_destroy

! Use CoCo here to make this buffer type generic for different types

! Let coco actually implement the buffer type given above for specific data
! types
  public :: tem_longbuffer_type

  !> process-wise buffer for data of type integer(kind=long_k)
  !!
  !! this datatype is used to describe the exchange with a specific process, in
  !! case of explicit buffers it provides the memory for them.
  type tem_longbuffer_type

    !> explicit buffer for data to be transferred
    integer(kind=long_k), pointer :: val(:) => null()

    !> explicit buffer in memory allocated by mpi
    type(c_ptr) :: mem_mpi

    !> position in the input vector from where to read the entries in
    !! val_long
    !!
    !! @note jz: in ateles we use this to specify the positions of the cell
    !!          states that have to be sent.
    integer, allocatable :: pos(:)

    !> number of values to exchange
    !!
    !! @note jz: in ateles this variable stores the number of coefficients we
    !!          transfer, i.e. number of cells to transfer times number of
    !!          degree of freedoms per cell times the number of scalar
    !!          variables.
    integer :: nvals

    !> handle for the mpi-datatype to describe the memory access,
    !! without explicit copying in the application.
    integer :: memindexed
  end type tem_longbuffer_type
  public :: tem_intbuffer_type

  !> process-wise buffer for data of type integer
  !!
  !! this datatype is used to describe the exchange with a specific process, in
  !! case of explicit buffers it provides the memory for them.
  type tem_intbuffer_type

    !> explicit buffer for data to be transferred
    integer, pointer :: val(:) => null()

    !> explicit buffer in memory allocated by mpi
    type(c_ptr) :: mem_mpi

    !> position in the input vector from where to read the entries in
    !! val_int
    !!
    !! @note jz: in ateles we use this to specify the positions of the cell
    !!          states that have to be sent.
    integer, allocatable :: pos(:)

    !> number of values to exchange
    !!
    !! @note jz: in ateles this variable stores the number of coefficients we
    !!          transfer, i.e. number of cells to transfer times number of
    !!          degree of freedoms per cell times the number of scalar
    !!          variables.
    integer :: nvals

    !> handle for the mpi-datatype to describe the memory access,
    !! without explicit copying in the application.
    integer :: memindexed
  end type tem_intbuffer_type
  public :: tem_realbuffer_type

  !> process-wise buffer for data of type real(kind=rk)
  !!
  !! this datatype is used to describe the exchange with a specific process, in
  !! case of explicit buffers it provides the memory for them.
  type tem_realbuffer_type

    !> explicit buffer for data to be transferred
    real(kind=rk), pointer :: val(:) => null()

    !> explicit buffer in memory allocated by mpi
    type(c_ptr) :: mem_mpi

    !> position in the input vector from where to read the entries in
    !! val_real
    !!
    !! @note jz: in ateles we use this to specify the positions of the cell
    !!          states that have to be sent.
    integer, allocatable :: pos(:)

    !> number of values to exchange
    !!
    !! @note jz: in ateles this variable stores the number of coefficients we
    !!          transfer, i.e. number of cells to transfer times number of
    !!          degree of freedoms per cell times the number of scalar
    !!          variables.
    integer :: nvals

    !> handle for the mpi-datatype to describe the memory access,
    !! without explicit copying in the application.
    integer :: memindexed
  end type tem_realbuffer_type

  ! ------------------------------------------------------------------------ !
  !> Description of communication data
  type tem_communication_type

    integer :: nProcs=0     !< amount of partitions to send to

    !> partition MPI rank
    integer,allocatable :: proc(:)

    !> How many data elements need to be exchanged with proc (per process).
    integer,allocatable :: nElemsProc(:)

    !> Request handle array
    integer,allocatable :: rqHandle(:)

    !> Data element positions in the actual arrays, used to built the pos
    !! information in the actual buffers (per process).
    type(grw_intArray_type), allocatable :: elemPos(:)

    !> declare communication buffers for each variable type

    type( tem_longbuffer_type ), allocatable :: buf_long(:)
    type( tem_intbuffer_type ), allocatable :: buf_int(:)
    type( tem_realbuffer_type ), allocatable :: buf_real(:)
  end type tem_communication_type
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> General description of the communication pattern to use.
  !!
  !! Depending on the chosen style, different exchange implementations are
  !! used. This data type provides the appropriate function pointers for
  !! initialization, finalization and exchange of the buffers.
  type tem_commPattern_type
    character(len=40) :: style


  procedure(tem_exchange_long), nopass, pointer :: exchange_long
  procedure(tem_commbuf_long_init), nopass, pointer :: initbuf_long
  procedure(tem_commbuf_long_fin), nopass, pointer :: finbuf_long
  procedure(tem_exchange_int), nopass, pointer :: exchange_int
  procedure(tem_commbuf_int_init), nopass, pointer :: initbuf_int
  procedure(tem_commbuf_int_fin), nopass, pointer :: finbuf_int
  procedure(tem_exchange_real), nopass, pointer :: exchange_real
  procedure(tem_commbuf_real_init), nopass, pointer :: initbuf_real
  procedure(tem_commbuf_real_fin), nopass, pointer :: finbuf_real
  end type tem_commPattern_type
  ! ------------------------------------------------------------------------ !

! Definition of the abstract interfaces for the function pointers in the
! TEM_commPattern_type

  abstract interface
    subroutine tem_exchange_long( send, recv, state, message_flag, &
      &                              send_state, comm                 )
      import :: rk, long_k, tem_communication_type
      type(tem_communication_type), intent(inout) :: send, recv
      integer(kind=long_k), intent(inout) :: state(*)
      integer, intent(in) :: message_flag
      integer(kind=long_k), intent(in), optional :: send_state(*)
      !> mpi communicator
      integer, intent(in) :: comm
    end subroutine tem_exchange_long

    subroutine tem_commbuf_long_init(me, pos, nvals)
      import :: tem_longbuffer_type
      type(tem_longbuffer_type), intent(inout) :: me
      integer, intent(in) :: nvals
      integer, intent(in) :: pos(nvals)
    end subroutine tem_commbuf_long_init

    subroutine tem_commbuf_long_fin(me)
      import :: tem_longbuffer_type
      type(tem_longbuffer_type), intent(inout) :: me
    end subroutine tem_commbuf_long_fin
  end interface

  abstract interface
    subroutine tem_exchange_int( send, recv, state, message_flag, &
      &                              send_state, comm                 )
      import :: rk, long_k, tem_communication_type
      type(tem_communication_type), intent(inout) :: send, recv
      integer, intent(inout) :: state(*)
      integer, intent(in) :: message_flag
      integer, intent(in), optional :: send_state(*)
      !> mpi communicator
      integer, intent(in) :: comm
    end subroutine tem_exchange_int

    subroutine tem_commbuf_int_init(me, pos, nvals)
      import :: tem_intbuffer_type
      type(tem_intbuffer_type), intent(inout) :: me
      integer, intent(in) :: nvals
      integer, intent(in) :: pos(nvals)
    end subroutine tem_commbuf_int_init

    subroutine tem_commbuf_int_fin(me)
      import :: tem_intbuffer_type
      type(tem_intbuffer_type), intent(inout) :: me
    end subroutine tem_commbuf_int_fin
  end interface

  abstract interface
    subroutine tem_exchange_real( send, recv, state, message_flag, &
      &                              send_state, comm                 )
      import :: rk, long_k, tem_communication_type
      type(tem_communication_type), intent(inout) :: send, recv
      real(kind=rk), intent(inout) :: state(*)
      integer, intent(in) :: message_flag
      real(kind=rk), intent(in), optional :: send_state(*)
      !> mpi communicator
      integer, intent(in) :: comm
    end subroutine tem_exchange_real

    subroutine tem_commbuf_real_init(me, pos, nvals)
      import :: tem_realbuffer_type
      type(tem_realbuffer_type), intent(inout) :: me
      integer, intent(in) :: nvals
      integer, intent(in) :: pos(nvals)
    end subroutine tem_commbuf_real_init

    subroutine tem_commbuf_real_fin(me)
      import :: tem_realbuffer_type
      type(tem_realbuffer_type), intent(inout) :: me
    end subroutine tem_commbuf_real_fin
  end interface



contains


  ! ************************************************************************ !
  !> This subroutine loads the communication pattern from a Lua script
  !! and sets the exchange routine to be used accordingly.
  !!
  !! The variable read from the script is "commpattern".
  !! Several patterns are available:
  !! * isend_irecv - Use explicit buffers, copy first the outgoing ones, then
  !!      post irecvs and isends, wait on all, and copy the incoming buffers to
  !!      their final locations (default).
  !! * isend_irecv_mpimem - Same as isend_irecv, but use memory that is
  !!      allocatd by MPI_Alloc_mem for the buffers.
  !! * isend_irecv_overlap - Similar to isend_irecv, but directly post sends,
  !!      after filling the outgoing buffer to each process, wait on any
  !!      incoming messages and only wait on sends, after everything is copied
  !!      into the final location.
  !! * overlap_mpimem - Same as isend_irecv_overlap, but use memory that is
  !!      allocated by MPI_Alloc_mem for the buffers.
  !! * typed_isend_irecv - Instead of copying the memory around, define a
  !!      indexed MPI datatype and use that in the exchange.
  !! * gathered_type - Similar to typed_isend_irecv, but with minimal number of
  !!      blocks in the indexed type, by tracking only contiguous blocks in the
  !!      memory layout.
  !!
  !! Instead of reading the style from a configuration script, it can also be
  !! directly set by the caller. If nothing is specified by the caller, the
  !! style will default to isend_irecv.
  !!
  !! Usage:
  !! ~~~~~~~~~~~~~~~~~~~~~{.lua}
  !! commpattern = 'isend_irecv'
  !! ~~~~~~~~~~~~~~~~~~~~~
  !!
  subroutine tem_load_commPattern( me, conf, style )
    ! -------------------------------------------------------------------- !
    !> commpattern to set
    type(tem_commPattern_type), intent(out) :: me
    !> handle to the Lua script
    type(flu_State), optional :: conf
    !> optional communication style
    character(len=*), intent(in), optional :: style
    ! -------------------------------------------------------------------- !
    integer :: iError
    ! -------------------------------------------------------------------- !

    write(logUnit(1),*)"Loading communication pattern:"

    if (present(conf)) then
      ! If a configuration is given, this trumps any other setting.
      ! Defaults to isend_irecv.
      call aot_get_val( L       = conf,          &
        &               key     = 'commpattern', &
        &               val     = me%style,      &
        &               ErrCode = iError,        &
        &               default = 'isend_irecv'  )
    else if (present(style)) then
      ! If a style is given directly by the caller, use that one.
      me%style = style
    else
      ! Default to isend_irecv if nothing provided by the caller.
      me%style = 'isend_irecv'
    end if


    select case(trim(me%style))
    case ('isend_irecv')
      me%exchange_long => comm_isend_irecv_long
      me%initbuf_long => tem_commbuf_long_fillpos
      me%finbuf_long => tem_commbuf_long_finpos
      me%exchange_int => comm_isend_irecv_int
      me%initbuf_int => tem_commbuf_int_fillpos
      me%finbuf_int => tem_commbuf_int_finpos
      me%exchange_real => comm_isend_irecv_real
      me%initbuf_real => tem_commbuf_real_fillpos
      me%finbuf_real => tem_commbuf_real_finpos

    case ('isend_irecv_overlap')
      me%exchange_long => comm_isend_irecv_overlap_long
      me%initbuf_long => tem_commbuf_long_fillpos
      me%finbuf_long => tem_commbuf_long_finpos
      me%exchange_int => comm_isend_irecv_overlap_int
      me%initbuf_int => tem_commbuf_int_fillpos
      me%finbuf_int => tem_commbuf_int_finpos
      me%exchange_real => comm_isend_irecv_overlap_real
      me%initbuf_real => tem_commbuf_real_fillpos
      me%finbuf_real => tem_commbuf_real_finpos

    case ('typed_isend_irecv')
      me%exchange_long => comm_typed_isend_irecv_long
      me%initbuf_long => tem_commbuf_long_fillindexed
      me%finbuf_long => tem_commbuf_long_fintyped
      me%exchange_int => comm_typed_isend_irecv_int
      me%initbuf_int => tem_commbuf_int_fillindexed
      me%finbuf_int => tem_commbuf_int_fintyped
      me%exchange_real => comm_typed_isend_irecv_real
      me%initbuf_real => tem_commbuf_real_fillindexed
      me%finbuf_real => tem_commbuf_real_fintyped

    case ('gathered_type')
      me%exchange_long => comm_typed_isend_irecv_long
      me%initbuf_long => tem_commbuf_long_gatherindexed
      me%finbuf_long => tem_commbuf_long_fintyped
      me%exchange_int => comm_typed_isend_irecv_int
      me%initbuf_int => tem_commbuf_int_gatherindexed
      me%finbuf_int => tem_commbuf_int_fintyped
      me%exchange_real => comm_typed_isend_irecv_real
      me%initbuf_real => tem_commbuf_real_gatherindexed
      me%finbuf_real => tem_commbuf_real_fintyped

    case ('isend_irecv_mpimem')
      me%exchange_long => comm_isend_irecv_long
      me%initbuf_long => tem_commbuf_long_fillmpimem
      me%finbuf_long => tem_commbuf_long_finmpimem
      me%exchange_int => comm_isend_irecv_int
      me%initbuf_int => tem_commbuf_int_fillmpimem
      me%finbuf_int => tem_commbuf_int_finmpimem
      me%exchange_real => comm_isend_irecv_real
      me%initbuf_real => tem_commbuf_real_fillmpimem
      me%finbuf_real => tem_commbuf_real_finmpimem

    case ('overlap_mpimem')
      me%exchange_long => comm_isend_irecv_overlap_long
      me%initbuf_long => tem_commbuf_long_fillmpimem
      me%finbuf_long => tem_commbuf_long_finmpimem
      me%exchange_int => comm_isend_irecv_overlap_int
      me%initbuf_int => tem_commbuf_int_fillmpimem
      me%finbuf_int => tem_commbuf_int_finmpimem
      me%exchange_real => comm_isend_irecv_overlap_real
      me%initbuf_real => tem_commbuf_real_fillmpimem
      me%finbuf_real => tem_commbuf_real_finmpimem

    case default
      write(logUnit(1),*) "ERROR, unknown commpattern: "//trim(me%style)
      write(logUnit(1),*) "available are: "
      write(logUnit(1),*) "* isend_irecv"
      write(logUnit(1),*) "* isend_irecv_overlap"
      write(logUnit(1),*) "* typed_isend_irecv"
      write(logUnit(1),*) "* gathered_type"
      write(logUnit(1),*) "* isend_irecv_mpimem"
      write(logUnit(1),*) "* overlap_mpimem"
      call tem_abort()

    end select

    write(logUnit(1),*) trim(me%style)

  end subroutine tem_load_commPattern
  ! ************************************************************************ !

! The following template defines the routines for various communication
! patterns, if you want to add a new pattern, you have to define an init,
! finalize and exchange routine in this template.
! To add a new data type, just add the appropriate CoCo copy line after this
! template.

  ! ************************************************************************ !
  !> fill the positions that describe how the data in the state
  !! vector relates to the entries in the buffer.
  !!
  subroutine tem_commbuf_long_fillpos( me, pos, nvals )
    ! -------------------------------------------------------------------- !
    type(tem_longbuffer_type), intent(inout) :: me
    integer, intent(in) :: nvals
    integer, intent(in) :: pos(nvals)
    ! -------------------------------------------------------------------- !

    me%nvals = nvals

    if ( allocated(me%pos) ) deallocate(me%pos)
    allocate(me%pos(nvals))
    me%pos = pos

    if ( associated(me%val) ) deallocate(me%val)
    allocate(me%val(nvals))

  end subroutine tem_commbuf_long_fillpos
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> free the communication buffer allocated by the fillpos routine again.
  !!
  subroutine tem_commbuf_long_finpos(me)
    ! -------------------------------------------------------------------- !
    type(tem_longbuffer_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !

    me%nvals = 0
    if ( allocated(me%pos) ) deallocate(me%pos)
    if ( associated(me%val) ) deallocate(me%val)

  end subroutine tem_commbuf_long_finpos
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> fill the positions that describe how the data in the state
  !! vector relates to the entries in the buffer use memory that
  !! is allocated by mpi for the buffer.
  !!
  subroutine tem_commbuf_long_fillmpimem( me, pos, nvals )
    ! -------------------------------------------------------------------- !
    type(tem_longbuffer_type), intent(inout) :: me
    integer, intent(in) :: nvals
    integer, intent(in) :: pos(nvals)
    ! -------------------------------------------------------------------- !
    integer(kind=long_k) :: typesample
    integer(kind=mpi_address_kind) :: typelen
    ! -------------------------------------------------------------------- !

    typelen = int(c_sizeof(typesample), kind=mpi_address_kind)

    me%nvals = nvals

    if ( allocated(me%pos) ) deallocate(me%pos)
    allocate(me%pos(nvals))
    me%pos = pos

    if ( associated(me%val) ) then
      nullify(me%val)
      call free_mpif_mem(me%mem_mpi)
    end if
    call alloc_mpif_mem( asize   = nvals*typelen, &
      &                  baseptr = me%mem_mpi     )
    call c_f_pointer(me%mem_mpi, me%val, [me%nvals])

  end subroutine tem_commbuf_long_fillmpimem
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> free the communication buffer allocated by the fillmpimem routine again.
  !!
  subroutine tem_commbuf_long_finmpimem(me)
    ! -------------------------------------------------------------------- !
    type(tem_longbuffer_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !

    me%nvals = 0
    if ( allocated(me%pos) ) deallocate(me%pos)
    if ( associated(me%val) ) then
      nullify(me%val)
      call free_mpif_mem(me%mem_mpi)
    end if

  end subroutine tem_commbuf_long_finmpimem
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> fill the indexed mpi datatype, which describes how the data in the state
  !! vector relates to the entries in the buffer.
  subroutine tem_commbuf_long_fillindexed(me, pos, nvals)
    ! -------------------------------------------------------------------- !
    type(tem_longbuffer_type), intent(inout) :: me
    integer, intent(in) :: nvals
    integer, intent(in) :: pos(nvals)
    ! -------------------------------------------------------------------- !
    integer :: ierror
    ! -------------------------------------------------------------------- !

    me%nvals = nvals
    ! call mpi_type_create_indexed_block(count, blocklength, &
    !   &                   array_of_displacements, &
    !   &                   oldtype, newtype, ierror)
    call mpi_type_create_indexed_block( nvals, 1, pos - 1,                   &
      &                                 long_k_mpi, me%memindexed, ierror )
    call check_mpi_error(ierror,'create indexed block in tem_commbuf_long_fillindexed')
    call mpi_type_commit(me%memindexed, ierror)
    call check_mpi_error(ierror,'commit memindexed in tem_commbuf_long_fillindexed')

  end subroutine tem_commbuf_long_fillindexed
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> free the communication buffer allocated by the fillindexed routine again.
  !!
  subroutine tem_commbuf_long_fintyped(me)
    ! -------------------------------------------------------------------- !
    type(tem_longbuffer_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !
    integer :: ierror
    ! -------------------------------------------------------------------- !

    me%nvals = 0
    call mpi_type_free(me%memindexed, ierror)
    call check_mpi_error(ierror,'free memindexed in tem_commbuf_long_fintyped')

  end subroutine tem_commbuf_long_fintyped
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> gather the indexed mpi datatype, which describes how the data in the state
  !! vector relates to the entries in the buffer.
  !! in contrast to the simple indexed type above, we try to minimize the number
  !! of blocks here, and gather contiguous blocks of memory together.
  !!
  subroutine tem_commbuf_long_gatherindexed( me, pos, nvals )
    ! -------------------------------------------------------------------- !
    type(tem_longbuffer_type), intent(inout) :: me
    integer, intent(in) :: nvals
    integer, intent(in) :: pos(nvals)
    ! -------------------------------------------------------------------- !
    type(grw_intarray_type) :: blocklength
    type(grw_intarray_type) :: displ
    integer :: ival, counter
    integer :: ierror
    ! -------------------------------------------------------------------- !

    me%nvals = nvals

    ! initialize growing arrays, a kb should be fine to start with...
    call init(blocklength, 256)
    call init(displ, 256)

    if (nvals > 0) then

      ! start with the displacement of the first entry in the list
      call append(displ, pos(1)-1)
      counter = 1

      do ival=2,nvals
        if (pos(ival) == pos(ival-1)+1) then
          ! contiguous memory location following the previous one, increase the
          ! the blocklength.
          counter = counter + 1
        else
          ! new block encountered, record the block found so far
          call append(blocklength, counter)

          ! start new block
          call append(displ, pos(ival)-1)
          counter = 1
        end if
      end do

      ! finish the last block, by recording its found length:
      call append(blocklength, counter)

    end if

    ! call mpi_type_indexed(count, array_of_blocklengths, &
    !   &                   array_of_displacements, oldtype, newtype, ierror)
    call mpi_type_indexed( displ%nvals, blocklength%val, displ%val, &
      &                    long_k_mpi, me%memindexed, ierror     )
    call check_mpi_error(ierror,'type indexed in tem_commbuf_long_gatherindexed')
    call mpi_type_commit(me%memindexed, ierror)
    call check_mpi_error(ierror,'commit memindexed in tem_commbuf_long_gatherindexed')

    call destroy(displ)
    call destroy(blocklength)

  end subroutine tem_commbuf_long_gatherindexed
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> exchange the communication buffers
  !! with a non-blocking mpi communication using preposted irecv and isend
  !! with a waitall
  !!
  subroutine comm_isend_irecv_long( send, recv, state, message_flag, &
    &                                  send_state, comm                 )
    ! -------------------------------------------------------------------- !
    type( tem_communication_type ), intent(inout) :: send
    type( tem_communication_type ), intent(inout) :: recv
    integer(kind=long_k), intent(inout) :: state(*)  !< state vector to update
    integer, intent(in) :: message_flag
    integer(kind=long_k), intent(in), optional :: send_state(*)  !< data to send
    !> mpi communicator
    integer, intent(in) :: comm
    ! -------------------------------------------------------------------- !
    ! request handle for messages
    ! @todo request handle array could exist during complete code runtime
    ! integer :: rq_handle( recv%nprocs + send%nprocs )
    integer :: status( mpi_status_size, max(recv%nprocs, send%nprocs) )
    integer :: ierr             !< error flag
    integer :: iproc, ival
    integer :: nsendvals, nrecvvals
    ! -------------------------------------------------------------------- !

    if (present(send_state)) then
      do iproc = 1, send%nprocs
        ! fill communication message
        nsendvals = send%buf_long( iproc )%nvals
        !$nec ivdep
        do ival = 1, nsendvals
          send%buf_long( iproc )%val( ival )                    &
            &  = send_state( send%buf_long( iproc )%pos( ival ) )
        end do
      end do
    else
      do iproc = 1, send%nprocs
        ! fill communication message
        nsendvals = send%buf_long( iproc )%nvals
        !$nec ivdep
        do ival = 1, nsendvals
          send%buf_long( iproc )%val( ival )               &
            &  = state( send%buf_long( iproc )%pos( ival ) )
        end do
      end do
    end if

    do iproc = 1, recv%nprocs
      ! start receive communications
      call mpi_irecv(                           &
       &      recv%buf_long( iproc )%val,    & ! me
       &      recv%buf_long( iproc )%nvals,  & ! me size
       &      long_k_mpi,                    & ! data type
       &      recv%proc(iproc),                 & ! target me
       &      message_flag,                     & ! flag
       &      comm,                             & ! communicator
       &      recv%rqhandle(iproc),             & ! handle
       &      ierr )                              ! error status
    enddo

    !  start the sending communications
    do iproc = 1, send%nprocs
      call mpi_isend(                           &
       &      send%buf_long( iproc )%val,    & ! buffer
       &      send%buf_long( iproc )%nvals,  & ! count
       &      long_k_mpi,                    & ! data type
       &      send%proc(iproc),                 & ! target
       &      message_flag,                     & ! tag
       &      comm,                             & ! communicator
       &      send%rqhandle( iproc ),           & ! handle
       &      ierr )                              ! error status
    enddo ! iproc

    ! wait for receive buffer to be ready
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall(recv%nprocs,               & ! count
        &              recv%rqhandle,             & ! request handles
        &              status,                    & ! mpi status
        &              ierr )                       ! error status
    end if

    ! now values from recv me can be copied to the actual state array
    do iproc = 1, recv%nprocs
      nrecvvals = recv%buf_long( iproc )%nvals
      !$nec ivdep
      do ival = 1, nrecvvals
        ! write the values from the recv me to the state array
        ! to positions given in recvme%pos
        state( recv%buf_long( iproc )%pos( ival ) ) &
          &  = recv%buf_long( iproc )%val( ival )
      end do
    end do

    ! wait for send buffer to be ready
    if ( send%nprocs /= 0 ) then
      call mpi_waitall(send%nprocs,   & ! count
        &              send%rqhandle, & ! request handles
        &              status,        & ! mpi status
        &              ierr )           ! error status
    end if

  end subroutine comm_isend_irecv_long
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine comm_isend_irecv_overlap_long( send, recv, state,             &
    &                                          message_flag, send_state, comm )
    ! -------------------------------------------------------------------- !
    type( tem_communication_type), intent(inout) :: send, recv
    integer(kind=long_k), intent(inout) :: state(*)
    integer, intent(in) :: message_flag
    integer(kind=long_k), intent(in), optional :: send_state(*)  !< data to send
    !> mpi communicator
    integer, intent(in) :: comm
    ! -------------------------------------------------------------------- !
    integer :: iproc, ival  ! counter for neigbor mees
    integer :: finproc
    ! integer :: rq_handle( recv%nprocs + send%nprocs )
    integer :: status( mpi_status_size, max(send%nprocs,1) )
    integer :: ierr ! error flag
    integer :: nsendvals, nrecvvals
    ! -------------------------------------------------------------------- !

    ! pre-post irecv
    do iproc = 1, recv%nprocs
      call mpi_irecv(                          &
       &      recv%buf_long( iproc )%val,   & ! buffer
       &      recv%buf_long( iproc )%nvals, & ! count
       &      long_k_mpi,                   & ! data type
       &      recv%proc(iproc),                & ! target me
       &      message_flag,                    & ! tag
       &      comm,                            & ! communicator
       &      recv%rqhandle(iproc),            & ! handle
       &      ierr                             ) ! error status
    end do

    !> fill send buffers and start sending
    do iproc = 1, send%nprocs
      nsendvals = send%buf_long( iproc )%nvals
      if (present(send_state)) then
        !$nec ivdep
        do ival = 1, nsendvals
          send%buf_long( iproc )%val( ival )                    &
            &  = send_state( send%buf_long( iproc )%pos( ival ) )
        end do
      else
        !$nec ivdep
        do ival = 1, nsendvals
          send%buf_long( iproc )%val( ival )               &
            &  = state( send%buf_long( iproc )%pos( ival ) )
        end do
      end if
      call mpi_isend(                               &
       &      send%buf_long( iproc )%val,        & ! buffer
       &      send%buf_long( iproc )%nvals,      & ! count
       &      long_k_mpi,                        & ! data type
       &      send%proc(iproc),                     & ! target
       &      message_flag,                         & ! tag
       &      comm,                                 & ! comm
       &      send%rqhandle( iproc ),               & ! handle
       &      ierr )                                  ! error status
    end do

    do iproc = 1, recv%nprocs
      ! wait for any of receive buffer to be ready
      call mpi_waitany(   &
        &    recv%nprocs, & ! count
        &    recv%rqhandle,   & ! request handles
        &    finproc,     & ! process that finished
        &    status(:,1), & ! mpi status
        &    ierr )         ! error status
      nrecvvals = recv%buf_long( finproc )%nvals
      !$nec ivdep
      do ival = 1, nrecvvals
        ! write the values from the recv me to the state array
        ! to positions given in recvme%pos
        state( recv%buf_long( finproc )%pos( ival ) ) &
          &  = recv%buf_long( finproc )%val( ival )
      end do
    end do

    if (send%nprocs > 0) then
      ! wait for send buffer to be ready
      call mpi_waitall(                  &
        &    send%nprocs,                & ! total number of comm.'s to wait for
        &    send%rqhandle,              & ! request handles
        &    status,                     & ! mpi status
        &    ierr )                        ! error status
    end if

  end subroutine comm_isend_irecv_overlap_long
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> exchange the communication mes
  !! with a non-blocking mpi communication using preposted irecv and isend
  !! with a waitall
  !!
  subroutine comm_typed_isend_irecv_long( send, recv, state,             &
    &                                        message_flag, send_state, comm )
    ! -------------------------------------------------------------------- !
    type( tem_communication_type ), intent(inout) :: send, recv
    integer(kind=long_k), intent(inout) :: state(*)  !< current state vector
    integer, intent(in) :: message_flag
    integer(kind=long_k), intent(in), optional :: send_state(*)  !< data to send
    !> mpi communicator
    integer, intent(in) :: comm
    ! -------------------------------------------------------------------- !
    ! @todo request handle array could exist during complete code runtime
    ! integer :: rq_handle( recv%nprocs + send%nprocs )
    integer :: status( mpi_status_size, max(recv%nprocs, send%nprocs) )
    integer :: ierr             ! error flag
    integer :: iproc
    ! -------------------------------------------------------------------- !

    !> values for send me must have been copied from the actual state array
    do iproc = 1, recv%nprocs
      !> start receive communications
      call mpi_irecv(                                     &
       &              state,                              & ! buffer
       &              1,                                  & ! count
       &              recv%buf_long(iproc)%memindexed, & ! type
       &              recv%proc(iproc),                   & ! source
       &              message_flag,                       & ! tag
       &              comm,                               & ! comm
       &              recv%rqhandle(iproc),               & ! request handle
       &              ierr )

    end do

    !>  start the sending communications
    if (present(send_state)) then
      do iproc = 1, send%nprocs
        call mpi_isend(                                     &
          &             send_state,                         & ! buffer
          &             1,                                  & ! count
          &             send%buf_long(iproc)%memindexed, & ! type
          &             send%proc(iproc),                   & ! target
          &             message_flag,                       & ! tag
          &             comm,                               & ! comm
          &             send%rqhandle( iproc ),             & ! handle
          &             ierr )
      end do !< iproc
    else
      do iproc = 1, send%nprocs
        call mpi_isend(                                     &
          &             state,                              & ! buffer
          &             1,                                  & ! count
          &             send%buf_long(iproc)%memindexed, & ! type
          &             send%proc(iproc),                   & ! target
          &             message_flag,                       & ! tag
          &             comm,                               & ! comm
          &             send%rqhandle( iproc ),             & ! handle
          &             ierr )
      end do !< iproc
    end if

    !> wait for above communications to complete
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall( recv%nprocs,     & ! count
        &               recv%rqhandle,   & ! request handles
        &               status,          & ! statuses
        &               ierr )
    end if
    if ( send%nprocs /= 0 ) then
      call mpi_waitall( send%nprocs,     & ! count
        &               send%rqhandle,   & ! request handles
        &               status,          & ! statuses
        &               ierr )
    end if

  end subroutine comm_typed_isend_irecv_long
  ! ************************************************************************ !
  ! ************************************************************************ !
  !> fill the positions that describe how the data in the state
  !! vector relates to the entries in the buffer.
  !!
  subroutine tem_commbuf_int_fillpos( me, pos, nvals )
    ! -------------------------------------------------------------------- !
    type(tem_intbuffer_type), intent(inout) :: me
    integer, intent(in) :: nvals
    integer, intent(in) :: pos(nvals)
    ! -------------------------------------------------------------------- !

    me%nvals = nvals

    if ( allocated(me%pos) ) deallocate(me%pos)
    allocate(me%pos(nvals))
    me%pos = pos

    if ( associated(me%val) ) deallocate(me%val)
    allocate(me%val(nvals))

  end subroutine tem_commbuf_int_fillpos
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> free the communication buffer allocated by the fillpos routine again.
  !!
  subroutine tem_commbuf_int_finpos(me)
    ! -------------------------------------------------------------------- !
    type(tem_intbuffer_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !

    me%nvals = 0
    if ( allocated(me%pos) ) deallocate(me%pos)
    if ( associated(me%val) ) deallocate(me%val)

  end subroutine tem_commbuf_int_finpos
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> fill the positions that describe how the data in the state
  !! vector relates to the entries in the buffer use memory that
  !! is allocated by mpi for the buffer.
  !!
  subroutine tem_commbuf_int_fillmpimem( me, pos, nvals )
    ! -------------------------------------------------------------------- !
    type(tem_intbuffer_type), intent(inout) :: me
    integer, intent(in) :: nvals
    integer, intent(in) :: pos(nvals)
    ! -------------------------------------------------------------------- !
    integer :: typesample
    integer(kind=mpi_address_kind) :: typelen
    ! -------------------------------------------------------------------- !

    typelen = int(c_sizeof(typesample), kind=mpi_address_kind)

    me%nvals = nvals

    if ( allocated(me%pos) ) deallocate(me%pos)
    allocate(me%pos(nvals))
    me%pos = pos

    if ( associated(me%val) ) then
      nullify(me%val)
      call free_mpif_mem(me%mem_mpi)
    end if
    call alloc_mpif_mem( asize   = nvals*typelen, &
      &                  baseptr = me%mem_mpi     )
    call c_f_pointer(me%mem_mpi, me%val, [me%nvals])

  end subroutine tem_commbuf_int_fillmpimem
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> free the communication buffer allocated by the fillmpimem routine again.
  !!
  subroutine tem_commbuf_int_finmpimem(me)
    ! -------------------------------------------------------------------- !
    type(tem_intbuffer_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !

    me%nvals = 0
    if ( allocated(me%pos) ) deallocate(me%pos)
    if ( associated(me%val) ) then
      nullify(me%val)
      call free_mpif_mem(me%mem_mpi)
    end if

  end subroutine tem_commbuf_int_finmpimem
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> fill the indexed mpi datatype, which describes how the data in the state
  !! vector relates to the entries in the buffer.
  subroutine tem_commbuf_int_fillindexed(me, pos, nvals)
    ! -------------------------------------------------------------------- !
    type(tem_intbuffer_type), intent(inout) :: me
    integer, intent(in) :: nvals
    integer, intent(in) :: pos(nvals)
    ! -------------------------------------------------------------------- !
    integer :: ierror
    ! -------------------------------------------------------------------- !

    me%nvals = nvals
    ! call mpi_type_create_indexed_block(count, blocklength, &
    !   &                   array_of_displacements, &
    !   &                   oldtype, newtype, ierror)
    call mpi_type_create_indexed_block( nvals, 1, pos - 1,                   &
      &                                 mpi_integer, me%memindexed, ierror )
    call check_mpi_error(ierror,'create indexed block in tem_commbuf_int_fillindexed')
    call mpi_type_commit(me%memindexed, ierror)
    call check_mpi_error(ierror,'commit memindexed in tem_commbuf_int_fillindexed')

  end subroutine tem_commbuf_int_fillindexed
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> free the communication buffer allocated by the fillindexed routine again.
  !!
  subroutine tem_commbuf_int_fintyped(me)
    ! -------------------------------------------------------------------- !
    type(tem_intbuffer_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !
    integer :: ierror
    ! -------------------------------------------------------------------- !

    me%nvals = 0
    call mpi_type_free(me%memindexed, ierror)
    call check_mpi_error(ierror,'free memindexed in tem_commbuf_int_fintyped')

  end subroutine tem_commbuf_int_fintyped
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> gather the indexed mpi datatype, which describes how the data in the state
  !! vector relates to the entries in the buffer.
  !! in contrast to the simple indexed type above, we try to minimize the number
  !! of blocks here, and gather contiguous blocks of memory together.
  !!
  subroutine tem_commbuf_int_gatherindexed( me, pos, nvals )
    ! -------------------------------------------------------------------- !
    type(tem_intbuffer_type), intent(inout) :: me
    integer, intent(in) :: nvals
    integer, intent(in) :: pos(nvals)
    ! -------------------------------------------------------------------- !
    type(grw_intarray_type) :: blocklength
    type(grw_intarray_type) :: displ
    integer :: ival, counter
    integer :: ierror
    ! -------------------------------------------------------------------- !

    me%nvals = nvals

    ! initialize growing arrays, a kb should be fine to start with...
    call init(blocklength, 256)
    call init(displ, 256)

    if (nvals > 0) then

      ! start with the displacement of the first entry in the list
      call append(displ, pos(1)-1)
      counter = 1

      do ival=2,nvals
        if (pos(ival) == pos(ival-1)+1) then
          ! contiguous memory location following the previous one, increase the
          ! the blocklength.
          counter = counter + 1
        else
          ! new block encountered, record the block found so far
          call append(blocklength, counter)

          ! start new block
          call append(displ, pos(ival)-1)
          counter = 1
        end if
      end do

      ! finish the last block, by recording its found length:
      call append(blocklength, counter)

    end if

    ! call mpi_type_indexed(count, array_of_blocklengths, &
    !   &                   array_of_displacements, oldtype, newtype, ierror)
    call mpi_type_indexed( displ%nvals, blocklength%val, displ%val, &
      &                    mpi_integer, me%memindexed, ierror     )
    call check_mpi_error(ierror,'type indexed in tem_commbuf_int_gatherindexed')
    call mpi_type_commit(me%memindexed, ierror)
    call check_mpi_error(ierror,'commit memindexed in tem_commbuf_int_gatherindexed')

    call destroy(displ)
    call destroy(blocklength)

  end subroutine tem_commbuf_int_gatherindexed
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> exchange the communication buffers
  !! with a non-blocking mpi communication using preposted irecv and isend
  !! with a waitall
  !!
  subroutine comm_isend_irecv_int( send, recv, state, message_flag, &
    &                                  send_state, comm                 )
    ! -------------------------------------------------------------------- !
    type( tem_communication_type ), intent(inout) :: send
    type( tem_communication_type ), intent(inout) :: recv
    integer, intent(inout) :: state(*)  !< state vector to update
    integer, intent(in) :: message_flag
    integer, intent(in), optional :: send_state(*)  !< data to send
    !> mpi communicator
    integer, intent(in) :: comm
    ! -------------------------------------------------------------------- !
    ! request handle for messages
    ! @todo request handle array could exist during complete code runtime
    ! integer :: rq_handle( recv%nprocs + send%nprocs )
    integer :: status( mpi_status_size, max(recv%nprocs, send%nprocs) )
    integer :: ierr             !< error flag
    integer :: iproc, ival
    integer :: nsendvals, nrecvvals
    ! -------------------------------------------------------------------- !

    if (present(send_state)) then
      do iproc = 1, send%nprocs
        ! fill communication message
        nsendvals = send%buf_int( iproc )%nvals
        !$nec ivdep
        do ival = 1, nsendvals
          send%buf_int( iproc )%val( ival )                    &
            &  = send_state( send%buf_int( iproc )%pos( ival ) )
        end do
      end do
    else
      do iproc = 1, send%nprocs
        ! fill communication message
        nsendvals = send%buf_int( iproc )%nvals
        !$nec ivdep
        do ival = 1, nsendvals
          send%buf_int( iproc )%val( ival )               &
            &  = state( send%buf_int( iproc )%pos( ival ) )
        end do
      end do
    end if

    do iproc = 1, recv%nprocs
      ! start receive communications
      call mpi_irecv(                           &
       &      recv%buf_int( iproc )%val,    & ! me
       &      recv%buf_int( iproc )%nvals,  & ! me size
       &      mpi_integer,                    & ! data type
       &      recv%proc(iproc),                 & ! target me
       &      message_flag,                     & ! flag
       &      comm,                             & ! communicator
       &      recv%rqhandle(iproc),             & ! handle
       &      ierr )                              ! error status
    enddo

    !  start the sending communications
    do iproc = 1, send%nprocs
      call mpi_isend(                           &
       &      send%buf_int( iproc )%val,    & ! buffer
       &      send%buf_int( iproc )%nvals,  & ! count
       &      mpi_integer,                    & ! data type
       &      send%proc(iproc),                 & ! target
       &      message_flag,                     & ! tag
       &      comm,                             & ! communicator
       &      send%rqhandle( iproc ),           & ! handle
       &      ierr )                              ! error status
    enddo ! iproc

    ! wait for receive buffer to be ready
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall(recv%nprocs,               & ! count
        &              recv%rqhandle,             & ! request handles
        &              status,                    & ! mpi status
        &              ierr )                       ! error status
    end if

    ! now values from recv me can be copied to the actual state array
    do iproc = 1, recv%nprocs
      nrecvvals = recv%buf_int( iproc )%nvals
      !$nec ivdep
      do ival = 1, nrecvvals
        ! write the values from the recv me to the state array
        ! to positions given in recvme%pos
        state( recv%buf_int( iproc )%pos( ival ) ) &
          &  = recv%buf_int( iproc )%val( ival )
      end do
    end do

    ! wait for send buffer to be ready
    if ( send%nprocs /= 0 ) then
      call mpi_waitall(send%nprocs,   & ! count
        &              send%rqhandle, & ! request handles
        &              status,        & ! mpi status
        &              ierr )           ! error status
    end if

  end subroutine comm_isend_irecv_int
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine comm_isend_irecv_overlap_int( send, recv, state,             &
    &                                          message_flag, send_state, comm )
    ! -------------------------------------------------------------------- !
    type( tem_communication_type), intent(inout) :: send, recv
    integer, intent(inout) :: state(*)
    integer, intent(in) :: message_flag
    integer, intent(in), optional :: send_state(*)  !< data to send
    !> mpi communicator
    integer, intent(in) :: comm
    ! -------------------------------------------------------------------- !
    integer :: iproc, ival  ! counter for neigbor mees
    integer :: finproc
    ! integer :: rq_handle( recv%nprocs + send%nprocs )
    integer :: status( mpi_status_size, max(send%nprocs,1) )
    integer :: ierr ! error flag
    integer :: nsendvals, nrecvvals
    ! -------------------------------------------------------------------- !

    ! pre-post irecv
    do iproc = 1, recv%nprocs
      call mpi_irecv(                          &
       &      recv%buf_int( iproc )%val,   & ! buffer
       &      recv%buf_int( iproc )%nvals, & ! count
       &      mpi_integer,                   & ! data type
       &      recv%proc(iproc),                & ! target me
       &      message_flag,                    & ! tag
       &      comm,                            & ! communicator
       &      recv%rqhandle(iproc),            & ! handle
       &      ierr                             ) ! error status
    end do

    !> fill send buffers and start sending
    do iproc = 1, send%nprocs
      nsendvals = send%buf_int( iproc )%nvals
      if (present(send_state)) then
        !$nec ivdep
        do ival = 1, nsendvals
          send%buf_int( iproc )%val( ival )                    &
            &  = send_state( send%buf_int( iproc )%pos( ival ) )
        end do
      else
        !$nec ivdep
        do ival = 1, nsendvals
          send%buf_int( iproc )%val( ival )               &
            &  = state( send%buf_int( iproc )%pos( ival ) )
        end do
      end if
      call mpi_isend(                               &
       &      send%buf_int( iproc )%val,        & ! buffer
       &      send%buf_int( iproc )%nvals,      & ! count
       &      mpi_integer,                        & ! data type
       &      send%proc(iproc),                     & ! target
       &      message_flag,                         & ! tag
       &      comm,                                 & ! comm
       &      send%rqhandle( iproc ),               & ! handle
       &      ierr )                                  ! error status
    end do

    do iproc = 1, recv%nprocs
      ! wait for any of receive buffer to be ready
      call mpi_waitany(   &
        &    recv%nprocs, & ! count
        &    recv%rqhandle,   & ! request handles
        &    finproc,     & ! process that finished
        &    status(:,1), & ! mpi status
        &    ierr )         ! error status
      nrecvvals = recv%buf_int( finproc )%nvals
      !$nec ivdep
      do ival = 1, nrecvvals
        ! write the values from the recv me to the state array
        ! to positions given in recvme%pos
        state( recv%buf_int( finproc )%pos( ival ) ) &
          &  = recv%buf_int( finproc )%val( ival )
      end do
    end do

    if (send%nprocs > 0) then
      ! wait for send buffer to be ready
      call mpi_waitall(                  &
        &    send%nprocs,                & ! total number of comm.'s to wait for
        &    send%rqhandle,              & ! request handles
        &    status,                     & ! mpi status
        &    ierr )                        ! error status
    end if

  end subroutine comm_isend_irecv_overlap_int
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> exchange the communication mes
  !! with a non-blocking mpi communication using preposted irecv and isend
  !! with a waitall
  !!
  subroutine comm_typed_isend_irecv_int( send, recv, state,             &
    &                                        message_flag, send_state, comm )
    ! -------------------------------------------------------------------- !
    type( tem_communication_type ), intent(inout) :: send, recv
    integer, intent(inout) :: state(*)  !< current state vector
    integer, intent(in) :: message_flag
    integer, intent(in), optional :: send_state(*)  !< data to send
    !> mpi communicator
    integer, intent(in) :: comm
    ! -------------------------------------------------------------------- !
    ! @todo request handle array could exist during complete code runtime
    ! integer :: rq_handle( recv%nprocs + send%nprocs )
    integer :: status( mpi_status_size, max(recv%nprocs, send%nprocs) )
    integer :: ierr             ! error flag
    integer :: iproc
    ! -------------------------------------------------------------------- !

    !> values for send me must have been copied from the actual state array
    do iproc = 1, recv%nprocs
      !> start receive communications
      call mpi_irecv(                                     &
       &              state,                              & ! buffer
       &              1,                                  & ! count
       &              recv%buf_int(iproc)%memindexed, & ! type
       &              recv%proc(iproc),                   & ! source
       &              message_flag,                       & ! tag
       &              comm,                               & ! comm
       &              recv%rqhandle(iproc),               & ! request handle
       &              ierr )

    end do

    !>  start the sending communications
    if (present(send_state)) then
      do iproc = 1, send%nprocs
        call mpi_isend(                                     &
          &             send_state,                         & ! buffer
          &             1,                                  & ! count
          &             send%buf_int(iproc)%memindexed, & ! type
          &             send%proc(iproc),                   & ! target
          &             message_flag,                       & ! tag
          &             comm,                               & ! comm
          &             send%rqhandle( iproc ),             & ! handle
          &             ierr )
      end do !< iproc
    else
      do iproc = 1, send%nprocs
        call mpi_isend(                                     &
          &             state,                              & ! buffer
          &             1,                                  & ! count
          &             send%buf_int(iproc)%memindexed, & ! type
          &             send%proc(iproc),                   & ! target
          &             message_flag,                       & ! tag
          &             comm,                               & ! comm
          &             send%rqhandle( iproc ),             & ! handle
          &             ierr )
      end do !< iproc
    end if

    !> wait for above communications to complete
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall( recv%nprocs,     & ! count
        &               recv%rqhandle,   & ! request handles
        &               status,          & ! statuses
        &               ierr )
    end if
    if ( send%nprocs /= 0 ) then
      call mpi_waitall( send%nprocs,     & ! count
        &               send%rqhandle,   & ! request handles
        &               status,          & ! statuses
        &               ierr )
    end if

  end subroutine comm_typed_isend_irecv_int
  ! ************************************************************************ !
  ! ************************************************************************ !
  !> fill the positions that describe how the data in the state
  !! vector relates to the entries in the buffer.
  !!
  subroutine tem_commbuf_real_fillpos( me, pos, nvals )
    ! -------------------------------------------------------------------- !
    type(tem_realbuffer_type), intent(inout) :: me
    integer, intent(in) :: nvals
    integer, intent(in) :: pos(nvals)
    ! -------------------------------------------------------------------- !

    me%nvals = nvals

    if ( allocated(me%pos) ) deallocate(me%pos)
    allocate(me%pos(nvals))
    me%pos = pos

    if ( associated(me%val) ) deallocate(me%val)
    allocate(me%val(nvals))

  end subroutine tem_commbuf_real_fillpos
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> free the communication buffer allocated by the fillpos routine again.
  !!
  subroutine tem_commbuf_real_finpos(me)
    ! -------------------------------------------------------------------- !
    type(tem_realbuffer_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !

    me%nvals = 0
    if ( allocated(me%pos) ) deallocate(me%pos)
    if ( associated(me%val) ) deallocate(me%val)

  end subroutine tem_commbuf_real_finpos
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> fill the positions that describe how the data in the state
  !! vector relates to the entries in the buffer use memory that
  !! is allocated by mpi for the buffer.
  !!
  subroutine tem_commbuf_real_fillmpimem( me, pos, nvals )
    ! -------------------------------------------------------------------- !
    type(tem_realbuffer_type), intent(inout) :: me
    integer, intent(in) :: nvals
    integer, intent(in) :: pos(nvals)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: typesample
    integer(kind=mpi_address_kind) :: typelen
    ! -------------------------------------------------------------------- !

    typelen = int(c_sizeof(typesample), kind=mpi_address_kind)

    me%nvals = nvals

    if ( allocated(me%pos) ) deallocate(me%pos)
    allocate(me%pos(nvals))
    me%pos = pos

    if ( associated(me%val) ) then
      nullify(me%val)
      call free_mpif_mem(me%mem_mpi)
    end if
    call alloc_mpif_mem( asize   = nvals*typelen, &
      &                  baseptr = me%mem_mpi     )
    call c_f_pointer(me%mem_mpi, me%val, [me%nvals])

  end subroutine tem_commbuf_real_fillmpimem
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> free the communication buffer allocated by the fillmpimem routine again.
  !!
  subroutine tem_commbuf_real_finmpimem(me)
    ! -------------------------------------------------------------------- !
    type(tem_realbuffer_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !

    me%nvals = 0
    if ( allocated(me%pos) ) deallocate(me%pos)
    if ( associated(me%val) ) then
      nullify(me%val)
      call free_mpif_mem(me%mem_mpi)
    end if

  end subroutine tem_commbuf_real_finmpimem
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> fill the indexed mpi datatype, which describes how the data in the state
  !! vector relates to the entries in the buffer.
  subroutine tem_commbuf_real_fillindexed(me, pos, nvals)
    ! -------------------------------------------------------------------- !
    type(tem_realbuffer_type), intent(inout) :: me
    integer, intent(in) :: nvals
    integer, intent(in) :: pos(nvals)
    ! -------------------------------------------------------------------- !
    integer :: ierror
    ! -------------------------------------------------------------------- !

    me%nvals = nvals
    ! call mpi_type_create_indexed_block(count, blocklength, &
    !   &                   array_of_displacements, &
    !   &                   oldtype, newtype, ierror)
    call mpi_type_create_indexed_block( nvals, 1, pos - 1,                   &
      &                                 rk_mpi, me%memindexed, ierror )
    call check_mpi_error(ierror,'create indexed block in tem_commbuf_real_fillindexed')
    call mpi_type_commit(me%memindexed, ierror)
    call check_mpi_error(ierror,'commit memindexed in tem_commbuf_real_fillindexed')

  end subroutine tem_commbuf_real_fillindexed
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> free the communication buffer allocated by the fillindexed routine again.
  !!
  subroutine tem_commbuf_real_fintyped(me)
    ! -------------------------------------------------------------------- !
    type(tem_realbuffer_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !
    integer :: ierror
    ! -------------------------------------------------------------------- !

    me%nvals = 0
    call mpi_type_free(me%memindexed, ierror)
    call check_mpi_error(ierror,'free memindexed in tem_commbuf_real_fintyped')

  end subroutine tem_commbuf_real_fintyped
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> gather the indexed mpi datatype, which describes how the data in the state
  !! vector relates to the entries in the buffer.
  !! in contrast to the simple indexed type above, we try to minimize the number
  !! of blocks here, and gather contiguous blocks of memory together.
  !!
  subroutine tem_commbuf_real_gatherindexed( me, pos, nvals )
    ! -------------------------------------------------------------------- !
    type(tem_realbuffer_type), intent(inout) :: me
    integer, intent(in) :: nvals
    integer, intent(in) :: pos(nvals)
    ! -------------------------------------------------------------------- !
    type(grw_intarray_type) :: blocklength
    type(grw_intarray_type) :: displ
    integer :: ival, counter
    integer :: ierror
    ! -------------------------------------------------------------------- !

    me%nvals = nvals

    ! initialize growing arrays, a kb should be fine to start with...
    call init(blocklength, 256)
    call init(displ, 256)

    if (nvals > 0) then

      ! start with the displacement of the first entry in the list
      call append(displ, pos(1)-1)
      counter = 1

      do ival=2,nvals
        if (pos(ival) == pos(ival-1)+1) then
          ! contiguous memory location following the previous one, increase the
          ! the blocklength.
          counter = counter + 1
        else
          ! new block encountered, record the block found so far
          call append(blocklength, counter)

          ! start new block
          call append(displ, pos(ival)-1)
          counter = 1
        end if
      end do

      ! finish the last block, by recording its found length:
      call append(blocklength, counter)

    end if

    ! call mpi_type_indexed(count, array_of_blocklengths, &
    !   &                   array_of_displacements, oldtype, newtype, ierror)
    call mpi_type_indexed( displ%nvals, blocklength%val, displ%val, &
      &                    rk_mpi, me%memindexed, ierror     )
    call check_mpi_error(ierror,'type indexed in tem_commbuf_real_gatherindexed')
    call mpi_type_commit(me%memindexed, ierror)
    call check_mpi_error(ierror,'commit memindexed in tem_commbuf_real_gatherindexed')

    call destroy(displ)
    call destroy(blocklength)

  end subroutine tem_commbuf_real_gatherindexed
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> exchange the communication buffers
  !! with a non-blocking mpi communication using preposted irecv and isend
  !! with a waitall
  !!
  subroutine comm_isend_irecv_real( send, recv, state, message_flag, &
    &                                  send_state, comm                 )
    ! -------------------------------------------------------------------- !
    type( tem_communication_type ), intent(inout) :: send
    type( tem_communication_type ), intent(inout) :: recv
    real(kind=rk), intent(inout) :: state(*)  !< state vector to update
    integer, intent(in) :: message_flag
    real(kind=rk), intent(in), optional :: send_state(*)  !< data to send
    !> mpi communicator
    integer, intent(in) :: comm
    ! -------------------------------------------------------------------- !
    ! request handle for messages
    ! @todo request handle array could exist during complete code runtime
    ! integer :: rq_handle( recv%nprocs + send%nprocs )
    integer :: status( mpi_status_size, max(recv%nprocs, send%nprocs) )
    integer :: ierr             !< error flag
    integer :: iproc, ival
    integer :: nsendvals, nrecvvals
    ! -------------------------------------------------------------------- !

    if (present(send_state)) then
      do iproc = 1, send%nprocs
        ! fill communication message
        nsendvals = send%buf_real( iproc )%nvals
        !$nec ivdep
        do ival = 1, nsendvals
          send%buf_real( iproc )%val( ival )                    &
            &  = send_state( send%buf_real( iproc )%pos( ival ) )
        end do
      end do
    else
      do iproc = 1, send%nprocs
        ! fill communication message
        nsendvals = send%buf_real( iproc )%nvals
        !$nec ivdep
        do ival = 1, nsendvals
          send%buf_real( iproc )%val( ival )               &
            &  = state( send%buf_real( iproc )%pos( ival ) )
        end do
      end do
    end if

    do iproc = 1, recv%nprocs
      ! start receive communications
      call mpi_irecv(                           &
       &      recv%buf_real( iproc )%val,    & ! me
       &      recv%buf_real( iproc )%nvals,  & ! me size
       &      rk_mpi,                    & ! data type
       &      recv%proc(iproc),                 & ! target me
       &      message_flag,                     & ! flag
       &      comm,                             & ! communicator
       &      recv%rqhandle(iproc),             & ! handle
       &      ierr )                              ! error status
    enddo

    !  start the sending communications
    do iproc = 1, send%nprocs
      call mpi_isend(                           &
       &      send%buf_real( iproc )%val,    & ! buffer
       &      send%buf_real( iproc )%nvals,  & ! count
       &      rk_mpi,                    & ! data type
       &      send%proc(iproc),                 & ! target
       &      message_flag,                     & ! tag
       &      comm,                             & ! communicator
       &      send%rqhandle( iproc ),           & ! handle
       &      ierr )                              ! error status
    enddo ! iproc

    ! wait for receive buffer to be ready
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall(recv%nprocs,               & ! count
        &              recv%rqhandle,             & ! request handles
        &              status,                    & ! mpi status
        &              ierr )                       ! error status
    end if

    ! now values from recv me can be copied to the actual state array
    do iproc = 1, recv%nprocs
      nrecvvals = recv%buf_real( iproc )%nvals
      !$nec ivdep
      do ival = 1, nrecvvals
        ! write the values from the recv me to the state array
        ! to positions given in recvme%pos
        state( recv%buf_real( iproc )%pos( ival ) ) &
          &  = recv%buf_real( iproc )%val( ival )
      end do
    end do

    ! wait for send buffer to be ready
    if ( send%nprocs /= 0 ) then
      call mpi_waitall(send%nprocs,   & ! count
        &              send%rqhandle, & ! request handles
        &              status,        & ! mpi status
        &              ierr )           ! error status
    end if

  end subroutine comm_isend_irecv_real
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine comm_isend_irecv_overlap_real( send, recv, state,             &
    &                                          message_flag, send_state, comm )
    ! -------------------------------------------------------------------- !
    type( tem_communication_type), intent(inout) :: send, recv
    real(kind=rk), intent(inout) :: state(*)
    integer, intent(in) :: message_flag
    real(kind=rk), intent(in), optional :: send_state(*)  !< data to send
    !> mpi communicator
    integer, intent(in) :: comm
    ! -------------------------------------------------------------------- !
    integer :: iproc, ival  ! counter for neigbor mees
    integer :: finproc
    ! integer :: rq_handle( recv%nprocs + send%nprocs )
    integer :: status( mpi_status_size, max(send%nprocs,1) )
    integer :: ierr ! error flag
    integer :: nsendvals, nrecvvals
    ! -------------------------------------------------------------------- !

    ! pre-post irecv
    do iproc = 1, recv%nprocs
      call mpi_irecv(                          &
       &      recv%buf_real( iproc )%val,   & ! buffer
       &      recv%buf_real( iproc )%nvals, & ! count
       &      rk_mpi,                   & ! data type
       &      recv%proc(iproc),                & ! target me
       &      message_flag,                    & ! tag
       &      comm,                            & ! communicator
       &      recv%rqhandle(iproc),            & ! handle
       &      ierr                             ) ! error status
    end do

    !> fill send buffers and start sending
    do iproc = 1, send%nprocs
      nsendvals = send%buf_real( iproc )%nvals
      if (present(send_state)) then
        !$nec ivdep
        do ival = 1, nsendvals
          send%buf_real( iproc )%val( ival )                    &
            &  = send_state( send%buf_real( iproc )%pos( ival ) )
        end do
      else
        !$nec ivdep
        do ival = 1, nsendvals
          send%buf_real( iproc )%val( ival )               &
            &  = state( send%buf_real( iproc )%pos( ival ) )
        end do
      end if
      call mpi_isend(                               &
       &      send%buf_real( iproc )%val,        & ! buffer
       &      send%buf_real( iproc )%nvals,      & ! count
       &      rk_mpi,                        & ! data type
       &      send%proc(iproc),                     & ! target
       &      message_flag,                         & ! tag
       &      comm,                                 & ! comm
       &      send%rqhandle( iproc ),               & ! handle
       &      ierr )                                  ! error status
    end do

    do iproc = 1, recv%nprocs
      ! wait for any of receive buffer to be ready
      call mpi_waitany(   &
        &    recv%nprocs, & ! count
        &    recv%rqhandle,   & ! request handles
        &    finproc,     & ! process that finished
        &    status(:,1), & ! mpi status
        &    ierr )         ! error status
      nrecvvals = recv%buf_real( finproc )%nvals
      !$nec ivdep
      do ival = 1, nrecvvals
        ! write the values from the recv me to the state array
        ! to positions given in recvme%pos
        state( recv%buf_real( finproc )%pos( ival ) ) &
          &  = recv%buf_real( finproc )%val( ival )
      end do
    end do

    if (send%nprocs > 0) then
      ! wait for send buffer to be ready
      call mpi_waitall(                  &
        &    send%nprocs,                & ! total number of comm.'s to wait for
        &    send%rqhandle,              & ! request handles
        &    status,                     & ! mpi status
        &    ierr )                        ! error status
    end if

  end subroutine comm_isend_irecv_overlap_real
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> exchange the communication mes
  !! with a non-blocking mpi communication using preposted irecv and isend
  !! with a waitall
  !!
  subroutine comm_typed_isend_irecv_real( send, recv, state,             &
    &                                        message_flag, send_state, comm )
    ! -------------------------------------------------------------------- !
    type( tem_communication_type ), intent(inout) :: send, recv
    real(kind=rk), intent(inout) :: state(*)  !< current state vector
    integer, intent(in) :: message_flag
    real(kind=rk), intent(in), optional :: send_state(*)  !< data to send
    !> mpi communicator
    integer, intent(in) :: comm
    ! -------------------------------------------------------------------- !
    ! @todo request handle array could exist during complete code runtime
    ! integer :: rq_handle( recv%nprocs + send%nprocs )
    integer :: status( mpi_status_size, max(recv%nprocs, send%nprocs) )
    integer :: ierr             ! error flag
    integer :: iproc
    ! -------------------------------------------------------------------- !

    !> values for send me must have been copied from the actual state array
    do iproc = 1, recv%nprocs
      !> start receive communications
      call mpi_irecv(                                     &
       &              state,                              & ! buffer
       &              1,                                  & ! count
       &              recv%buf_real(iproc)%memindexed, & ! type
       &              recv%proc(iproc),                   & ! source
       &              message_flag,                       & ! tag
       &              comm,                               & ! comm
       &              recv%rqhandle(iproc),               & ! request handle
       &              ierr )

    end do

    !>  start the sending communications
    if (present(send_state)) then
      do iproc = 1, send%nprocs
        call mpi_isend(                                     &
          &             send_state,                         & ! buffer
          &             1,                                  & ! count
          &             send%buf_real(iproc)%memindexed, & ! type
          &             send%proc(iproc),                   & ! target
          &             message_flag,                       & ! tag
          &             comm,                               & ! comm
          &             send%rqhandle( iproc ),             & ! handle
          &             ierr )
      end do !< iproc
    else
      do iproc = 1, send%nprocs
        call mpi_isend(                                     &
          &             state,                              & ! buffer
          &             1,                                  & ! count
          &             send%buf_real(iproc)%memindexed, & ! type
          &             send%proc(iproc),                   & ! target
          &             message_flag,                       & ! tag
          &             comm,                               & ! comm
          &             send%rqhandle( iproc ),             & ! handle
          &             ierr )
      end do !< iproc
    end if

    !> wait for above communications to complete
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall( recv%nprocs,     & ! count
        &               recv%rqhandle,   & ! request handles
        &               status,          & ! statuses
        &               ierr )
    end if
    if ( send%nprocs /= 0 ) then
      call mpi_waitall( send%nprocs,     & ! count
        &               send%rqhandle,   & ! request handles
        &               status,          & ! statuses
        &               ierr )
    end if

  end subroutine comm_typed_isend_irecv_real
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Write communication type data to nUnit (debugging routine)
  !!
  subroutine tem_comm_dumpType(me, nUnit)
    ! -------------------------------------------------------------------- !
    type( tem_communication_type ), intent(in) :: me
    integer, intent(in) :: nUnit
    ! -------------------------------------------------------------------- !

    write(nUnit,*)        '----   Communication type    -----------------------'
    write(nUnit,"(A,I0)") '       nProcs: ', me%nProcs
    write(nUnit,"(A   )") '         proc: '//trim(tem_toStr( me%proc, ',' ))
    write(nUnit,"(A   )") '   nElemsProc: ' &
      &                   // trim(tem_toStr( me%nElemsProc, ',' ))
    write(nUnit,*)        '---------------------------------------------'

  end subroutine tem_comm_dumpType
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Allocate tem_communication_type and its variables
  !!
  subroutine tem_comm_init( me, nProcs )
    ! -------------------------------------------------------------------- !
    type( tem_communication_type ), intent(inout) :: me
    integer, intent(in) :: nProcs
    ! -------------------------------------------------------------------- !

    if ( allocated(me%proc) )       deallocate( me%proc )
    if ( allocated(me%nElemsProc) ) deallocate( me%nElemsProc )
    if ( allocated(me%elemPos) )    deallocate( me%elemPos )
    if ( allocated(me%rqHandle) )   deallocate( me%rqHandle )
    if ( allocated(me%buf_long) )   deallocate( me%buf_long )
    if ( allocated(me%buf_int) )    deallocate( me%buf_int )
    if ( allocated(me%buf_real) )   deallocate( me%buf_real )

    me%nProcs = nProcs
    allocate( me%proc      ( nProcs ) )
    allocate( me%nElemsProc( nProcs ) )
    allocate( me%elemPos   ( nProcs ) )
    allocate( me%rqHandle  ( nProcs ) )

    allocate( me%buf_long  ( nProcs ) )
    allocate( me%buf_int   ( nProcs ) )
    allocate( me%buf_real  ( nProcs ) )

  end subroutine tem_comm_init
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Allocate tem_communication_type and its variables
  !!
  subroutine tem_comm_count( me, comm_size, nHalos )
    ! -------------------------------------------------------------------- !
    type( tem_communication_type ), intent(inout) :: me
    !> communicator size
    integer, intent(in) :: comm_size
    !> number of halos for each other processes
    integer, intent(in) :: nHalos( comm_size )
    ! -------------------------------------------------------------------- !
    integer :: iPartner, iProc
    ! -------------------------------------------------------------------- !

    iPartner = 0
    do iProc = 1, comm_size
      if( nHalos( iProc ) > 0 ) then
        ! Store the processes numbers to receive from
        iPartner = iPartner + 1
        me%nElemsProc( iPartner ) = nHalos( iProc )
        me%proc( iPartner )       = iProc - 1
      end if
    end do

  end subroutine tem_comm_count
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Routine to build communication buffer using elemRanks.
  !! This routine can be used only if all elements need to be communicated
  !! but they need process-wise seperation.
  !! Uses nScalars to get position in the value array to communicate.
  !! For send buffer: elemRanks contains target ranks to send data to
  !! For recv buffer: elemRanks contains source ranks to recv data from
  subroutine tem_comm_createBuffer( commBuffer, nScalars, nElems, elemRanks )
    ! -------------------------------------------------------------------------!
    !> send or recv communication buffer to be created
    type(tem_communication_type), intent(out) :: commBuffer
    !> Number of scalars per element
    integer, intent(in) :: nScalars
    !> Total number of elements or points to communicate
    integer, intent(in) :: nElems
    !> Target or source rank for each element or point
    integer, intent(in) :: elemRanks(nElems)
    ! -------------------------------------------------------------------------!
    integer :: iElem, iProc, iVar, counter, pntPos
    type(dyn_intArray_type) :: partnerProc
    integer, allocatable :: pos(:)
    ! -------------------------------------------------------------------------!
    ! Create dynamic array of rank ids to communicate to initialize commBuffer
    do iElem = 1, nElems
      call append( me  = partnerProc,    &
        &          val = elemRanks(iElem) )
    end do

    ! Initialize commBuffer
    call tem_comm_init( me     = commBuffer,       &
      &                 nProcs = partnerProc%nVals )
    ! Store rank id to communicate
    do iProc = 1, partnerProc%nVals
      commBuffer%proc(iProc) = partnerProc%val(iProc)
    end do

    ! Create map from commBuffer to elem array
    do iElem = 1, nElems
      iProc = PositionOfVal( me  = partnerProc,     &
        &                    val = elemRanks(iElem) )
      call append( me  = commBuffer%elemPos(iProc), &
        &          val = iElem                       )
    end do
    call destroy(partnerProc)

    commBuffer%nElemsProc(:) = commBuffer%elemPos(:)%nVals

    ! Get position in state array to send or recv data
    allocate( pos( nScalars * maxVal( commBuffer%nElemsProc(:) ) ) )

    ! Assign comm buffer positions
    do iProc = 1, commBuffer%nProcs
      counter = 0

      ! loop of nElems per proc
      do iElem = 1, commBuffer%nElemsProc(iProc)
        ! position of this proc point in the point array
        pntPos = commBuffer%elemPos(iProc)%val(iElem)
        do iVar = 1, nScalars
          counter = counter + 1
          ! position in evalVal array which has size: nElems*nScalars
          pos(counter) = (pntPos-1)*nScalars + iVar
        end do !iVar
      end do !iElem
      ! copy position array to me%pos, allocate me%val array
      call tem_commbuf_real_fillPos( me    = commBuffer%buf_real(iProc), &
        &                            pos   = pos,                        &
        &                            nVals = counter                     )
    end do !iProc
    deallocate(pos)

  end subroutine tem_comm_createBuffer
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> All to all exchange of a single integer.
  !!
  !! This is a wrapper around the sparse alltoall implementation and overcome
  !! the lack of non-blocking collectives on some systems.
  subroutine tem_comm_alltoall_int( targets, send_buffer, &
    &                               sources, recv_buffer, &
    &                               comm, tag             )

    !> List of target ranks to send an integer to.
    integer, intent(in) :: targets(:)

    !> Data to send to the respective target ranks. This array has to have the
    !! same ordering as targets.
    integer, intent(in) :: send_buffer(:)

    !> List of ranks we received data from (source ranks).
    !! The array will be allocated with a size according to the number of
    !! processes that send a request to this process.
    integer, intent(out), allocatable :: sources(:)

    !> Received data from the sources. The array has the same size and ordering
    !! as the sources array.
    integer, intent(out), allocatable :: recv_buffer(:)

    !> MPI Communicator to use for this data exchange.
    integer, intent(in) :: comm

    !> Tag to use in the communications. Defaults to 22.
    integer, intent(in), optional :: tag

    integer :: nProcs
    integer :: nSources
    integer :: iProc, iSource
    integer :: iError
    integer, allocatable :: buf(:)

    if (use_sparse_alltoall) then

      call tem_sparse_alltoall_int( targets     = targets,     &
        &                           send_buffer = send_buffer, &
        &                           sources     = sources,     &
        &                           recv_buffer = recv_buffer, &
        &                           comm        = comm,        &
        &                           tag         = tag          )

    else

      call MPI_Comm_size(comm, nProcs, iError)
      allocate(buf(0:nProcs-1))
      buf = 0
      buf(targets(:)) = send_buffer
      call MPI_Alltoall( MPI_IN_PLACE, 1, MPI_INTEGER,     &
        &                buf, 1, MPI_INTEGER, comm, iError )
      nSources = count(buf/=0)
      allocate(sources(nSources))
      allocate(recv_buffer(nSources))
      recv_buffer = 0
      iSource = 1
      do iProc=0,nProcs-1
        if (buf(iProc) /= 0) then
          sources(iSource) = iProc
          recv_buffer(iSource) = buf(iProc)
          iSource = iSource + 1
        end if
      end do
      deallocate(buf)
    end if

  end subroutine tem_comm_alltoall_int
  ! *************************************************************************** !

  ! *************************************************************************** !
  subroutine tem_comm_destroy( me, commPattern )
    ! -------------------------------------------------------------------------!
    !> communication type to be destroyed
    type(tem_communication_type), intent(inout) :: me
    !> Communication pattern
    type(tem_commPattern_type), intent(in) :: commPattern
    ! -------------------------------------------------------------------------!
    integer :: iProc
    ! -------------------------------------------------------------------------!

    do iProc = 1, me%nProcs
      call commPattern%finbuf_real( me%buf_real(iProc) )
      call commPattern%finbuf_long( me%buf_long(iProc) )
      call commPattern%finbuf_int(  me%buf_int(iProc)  )
    end do

    deallocate( me%buf_real )
    deallocate( me%buf_long )
    deallocate( me%buf_int  )

    deallocate( me%proc )
    deallocate( me%nElemsProc )
    do iProc = 1, me%nProcs
      call destroy( me%elemPos(iProc) )
    end do
    deallocate( me%elemPos )
    deallocate( me%rqHandle )

  end subroutine tem_comm_destroy
  ! *************************************************************************** !

end module tem_comm_module
! **************************************************************************** !
