! Copyright (c) 2011-2013, 2017 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2012 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2012-2013 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
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
!> Declaration of the element property header data, which is globally defined
!! for all processes.
!!
!! Each property is described by a label, the bit position it is encoded in
!! in the `elempropertybits` of the mesh and the total number of elements with
!! this property found in the mesh.
module tem_prophead_module

  ! include treelm modules
  use mpi
  use env_module, only: long_k

  ! include aotus modules
  use aot_out_module,   only: aot_out_type, aot_out_val,                 &
    &                         aot_out_open_table, aot_out_close_table
  use aot_table_module, only: aot_table_open, aot_table_close, aot_get_val
  use aotus_module,     only: flu_State

  implicit none

  private

  public :: tem_prophead_type, tem_prop_countelems
  public :: load_tem_prophead, dump_tem_prophead

  type tem_prophead_type
    !> A label to describe this property
    character(len=64) :: label

    !> The position of the bit, which flags the presence
    !! of this property for an element.
    integer :: bitpos

    !> Total number of elements with this property
    integer(kind=long_k) :: nElems
  end type tem_prophead_type


contains


  ! ------------------------------------------------------------------------ !
  !> Load the property header from the mesh file.
  !!
  subroutine load_tem_prophead(me, myPart, comm, conf, root)
    ! -------------------------------------------------------------------- !
    !> Property to load
    type(tem_prophead_type), intent(out) :: me(:)

    !> The process within comm, calling the routine
    integer, intent(in) :: myPart

    !> MPI Communicator within which the mesh is loaded.
    integer, intent(in) :: comm

    !> Lua state containing the mesh description.
    type( flu_State ) :: conf

    !> The root process doing the IO, defaults to 0.
    integer, optional, intent(in) :: root
    ! -------------------------------------------------------------------- !
    integer :: locroot, iprop, sub_handle
    integer :: iError
    integer :: thandle ! handle of property table
    ! -------------------------------------------------------------------- !

    if (present(root)) then
      locroot = root
    else
      locroot = 0
    end if

    if (myPart == locroot) then
      !open properties table
      call aot_table_open( L       = conf,      &
        &                  thandle = thandle,   &
        &                  key     = 'property' )
      do iprop = 1, size(me)
        call aot_table_open( L       = conf,       &
          &                  parent  = thandle,    &
          &                  thandle = sub_handle, &
          &                  pos     = iprop       )
        call aot_get_val( L       = conf,            &
          &               thandle = sub_handle,      &
          &               key     = 'label',         &
          &               val     = me(iprop)%label, &
          &               ErrCode = iError           )
        call aot_get_val( L       = conf,             &
          &               thandle = sub_handle,       &
          &               key     = 'bitpos',         &
          &               val     = me(iprop)%bitpos, &
          &               ErrCode = iError            )
        call aot_get_val( L       = conf,             &
          &               thandle = sub_handle,       &
          &               key     = 'nElems',         &
          &               val     = me(iprop)%nElems, &
          &               ErrCode = iError            )
        call aot_table_close( L = conf, thandle = sub_handle )
      end do
      call aot_table_close( L = conf, thandle = thandle )
    end if

    do iprop = 1, size(me)
      call MPI_Bcast( me(iprop)%nElems, 1, MPI_INTEGER8, locroot, comm, iError)
      call MPI_Bcast( me(iprop)%bitpos, 1, MPI_INTEGER, locroot, comm, iError)
      call MPI_Bcast( me(iprop)%label, 64, MPI_CHARACTER, locroot, comm, iError)
    end do

  end subroutine load_tem_prophead
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Subroutine to dump property header informations into a Lua file.
  subroutine dump_tem_prophead(me, conf)
    ! -------------------------------------------------------------------- !
    !> Property to write
    type(tem_prophead_type), intent(in) :: me(:)
    !> aotus lua state to write output
    type(aot_out_type) :: conf
    ! -------------------------------------------------------------------- !
    integer :: iprop
    ! -------------------------------------------------------------------- !

    ! open property table
    call aot_out_open_table( conf, 'property' )
    do iprop=1, size(me)
      call aot_out_open_table( conf )
      call aot_out_val( put_conf = conf,                 &
        &               vname    = 'label',              &
        &               val      = trim(me(iprop)%label) )
      call aot_out_val( put_conf = conf,            &
        &               vname    = 'bitpos',        &
        &               val      = me(iprop)%bitpos )
      call aot_out_val( put_conf = conf,            &
        &               vname    = 'nElems',        &
        &               val      = me(iprop)%nElems )
      call aot_out_close_table(conf)
    end do
    call aot_out_close_table(conf)

  end subroutine dump_tem_prophead
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Subroutine to (re-)count the global number of elements the property.
  !!
  !! All elements with the given property are counted across all elements and
  !! the nElems setting is updated accordingly.
  subroutine tem_prop_countelems(me, elempropertybits, comm)
    ! -------------------------------------------------------------------- !
    !> Property to count and set the number of elements for.
    type(tem_prophead_type), intent(inout) :: me

    !> Propertybits of all elements.
    integer(kind=long_k), intent(in) :: elempropertybits(:)

    !> MPI communicator within which the mesh is distributed.
    integer, intent(in) :: comm
    ! -------------------------------------------------------------------- !
    integer(kind=long_k) :: localElems
    integer :: iError
    ! -------------------------------------------------------------------- !

    localElems = count(btest(elempropertybits, me%BitPos))

    ! MPI_ALLREDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)
    call MPI_Allreduce( localElems, me%nElems, 1, MPI_INTEGER8, MPI_SUM, &
      &                 comm, iError                                     )

  end subroutine tem_prop_countelems
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !

end module tem_prophead_module
! ---------------------------------------------------------------------------- !
