! Copyright (c) 2014-2015, 2021 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2014 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
!> This module provides the description of colors for elements.
!!
!! Colors might be used to describe different areas in the domain, for example
!! varying material parameters.
module tem_color_prop_module
  use mpi

  use env_module,          only: labelLen, PathLen, long_k, rk, rk_mpi
  use treelmesh_module,    only: treelmesh_type
  use tem_aux_module,      only: tem_open
  use tem_prophead_module, only: tem_prophead_type
  use tem_property_module, only: tem_property_type, prp_isColored

  use aotus_module,     only: flu_state, aot_get_val, &
    &                         open_config_file, close_config
  use aot_table_module, only: aot_table_open, aot_table_close
  use aot_out_module,   only: aot_out_type, aot_out_val, aot_out_open, &
    &                         aot_out_close, aot_out_open_table,       &
    &                         aot_out_close_table

  implicit none

  private

  public :: tem_color_prop_type
  public :: tem_color_prop_load
  public :: tem_color_prop_out

  !> Number of colors to store per character.
  !!
  !! We use the ASCII character set to encode the colors.
  !! Thus, 7 colors can be stored in each of them.
  !! This leaves us with unused space of 1/8 in the files.
  !! However, it is pretty likely that we will use less then
  !! 7 colors in common settings and in this case using a single
  !! byte for all colors is the smallest packing we can use
  !! conveniently.
  integer, parameter, public :: colors_per_char = 7

  type tem_color_prop_type
    !> Pointer to treelmesh_type%global%property
    type(tem_prophead_type),  pointer :: header => null()

    !> Number of colors present in the mesh.
    integer :: nColors

    !> Number of characters required to store a bit for each color.
    integer :: nChars

    !> Array of labels identifying each of the colors.
    !! This array has a length of nColors
    character(len=LabelLen), allocatable :: color_label(:)

    !> Value to use, where the color is present.
    real(kind=rk), allocatable :: color_fill(:)

    !> Value to use, where the color is not present.
    real(kind=rk), allocatable :: color_void(:)

    !> Actual color identification.
    !!
    !! For every element with this property a bitfield is stored, to
    !! indicate which colors it has.
    !! The first index has length nChars, and the second runs over all
    !! elements with an attached color.
    !! Using characters here, to minimize the required space for color
    !! encoding of few colors.
    character, allocatable :: colored_bit(:,:)

    !> Pointer to treelmesh_type%property
    type(tem_property_type),  pointer :: property => null()

  end type tem_color_prop_type


contains


  ! **************************************************************************** !
  !> Load the color property from disk.
  subroutine tem_color_prop_load( me, tree, myPart, comm )
    ! --------------------------------------------------------------------------!
    !> Color definitions to load.
    type(tem_color_prop_type), intent(out) :: me

    !> Tree to build the polynomial subresolution information for
    type(treelmesh_type), intent(in) :: tree

    !> Partition to load
    integer, intent(in) :: myPart

    !> Communicator to use
    integer, intent(in) :: comm
    ! --------------------------------------------------------------------------!
    integer, parameter :: root = 0
    type(flu_State) :: conf
    integer :: iError
    integer :: rl
    integer :: thandle
    integer :: fUnit
    integer :: i
    integer :: iProp
    character(len=pathLen) :: headerfile
    character(len=pathLen) :: datafile
    ! --------------------------------------------------------------------------!

    me%nColors = 0

    prp_loop: do iprop=1, tree%global%nProperties
      if (tree%global%Property(iprop)%bitpos == prp_isColored) then
        me%header => tree%global%Property(iprop)
        me%property => tree%property(iprop)

        headerfile = trim(tree%global%dirname)//'colors.lua'
        datafile   = trim(tree%global%dirname)//'colors.ascii'

        if (myPart == root) then
          ! Read the header only on the root process, broadcast to all others
          ! open mesh header file
          call open_config_file( L = conf, filename = headerfile )
          call aot_get_val( L       = conf,       &
            &               key     = 'nColors',  &
            &               val     = me%nColors, &
            &               ErrCode = iError      )
        end if

        call MPI_Bcast(me%nColors, 1, MPI_INTEGER, root, comm, iError)

        ! The number of colors that can be stored per character is fixed, thus
        ! the number of characters required by a given number of colors is
        ! immediatly known.
        me%nChars = ceiling(real(me%nColors)/real(colors_per_char))

        allocate(me%color_label(me%nColors))
        allocate(me%color_fill(me%nColors))
        allocate(me%color_void(me%nColors))
        allocate(me%colored_bit(me%nChars, me%property%nElems))

        if (myPart == root) then
          ! Now read the color labels on the root process.
          call aot_table_open( L = conf, thandle = thandle, &
            &                  key = 'color_label' )
          do i=1,me%nColors
            call aot_get_val( L       = conf,              &
              &               thandle = thandle,           &
              &               pos     = i,                 &
              &               val     = me%color_label(i), &
              &               ErrCode = iError             )
          end do
          call aot_table_close( L = conf, thandle = thandle )

          ! Now read the color fill values on the root process.
          call aot_table_open( L = conf, thandle = thandle, &
            &                  key = 'color_fill' )
          do i=1,me%nColors
            call aot_get_val( L       = conf,             &
              &               thandle = thandle,          &
              &               pos     = i,                &
              &               val     = me%color_fill(i), &
              &               ErrCode = iError            )
          end do
          call aot_table_close( L = conf, thandle = thandle )

          ! Now read the color void values on the root process.
          call aot_table_open( L = conf, thandle = thandle, &
            &                  key = 'color_void' )
          do i=1,me%nColors
            call aot_get_val( L       = conf,             &
              &               thandle = thandle,          &
              &               pos     = i,                &
              &               val     = me%color_void(i), &
              &               ErrCode = iError            )
          end do
          call aot_table_close( L = conf, thandle = thandle )

          call close_config(conf)
        end if

        call MPI_Bcast( me%color_label, LabelLen*me%nColors, MPI_CHARACTER, &
          &             root, comm, iError                                  )

        call MPI_Bcast( me%color_fill, me%nColors, rk_mpi, &
          &             root, comm, iError                 )
        call MPI_Bcast( me%color_void, me%nColors, rk_mpi, &
          &             root, comm, iError                 )

        ! If there are actually colored elements on the local process,
        ! read them now.
        if (me%property%nElems > 0) then

          inquire(iolength=rl) me%colored_bit(:,1)
          call tem_open( newunit = fUnit,         &
            &            file    = datafile,      &
            &            action  = 'read',        &
            &            access  = 'stream',      &
            &            form    = 'unformatted', &
            &            status  = 'old'          )
          read(fUnit, pos=me%property%offset+1) me%colored_bit
          close(fUnit)

        end if

        EXIT prp_loop

      end if
    end do prp_loop

  end subroutine tem_color_prop_load
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Write the color property to disk.
  subroutine tem_color_prop_out( me, offset, nElems, basename, myPart, comm )
    ! --------------------------------------------------------------------------!
    !> Color definitions to load.
    type(tem_color_prop_type), intent(out) :: me

    !> Offset of the local set of elements in the global list
    integer(kind=long_k), intent(in) :: offset

    !> Local number of elements with this property
    integer, intent(in) :: nElems

    !> Name of the file, the data is stored in, will be appended with
    !! ".lua" for the header information and ".lsb" or ".msb" for the
    !! binary data. An ".ascii" extension will be used for ASCII data.
    character(len=*), intent(in) :: basename

    !> Partition to load
    integer, intent(in) :: myPart

    !> Communicator to use
    integer, intent(in) :: comm
    ! --------------------------------------------------------------------------!
    type(aot_out_type) :: conf
    integer :: locomm
    integer :: root
    integer :: rl
    integer :: fUnit
    integer :: i
    character(len=pathLen) :: headerfile
    character(len=pathLen) :: datafile
    ! --------------------------------------------------------------------------!

    root = 0

    locomm = comm

    headerfile = trim(basename)//'.lua'
    datafile   = trim(basename)//'.ascii'

    if (myPart == root) then
      ! Write the header only on the root process
      call aot_out_open( put_conf = conf, filename = headerfile )
      call aot_out_val( put_conf = conf, vname='nColors', val = me%nColors )
      call aot_out_open_table(conf, 'color_label')
      do i=1,me%nColors
        call aot_out_val( put_conf = conf, val = me%color_label(i) )
      end do
      call aot_out_close_table(conf)
      call aot_out_close(conf)
    end if

    ! If there are actually colored elements on the local process,
    ! write them now.
    ! @todo HK: replace this color output by MPI-IO.
    if (nElems > 0) then

      inquire(iolength=rl) me%colored_bit(:,1)
      call tem_open( newunit = fUnit,         &
        &            file    = datafile,      &
        &            action  = 'write',       &
        &            access  = 'direct',      &
        &            form    = 'unformatted', &
        &            recl    = rl,            &
        &            status  = 'replace'      )
      write(fUnit, rec=offset+1) me%colored_bit
      close(fUnit)

    end if

  end subroutine tem_color_prop_out
  ! **************************************************************************** !

end module tem_color_prop_module
