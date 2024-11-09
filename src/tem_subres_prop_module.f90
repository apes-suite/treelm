! Copyright (c) 2014, 2019, 2021 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
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
! **************************************************************************** !
!> This module provides the description of subresolution for elements.
!!
!! Subresolutions provide polynomial respresentations of boundaries within
!! elements. This property depends on the color property. Only colored elements
!! can provide subresolution data. For each color there can be an independent
!! subresolved representation of the according boundaries.
module tem_subres_prop_module
  use mpi
  use env_module, only: pathlen, long_k

  use tem_aux_module, only: tem_open
  use treelmesh_module, only: treelmesh_type
  use tem_prophead_module, only: tem_prophead_type
  use tem_property_module, only: tem_property_type, prp_hasPolynomial
  use tem_color_prop_module, only: tem_color_prop_type, colors_per_char

  implicit none

  type elemid_list_type
    integer, allocatable :: id(:)
  end type elemid_list_type

  type tem_subres_prop_type
    !> Pointer to treelmesh_type%global%property
    type(tem_prophead_type),  pointer :: header => null()

    !> Indication which of the colors contain subresolution information.
    !!
    !! The first index has length nChars from the color property.
    !! The second runs through all elements with with subresolution property.
    character, allocatable :: subresolved_colors(:,:)

    !> Pointer to treelmesh_type%property
    type(tem_property_type),  pointer :: property => null()

    !> Number of subresolved elements in each color.
    integer, allocatable :: nElems(:)

    !> List of indices of elements with subresolution for each color.
    type(elemid_list_type), allocatable :: elem(:)

    !> Offset for the subresolved elements on this partition for each color.
    integer(kind=long_k), allocatable :: offset(:)

  end type tem_subres_prop_type


contains


  ! ************************************************************************** !
  !> Load the subresolution property from disk.
  !!
  !! Before this can be done, the coloring information has to have been loaded.
  subroutine tem_subres_prop_load( me, tree, coloring )
    ! -----------------------------------------------------------------------!
    !> Color definitions to load.
    type(tem_subres_prop_type), intent(out) :: me

    !> Tree to build the polynomial subresolution information for
    type(treelmesh_type), intent(in) :: tree

    !> Information on the colors in the mesh.
    type(tem_color_prop_type), intent(in) :: coloring
    ! -----------------------------------------------------------------------!
    integer :: rl
    integer :: fUnit
    character(len=pathLen) :: datafile
    integer :: iColor
    integer :: iProp
    integer :: colChar, colBit
    integer(kind=long_k), allocatable :: long_counts(:)
    integer :: iError
    integer :: iElem
    integer :: ice
    ! -----------------------------------------------------------------------!

    prp_loop: do iprop=1, tree%global%nProperties
      if (tree%global%Property(iprop)%bitpos == prp_hasPolynomial) then
        me%header => tree%global%Property(iprop)
        me%property => tree%property(iprop)

        datafile   = trim(tree%global%dirname)//'subres.ascii'

        allocate(me%subresolved_colors(coloring%nChars, me%property%nElems))
        allocate(me%nElems(coloring%nColors))
        allocate(me%Offset(coloring%nColors))
        me%Offset = 0_long_k

        ! If there are actually subresolved elements on the local process,
        ! read them now.
        if (me%property%nElems > 0) then

          allocate(me%elem(coloring%nColors))
          inquire(iolength=rl) me%subresolved_colors(:,1)
          call tem_open( newunit = fUnit,         &
            &            file    = datafile,      &
            &            action  = 'read',        &
            &            access  = 'stream',      &
            &            form    = 'unformatted', &
            &            status  = 'old'          )
          read(fUnit, pos=me%property%offset+1) me%subresolved_colors
          close(fUnit)

          do iColor=1,coloring%nColors
            colChar = (iColor-1)/colors_per_char + 1
            colBit = mod(iColor-1, colors_per_char)
            me%nElems(iColor)                                            &
              &  = count( btest(ichar(me%subresolved_colors(ColChar,:)), &
              &                 ColBit)                                  )
            allocate( me%elem(iColor)%id(me%nElems(iColor)) )

            ! Store the element link for each color.
            if (me%nElems(iColor) > 0) then
              ice = 0
              do iElem=1,me%property%nElems
                if ( btest(ichar(me%subresolved_colors(ColChar,iElem)), &
                  &        colBit) ) then
                  ice = ice + 1
                  me%elem(iColor)%id(ice) = me%property%elemid(iElem)
                end if
              end do
            end if

          end do

        else
          ! No subresolved elements on the local partition:
          me%nElems = 0
        end if

        allocate( long_counts(coloring%nColors) )
        long_counts = me%nElems
        call MPI_Exscan( long_counts, me%Offset, coloring%nColors,       &
          &              MPI_INTEGER8, MPI_SUM, tree%global%comm, iError )

        EXIT prp_loop

      end if
    end do prp_loop

  end subroutine tem_subres_prop_load

end module tem_subres_prop_module
