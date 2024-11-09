! Copyright (c) 2011-2014, 2019, 2023 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2011-2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2011-2012, 2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2014 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Raphael Haupt <Raphael.Haupt@student.uni-siegen.de>
! Copyright (c) 2018 Daniel Fleischer <daniel.fleischer@student.uni-siegen.de>
! Copyright (c) 2021 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
!> Module to describe additional elmental properties, like boundary conditions.
!!
!! Every element in the octree structure can be
!! assigned several properties.
!! The description of the properties is kept abstract, such that very different
!! properties can be assigned to elements in a similar manner.
!! As 64 bits are available to describe additional properties, we could attach
!! 64 different kinds of data to each element and each one combinable with all
!! other properties.
!! The most important data that has to be attached to elements in nearly every
!! mesh are the boundary conditions in [[tem_bc_module]], others are for
!! example vertex displacements in [[tem_vrtx_module]]
!! and regional influences.
!! As the boundary conditions are such a basic part of the general mesh
!! description, we present their property implementation exemplary in more
!! detail.
!! Each property is indicated by a bit in the property bit-mask.
!! Which bit is attached to what kind of property is defined in the header file.
!! The property might then be linked to further information stored in additional
!! files.
!! See [[tem_bc_module]] for details of the boundary implementation.
!!
!!
!! The amount and kind of additional data that might be attached to the elements
!! by the properties is arbitrary and could involve several more
!! indirections.
!! For example an additional non-elemental file is typically required to
!! describe vertex displacements of elements.
!! Thus, for elements, where the vertices should be moved, there is an elemental
!! list of integers indicating which vertices to use for the eight corners of
!! the element.
!! The actual coordinates are stored in a separate file with the coordinates for
!! each vertex.
!!
module tem_property_module

  ! include treelm modules
  use mpi
  use env_module,          only: long_k
  use tem_prophead_module, only: tem_prophead_type

  implicit none

  private

  public :: tem_property_type
  public :: init_propertylist
  public :: destroy_propertylist
  public :: gather_property

  !> An auxilary data type to describe modifications and additional features
  !! for some elements.
  type tem_property_type
    !> Number of local elements with this property
    integer :: nElems

    !> Offset of the local chunk of elements with this property in the list of
    !! all elements with that properties on disk
    integer(kind=long_k) :: Offset

    !> The indices of elements in the local partition, that have this
    !! property.
    integer, allocatable :: ElemID(:)
  end type

  !! Element properties:
  !! each element might have up to 64 properties attached to it
  !! the bitpos indicates just the position of the bit in the field
  !! indicating, which should flag the corresponding property
  !! A 0 indicates the absence of the property and a 1 the presence
  !! Thus, an element with all bits set to 0 has no additional properties
  !! Must be smaller than 64, because the respective bits will be set in
  !! a 64-bit integer.

  !> fluid
  integer,parameter,public :: prp_fluid    = 1
  !> solid
  integer,parameter,public :: prp_solid    = 2
  !> The element has boundaries. closer defined in the property module
  integer,parameter,public :: prp_hasBnd   = 3
  !> deformed Elements
  integer,parameter,public :: prp_defElems = 4
  !> changed Elements
  !> This is now used for adaptive refinement
  integer,parameter,public :: prp_chgElems = 5
  !> Fluidifiable Elements
  integer,parameter,public :: prp_fluidify = 6
  !> vertex information
  integer,parameter,public :: prp_vrtx     = 7
  !> has qVal elements
  integer,parameter,public :: prp_hasQVal  = 8
  !> the element is a parent element of one or more surface data points
  integer,parameter,public :: prp_hasIBM   = 9
  !> element has a remote neighbor
  integer,parameter,public :: prp_hasRemoteNgh = 10
  !> is this a fluid element which needs to be send to an attached proc
  integer,parameter,public :: prp_sendHalo = 12
  !> Solidification is not allowed for this element.
  integer,parameter,public :: prp_noSolidification = 13
  !> There is a boundary normal stored for this element
  integer,parameter,public :: prp_hasNormal = 15
  !> the complete neighborhood is required for this element
  integer,parameter,public :: prp_requireFullNeighborhood = 20

  !> indicates whether this fine ghost cell is the closest to the fluid cell
  integer,parameter,public :: prp_fineGhostClosestToFluid = 25

  !> Indicate if the element has some coloring attached to it.
  integer, parameter, public :: prp_isColored = 30

  !> This bit indicates wether the color in the element is further resolved by
  !! polynomial information.
  integer, parameter, public :: prp_hasPolynomial = 31

  !-- TV: PARTICLE PROPERTY --!
  integer, parameter, public :: prp_particle = 32


contains


  ! ************************************************************************** !
  !> Defines how many properties there are.
  !!
  !! Allocate the property and header lists accordingly.
  !!
  subroutine init_propertylist(header, property, nproperties)
    ! ---------------------------------------------------------------------- !
    !> Pointer to the list of headers for the properties to be allocated
    type(tem_prophead_type), pointer :: header(:)
    !> Pointer to the list of properties to be allocated
    type(tem_property_type), pointer :: property(:)
    !> Number of properties to allocate in the list
    integer, intent(in) :: nproperties
    ! ---------------------------------------------------------------------- !

    allocate(property(nProperties), header(nProperties))

  end subroutine init_propertylist
  ! ************************************************************************** !


  ! ************************************************************************** !
  !! deallocate the property and header lists accordingly.
  !!
  subroutine destroy_propertylist(header, property)
    ! ---------------------------------------------------------------------- !
    !> Pointer to list of property headers to destroy
    type(tem_prophead_type), pointer :: header(:)
    !> Pointer to list of properties to destroy
    type(tem_property_type), pointer :: property(:)
    ! ---------------------------------------------------------------------- !

    deallocate(property)
    deallocate(header)
  end subroutine destroy_propertylist
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Gather the information on a property from the bit fields of all elements
  !!
  subroutine gather_property(Property, Header, BitField, comm)
    ! ---------------------------------------------------------------------- !
    !> Property to gather
    type(tem_property_type), intent(out) :: Property
    !> Header for this property
    type(tem_prophead_type), intent(in) :: Header
    !> The BitField for the properties of all local elements
    integer(kind=long_k), intent(in) :: BitField(:)
    !> Communicator to act on
    integer, intent(in) :: comm
    ! ---------------------------------------------------------------------- !
    integer :: nElems
    integer :: iElem, PropElem
    integer(kind=long_k) :: myElems
    integer :: iError
    ! ---------------------------------------------------------------------- !

    nElems = size(BitField)

    ! First count the number of local elements with the given property.
    Property%nElems = count(btest(BitField, Header%BitPos))
    myElems = Property%nElems

    Property%Offset = 0
    ! Calculate offset on each process, by summing the number of elements
    ! on all lower ranks.
    call MPI_Exscan( myElems, Property%Offset, 1, MPI_INTEGER8, MPI_SUM, &
      &              comm, iError)

    ! Allocate an array to store the link from the list of elements with this
    ! property to the list of all elements.
    ! (Property%ElemID -> tree%treeID)
    allocate(Property%ElemID(Property%nElems))

    PropElem = 0
    do iElem=1,nElems
      ! Run over all elements.
      if (btest(BitField(iElem), Header%BitPos)) then
        ! If the element has the property, increase the counter for elements
        ! with this property and store this position for later lookups.
        PropElem = PropElem + 1
        Property%ElemID(PropElem) = iElem
      end if
    end do

  end subroutine gather_property
  ! ************************************************************************** !

end module tem_property_module
