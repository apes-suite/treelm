! Copyright (c) 2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2013, 2015-2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2014-2015 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
! ****************************************************************************** !
!> This module  provides descriptions for elements
module tem_halo_module

  ! include treelm modules
  use tem_dyn_array_module, only: dyn_intArray_type, init, append, &
    &                             destroy
  use tem_arrayofarrays_module, only: grw_dynintArray_type, init, append, &
    &                                 destroy

  implicit none

  private

  public :: tem_haloList_type
  public :: tem_halo_init
  public :: tem_halo_destroy
  public :: tem_halo_append

  !> Data structure to store the halo elements according to the
  !! partner process with which they need to be exchanged.
  !!
  !! Partner processes are identified by their number (rank+1), and
  !! stored in a dynamic array of integers (partnerProc).
  !! For each partnerProc there is a dynamic array of halo element ids,
  !! stored in a growing array of dynamic integer arrays (halos).
  type tem_haloList_type
    !> Process numbers (rank+1) of partner processes
    type(dyn_intArray_type) :: partnerProc

    !> List of my halo elements, which I request from partnerProc.
    !!
    !! The growing array follows the same ordering as the partnerProc list.
    type(grw_dynintArray_type) :: halos
  end type tem_haloList_type


contains


  !> Initialize the list of halo elements.
  subroutine tem_halo_init(me)
    !> The list of halo elements to initialize
    type(tem_haloList_type), intent(out) :: me

    call init(me%partnerProc)
    call init(me%halos)
  end subroutine tem_halo_init


  !> Destroy the list of halo elements.
  subroutine tem_halo_destroy(me)
    !> The list of halo elements to destroy.
    type(tem_haloList_type), intent(inout) :: me

    call destroy(me%partnerProc)
    call destroy(me%halos)
  end subroutine tem_halo_destroy


  !> Append an element to the list of halo elements.
  !!
  !! Each entry needs to be identified by the process (proc), the halo has to be
  !! exchanged with, and the position of the halo element in the local array
  !! of elements (elemPos).
  subroutine tem_halo_append(me, proc, elemPos, wasAdded)
    !> List of halo elements, this entry has to be appended to.
    type(tem_haloList_type), intent(inout) :: me

    !> Process this element is exchanged with.
    integer, intent(in) :: proc

    !> Local position of the halo element.
    integer, intent(in) :: elemPos

    !> Flag, wether this halo element was newly added, or already there.
    logical, optional, intent(out) :: wasAdded

    type(dyn_intArray_type) :: newHalos
    integer :: partnerPos
    integer :: haloPos
    logical :: newProc

    ! First attempt to append the process, this element is to be exchanged with.
    ! The call will return, where the process is to be found, and wether it was
    ! a new process or not.
    call append( me       = me%partnerProc, &
      &          val      = proc,           &
      &          pos      = partnerPos,     &
      &          wasAdded = newProc         )

    if (newProc) then
      ! If the process was newly added, we need to create a new list for halo
      ! elements. We do this by setting up a temporary dynamic array of integers
      ! and append this then to the halos component (which is a growing array
      ! of dynamic arrays).
      call append( me  = me%halos, &
        &          val = newHalos  )
    end if

    ! Now, the necessary arrays for this process should be available, and we
    ! can add the actual halo element to it.
    call append( me       = me%halos%val(partnerPos), &
      &          val      = elemPos,                  &
      &          pos      = haloPos,                  &
      &          wasAdded = wasAdded                  )

  end subroutine tem_halo_append

end module tem_halo_module
! ****************************************************************************** !
