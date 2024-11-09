! Copyright (c) 2011 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011-2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Manuel Hasert <m.hasert@grs-sim.de>
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
! ****************************************************************************** !
!> author: Jens Zudrop, Harald Klimach
!! this module is intended to collect all data types and informations for
!! interpolation of data for cells of different refinement levels (so
!! ghost or halo cells).
!!
module tem_interpolation_module

  ! include treelm modules
  use tem_construction_module, only: tem_levelDesc_type
  use treelmesh_module,        only: treelmesh_type

  implicit none

  private

! ****************************************************************************** !
  !> implicit interface which can be used to create interpolation routines
  !! for ghost and halo cells within the solvers and to pass them as an
  !! arguement to the interpolation subroutine interpolate.
  abstract interface
    subroutine solver_interpolation_down( target_cell_data,                    &
      &                                   target_cell_variables,               &
      &                                   target_cell_level,                   &
      &                                   source_cell_data,                    &
      &                                   source_cell_variables,               &
      &                                   source_child_numbers )
      !> The data of the cell you want to reconstruct. The subroutine
      !! implementing this interface does not allocate memory, so the calling
      !! subroutine has to allocate memory for this subroutine. The length
      !! of this array is the number of degrees of freedom (in case of PNPM
      !! or FD) or number of links (in case of LBM) you want to reconstruct.
      real, allocatable, intent(inout) :: target_cell_data(:)

      !> the level of the cell you want to reconstruct inside the octree.
      !! Since we do a downwards interpolation (from the higher refinement
      !! level to the lower one) we implicitly assume that the source cell
      !! levels are all equal to target_cell_level - 1.
      integer,intent(in) :: target_cell_level

      !> the number of variables you want to interpolate for the target cell.
      !! In case of PNPM it is the number of degrees of freedom (e.g. see
      !! ATELES: scheme_type%nDoFs ). In case of LBM it is the number of
      !! links (e.g. in case of D3Q19 it is 19).
      integer,intent(in) :: target_cell_variables

      !> Data of the source cells, which are necessary to interpolate the
      !! target cell on a lower level. Since we are working on an octree
      !! structure we assume to get everytime the data of all 8 children.
      !! So the dimensions of this array are:
      !! First dimension is 8, second dimension is source_cell_variables.
      !! If one source cell is a solid cell put a negative integer number into
      !! source_child_numbers at the right position to indicate that this child
      !! is not valid.
      real, allocatable, intent(in) :: source_cell_data(:,:)

      !> array of integer numbers of length 8 indicating the position of the
      !! children within the space filling curve. If one child is a solid cell
      !! and therefore not part of the interpolation, put a negative integer
      !! inside this array. This tells the subroutine, that this child will not
      !! be part of the interpolation step.
      integer, allocatable, intent(in) :: source_child_numbers(:)

      !> the number of variables of the source cell.
      !! In case of PNPM it is the number of degrees of freedom (e.g. see
      !! ATELES: scheme_type%nDoFs ) of the source cells.
      !! In case of LBM it is the number of links (e.g. in case of D3Q19
      !! it is 19).
      integer,intent(in) :: source_cell_variables
    end subroutine
  end interface
! ****************************************************************************** !

!> @todo Actually implement something in this interpolation module.
!HK!contains
!HK!
!HK!!> subroutine to interpolate data for ghost and halo cells (which are not
!HK!!! existing on this refinement level and therefore have to be reconstructed
!HK!!! from cells on other refinement levels).
!HK!subroutine interpolate(tree, levelDesc, numRefineLevels, f)
!HK!  ! ----------------------------------------------------------------------------
!HK!  !> procedure pointer to the solver specific interpolation routine. If you
!HK!  !! want to use the interpolate routine inside your own solver you should
!HK!  !! have a look at the Ateles implementation of an interpolation routine to
!HK!  !! see how it can be implemented.
!HK!  procedure(solver_interpolation_down), pointer, intent(in) :: f
!HK!  !> array of descriptors describing the different levels of the mesh.
!HK!  !! These descriptors should be generated by the buildDependencyList of the
!HK!  !! tem_construction_module. The dimension of this array has to be
!HK!  !! numRefineLevels.
!HK!  type(tem_levelDesc_type), intent(in), allocatable :: levelDesc(:)
!HK!  !> tree containing the informations of this partition of the mesh.
!HK!  type(treelmesh_type) , intent(in) :: tree
!HK!  !> the number of refinement levels, i.e. the size of the levelDesc
!HK!  !! array.
!HK!  integer, intent(in) :: numRefineLevels
!HK!
!HK!!> @todo we need the data array of the solver for this!
!HK!
!HK!  ! ----------------------------------------------------------------------------
!HK!  ! ----------------------------------------------------------------------------
!HK!
!HK!
!HK!!  call f(levelDesc)
!HK!
!HK!  !> @todo implement this subroutine
!HK!
!HK!
!HK!
!HK!end subroutine interpolate


end module tem_interpolation_module
! ****************************************************************************** !
