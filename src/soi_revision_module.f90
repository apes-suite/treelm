!> SOIL module for holding the revision and compilation information
!! of the executable.
!!
!! This information will be written by soi_world_init.
!! This source file is generated during compilation!
! *************************************************************************** !
! WARNING: Do NOT change this file, as it will be overwritten during
!          compilation.
!          See bin/revision_module.py for the generating script.
!          (in the apes build infrastructure project)
! *************************************************************************** !

module soi_revision_module
  implicit none
  !> The HG revision of the application used for this executable.
  character(len=13), parameter :: soi_solver_revision &
    &                            = '75d5a154054a'

  !> Name of the compiler.
  character(len=32), parameter :: soi_FC_name &
    &                            = 'GFORTRAN'

  !> The compilation command that was used to build this executable.
  character(len=32), parameter :: soi_FC_command &
    &                            = 'mpif90'

  !> The version of the Fortran compiler used in the compilation of this
  !! executable.
  character(len=32), parameter :: soi_FC_version &
    &                            = '13.3.0'

  !> Number of lines needed to represent the compiler flags
  integer, parameter :: soi_FC_nFlagLines = 2

  !> The Fortran compiler flags used to compile this executable.
  character(len=72), parameter :: soi_FC_flags(soi_FC_nFlagLines) &
    & = [ '-O3 -march=native -Wall -Wno-maybe-uninitialized -Werror -std=f2008 -ped', &
    &     'antic-errors -O3 -march=native -Wall -Wno-maybe-uninitialized -Werror   ' ]

  !> The date when this executable was built.
  character(len=10), parameter :: soi_build_date &
    &                            = '2025-04-28'
end module soi_revision_module
