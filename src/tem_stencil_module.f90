! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011-2015,2019,2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2017, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2014-2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2015-2016, 2019-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2020 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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
!> Stencil definitions for the elements
!!
!! A stencil is basically a set of element-offsets \( (s_x, s_y, s_z) \),
!! describing the relative positions of all required elements for a given
!! element.
!!
!! There is an example in the Documentation named
!! 'Example for the stencil construction'.
!!
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012, 2015-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013 Daniel Harlacher <d.harlacher@grs-sim.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2015-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2015-2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Daniel Fleischer <daniel.fleischer@student.uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
!
! Parts of this file were written by Harald Klimach, Simon Zimny and Manuel
! Hasert for German Research School for Simulation Sciences GmbH.
!
! Parts of this file were written by Harald Klimach, Kannan Masilamani,
! Daniel Harlacher, Kartik Jain, Verena Krupp, Jiaxing Qi, Peter Vitt,
! Daniel Fleischer, Tobias Girresser and Daniel Petró for University Siegen.
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

! This file contains the source code for growing and dynamic arrays.
! This is used for arrays of primitives (int, long_int, real, ...) as well as
! for arrays of derived datatypes (tem_variable_type,...).
!
! To use these macros include the following to your source file.
!
! Smart growing array (GA) for ?tstring?
! Growing Arrays:
!
! declaration
!
!
! implementation
!

! -----------------------------------------------------------------
! 2d Array, which can grow in second dimension only (GA2d)
! tname ... indicates type of dynamic array (long, int, real, ...)


!
!------------------------------------------------------------------------------
!
! dynamic Arrays:
!
! declaration
!
!
! implementation
!
!!
module tem_stencil_module

  ! include treelm modules
  use env_module,            only: rk, long_k, zeroLength, minLength, labelLen
  use tem_aux_module,        only: tem_abort
  use tem_param_module,      only: qQQQ, qOffset, qDirName
  use tem_comm_module,       only: tem_communication_type, tem_commPattern_type
  use tem_comm_env_module,   only: tem_comm_env_type
  use tem_grow_array_module, only: grw_intArray_type, init, destroy
  use tem_dyn_array_module,  only: dyn_longArray_type
  use tem_logging_module,    only: logUnit, tem_toStr
  use tem_geometry_module,   only: tem_determine_discreteVector
  use tem_tools_module,      only: tem_printArray

  ! include aotus modules
  use aotus_module,     only: flu_State
  use aot_table_module, only: aot_table_open, aot_table_close, &
    &                         aot_table_length, aot_get_val

  implicit none

  private

  public :: tem_stencilHeader_type
  public :: tem_stencilElement_type
  public :: tem_stencil_dump
  public :: tem_stencil_getHeaderPos
  public :: tem_loadStencil
  public :: tem_identify_prevailDirections
  public :: tem_identify_inverseDirections
  public :: tem_stencil_findIndexOfDir
  public :: tem_stencil_map_toTreelmDef
  public :: tem_treelmDef_map_toStencil
  public :: tem_stencil_zeroPos
  public :: tem_stencil_getLabelForcxDir
  public :: d3q125_cxDir, d3q81_cxDir
  public :: d3q19_cxDir, d3q27_cxDir
  public :: d3q7_cxDir, d2q9_cxDir, d1q3_cxDir, d2q5_cxDir
  public :: d3q6_cxDir, d2q4_cxDir, d1q2_cxDir
  public :: tem_create_stencil
  public :: grw_stencilHeaderArray_type
  public :: grw_stencilElementArray_type
  public :: init, destroy, append, expand, empty, truncate, placeAt

  ! cartesian directions of discrete velocities
  integer, parameter :: d3q19_cxDir(3,19)          &
    !                              axis parallel offsets
    &                   = reshape( [ -1,  0,  0,   & ! West(W)
    &                                 0, -1,  0,   & ! South(S)
    &                                 0,  0, -1,   & ! Bottom(B)
    &                                 1,  0,  0,   & ! East(E)
    &                                 0,  1,  0,   & ! North(N)
    &                                 0,  0,  1,   & ! Top(T)
    !                              diagonal offsets
    &                                 0, -1, -1,   & ! Bottom South(BS)
    &                                 0, -1,  1,   & ! Top South(TS)
    &                                 0,  1, -1,   & ! Bottom North(BN)
    &                                 0,  1,  1,   & ! Top North(TN)
    &                                -1,  0, -1,   & ! Bottom West(BW)
    &                                 1,  0, -1,   & ! Bottom East(BE)
    &                                -1,  0,  1,   & ! Top West(TW)
    &                                 1,  0,  1,   & ! Top East(TE)
    &                                -1, -1,  0,   & ! South West(SW)
    &                                -1,  1,  0,   & ! North West(NW)
    &                                 1, -1,  0,   & ! South East(SE)
    &                                 1,  1,  0,   & ! North East(NE)
    !                              central offset
    &                                 0,  0,  0 ], & ! Center(C)
    &                              [ 3, 19 ]       )

  ! cartesian directions of discrete velocities
  integer, parameter :: d3q13_cxDir(3,13)          &
    !                              diagonal offsets
    &                   = reshape( [ 0, -1, -1,    & ! Bottom South(BS)
    &                                0, -1,  1,    & ! Top South(TS)
    &                                0,  1, -1,    & ! Bottom North(BN)
    &                                0,  1,  1,    & ! Top North(TN)
    &                               -1,  0, -1,    & ! Bottom West(BW)
    &                                1,  0, -1,    & ! Bottom East(BE)
    &                               -1,  0,  1,    & ! Top West(TW)
    &                                1,  0,  1,    & ! Top East(TE)
    &                               -1, -1,  0,    & ! South West(SW)
    &                               -1,  1,  0,    & ! North West(NW)
    &                                1, -1,  0,    & ! South East(SE)
    &                                1,  1,  0,    & ! North East(NE)
    !                              central offset
    &                                 0,  0,  0 ], & ! Center(C)
    &                              [ 3, 13 ]       )

  ! cartesian directions of discrete velocities
  integer, parameter :: d3q27_cxDir(3,27)          &
    !                              axis parallel offsets
    &                   = reshape( [ -1,  0,  0,   & ! West(W)
    &                                 0, -1,  0,   & ! South(S)
    &                                 0,  0, -1,   & ! Bottom(B)
    &                                 1,  0,  0,   & ! East(E)
    &                                 0,  1,  0,   & ! North(N)
    &                                 0,  0,  1,   & ! Top(T)
    !                              diagonal offsets
    &                                 0, -1, -1,   & ! Bottom South(BS)
    &                                 0, -1,  1,   & ! Top South(TS)
    &                                 0,  1, -1,   & ! Bottom North(BN)
    &                                 0,  1,  1,   & ! Top North(TN)
    &                                -1,  0, -1,   & ! Bottom West(BW)
    &                                 1,  0, -1,   & ! Bottom East(BE)
    &                                -1,  0,  1,   & ! Top West(TW)
    &                                 1,  0,  1,   & ! Top East(TE)
    &                                -1, -1,  0,   & ! South West(SW)
    &                                -1,  1,  0,   & ! North West(NW)
    &                                 1, -1,  0,   & ! South East(SE)
    &                                 1,  1,  0,   & ! North East(NE)
    !                              corner offsets
    &                                -1, -1, -1,   & ! Bottom South West(BSW)
    &                                -1, -1,  1,   & ! Top South West(TSW)
    &                                -1,  1, -1,   & ! Bottom North West(BNW)
    &                                -1,  1,  1,   & ! Top North West(TNW)
    &                                 1, -1, -1,   & ! Bottom South East(BSE)
    &                                 1, -1,  1,   & ! Top South East(TSE)
    &                                 1,  1, -1,   & ! Bottom North East(BNE)
    &                                 1,  1,  1,   & ! Top North East(TNE)
    !                              central offset
    &                                 0,  0,  0 ], & ! Center(C)
    &                              [ 3, 27 ]                            )

  integer, parameter :: d3q7_cxDir(3,7)            &
    !                              axis parallel offsets
    &                   = reshape( [ -1,  0,  0,   & ! West(W)
    &                                 0, -1,  0,   & ! South(S)
    &                                 0,  0, -1,   & ! Bottom(B)
    &                                 1,  0,  0,   & ! East(E)
    &                                 0,  1,  0,   & ! North(N)
    &                                 0,  0,  1,   & ! Top(T)
    !                              central offset
    &                                 0,  0,  0 ], & ! Center(C)
    &                              [ 3, 7 ]        )

  integer, parameter :: d3q6_cxDir(3,6)            &
    !                              axis parallel offsets
    &                   = reshape( [ -1,  0,  0,   & ! West(W)
    &                                 0, -1,  0,   & ! South(S)
    &                                 0,  0, -1,   & ! Bottom(B)
    &                                 1,  0,  0,   & ! East(E)
    &                                 0,  1,  0,   & ! North(N)
    &                                 0,  0,  1 ], & ! Top(T)
    &                              [ 3, 6 ]        )

  ! cartesian directions of discrete velocities
  integer, parameter :: d2q9_cxDir(3,9)            &
    !                              axis parallel offsets
    &                   = reshape( [ -1,  0,  0,   & ! West(W)
    &                                 0, -1,  0,   & ! South(S)
    &                                 1,  0,  0,   & ! East(E)
    &                                 0,  1,  0,   & ! North(N)
    !                              diagonal offsets
    &                                -1, -1,  0,   & ! South West(SW)
    &                                -1,  1,  0,   & ! North West(NW)
    &                                 1, -1,  0,   & ! South East(SE)
    &                                 1,  1,  0,   & ! North East(NE)
    !                              central offset
    &                                 0,  0,  0 ], & ! Center(C)
    &                              [ 3, 9 ]        )

  ! cartesian directions of discrete velocities
  integer, parameter :: d2q5_cxDir(3,5)            &
    !                              axis parallel offsets
    &                   = reshape( [ -1,  0,  0,   & ! West(W)
    &                                 0, -1,  0,   & ! South(S)
    &                                 1,  0,  0,   & ! East(E)
    &                                 0,  1,  0,   & ! North(N)
    !                              central offset
    &                                 0,  0,  0 ], & ! Center(C)
    &                              [ 3, 5 ]        )

  ! cartesian directions of discrete velocities
  integer, parameter  :: d2q4_cxDir(3,4)            &
    !                               axis parallel offsets
    &                    = reshape( [ -1,  0,  0,   & ! West(W)
    &                                  0, -1,  0,   & ! South(S)
    &                                  1,  0,  0,   & ! East(E)
    &                                  0,  1,  0 ], & ! North(N)
    &                               [ 3,4 ]         )


  ! cartesian directions of discrete velocities
  integer, parameter  :: d1q3_cxDir(3,3)            &
    !                               axis parallel offsets
    &                    = reshape( [ -1,  0,  0,   & ! West(W)
    &                                  1,  0,  0,   & ! East(E)
    &                                  0,  0,  0 ], & ! Center(C)
    &                               [ 3, 3 ]        )

  ! cartesian directions of discrete velocities
  integer, parameter  :: d1q2_cxDir(3,2)            &
    !                               axis parallel offsets
    &                    = reshape( [ -1,  0,  0,   & ! West(W)
    &                                  1,  0,  0 ], & ! East(E)
    &                               [ 3,2 ]         )

  !> Stencil definitions
  !!
  !! A stencil is basically a set of element-offsets $(s_x, s_y, s_z)$,
  !! describing the relative positions of all required elements for
  !! a given element.
  type tem_stencilHeader_type

    !> a stencil label
    character(len=labelLen) :: label

    !> number of directions
    integer :: QQ

    !>@todo HK: do we really need this? It is mostly confusing, and solver
    !!          specific!
    !!
    !! number of directions excluding the central (0,0,0)
    integer :: QQN
    integer :: nDims = 0

    !> cartesian directions of discrete velocities.
    !! Integer number version.
    !! Size is (3,nDir)
    integer, allocatable :: cxDir(:,:)

    !> Cartesian directions of discrete velocities.
    !! Real number version.
    !! Size is (3,nDir)
    real(kind=rk), allocatable :: cxDirRK(:,:)

    !> inverted cartesian direction indices.
    !! this is well defined for symmetric stencils only.
    integer, allocatable :: cxDirInv(:)

    !> 2nd order tensor of discrete velocities
    !! Size is (6,QQ) for 3D: xx, yy, zz, xy, yz, xz
    !!         (3,QQ) for 2D: xx, yy, xy
    !!         (1,QQ) for 1D: xx
    real(kind=rk), allocatable :: cxcx(:,:)

    !> mapping of stencil entries to treelm definition, if possible
    !! entry is zero if nothing defined,
    !! @todo: SZ: is this really the case?
    !! -1 if a stencil entry with entries further than neighbor are encountered
    integer, allocatable :: map(:)

    !> mapping of treelm definition to stencil entries, if possible
    !! entry is zero if nothing defined
    integer, allocatable :: map2treeDef(:)

    !> Rest-density position in stencil
    integer :: restPosition

    !> the stencil on which the current one depends
    !! this is required for boundary stencils which elements actually require
    !! the neighbors of the compute stencil
    integer :: depStencil

    !> use this stencil for all elements
    logical :: useAll

    !> requires valid neighbors of the stencil neighbors
    logical :: requireNeighNeigh = .false.

    !> requires valid neighbors of the stencil neighbors
    logical :: requireAll = .false.

    !> the number of elements using this stencil
    integer :: nElems

    ! @todo: SZ: elemLvl and elem host more or less the same data maybe one
    !            of them could be removed.
    !            - elem is used in tem_construction (not in musubi, as far as
    !              I saw)
    !            - elemLvl is used in the tem_construction (not in musubi, as
    !              far as I saw)

    !> list of elements on which this stencil should be applied
    !! Both elemLvl and elem array are used for stencil
    !! other than fluid stencil
    type(grw_intArray_type), allocatable :: elemLvl(:)

    !> list of elements on which this stencil should be applied
    !! In build_BCstencil they point to original treeID list
    !! and used in tem_initelemLevels.
    !! Later this position is updated such that they point to
    !! the level wise total list in update_elemPosToTotalPos.
    ! KM: changed from allocatable array to growing array since
    ! musubi needs to create growing array of stencil headers
    ! and element list for each stencil header
    type(grw_intArray_type) :: elem
  end type tem_stencilHeader_type


  !> Element stencil definition
  type tem_stencilElement_type

    !> number of entries in pos
    integer :: QQN

    !> position in the tem_element_type::neighID%val( elemPos )
    !! Array size: QQN
    integer, allocatable :: tIDpos(:)

    !> the stencil on which the current one depends
    !! this is required for boundary stencils which elements actually require
    !! the neighbors of the compute stencil
    integer :: headerPos

    !> pointer to tem_element_type::tID( )
    !! For neighbor of an element in every stencil direction.
    !! It also inclues neighbor which is either fluid/halo/ghost.
    !! Array size: QQN
    !! Set in routine: build_levelElements, identify_stencilNeigh,
    !!                 identify_additionalNeigh
    integer, allocatable :: totalPos(:)

  end type tem_stencilElement_type


  ! declare the growing array of stencil type
  !> growing array type for type(tem_stencilelement_type)
  type grw_stencilelementarray_type
    integer :: nvals = 0
    integer :: containersize = 0
    type(tem_stencilelement_type), allocatable :: val(:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga_stencilelement
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_ga_stencilelement
  end interface

  !> empty the entries  without changing arrays
  interface empty
    module procedure empty_ga_stencilelement
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga_stencilelement
  end interface

  !> insert an element at a given position
  interface placeat
    module procedure placeat_ga_stencilelement
    module procedure placeat_ga_stencilelement_vec
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_ga_stencilelement
    module procedure append_ga_stencilelement_vec
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga_stencilelement
  end interface

  !> growing array type for type(tem_stencilheader_type)
  type grw_stencilheaderarray_type
    integer :: nvals = 0
    integer :: containersize = 0
    type(tem_stencilheader_type), allocatable :: val(:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga_stencilheader
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_ga_stencilheader
  end interface

  !> empty the entries  without changing arrays
  interface empty
    module procedure empty_ga_stencilheader
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga_stencilheader
  end interface

  !> insert an element at a given position
  interface placeat
    module procedure placeat_ga_stencilheader
    module procedure placeat_ga_stencilheader_vec
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_ga_stencilheader
    module procedure append_ga_stencilheader_vec
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga_stencilheader
  end interface


  interface init
    module procedure init_stencilHeader
    module procedure init_stencilElement
  end interface

  interface destroy
    module procedure empty_stencil
    module procedure destroy_stencilElement
  end interface

  interface empty
    module procedure empty_stencil
  end interface

  interface tem_stencil_dump
    module procedure tem_stencilHeader_dump
    module procedure tem_stencilElement_dump
  end interface

  interface assignment(=)
    module procedure copy_stencilElement
    module procedure copy_stencilHeader
  end interface


contains


! **************************************************************************** !
  !> Include the subroutines for the growing array.
  subroutine init_ga_stencilelement(me, length)
    type(grw_stencilelementarray_type), intent(out) :: me !< dynamic array to init
    integer, intent(in), optional :: length !< initial length of the container

    if (present(length)) then
      me%containersize = length
    else
      me%containersize = zerolength
    end if
    ! deallocate ...
    if( allocated( me%val ))     &
      deallocate(me%val)
    ! ... and reallocate
    allocate(me%val(me%containersize))
    me%nvals = 0

  end subroutine init_ga_stencilelement

  subroutine destroy_ga_stencilelement(me)
    type(grw_stencilelementarray_type), intent(inout) :: me !< dynamic array to destroy

    me%containersize = 0
    me%nvals = 0
    if( allocated( me%val ) ) deallocate(me%val)

  end subroutine destroy_ga_stencilelement


  subroutine truncate_ga_stencilelement(me)
    !------------------------------------------------------------------------
    type(grw_stencilelementarray_type) :: me !< array to truncate
    !------------------------------------------------------------------------
    type(tem_stencilelement_type), allocatable :: tarray(:)
    !------------------------------------------------------------------------
    integer :: ii
    !------------------------------------------------------------------------

    ! nothing to do if container size is not larger than the number of values
    ! in the array.
    if (me%containersize > me%nvals) then
      allocate(tarray(me%nvals))
      do ii = 1, me%nvals
        tarray(ii) = me%val(ii)
      end do
      call move_alloc(tarray, me%val)
      me%containersize = me%nvals
    end if

  end subroutine truncate_ga_stencilelement


  subroutine empty_ga_stencilelement(me)
    !------------------------------------------------------------------------
    type(grw_stencilelementarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------

    me%nvals = 0

  end subroutine empty_ga_stencilelement

  !> adds the value to a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! element at the requested position will be replaced.
  subroutine placeat_ga_stencilelement(me, val, pos, length)
    type(grw_stencilelementarray_type) :: me !< array to place the value into
    type(tem_stencilelement_type), intent(in) :: val !< value to place at the given position
    integer, intent(in) :: pos !< predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length


    ! value to append is larger than all existing ones,
    ! just put it to the end of the list, this captures
    ! also the case of empty lists.
    ! in this case foundpos = me%nvals + 1 holds.
    if (pos > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand(me = me, pos = pos, length = length)
    end if

    me%nvals = max( pos, me%nvals )
    me%val(pos) = val

  end subroutine placeat_ga_stencilelement


  !> adds the values starting from a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! elements starting from the requested position will be replaced up to
  !! the element at position `pos + size(val) - 1`.
  subroutine placeat_ga_stencilelement_vec(me, val, pos, length)
    type(grw_stencilelementarray_type) :: me !< array to append the value to
    type(tem_stencilelement_type), intent(in) :: val(:) !< values to append
    integer, intent(in) :: pos !< predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    integer :: ub, ii

    if (me%nvals == huge(me%nvals)) then
      write(*,*) "reached end of integer range for growing array!"
      write(*,*) "aborting!!"
      stop
    end if

    ub = pos + size(val) - 1

    if (ub > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand(me = me, pos = ub, length = length)
    end if

    me%nvals = max( ub, me%nvals )
    do ii = pos, ub
      me%val(ii) = val(1+ii-pos)
    end do

  end subroutine placeat_ga_stencilelement_vec


  subroutine append_ga_stencilelement(me, val, length)
    type(grw_stencilelementarray_type) :: me !< array to append the value to
    type(tem_stencilelement_type), intent(in) :: val !< value to append
    !> optional length to expand the array
    integer, intent(in), optional :: length

    ! value to append is larger than all existing ones,
    ! just put it to the end of the list, this captures
    ! also the case of empty lists.
    ! in this case foundpos = me%nvals + 1 holds.
    if (me%nvals+1 > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand(me = me, length = length)
    end if

    me%nvals = me%nvals+1
    me%val(me%nvals) = val

  end subroutine append_ga_stencilelement

  subroutine append_ga_stencilelement_vec(me, val, length)
    type(grw_stencilelementarray_type) :: me !< array to append the value to
    type(tem_stencilelement_type), intent(in) :: val(:) !< values to append
    !> optional length to expand the array
    integer, intent(in), optional :: length

    integer :: lb, ub, ii

    if (me%nvals == huge(me%nvals)) then
      write(*,*) "reached end of integer range for growing array!"
      write(*,*) "aborting!!"
      stop
    end if

    lb = me%nvals + 1
    ub = lb + size(val) - 1

    if (ub > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand(me = me, pos = ub, length = length)
    end if

    me%nvals = max( ub, me%nvals )
    do ii = lb, ub
      me%val(ii) = val(1+ii-lb)
    end do

  end subroutine append_ga_stencilelement_vec


  subroutine expand_ga_stencilelement(me, pos, length)
    type(grw_stencilelementarray_type) :: me !< array to resize
    integer, intent(in), optional :: pos !< optional predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    type(tem_stencilelement_type), allocatable :: swpval(:)
    integer :: explen, ii

    explen = 0
    ! increase the container by the requested length of double it
    if( present(length) ) then
      explen = max( length, minlength )
    else
      ! set the global minimum length, if doubling would be smaller than that
      explen = max(me%containersize, minlength)
    end if

    ! if a position is given, increase the container to at least the size to
    ! fit the position.
    if( present(pos) ) explen = max(explen, pos-me%containersize)

    ! if the current size plus explen exceeds the max container size,
    ! reduce the size to the max container size.
    if( (huge(me%containersize) - explen) <= me%containersize) then
      ! set max container size
      me%containersize = huge(me%containersize)
    else
      ! set the new container size
      me%containersize = me%containersize + explen
    end if

    if ( me%nvals > 0 ) then
      allocate(swpval(me%containersize))
      do ii = 1, me%nvals
        swpval(ii) = me%val(ii)
      end do
      call move_alloc( swpval, me%val )
    else ! me%nvals == 0
      if ( allocated(me%val) ) deallocate( me%val )
      allocate( me%val(me%containersize) )
    end if

  end subroutine expand_ga_stencilelement

  subroutine init_ga_stencilheader(me, length)
    type(grw_stencilheaderarray_type), intent(out) :: me !< dynamic array to init
    integer, intent(in), optional :: length !< initial length of the container

    if (present(length)) then
      me%containersize = length
    else
      me%containersize = zerolength
    end if
    ! deallocate ...
    if( allocated( me%val ))     &
      deallocate(me%val)
    ! ... and reallocate
    allocate(me%val(me%containersize))
    me%nvals = 0

  end subroutine init_ga_stencilheader

  subroutine destroy_ga_stencilheader(me)
    type(grw_stencilheaderarray_type), intent(inout) :: me !< dynamic array to destroy

    me%containersize = 0
    me%nvals = 0
    if( allocated( me%val ) ) deallocate(me%val)

  end subroutine destroy_ga_stencilheader


  subroutine truncate_ga_stencilheader(me)
    !------------------------------------------------------------------------
    type(grw_stencilheaderarray_type) :: me !< array to truncate
    !------------------------------------------------------------------------
    type(tem_stencilheader_type), allocatable :: tarray(:)
    !------------------------------------------------------------------------
    integer :: ii
    !------------------------------------------------------------------------

    ! nothing to do if container size is not larger than the number of values
    ! in the array.
    if (me%containersize > me%nvals) then
      allocate(tarray(me%nvals))
      do ii = 1, me%nvals
        tarray(ii) = me%val(ii)
      end do
      call move_alloc(tarray, me%val)
      me%containersize = me%nvals
    end if

  end subroutine truncate_ga_stencilheader


  subroutine empty_ga_stencilheader(me)
    !------------------------------------------------------------------------
    type(grw_stencilheaderarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------

    me%nvals = 0

  end subroutine empty_ga_stencilheader

  !> adds the value to a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! element at the requested position will be replaced.
  subroutine placeat_ga_stencilheader(me, val, pos, length)
    type(grw_stencilheaderarray_type) :: me !< array to place the value into
    type(tem_stencilheader_type), intent(in) :: val !< value to place at the given position
    integer, intent(in) :: pos !< predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length


    ! value to append is larger than all existing ones,
    ! just put it to the end of the list, this captures
    ! also the case of empty lists.
    ! in this case foundpos = me%nvals + 1 holds.
    if (pos > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand(me = me, pos = pos, length = length)
    end if

    me%nvals = max( pos, me%nvals )
    me%val(pos) = val

  end subroutine placeat_ga_stencilheader


  !> adds the values starting from a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! elements starting from the requested position will be replaced up to
  !! the element at position `pos + size(val) - 1`.
  subroutine placeat_ga_stencilheader_vec(me, val, pos, length)
    type(grw_stencilheaderarray_type) :: me !< array to append the value to
    type(tem_stencilheader_type), intent(in) :: val(:) !< values to append
    integer, intent(in) :: pos !< predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    integer :: ub, ii

    if (me%nvals == huge(me%nvals)) then
      write(*,*) "reached end of integer range for growing array!"
      write(*,*) "aborting!!"
      stop
    end if

    ub = pos + size(val) - 1

    if (ub > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand(me = me, pos = ub, length = length)
    end if

    me%nvals = max( ub, me%nvals )
    do ii = pos, ub
      me%val(ii) = val(1+ii-pos)
    end do

  end subroutine placeat_ga_stencilheader_vec


  subroutine append_ga_stencilheader(me, val, length)
    type(grw_stencilheaderarray_type) :: me !< array to append the value to
    type(tem_stencilheader_type), intent(in) :: val !< value to append
    !> optional length to expand the array
    integer, intent(in), optional :: length

    ! value to append is larger than all existing ones,
    ! just put it to the end of the list, this captures
    ! also the case of empty lists.
    ! in this case foundpos = me%nvals + 1 holds.
    if (me%nvals+1 > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand(me = me, length = length)
    end if

    me%nvals = me%nvals+1
    me%val(me%nvals) = val

  end subroutine append_ga_stencilheader

  subroutine append_ga_stencilheader_vec(me, val, length)
    type(grw_stencilheaderarray_type) :: me !< array to append the value to
    type(tem_stencilheader_type), intent(in) :: val(:) !< values to append
    !> optional length to expand the array
    integer, intent(in), optional :: length

    integer :: lb, ub, ii

    if (me%nvals == huge(me%nvals)) then
      write(*,*) "reached end of integer range for growing array!"
      write(*,*) "aborting!!"
      stop
    end if

    lb = me%nvals + 1
    ub = lb + size(val) - 1

    if (ub > me%containersize) then
      ! expand the array, if its boundary is reached
      call expand(me = me, pos = ub, length = length)
    end if

    me%nvals = max( ub, me%nvals )
    do ii = lb, ub
      me%val(ii) = val(1+ii-lb)
    end do

  end subroutine append_ga_stencilheader_vec


  subroutine expand_ga_stencilheader(me, pos, length)
    type(grw_stencilheaderarray_type) :: me !< array to resize
    integer, intent(in), optional :: pos !< optional predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    type(tem_stencilheader_type), allocatable :: swpval(:)
    integer :: explen, ii

    explen = 0
    ! increase the container by the requested length of double it
    if( present(length) ) then
      explen = max( length, minlength )
    else
      ! set the global minimum length, if doubling would be smaller than that
      explen = max(me%containersize, minlength)
    end if

    ! if a position is given, increase the container to at least the size to
    ! fit the position.
    if( present(pos) ) explen = max(explen, pos-me%containersize)

    ! if the current size plus explen exceeds the max container size,
    ! reduce the size to the max container size.
    if( (huge(me%containersize) - explen) <= me%containersize) then
      ! set max container size
      me%containersize = huge(me%containersize)
    else
      ! set the new container size
      me%containersize = me%containersize + explen
    end if

    if ( me%nvals > 0 ) then
      allocate(swpval(me%containersize))
      do ii = 1, me%nvals
        swpval(ii) = me%val(ii)
      end do
      call move_alloc( swpval, me%val )
    else ! me%nvals == 0
      if ( allocated(me%val) ) deallocate( me%val )
      allocate( me%val(me%containersize) )
    end if

  end subroutine expand_ga_stencilheader


  ! ------------------------------------------------------------------------ !
  !> Load the stencil configuration from the lua file
  subroutine tem_loadStencil( stencil, parent_handle, conf )
    ! -------------------------------------------------------------------- !
    !> stencil type to fill
    type(tem_stencilHeader_type) :: stencil
    !> handle of the parent table
    integer, intent(in)    :: parent_handle
    !> lua state type
    type(flu_State) :: conf
    ! -------------------------------------------------------------------- !
    integer :: nCartDiscVel, iCartDiscVel
    integer :: nDiscreteVel, iDiscreteVel
    integer :: discVel_handle
    integer :: discVel_cart_handle
    integer :: iError
    ! -------------------------------------------------------------------- !
    ! Usually, the stencil to read applies for all elements
    ! (Compute stencil)
    stencil%useAll = .true.

    ! read the number of discrete velocities QQ
    ! @todo{ MH: Do we really need to read QQ explicitly?
    ! Isn't it enough to just count the number of entries below?}
    call aot_get_val( L       = conf,          &
      &               thandle = parent_handle, &
      &               val     = stencil%QQ,    &
      &               ErrCode = iError,        &
      &               key     = 'QQ'           )

    ! read the discrete velocities
    call aot_table_open( L       = conf,           &
      &                  parent  = parent_handle,  &
      &                  thandle = discVel_handle, &
      &                  key     = 'disc_vel'      )

    ! Check the length of the discrete velocity table
    nDiscreteVel = aot_table_length( L       = conf,          &
      &                              thandle = discVel_handle )

    if (nDiscreteVel == stencil%QQ ) then
      ! If the number of discrete velocities matches the given stencil QQ
      call init_stencilHeader( me     = stencil,     &
        &                      QQN    = stencil%QQN, &
        &                      QQ     = stencil%QQ,  &
        &                      useAll = .true.       )

      do iDiscreteVel = 1, nDiscreteVel

        call aot_table_open( L       = conf,                &
          &                  parent  = discVel_handle,      &
          &                  thandle = discVel_cart_handle, &
          &                  pos     = iDiscreteVel         )

        nCartDiscVel = aot_table_length( L       = conf,               &
          &                              thandle = discVel_cart_handle )

        if (nCartDiscVel == 3) then

          do iCartDiscVel = 1, nCartDiscVel
            call aot_get_val(                                                  &
              &         L       = conf,                                        &
              &         thandle = discVel_cart_handle,                         &
              &         val     = stencil%cxDir( iCartDiscVel, iDiscreteVel ), &
              &         ErrCode = iError,                                      &
              &         pos     = iCartDiscVel                                 )

            stencil%cxDirRK(iCartDiscVel, iDiscreteVel) &
              & = real( stencil%cxDir(iCartDiscVel, iDiscreteVel), rk )
          end do ! iCartDiscVel

        else

          write(logUnit(1),*) "Not every velocity has three cartesian " &
            &                 //"coordinates!!"

        end if

      end do ! iDiscreteVel

      call aot_table_close( L       = conf,               &
        &                   thandle = discVel_cart_handle )

    else
      write(logUnit(1),*) "The number of discrete velocities for the " &
        &                 //"new stencil does not match QQ!!"
      STOP
    end if

  end subroutine tem_loadStencil
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Delete all entries in the stencil
  subroutine empty_stencil( me )
    ! -------------------------------------------------------------------- !
    !> stencil header type to empty
    type( tem_stencilHeader_type ), intent(out) :: me
    ! -------------------------------------------------------------------- !
    me%QQ  = 0
    me%QQN = 0
    me%nDims = 0
    if( allocated( me%cxDir )) deallocate( me%cxDir )
    if( allocated( me%cxDirRK )) deallocate( me%cxDirRK )
    if( allocated( me%cxDirInv )) deallocate( me%cxDirInv )
    if( allocated( me%cxcx )) deallocate( me%cxcx )
    if( allocated( me%map )) deallocate( me%map )
    if ( me%elem%nVals > 0 ) call destroy( me%elem )
    if( allocated( me%elemLvl )) deallocate( me%elemLvl )

  end subroutine empty_stencil
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This subroutine fills the array of inverse directions according to the
  !! given array of directions.
  !!
  subroutine tem_identify_inverseDirections( me, cxDir )
    ! -------------------------------------------------------------------- !
    !> array of positions of the inverse offsets in cxDir
    integer, intent(out) :: me(:)
    !> array of given offsets
    integer, intent(in)  :: cxDir(:,:)
    ! -------------------------------------------------------------------- !
    ! counter variable
    integer :: iOffset, iCompare, QQ
    ! -------------------------------------------------------------------- !

    QQ = size(me)

    do iOffset=1,QQ
      do iCompare=1,QQ

        if ( (( -cxDir(1, iOffset) ) == ( cxDir(1, iCompare) ))       &
          &  .and. (( -cxDir(2, iOffset) ) == ( cxDir(2, iCompare) )) &
          &  .and. (( -cxDir(3, iOffset) ) == ( cxDir(3, iCompare) )) ) then
          me( iOffset ) = iCompare
          EXIT
        end if

      end do ! iCompare
    end do ! iOffset

  end subroutine tem_identify_inverseDirections
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This subroutine fills the array of prevailing directions according to the
  !! given array of directions.
  !!
  subroutine tem_identify_prevailDirections( me, cxDir )
    ! -------------------------------------------------------------------- !
    !> growing array of prevailing directions
    real( kind=rk ), allocatable, intent(out) :: me(:,:)
    !> array of given offsets
    integer, intent(in)  :: cxDir(:,:)
    ! -------------------------------------------------------------------- !
    real(kind=rk),allocatable :: tmpPrevailDir(:,:)
    real(kind=rk) :: length_rk
    integer :: length
    ! counters
    integer :: iOffset
    integer :: nEntries ! size of array
    integer :: counter
    ! -------------------------------------------------------------------- !

    nEntries = size(cxDir, 2)
    allocate( tmpPrevailDir(3, nEntries) )
    counter = 0

    do iOffset=1,nEntries

      length =   cxDir(1,iOffset)**2 &
        &      + cxDir(2,iOffset)**2 &
        &      + cxDir(3,iOffset)**2

      if (length /= 0) then
        counter = counter + 1
        length_rk = sqrt(real(length, kind=rk))
        tmpPrevailDir(1:3, counter) = real(cxDir(1:3, iOffset), kind=rk) &
          &                          / length_rk
      end if

    end do

    allocate( me( 3, counter ))

    do iOffset = 1, counter
      me(:, iOffset) = tmpPrevailDir( :, iOffset)
    end do

  end subroutine tem_identify_prevailDirections
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Find the index of the given direction in the list of all directions.
  !! -1 if not found.
  pure function tem_stencil_findIndexOfDir(findDir, cxDir) result(idx)
    ! -------------------------------------------------------------------- !
    !> vector index to find in cxDir
    integer, intent(in) :: findDir(3)
    !> array of vectors
    integer, intent(in)  :: cxDir(:,:)
    integer :: idx
    ! -------------------------------------------------------------------- !
    integer :: iDir
    ! -------------------------------------------------------------------- !

    idx = -1

    do iDir = 1, size(cxDir, 2)

      if ( cxDir(1, iDir) == findDir(1)       &
        &  .and. cxDir(2, iDir) == findDir(2) &
        &  .and. cxDir(3, iDir) == findDir(3) ) then
        idx = iDir
      end if

    end do

  end function tem_stencil_findIndexOfDir
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Map the stencil offsets to the internally defined treelm offsets.
  !!
  subroutine tem_stencil_map_toTreelmDef( me )
    ! -------------------------------------------------------------------- !
    !> stencil to map
    type(tem_stencilHeader_type) :: me
    ! -------------------------------------------------------------------- !
    ! counters
    integer :: iQQQ, iDir
    ! -------------------------------------------------------------------- !

    if ( allocated(me%map) ) deallocate(me%map)

    allocate( me%map(me%QQ) )
    me%map = 0

    do iDir = 1,me%QQ

      if ( maxval( abs( me%cxDir(:,iDir) ) ) <= 1 ) then

        ! Now match the stencil offset with the treelm definition offsets
        qqqLoop: do iQQQ=1,qQQQ
          if ( me%cxDir(1, iDir) == qOffset(iQQQ, 1)      &
            & .and. me%cxDir(2, iDir) == qOffset(iQQQ, 2) &
            & .and. me%cxDir(3, iDir) == qOffset(iQQQ, 3) ) then
            me%map( iDir ) = iQQQ
            exit qqqLoop
          end if
        end do qqqLoop

      end if

    end do ! iDir = 1, QQ

  end subroutine tem_stencil_map_toTreelmDef
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Map the internally defined treelm offsets to the stencil offsets.
  !!
  subroutine tem_treelmDef_map_toStencil( me )
    ! -------------------------------------------------------------------- !
    !> stencil to map
    type(tem_stencilHeader_type) :: me
    ! -------------------------------------------------------------------- !
    ! counters
    integer :: iDir, iQQQ
    ! -------------------------------------------------------------------- !

    ! allocate the map2treeDef with all 26 direct neighbors
    if ( allocated(me%map2treeDef) ) deallocate(me%map2treeDef)

    allocate( me%map2treeDef(qQQQ) )
    me%map2treeDef = 0

    ! loop over all 26 neighbors
    do iQQQ=1,qQQQ

      ! search for the correct entry in the stencil%cxDir array
      do iDir = 1, me%QQ

        if ( me%cxDir(1, iDir) == qOffset(iQQQ, 1)       &
          &  .and. me%cxDir(2, iDir) == qOffset(iQQQ, 2) &
          &  .and. me%cxDir(3, iDir) == qOffset(iQQQ, 3) ) then
          ! assign the correct position
          me%map2treeDef(iQQQ) = iDir
          EXIT
        end if

      end do

    end do

  end subroutine tem_treelmDef_map_toStencil
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Write element information to disk
  !!
  subroutine tem_stencilHeader_dump( me, nUnit )
    ! -------------------------------------------------------------------- !
    !> stencil header to be dumped
    type( tem_stencilHeader_type ), intent(in) :: me
    integer, optional, intent(in) :: nUnit
    ! -------------------------------------------------------------------- !
    integer :: iDir
    integer :: dmpU
    ! -------------------------------------------------------------------- !

    ! Default to write to logUnit
    dmpU = logUnit(6)

    if (present(nUnit)) then
      dmpU = nUnit
    end if

    write(dmpU, *) '   Stencil  '
    write(dmpU, *) '     QQ   : ', me%QQ
    write(dmpU, *) '     QQN  : ', me%QQN

    if ( allocated(me%cxDir) ) then
      write(dmpU, *) '     Offset directions:'

      do iDir=1,me%QQN
        if ( allocated(me%cxDir) ) then
          write(dmpU, *) '       cxDir(:,', iDir, ') = ', &
            &            me%cxDir(:,iDir)
        end if
      end do

    end if

    write(dmpU, *) '     useAll: ', me%useAll

    if ( .not. me%useAll ) then
      write(dmpU, *) '     nElems: ', me%nElems

      if ( me%elem%nVals > 0 )                              &
        write(dmpU, *) '     Elems : '                      &
          &            //trim( tem_toStr(val = me%elem%val, &
          &                              sep = ','        ) )
    end if

  end subroutine tem_stencilHeader_dump
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Write element information to disk
  !!
  subroutine tem_stencilElement_dump( me, nUnit, neighID, tIDonly )
    ! -------------------------------------------------------------------- !
    !> stencil element to be dumped
    type( tem_stencilElement_type ), intent(in) :: me
    !> unit to dump to
    integer, intent(in) :: nUnit
    !> neighbor ID
    integer(kind=long_k), allocatable, intent(in), optional :: neighID(:)
    !> only use tree IDs
    logical, intent(in), optional :: tIDonly
    ! -------------------------------------------------------------------- !
    integer :: iVal
    logical :: tIDonly_loc
    integer(kind=long_k) :: treeIDs( me%QQN )
    ! -------------------------------------------------------------------- !

    if( present( tIDonly )) then
      tIDonly_loc = tIDonly
    else
      tIDonly_loc = .false.
    endif

    write(nUnit, "(' QQN: ', I2, ', headerPos: ', I0 )" ) me%QQN, me%headerPos

    if ( present(neighID) ) then

      if ( allocated(neighID) .and. allocated(me%tIDpos) ) then
        ! get neighIDs
        do iVal=1,size(me%tIDpos)
          if ( me%tIDpos(iVal) == 0 ) then
            treeIDs(iVal) = 0
          else
            treeIDs(iVal) = neighID( me%tIDpos(iVal) )
          end if
        end do
        call tem_printArray( treeIDs, 8, 'neighIDs', 120, nUnit)
      end if

    end if

  end subroutine tem_stencilElement_dump
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Definition of the d3q125 neighborhood
  pure function d3q125_cxDir() result(cxdir)
    ! -------------------------------------------------------------------- !
    integer :: cxDir(3,125)
    ! -------------------------------------------------------------------- !
    integer :: iX, iY, iZ
    integer :: dirPos
    integer :: length
    ! -------------------------------------------------------------------- !

    ! save the cartesian directions of the discrete velocities
    dirPos = 0
    do iX=-2,2
      do iY=-2,2
        do iZ=-2,2
          ! calculate the square of the length
          length = iX*iX + iY*iY + iZ*iZ
          ! if the actual length is smaller than 2.5
          ! this is the same as length^2 .lt. 8 for integer stencil directions
          if (length > 0) then
            dirPos = dirPos + 1
            cxDir(:, dirPos) = [ iX, iY, iZ ]
          end if
        end do
      end do
    end do

    ! add the rest direction at the last position
    cxDir(:, 125) = [0, 0, 0]

  end function d3q125_cxDir
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Definition of the d3q81 neighborhood.
  pure function d3q81_cxDir() result(cxdir)
    ! -------------------------------------------------------------------- !
    integer :: cxDir(3,81)
    ! -------------------------------------------------------------------- !
    integer :: iX, iY, iZ
    integer :: dirPos
    integer :: length
    ! -------------------------------------------------------------------- !

    ! save the cartesian directions of the discrete velocities
    dirPos = 0

    do iX=-2,2
      do iY=-2,2
        do iZ=-2,2

          ! calculate the square of the length
          length = iX*iX + iY*iY + iZ*iZ

          ! if the actual length is smaller than 2.5
          ! this is the same as length^2 .lt. 8 for integer stencil directions
          if ((length > 0) .and. (length < 8)) then
            dirPos = dirPos + 1
            cxDir(:, dirPos) = [ iX, iY, iZ ]
          end if

        end do
      end do
    end do

    ! add the rest direction at the last position
    cxDir(:, 81) = [0, 0, 0]

  end function d3q81_cxDir
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This subroutine creates the required stencil.
  !!
  subroutine tem_create_stencil( stencil, stencilKind )
    ! -------------------------------------------------------------------- !
    !> stencil type to be defined
    type( tem_stencilHeader_type), intent(out) :: stencil
    !> stencil kind to decide which create function to call
    character(len=*) :: stencilKind
    ! -------------------------------------------------------------------- !

    write(logUnit(1),*) 'Creating stencil kind: '//trim(stencilKind)

    select case (trim(stencilKind))
      case ('d3q19')
        call init( me     = stencil,    &
          &        QQN    = 18,         &
          &        QQ     = 19,         &
          &        useAll = .true.,     &
          &        nDims  = 3,          &
          &        label  = 'd3q19',    &
          &        cxDir  = d3q19_cxDir )
      case ('d3q13')
        call init( me     = stencil,    &
          &        QQN    = 12,         &
          &        QQ     = 13,         &
          &        useAll = .true.,     &
          &        nDims  = 3,          &
          &        label  = 'd3q13',    &
          &        cxDir  = d3q13_cxDir )
      case ('d3q27')
        call init( me     = stencil,    &
          &        QQN    = 26,         &
          &        QQ     = 27,         &
          &        useAll = .true.,     &
          &        nDims  = 3,          &
          &        label  = 'd3q27',    &
          &        cxDir  = d3q27_cxDir )
      case ('d3q7')
        call init( me     = stencil,   &
          &        QQN    = 6,         &
          &        QQ     = 7,         &
          &        useAll = .true.,    &
          &        nDims  = 3,         &
          &        label  = 'd3q7',    &
          &        cxDir  = d3q7_cxDir )
      case ('d3q6', 'flekkoy')
        call init( me     = stencil,   &
          &        QQN    = 6,         &
          &        QQ     = 6,         &
          &        useAll = .false.,   &
          &        nDims  = 3,         &
          &        label  = 'd3q6',    &
          &        cxDir  = d3q6_cxDir )
      case ('d2q9')
        call init( me     = stencil,   &
          &        QQN    = 8,         &
          &        QQ     = 9,         &
          &        useAll = .true.,    &
          &        nDims  = 2,         &
          &        label  = 'd2q9',    &
          &        cxDir  = d2q9_cxDir )
      case ('d2q5')
        call init( me     = stencil,   &
          &        QQN    = 4,         &
          &        QQ     = 5,         &
          &        useAll = .true.,    &
          &        nDims  = 2,         &
          &        label  = 'd2q5',    &
          &        cxDir  = d2q5_cxDir )
      case ('d1q3')
        call init( me     = stencil,   &
          &        QQN    = 2,         &
          &        QQ     = 3,         &
          &        useAll = .true.,    &
          &        nDims  = 1,         &
          &        label  = 'd1q3',    &
          &        cxDir  = d1q3_cxDir )
      case ('default')
        write(logUnit(1),*) "Stencil kind Not defined yet!!"
        call tem_abort()
    end select

  end subroutine tem_create_stencil
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Identify the zero-position in the stencil
  !! Return -1 if not found in stencil
  !!
  function tem_stencil_zeroPos( me ) result( pos )
    ! -------------------------------------------------------------------- !
    !> stencil to map
    type(tem_stencilHeader_type),intent(in) :: me
    !> Position of zero-entry in the stencil
    integer :: pos
    ! -------------------------------------------------------------------- !
    ! counter
    integer :: iQQ
    ! -------------------------------------------------------------------- !
    pos = -1
    do iQQ = 1, me%QQ
      if( me%cxDir( 1, iQQ ) == 0 .and.  &
        & me%cxDir( 2, iQQ ) == 0 .and.  &
        & me%cxDir( 3, iQQ ) == 0 ) then
        pos = iQQ
        exit
      end if
    end do
  end function tem_stencil_zeroPos
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> initialize stencil
  !!
  subroutine init_stencilHeader( me, QQN, QQ, nElems, useAll, &
    &                            nDims, label, cxDir          )
    ! -------------------------------------------------------------------- !
    !> stencil header to be initialized
    type( tem_stencilHeader_type ), intent(out) :: me
    !> number of discrete velocities in the model (without the center one)
    integer, intent(in) :: QQN
    !> number of discrete velocities in the model (incl. the center one)
    integer, intent(in), optional :: QQ
    !> The number of elements to use this stencil for
    integer, intent(in), optional :: nElems
    !> use this stencil for all elements?
    logical, intent(in), optional :: useAll

    ! Number of dimensions considered by the stencil
    integer, intent(in), optional :: nDims

    ! A label to describe the stencil
    character(len=*), intent(in), optional :: label

    ! Directions describing the neighboring elements in the stencil
    integer, intent(in), optional :: cxDir(:,:)
    ! -------------------------------------------------------------------- !
    integer :: iElem
    ! -------------------------------------------------------------------- !

    me%QQN = QQN
    if( present( QQ )) me%QQ  = QQ
    if( present( useAll )) then
      me%useAll = useAll
      me%nElems = 0
    else if( present( nElems ) ) then
      me%useAll = .false.
      me%nElems = nElems
      call init( me%elem, nElems )
      allocate( me%elemLvl( nElems ) )
      do iElem = 1, nElems
        call init( me%elemLvl(iElem) )
      end do
    else
      me%useAll = .false.
      me%nElems = 0
    end if

    if ( .not. allocated(me%cxDir) )    allocate( me%cxDir(   3, me%QQ ) )
    if ( .not. allocated(me%cxDirRK) )  allocate( me%cxDirRK( 3, me%QQ ) )
    if ( .not. allocated(me%cxDirInv) ) allocate( me%cxDirInv(   me%QQ ) )
    if ( .not. allocated(me%cxcx) )     allocate( me%cxcx(    6, me%QQ ) )

    if (present(nDims)) me%nDims = nDims
    if (present(label)) me%label = label

    if (present(cxDir)) then
      me%cxDir = cxDir
      me%cxDirRK(:,:) = real( cxDir(:,:), rk )
      call tem_identify_inverseDirections( me%cxDirInv, me%cxDir )
      call tem_stencil_createCxcx( me )
    else
      me%cxDir    = 0
      me%cxDirRK  = 0.0_rk
      me%cxDirInv = 0
      me%cxcx     = 0.0_rk
    end if

  end subroutine init_stencilHeader
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> initialize stencil
  !!
  subroutine init_stencilElement( me, QQN, headerPos, tIDpos )
    ! -------------------------------------------------------------------- !
    !> stencil element type to be initialized
    type( tem_stencilElement_type ), intent(out) :: me
    !> number of discrete velocities in the model (without the center one)
    integer, intent(in) :: QQN
    !>
    integer, intent(in), optional :: headerPos
    !>
    integer, intent(in), optional :: tIDpos(:)
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    me%QQN = QQN
    if ( allocated(me%tIDpos) ) deallocate(me%tIDpos)
    allocate( me%tIDpos(me%QQN) )

    if ( present(tIDpos) ) then
      me%tIDpos  = tIDpos
    else
      me%tIDpos  = 0
    end if

    me%headerPos = 0

    if ( present(headerPos) ) me%headerPos = headerPos

    allocate( me%totalPos(me%QQN) )
    me%totalPos = 0

  end subroutine init_stencilElement
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> destroy stencil
  !!
  subroutine destroy_stencilElement( me )
    ! -------------------------------------------------------------------- !
    !> stencil element type to be destroyed
    type( tem_stencilElement_type ), intent(out) :: me
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    me%QQN = 0
    if ( allocated(me%tIDpos) ) deallocate(me%tIDpos)

    me%headerPos = 0

    if( allocated( me%totalPos ) ) deallocate( me%totalPos )

  end subroutine destroy_stencilElement
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> find position stencil
  !!
  function tem_stencil_getHeaderPos( me, val ) result( headerPos )
    ! -------------------------------------------------------------------- !
    !> List of element stencil definitions
    type( grw_stencilElementArray_type ), intent(in) :: me
    !> Stencil to search for
    integer, intent(in) :: val
    !> Position of the header in the list, 0 if not found
    integer :: headerPos
    ! -------------------------------------------------------------------- !
    integer :: iStencil
    ! -------------------------------------------------------------------- !

    headerPos = 0

    ! Run over all the stencils
    do iStencil = 1, me%nVals
      if ( me%val(iStencil)%headerPos == val ) then
        headerPos = iStencil
      end if
    end do

  end function tem_stencil_getHeaderPos
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function returns a unique label for given stencil cxDir
  function tem_stencil_getLabelForcxDir( me, prevailDir ) result (uLabel)
    ! -------------------------------------------------------------------- !
    !> stencil
    type(tem_stencilHeader_type), intent(in) :: me
    !> prevail directions
    real(kind=rk), intent(in) :: prevailDir(:,:)
    !> return unique label
    character(len=labelLen) :: uLabel
    ! -------------------------------------------------------------------- !
    integer :: vector(3)
    integer :: iDir, iDir_tem, length
    character(len=4) :: buffer
    ! -------------------------------------------------------------------- !

    ! if stencil label is one of below then set label to stencil label
    select case(trim(me%label))
    case ( 'd3q19', 'd3q27', 'd3q13', 'd3q7', 'd3q6', 'flekkoy', &
      &    'd2q9', 'd2q5',                                       &
      &    'd1q3',                                               &
      &    'new_stencil'                                         )
      uLabel = trim(me%label)

    case default
      uLabel = ''
      ! loop over each direction of stencil
      do iDir = 1, me%QQN

        vector = me%cxDir(:, iDir)

        ! if length of vector is > 1 append length to label
        length = ceiling( sqrt( dot_product( real(vector, kind=rk), &
          &                                  real(vector, kind=rk)  &
          &                                )                        &
          &                    )                                    )

        ! determine discrete vector from prevailDir
        call tem_determine_discreteVector(vector, prevailDir)
        ! find direction label
        do iDir_tem=1,qQQQ
          if ( all(vector == qOffset(iDir_tem,:)) ) then
            uLabel = trim(uLabel)//trim(adjustl(qDirName(iDir_tem)))
          end if
        end do

        ! append length
        if (length > 1) then
          write(buffer,'(i4)') length
          uLabel = trim(uLabel)//trim(adjustl(buffer))
        end if

        ! if cxDir is 0 then set label to restPosition
        if ( all(me%cxDir(:, iDir) == 0 ) ) then
          uLabel = 'restPosition'
        end if
      end do

      write(buffer,'(i4)') me%QQN
      uLabel = trim(uLabel)//'_QQN'//trim(adjustl(buffer))

    end select

  end function tem_stencil_getLabelForcxDir
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Create cxcx for a given stencil
  !!
  subroutine tem_stencil_createCxcx( me )
    ! -------------------------------------------------------------------- !
    !> stencil
    type(tem_stencilHeader_type) :: me
    ! -------------------------------------------------------------------- !
    me%cxcx = 0.0_rk

    ! allocate cxcx according to dimension and compute cxcx using cxDir
    select case(me%nDims)
    case (3)
      ! create cxcx for 3D, xx, yy, zz, xy, yz, xz
      me%cxcx(1,:) = real(me%cxDir(1,:) * me%cxDir(1,:), rk)
      me%cxcx(2,:) = real(me%cxDir(2,:) * me%cxDir(2,:), rk)
      me%cxcx(3,:) = real(me%cxDir(3,:) * me%cxDir(3,:), rk)
      me%cxcx(4,:) = real(me%cxDir(1,:) * me%cxDir(2,:), rk)
      me%cxcx(5,:) = real(me%cxDir(2,:) * me%cxDir(3,:), rk)
      me%cxcx(6,:) = real(me%cxDir(1,:) * me%cxDir(3,:), rk)

    case (2)
      ! create cxcx for 2D, xx, yy, xy
      me%cxcx(1,:) = real(me%cxDir(1,:) * me%cxDir(1,:), rk)
      me%cxcx(2,:) = real(me%cxDir(2,:) * me%cxDir(2,:), rk)
      me%cxcx(3,:) = real(me%cxDir(1,:) * me%cxDir(2,:), rk)

    case (1)
      ! create cxcx for 1D, xx
      me%cxcx(1,:) = real(me%cxDir(1,:) * me%cxDir(1,:), rk)

    case default
      write( logUnit(1), "(A)" ) " Can't create cxcx for stencil " &
        &                        // trim(me%label)                 &
        &                        // ", as dimension not defined!"

    end select

  end subroutine tem_stencil_createCxcx
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function provides copy assigment for tem_stencilHeader_type
  subroutine copy_stencilHeader(left, right)
    ! -------------------------------------------------------------------- !
    !> tem_stencilHeader to copy to
    type(tem_stencilHeader_type), intent(out) :: left
    !> tem_stencilHeader to copy from
    type(tem_stencilHeader_type), intent(in) :: right
    ! -------------------------------------------------------------------- !

    left%label = right%label
    left%QQ = right%QQ
    left%QQN = right%QQN
    left%nDims = right%nDims

    if ( allocated(right%cxDir) ) then
      allocate(left%cxDir(size(right%cxDir,1),size(right%cxDir,2)))
      left%cxDir = right%cxDir
    end if

    if ( allocated(right%cxDirRK) ) then
      allocate(left%cxDirRK(size(right%cxDirRK,1),size(right%cxDirRK,2)))
      left%cxDirRK = right%cxDirRK
    end if

    if ( allocated(right%cxDirInv) ) then
      allocate(left%cxDirInv(size(right%cxDirInv)))
      left%cxDirInv = right%cxDirInv
    end if

    if ( allocated(right%cxcx) ) then
      allocate(left%cxcx(size(right%cxcx,1),size(right%cxcx,2)))
      left%cxcx = right%cxcx
    end if

    if ( allocated(right%map) ) then
      allocate(left%map(size(right%map)))
      left%map = right%map
    end if

    if ( allocated(right%map2treeDef) ) then
      allocate(left%map2treeDef(size(right%map2treeDef)))
      left%map2treeDef = right%map2treeDef
    end if

    left%restPosition      = right%restPosition
    left%depStencil        = right%depStencil
    left%useAll            = right%useAll
    left%requireNeighNeigh = right%requireNeighNeigh
    left%requireAll        = right%requireAll

    left%nElems = right%nElems

    if ( allocated(right%elemLvl) ) then
      allocate( left%elemLvl(lbound(right%elemLvl,1):ubound(right%elemLvl,1)) )
      left%elemLvl = right%elemLvl
    end if

    left%elem = right%elem
  end subroutine copy_stencilHeader
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function provides copy assigment for tem_stencilElement_type
  subroutine copy_stencilElement(left, right)
    ! -------------------------------------------------------------------- !
    !> tem_stencilElement to copy to
    type(tem_stencilElement_type), intent(out) :: left
    !> tem_stencilElement to copy from
    type(tem_stencilElement_type), intent(in) :: right
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    left%QQN = right%QQN
    if ( allocated(right%tIDpos) ) then
      allocate(left%tIDpos(size(right%tIDpos)))
      left%tIDpos = right%tIDpos
    end if

    left%headerPos = right%headerPos

    if ( allocated(right%totalPos) ) then
      allocate(left%totalPos(size(right%totalPos)))
      left%totalPos = right%totalPos
    end if

  end subroutine copy_stencilElement
  ! ------------------------------------------------------------------------ !


end module tem_stencil_module
! **************************************************************************** !
