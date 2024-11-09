! Copyright (c) 2016, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016-2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2019 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> This module contains data type and subroutines required to store levelWise
!! points and values for variable access
!!
!! \author Kannan Masilamani
module tem_pointData_module

  ! include treelm modules
  use env_module,            only: rk, long_k, globalMaxLevels
  use treelmesh_module,      only: treelmesh_type
  use tem_grow_array_module, only: grw_realArray_type,    &
    &                              grw_intArray_type,     &
    &                              grw_charArray_type,    &
    &                              init, append,          &
    &                              truncate, destroy
  use tem_dyn_array_module,  only: dyn_longArray_type, append
  use tem_logging_module,    only: logUnit
  use tem_aux_module,        only: tem_abort
  use tem_geometry_module,   only: tem_CoordOfReal
  use tem_topology_module,   only: tem_IdOfCoord

  implicit none

  private

  public :: tem_pointData_type, tem_pointData_list_type
  public :: init, append, truncate, destroy
  public :: tem_grwPoints_type
  public :: tem_sourceElems_type

  !> Data type contain 1D growing array of points for each dimension.
  !! In general, they are space coordinate in the treelmesh but
  !! for Ateles solver variables, they are used to store local coordinate
  !! with in a reference element.
  !!
  !! Points are stored in growing array since same spacetime function or
  !! state variables can be used by multiple boundaries or sources
  type tem_grwPoints_type
    !> X-coordinate points
    type(grw_realArray_type) :: coordX
    !> Y-coordinate points
    type(grw_realArray_type) :: coordY
    !> Z-coordinate points
    type(grw_realArray_type) :: coordZ
  end type tem_grwPoints_type


  !> Contains source elements required for interpolation to derive
  !! solver variables
  type tem_sourceElems_type
    !> First position of source elements in srcElemPos per point
    !! size: nPnts
    type(grw_intArray_type) :: first

    !> last position of source elements in srcElemPos per point
    !! size: nPnts
    type(grw_intArray_type) :: last

    !> Position of source elements in levelwise list for solver variable.
    !! position in leveldesc is used so interpolation can be done level wise
    !! using ghost elements.
    !! In Musubi, the genertic get_valOfIndex routine and get_valOfIndex routine
    !! of state variables does not include halo elements for interpolation
    type(grw_intArray_type) :: elemPos

    !> Interpolation weight for each source elements
    type(grw_realArray_type) :: weight
  end type tem_sourceElems_type

  !> Data type contains growing array of points, evaluated variable value
  !! on those points per level
  !! For solver variables, additional information like elemPos and local coord
  !! are stored for every point.
  type tem_pointData_type
    !> Number of points in the growing array
    integer :: nPnts = 0

    !> growing array of points either global space coorinate in treelmmesh or
    !! local coordinate with on a reference element
    type(tem_grwPoints_type) :: grwPnt

    !> Unique treeID in finest level to create unique list of grwPnt
    type(dyn_longArray_type) :: treeID

    !> Offset bit encodes direction of boundary for surface coupling.
    !! used to translate space coordinate in the offset direction
    !! to determine the treeID in remote domain
    type(grw_charArray_type) :: offset_bit

    !> elemPos refer to position in linearized tree.
    !! size: number of points
    type(grw_intArray_type) :: elemPos

    !> For solver variables, pntlevel refer to the local level of the point
    !! In case of coupling, the requested level could differ from local
    !! level, hence we stored the local level explicitly
    !! size: number of points
    type(grw_intArray_type) :: pntlevel

    !> Contains source elements required for interpolation to derive
    !! solver variables for every point
    type(tem_sourceElems_type) :: srcElem
  end type tem_pointData_type

  !> Data type contains pointData information for all levels.
  type tem_pointData_list_type
    !> pointData for all levels
    type(tem_pointData_type) :: pntLvl(globalMaxLevels)

    !> Contains source elements position in global tree for interpolation to derive
    !! solver variables for every point. srcElem in tem_pointData_type contains
    !! source element position in level-wise list which is used in get_valOfIndex.
    !! This routine is used in solver get_point routine which interpolate a
    !! variable to a point from elements in fluid tree excluding ghosts.
    !! Halo elements are considered as source only for auxField variable.
    type(tem_pointData_type) :: mapToTree
  end type tem_pointData_list_type

  !> Interface to initialize growing array of points
  interface init
    module procedure init_grwPoints
  end interface init

  !> Interface to append a single point or array of points to growing
  !! array of points
  interface append
    module procedure append_singlePnt2grwPoints
    module procedure append_vectorPnt2grwPoints
  end interface append

  !> Interface to truncate growing array of points
  interface truncate
    module procedure truncate_grwPoints
  end interface truncate

  !> Interface to destroy growing array of points
  interface destroy
    module procedure destroy_grwPoints
  end interface destroy

  !> Interface to append point, offset_bit and elemPos to pointData
  interface append
    module procedure append_pointData
  end interface append

contains


  ! **************************************************************************** !
  !> Routine to append point Datas like points, offset_bit and elemPos
  !! Append point datas only if treeID of a point in max level
  !! is newly added
  subroutine append_pointData(me, point, storePnt, offset_bit, storeOffsetBit, &
    &                         elemPos, tree, pos, wasAdded )
    !---------------------------------------------------------------------------
    !> Point data type to be filled
    type(tem_pointData_type), intent(inout) :: me
    !> space coordinate to append
    real(kind=rk), intent(in) :: point(1:3)
    !> logical to store point into me%grwPnt
    logical, intent(in) :: storePnt
    !> offset bit to append
    character, intent(in) :: offset_bit
    !> logical to store offset bit into me%offset_bit
    logical, intent(in) :: storeOffsetBit
    !> Position of element which contains given point in global tree%treeID
    integer, intent(in) :: elemPos
    !> global tree
    type(treelmesh_type), intent(in) :: tree
    !> return position of treeID of a point in maxLevel in me%treeID
    integer, intent(out) :: pos
    !> If point is new and added to pointData
    logical, intent(out) :: wasAdded
    !---------------------------------------------------------------------------
    integer(kind=long_k) :: treeID
    !---------------------------------------------------------------------------
    treeID = tem_IdOfCoord(                            &
      &                 tem_CoordOfReal(tree, point) )
    call append(me       = me%treeID, &
      &         val      = treeID,    &
      &         pos      = pos,       &
      &         wasAdded = wasAdded   )

    if (wasAdded) then
      if (storePnt) then
        call append(me  = me%grwPnt, &
          &         val = point(1:3) )
      end if !storePnt

      if (storeOffsetBit) then
        call append(me  = me%offset_bit, &
          &         val = offset_bit     )
      end if !storeOffsetBit

      call append( me  = me%elemPos, &
        &          val = elemPos     )

      ! number of points in pointdata
      me%nPnts = me%treeID%nVals
    end if ! new point

  end subroutine append_pointData
  ! **************************************************************************** !

  ! **************************************************************************** !
  !> This routine initialize growing array of points
  subroutine init_grwPoints(me, length)
    !---------------------------------------------------------------------------
    !> Growing array of points in each dimension
    type(tem_grwPoints_type), intent(out) :: me
    !> Initial length of the container
    integer, optional, intent(in) :: length
    !---------------------------------------------------------------------------
    call init(me = me%coordX, length = length)
    call init(me = me%coordY, length = length)
    call init(me = me%coordZ, length = length)
  end subroutine init_grwPoints
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> This routine append a single point to growing array of points
  subroutine append_singlePnt2grwPoints(me, val)
    !---------------------------------------------------------------------------
    !> Growing array of points in each dimension
    type(tem_grwPoints_type), intent(inout) :: me
    !> single point to append
    real(kind=rk), intent(in) :: val(3)
    !---------------------------------------------------------------------------
    call append(me = me%coordX, val = val(1))
    call append(me = me%coordY, val = val(2))
    call append(me = me%coordZ, val = val(3))
  end subroutine append_singlePnt2grwPoints
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> This routine append a array of points to growing array of points
  subroutine append_vectorPnt2grwPoints(me, val)
    !---------------------------------------------------------------------------
    !> Growing array of points in each dimension
    type(tem_grwPoints_type), intent(inout) :: me
    !> Array of points to append
    real(kind=rk), intent(in) :: val(:,:)
    !---------------------------------------------------------------------------
    ! Expect points in shape(n,3) if not terminate
    if (size(val,dim=2)/=3) then
      write(logUnit(1),*) 'Error: Appending array of points to 1D growing'
      write(logUnit(1),*) '       array. Expects shape (n,3)'
      call tem_abort()
    end if

    call append(me = me%coordX, val = val(:,1))
    call append(me = me%coordY, val = val(:,2))
    call append(me = me%coordZ, val = val(:,3))
  end subroutine append_vectorPnt2grwPoints
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> This routine truncates growing array of points
  subroutine truncate_grwPoints(me)
    !---------------------------------------------------------------------------
    !> Growing array of points in each dimension
    type(tem_grwPoints_type), intent(inout) :: me
    !---------------------------------------------------------------------------
    call truncate(me = me%coordX)
    call truncate(me = me%coordY)
    call truncate(me = me%coordZ)
  end subroutine truncate_grwPoints
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> This routine destroys growing array of points
  subroutine destroy_grwPoints(me)
    !---------------------------------------------------------------------------
    !> Growing array of points in each dimension
    type(tem_grwPoints_type), intent(inout) :: me
    !---------------------------------------------------------------------------
    call destroy(me = me%coordX)
    call destroy(me = me%coordY)
    call destroy(me = me%coordZ)
  end subroutine destroy_grwPoints
  ! **************************************************************************** !

end module tem_pointData_module
