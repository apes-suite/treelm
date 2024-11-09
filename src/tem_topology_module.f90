! Copyright (c) 2011-2013, 2015 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2012, 2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2014-2015 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
!> author: Harald Klimach
!! This module provides informations and operations on
!! the topology of the complete linearized octree.
!! It is completely independent of the actual sparse
!! octree and therefore this module is independent of the Treelmesh_module.
!!
!! Node identification in the octree
!! ---
!! Due to the breadth first numbering of the octree
!! it is straight forward to move in vertical direction through the tree.
!!
!! The parent \(i_p\) of a node with a certain
!! [[treelmesh_type:treeID]] \(i_n\) is given by:
!! \[  i_p = \left\lfloor \frac{i_n - 1}{8} \right\rfloor \]
!! Where the bracket \(\left\lfloor\cdot\right\rfloor\) denotes the integer floor
!! rounding to the next smaller integer number.
!! While the children \(i_{c1..8}\) of a given node \(i_n\) can be obtained by:
!! \[   i_{cj} = 8 \cdot i_n + j \]
!! The numbering scheme begins at the root node with a
!! [[treelmesh_type:treeID]] of 0 for
!! this element, that covers the complete enclosed cubic domain.
!! All further [[treelmesh_type:treeID]]s are then
!! given by the rules stated above.
!! What remains to be defined is the actual spatial placement of the nodes in
!! the tree as elements in the mesh.
!! This is achieved by a space filling curve, that describes the numbering
!! scheme of the 8 children for a given node.
!!
!! Our choice here is the Morton ordering [6] (bibliography of treelm) , that
!! allows an efficient computation of the index from the coordinates with the
!! same orientation on each refinement level.
!! For this ordering the coordinates of the children within their parent element
!! are first noted as a three dimensional tuple, where each index might be 0 or
!! 1.
!! Now this tuple is interpreted as a single number of the basis 2 with three
!! digits.
!! The recursive in depth traversal results then in a concatenation of these
!! numbers visited during the traversal.
!! This concatenated number provides the sorting key for all elements on a given
!! level according to the space filling curve.
!!
module tem_topology_module

  ! include treelm modules
  use env_module,       only: long_k, globalMaxLevels

  implicit none

  private

  public :: tem_ParentOf
  public :: tem_childNumber
  public :: tem_SiblingsOfId
  public :: tem_PathOf
  public :: tem_LevelOf
  public :: tem_FirstIdAtLevel, tem_LastIdAtLevel
  public :: tem_directChildren
  public :: tem_TreeIDComparison
  public :: tem_PathComparison
  public :: tem_path_type
  public :: tem_coordOfID, tem_IDofCoord

  integer(kind=long_k), parameter, public :: firstIdAtLevel(20) =  [&
    &                         1_long_k, & ! level =  1
    &                         9_long_k, & ! level =  2
    &                        73_long_k, & ! level =  3
    &                       585_long_k, & ! level =  4
    &                      4681_long_k, & ! level =  5
    &                     37449_long_k, & ! level =  6
    &                    299593_long_k, & ! level =  7
    &                   2396745_long_k, & ! level =  8
    &                  19173961_long_k, & ! level =  9
    &                 153391689_long_k, & ! level = 10
    &                1227133513_long_k, & ! level = 11
    &                9817068105_long_k, & ! level = 12
    &               78536544841_long_k, & ! level = 13
    &              628292358729_long_k, & ! level = 14
    &             5026338869833_long_k, & ! level = 15
    &            40210710958665_long_k, & ! level = 16
    &           321685687669321_long_k, & ! level = 17
    &          2573485501354569_long_k, & ! level = 18
    &         20587884010836553_long_k, & ! level = 19
    &        164703072086692425_long_k  ] ! level = 20

  !> Paths of elements from the root node to themselves going through the
  !! hierarchy of the tree. Is used to compare to elements
  !!
  !! For all
  !! [[tem_construction_module::identify_local_element]] "local elements",
  !! their actual position in the sparse mesh has to be identified for a given
  !! \ref treelmesh_module::treelmesh_type::treeid "treeID".
  !! Due to the sparsity of the mesh, the position of a certain element in the
  !! total list of elements can not be directly deduced from its
  !! \ref treelmesh_module::treelmesh_type::treeid "treeID".
  !! As already pointed out, we can rely on binary searches in the sorted list
  !! of elements for this task, as long as a comparison between any two elements
  !! is possible.
  !! A comparison of two nodes to decide their relative position in the ordered
  !! list of elements has to take into account the hierarchy of the mesh.
  !! The comparison operator is therefore defined as follows.
  !!
  !! - Build the path from each leaf to the root of the tree, by repeatedly
  !!   computing the parents with the equation used in the
  !!   [tem_ParentOf function] (@ref tem_parentatlevel), until the root with
  !!   \ref treelmesh_module::treelmesh_type::treeid "treeID" \f$0\f$ is reached.
  !! - After the path through the tree is known for both leaves,
  !!   comparison can be done on the greatest common level by a simple integer
  !!   subtraction of the treeIDs in the path on this level.
  !! - The result of this subtraction indicates the relation between the two
  !!   paths, if it is zero they are considered to be equal.
  !!   Otherwise the sign indicates the ordering of the two compared leaves.
  !!
  !! In a tree with leaves on different levels, this results in children to be
  !! considered to be equal to their parents.
  type tem_path_type
    !> Levels counted from 1
    !! This level is 1 higher than the actual level.
    integer :: Level
    !> treeIDs from current leaf to root
    integer(kind=long_k) :: Node(globalMaxLevels+1)
  end type

  ! interface tem_LevelOf
  !   module procedure tem_levelOf_singleTreeID
  !   module procedure tem_levelOf_arrayOfTreeIDs
  ! end interface

  interface tem_ParentOf
    module procedure tem_directParent
    module procedure tem_ParentAtLevel
  end interface

  contains

! ****************************************************************************** !
  !> This function delivers the refinement level of a TreeID
  !!
  !! The refinement level of a
  !! \ref treelmesh_module::treelmesh_type::treeid "treeID" is found by
  !! repeatedly computing the parents with the equation used in the
  !! [tem_ParentOf function] (@ref tem_parentatlevel) and
  !! counting the iterations, until the universal root is reached.
  !!
  ! elemental function tem_LevelOf_singleTreeID( TreeID ) result(res)
  elemental function tem_LevelOf( TreeID ) result(res)
    ! ---------------------------------------------------------------------------
    !> ID in the complete tree to look up
    integer(kind=long_k), intent(in) :: TreeID
    !> Level of the given TreeID
    integer :: res
    ! ---------------------------------------------------------------------------
    integer(kind=long_k) :: tElem
    integer :: level
    ! ---------------------------------------------------------------------------

    level = 0
    tElem = TreeID
    do while (tElem /= 0)
      ! bit operation is more efficient than division
      ! tElem = (tElem-1) / 8
      tElem = ishft( (tElem-1), -3 )
      level = level + 1
    end do
    res = level

  ! end function tem_LevelOf_singleTreeID
  end function tem_LevelOf
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function delivers the refinement level of a TreeID
  !!
  !! The refinement level of a
  !! \ref treelmesh_module::treelmesh_type::treeid "treeID" is found by
  !! repeatedly computing the parents with the equation used in the
  !! [tem_ParentOf function] (@ref tem_parentatlevel) and
  !! counting the iterations, until the universal root is reached.
  !!
  ! function tem_LevelOf_arrayOfTreeIDs(TreeID, nElems) result(res)
  !   ! ---------------------------------------------------------------------------
  !   !> number of elements in the array
  !   integer, intent(in) :: nElems
  !   !> Levels of the given TreeIDs
  !   integer :: res(nElems)
  !   !> ID in the complete tree to look up
  !   integer(kind=long_k), intent(in) :: TreeID(nElems)
  !   ! ---------------------------------------------------------------------------
  !   integer(kind=long_k) :: tElem
  !   integer :: level, iElem
  !   ! ---------------------------------------------------------------------------

  !   do iElem = 1, size(treeID)
  !     level = 0
  !     tElem = TreeID(iElem)
  !     do while (tElem /= 0)
  !       ! bit operation is more efficient than division
  !       ! division by 8 is equvalent to right shift 3 bits
  !       ! tElem = (tElem-1) / 8
  !       tElem = ishft( (tElem-1), -3 )
  !       level = level + 1
  !     end do
  !     res(iElem ) = level
  !   end do

  ! end function tem_LevelOf_arrayOfTreeIDs
! ****************************************************************************** !


! ****************************************************************************** !
  !> First TreeID in the complete tree on a given level
  !!
  !! e.g. level 1:  1 ..  8   -> offset = 1
  !!      level 2:  9 .. 72   -> offset = 9
  !!      level 3: 73 ..      -> offset = 73
  !!
  !! A level offset is calculated by summing all the nodes of coarser levels.
  !! With a level count starting at \f$0\f$ on the coarsest level (the universal
  !! root node), this offset \f$s\f$ for a given level \f$L\f$ can be computed by
  !! the follwing equation:
  !! \f[
  !!   s = \sum\limits_{i=0}^{L-1} 8^i = \frac{8^{L} - 1}{7}
  !! \f]
  !!
  !! @todo: calling to this function can be replaced by accessing array
  !!        firstIdAtLevel which explicitly stores the results
  elemental function tem_FirstIdAtLevel( level ) result(res)
    ! ---------------------------------------------------------------------------
    !> level to check
    integer, intent(in) :: level
    !> resulting first treeID
    integer(kind=long_k) :: res
    ! ---------------------------------------------------------------------------

    res = ( 8_long_k**level - 1_long_k ) / 7_long_k

  end function tem_FirstIdAtLevel
! ****************************************************************************** !


! ****************************************************************************** !
  !> Last ID in the complete tree on a given level
  !!
  elemental function tem_LastIdAtLevel( level ) result(res)
    ! ---------------------------------------------------------------------------
    !> level to check
    integer, intent(in) :: level
    !> resulting last treeID
    integer(kind=long_k) :: res
    ! ---------------------------------------------------------------------------
    integer :: i
    integer(kind=long_k) :: monomonial
    ! ---------------------------------------------------------------------------

    res = -1
    monomonial = 1
    do i = 0, level
      res = res + monomonial
      monomonial = monomonial * 8
    end do

  end function tem_LastIdAtLevel
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function delivers all sibling treeIDs of TreeID in an array
  !!
  function tem_SiblingsOfId( TreeID ) result(siblings)
    ! ---------------------------------------------------------------------------
    !> current treeID
    integer(kind=long_k), intent(in) :: TreeID
    !> treeIDs siblings
    integer(kind=long_k) :: siblings(7)
    ! ---------------------------------------------------------------------------
    integer :: myNumber, iElem, i
    integer(kind=long_k) :: offset
    ! ---------------------------------------------------------------------------

    offset   = ishft((TreeID-1_long_k), -3) * 8_long_k
    myNumber = int( TreeID - offset )
    iElem = 0
    do i = 1,8
      if( i /= myNumber ) then
        iElem = iElem + 1
        siblings( iElem ) = offset + int(i,long_k)
      endif
    enddo

  end function tem_SiblingsOfId
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function delivers of TreeID, which child number it is from its parent
  !!
  elemental function tem_childNumber( TreeID ) result(res)
    ! ---------------------------------------------------------------------------
    !> current treeID
    integer(kind=long_k), intent(in) :: TreeID
    !> result of function containing child number
    integer :: res
    ! ---------------------------------------------------------------------------
    integer(kind=long_k) :: offset
    ! ---------------------------------------------------------------------------

    ! offset = ((treeID-1)/8) * 8
    offset = ishft((TreeID-1_long_k), -3) * 8_long_k
    res    = int( TreeID - offset )

  end function tem_childNumber
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function delivers the parent ID of a given TreeID
  !!
  elemental function tem_directParent(TreeID) result(res)
    ! ---------------------------------------------------------------------------
    !> current treeID
    integer(kind=long_k), intent(in) :: TreeID
    !> result of function containing parent ID
    integer(kind=long_k) :: res
    ! ---------------------------------------------------------------------------

    ! res = ( treeID - 1 ) / 8
    res = ishft( (treeID-1_long_k), -3 )

  end function tem_directParent
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function returns the complete path through the tree from
  !! a given treeID to the root (all parents).
  !!
  elemental function tem_PathOf( TreeID ) result(res)
    ! ---------------------------------------------------------------------------
    !> current treeID
    integer(kind=long_k), intent(in) :: TreeID
    !> resulting path
    type(tem_path_type) :: res
    ! ---------------------------------------------------------------------------

    res%level = 1
    res%Node(res%Level) = treeID
    do while (res%Node(res%Level) > 0)
      res%Level = res%Level + 1
      ! res%Node(res%Level) = (res%Node(res%Level-1)-1) / 8
      res%Node(res%Level) = ishft( (res%Node(res%Level-1)-1), -3 )
    end do
  end function tem_PathOf
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function provides the parent ID of a given tree ID on a given level.
  !!
  !! The parent ID on a given level is calculated using the equation
  !! \f[
  !!    id_{l-1}=\frac{id_{l}-1}{8}
  !! \f]
  !! in a loop.
  !!
  function tem_ParentAtLevel(TreeID, level) result(parentID)
    ! ---------------------------------------------------------------------------
    !> treeID of which the pID is requested
    integer(kind=long_k),intent(in) :: TreeID
    !> the level on which the pID is requested
    integer,intent(in) :: level
    !> resulting parent ID
    integer(kind=long_k) :: parentID
    ! ---------------------------------------------------------------------------
    integer :: cLevel ! the current level of the tree ID
    integer :: iLevel ! level iterator
    ! ---------------------------------------------------------------------------

    parentID = TreeID
    cLevel = tem_LevelOf(TreeID)

    do iLevel = 1, cLevel-level
      ! parentID = (parentID - 1) / 8
      parentID = ishft( (parentID-1_long_k), -3 )
    end do

  end function tem_ParentAtLevel
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function delivers the direct children in the full tree for
  !! a given tree ID
  !!
  function tem_directChildren( TreeID ) result(childrenIDs)
    ! ---------------------------------------------------------------------------
    !> given treeID
    integer(kind=long_k), intent(in) :: TreeID
    !> Array for the treeIDs of the 8 direct children
    integer(kind=long_k) :: childrenIDs(8)
    ! ---------------------------------------------------------------------------
    integer :: childCount
    integer(kind=long_k) :: off
    ! ---------------------------------------------------------------------------

    childrenIDs(:) = 0

    off = treeID * 8_long_k
    do childCount = 1,8
      childrenIDs(childCount) = off + int(childCount, kind=long_k)
    end do

  end function tem_directChildren
! ****************************************************************************** !


! ****************************************************************************** !
  !> Compare Position of two treeIDs in the linearized tree
  !! This is relatively expansive, therefore the result should
  !! be stored, if more than one comparison of the same treeIDs
  !! has to be done.
  !! Result:
  !! -1: left is lower  than right
  !!  0: left is same   than right
  !!  1: left is higher than right
  !!
  elemental function tem_TreeIDComparison( left, right ) result(relation)
    ! ---------------------------------------------------------------------------
    !> candidate treeID
    integer(kind=long_k), intent(in) :: left
    !> candidate treeID
    integer(kind=long_k), intent(in) :: right
    !> relation between the treeIDs
    integer :: relation
    ! ---------------------------------------------------------------------------
    type(tem_path_type) :: left_path
    type(tem_path_type) :: right_path
    ! ---------------------------------------------------------------------------

    ! Build the path to the left from the leaf to the root (0).
    left_path = tem_PathOf(left)

    ! Build the path to the right from the leaf to the root (0).
    right_path = tem_PathOf(right)

    relation = tem_PathComparison(left_path, right_path)

  end function tem_TreeIDComparison
! ****************************************************************************** !


! ****************************************************************************** !
  !> Compare two Paths in the linearized tree
  !! Result:
  !! -1: left is lower  than right
  !!  0: left is same   than right
  !!  1: left is higher than right
  !!
  elemental function tem_PathComparison( left, right ) result(relation)
    ! ---------------------------------------------------------------------------
    !> candidate path
    type(tem_path_type), intent(in) :: left
    !> candidate path
    type(tem_path_type), intent(in) :: right
    !> relation between the paths
    integer :: relation
    ! ---------------------------------------------------------------------------
    integer :: maxPathLen
    integer(kind=long_k) :: diff
    ! ---------------------------------------------------------------------------

    ! Find the greatest common level.
    maxPathLen = min(left%level-1, right%level-1)

    ! Difference of treeIDs on that level indicates the relation of the two
    ! paths.
    diff =  left%Node( left%Level - maxPathLen)                                &
      &  - right%Node(right%Level - maxPathLen)

    ! Normalize the relation.
    relation = int(diff / max(1_long_k, abs(diff)))

  end function tem_PathComparison
! ****************************************************************************** !


! ****************************************************************************** !
  !> The following function provides the spatial index triple for a given ID in
  !! the complete mesh, on its refinement level.
  !! The returned coordinates start at 0.
  !! The fourth entry is the level on which the coordinates are located in the
  !! tree.
  !! The steps are following:
  !! 1. calculate the refinement level, store to coord(4).
  !! 2. Calcualte the level offset.
  !! 3. the Morton index is obtained by subtracting the offset from treeID.
  !! 4. Apply bit interleaving rule to generate coordinates.
  !!
  pure function tem_CoordOfId(TreeID, offset) result(coord)
    ! ---------------------------------------------------------------------------
    !> input element ID
    integer(kind=long_k), intent(in) :: TreeID
    !> First Elem of current level, if known, to avoid recomputations
    integer(kind=long_k), intent(in), optional :: offset
    !> coordinate of element return value
    integer :: coord(4)
    ! ---------------------------------------------------------------------------
    integer(kind=long_k) :: tElem ! temporary element
    integer(kind=long_k) :: fak(3)
    integer :: bitlevel
    integer :: i
    ! ---------------------------------------------------------------------------

    coord(:) = 0
    coord(4) = tem_LevelOf(treeID)
    fak(1) = 1
    fak(2) = 2
    fak(3) = 4
    bitlevel = 1

    ! offset has to be subtracted to get the "real" TreeID
    if (present(offset)) then
      tElem = TreeID - offset
    else
      tElem = TreeID - FirstIdAtLevel(coord(4))
    end if

    ! get x coordinate from
    !
    ! x = sum(iLevel=0, inf) { 2**iLevel * mod(tElem / (8**iLevel) ),2)   }
    ! y = sum(iLevel=0, inf) { 2**iLevel * mod(tElem /(2 * (8**iLevel) ),2)   }
    ! z = sum(iLevel=0, inf) { 2**iLevel * mod(tElem /(4 * (8**iLevel) ),2)   }
    !
    do while ((tElem / fak(1)) > 0)
      do i=1,3
        coord(i) = coord(i) + bitlevel * int(mod(tElem / fak(i), 2_long_k))
      end do
      bitlevel = bitlevel * 2
      fak = fak * 8
    end do
  end function tem_CoordOfId
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function calculates the tree ID for a given x,y and z on the given
  !! refinement level. This ID in the complete tree does not have to be in the
  !! list of leafs (treeIDlist)
  !! An example of this procedure is following:
  !! 1. Convert coordinates into binary representation:
  !!   (x,y,z) = (5,9,1) = (0101,1001,0001)
  !! 2. Interleaving the bits by 3 digits for each direction:
  !!        x(0101):   **0 **1 **0 **1
  !!        y(1001):   *1* *0* *0* *1*
  !!        z(0001):   0** 0** 0** 1**
  !!    Combining these bits results in a binary number: 010 001 000 111,
  !!    which is 1095 when seen as a 10-base number.
  !! 3. This number is the cell position in a single level sorted element list.
  !!    Its treeID can be obtained by adding the correspoding level offset.
  !!
  pure function tem_IdOfCoord(coord, offset) result(nTreeID)
    ! ---------------------------------------------------------------------------
    !> 3D Coordinates to translate
    integer, intent(in) :: coord(4)
    !>
    integer(kind=long_k), intent(in), optional :: offset
    !> resulting treeID
    integer(kind=long_k) :: nTreeID
    ! ---------------------------------------------------------------------------
    integer :: iLevel    ! level iterator
    ! this variable defines the upper threshold for x, y and z coordinates
    ! for a given refinement level e.g. for level = 1 only x, y and z
    ! from 0 to 1 are allowed, otherwise this neighbor is outside of the
    ! geometry
    integer :: upperBoundary
    integer :: myCoord(3)
    integer :: fak2
    integer(kind=long_k) :: fak8
    ! ---------------------------------------------------------------------------

    ! the offset has to be added to the level summation
    ! @todo: level offset is hard coded in array firstIdAtLevel,
    !        thus no need to pass offset as an argument.
    if (present(offset)) then
      nTreeID = offset
    else
      nTreeID = firstIdAtLevel(coord(4))
    end if

    !! Periodic displacement for requested elements outside the cube
    upperBoundary = 2**coord(4)
    myCoord = mod(coord(1:3) + upperBoundary, upperBoundary)

    ! To compute the tree ID of a neighbor on the same refinement level the
    ! following formula is used:
    ! newTreeID = sum(i=0, level) {8**i * (mod( (x / (2**level) ) , 2)
    !                                      + 2 * mod( ( y / (2**level) ) , 2)
    !                                      + 4 * mod( (z / (2**level) ) , 2) ) }
    fak8 = 1
    fak2 = 1
    do iLevel = 0, coord(4) - 1
      nTreeID = nTreeID + fak8 * (      mod( (myCoord(1) / fak2), 2) &
        &                         + 2 * mod( (myCoord(2) / fak2), 2) &
        &                         + 4 * mod( (myCoord(3) / fak2), 2) )
      fak2 = fak2 * 2
      fak8 = fak8 * 8
    end do

  end function tem_IdOfCoord
! ****************************************************************************** !


end module tem_topology_module
! ****************************************************************************** !
