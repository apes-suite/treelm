! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2013 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011-2013, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011 Aravindh Krishnamoorthy <aravindh28.4@gmail.com>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2011 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2014, 2016, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
!> Additional methods for the TreElM module
!!
module tem_tools_module

  ! include treelm modules
  use env_module,          only: rk, int_k, long_k, single_k, stdOutUnit, PathLen

  implicit none

  private

  public :: append, appendUnique
  public :: resize
  public :: destroy
  public :: upper_to_lower
  public :: tem_intList, tem_longList
  public :: tem_horizontalSpacer
  public :: tem_PositionInSorted
  public :: tem_printArray
  public :: tem_file_to_string
  public :: tem_getOptValOrDef

  !> linked list of integers
  type tem_intList
    integer :: elem !< Position of Element in TreeIDlist
    type(tem_intList), pointer :: next => null()
  endType tem_intList

  !> linked list of long integers
  type tem_longList
    integer(kind=long_k) :: elem !< Position of Element in TreeIDlist
    type(tem_longList), pointer :: next => null()
  endType tem_longList

  interface tem_positioninsorted
    module procedure tem_positioninsorted_long
  end interface
  interface tem_positioninsorted
    module procedure tem_positioninsorted_int
  end interface

  !>@todo append position should be only int type. Only the value
  !! should be long_k
  interface append
    module procedure tem_appendIntList
    module procedure tem_appendLongList
    module procedure tem_appendInt1dArray
    module procedure tem_appendInt2dArray
    module procedure tem_appendIntLong1dArray
    module procedure tem_appendIntLong2dArray
    module procedure tem_appendSp1dArray
    module procedure tem_appendSp2dArray
    module procedure tem_appendDp1dArray
    module procedure tem_appendDp2dArray
    module procedure tem_appendIntLongArrayTo1dArray
  end interface

  interface appendUnique
    module procedure tem_appendIntListUnique
    module procedure tem_appendLongListUnique
    module procedure tem_appendInt1dArrayUnique
    module procedure tem_appendIntLong1dArrayUnique
  end interface

  interface resize
    module procedure tem_resizeInt1dArray
    module procedure tem_resizeInt2dArray
    module procedure tem_resizeIntLong1dArray
    module procedure tem_resizeIntLong2dArray
    module procedure tem_resizeDp1dArray
    module procedure tem_resizeDp2dArray
  end interface

  interface destroy
    module procedure tem_destroyIntList
    module procedure tem_destroyLongList
  end interface

  !> initial size of arrays, if size was 0
  integer, parameter :: InitialSize = 8



  interface tem_getOptValOrDef
    module procedure tem_getOptValOrDef_logical
    module procedure tem_getOptValOrDef_int
    module procedure tem_getOptValOrDef_long
    module procedure tem_getOptValOrDef_real
    module procedure tem_getOptValOrDef_char
  end interface tem_getOptValOrDef

contains

! ------------------------------------------------------------------------------
  !> returns the optional value, if present, or the default.
  !!
  !! this convenience routine encapsulates the check for an optional argument.
  function tem_getoptvalordef_logical( value, default ) result(res)
    ! --------------------------------------------------------------------------
    !> the optional value to check for
    logical, optional :: value
    !> the default to use when the optional value is not present.
    logical           :: default
    !> the result
    logical           :: res
    ! --------------------------------------------------------------------------
    if( present( value ) ) then
      res = value
    else
      res = default
    end if
  end function
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
  !> returns the optional value, if present, or the default.
  !!
  !! this convenience routine encapsulates the check for an optional argument.
  function tem_getoptvalordef_int( value, default ) result(res)
    ! --------------------------------------------------------------------------
    !> the optional value to check for
    integer, optional :: value
    !> the default to use when the optional value is not present.
    integer           :: default
    !> the result
    integer           :: res
    ! --------------------------------------------------------------------------
    if( present( value ) ) then
      res = value
    else
      res = default
    end if
  end function
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
  !> returns the optional value, if present, or the default.
  !!
  !! this convenience routine encapsulates the check for an optional argument.
  function tem_getoptvalordef_long( value, default ) result(res)
    ! --------------------------------------------------------------------------
    !> the optional value to check for
    integer(kind=long_k), optional :: value
    !> the default to use when the optional value is not present.
    integer(kind=long_k)           :: default
    !> the result
    integer(kind=long_k)           :: res
    ! --------------------------------------------------------------------------
    if( present( value ) ) then
      res = value
    else
      res = default
    end if
  end function
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
  !> returns the optional value, if present, or the default.
  !!
  !! this convenience routine encapsulates the check for an optional argument.
  function tem_getoptvalordef_real( value, default ) result(res)
    ! --------------------------------------------------------------------------
    !> the optional value to check for
    real(kind=rk), optional :: value
    !> the default to use when the optional value is not present.
    real(kind=rk)           :: default
    !> the result
    real(kind=rk)           :: res
    ! --------------------------------------------------------------------------
    if( present( value ) ) then
      res = value
    else
      res = default
    end if
  end function
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
  !> returns the optional value, if present, or the default.
  !!
  !! this convenience routine encapsulates the check for an optional argument.
  function tem_getoptvalordef_char( value, default ) result(res)
    ! --------------------------------------------------------------------------
    !> the optional value to check for
    character, optional :: value
    !> the default to use when the optional value is not present.
    character           :: default
    !> the result
    character           :: res
    ! --------------------------------------------------------------------------
    if( present( value ) ) then
      res = value
    else
      res = default
    end if
  end function
! ------------------------------------------------------------------------------

! ****************************************************************************** !
  !> Function to turn all upper case characters to lower case.
  !!
  !! The resulting string returned has the same length as the
  !! input string, and all upper case characters turned into
  !! lower case.
  !!
  function upper_to_lower(string) result(result_string)
    ! ---------------------------------------------------------------------------
    !> string to be converted
    character(len=*)    :: string
    !> converted string
    character(len=len(string))    :: result_string
    ! ---------------------------------------------------------------------------
    integer :: mm
    ! ---------------------------------------------------------------------------

    result_string=''
    do mm=1,len_trim(string)
      if( ('A' <= string(mm:mm)) .and. (string(mm:mm) <= 'Z') ) then
        result_string(mm:mm) = char(iachar(string(mm:mm))+32)
      else
        result_string(mm:mm)=string(mm:mm)
      end if
    end do
  end function upper_to_lower
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry at the end of the integer list
  !! If the first entry is zero, write into that one
  !! Check, if the current entry already exists
  !! Count, how many elements there are in the list
  subroutine tem_appendLongListUnique( firstEntry, entryPos, nItems, added )
    ! ---------------------------------------------------------------------------
    !> linked list of resulting elements building the neighbor
    type(tem_longList),pointer :: firstEntry
    !> Add that element
    integer(kind=long_k),intent(in) :: entryPos
    !> how many items are in list
    integer,intent(inout) :: nItems
    !> has the current item been added?
    logical,intent(out),optional :: added
    ! ---------------------------------------------------------------------------
    type(tem_longList),pointer :: currPos    ! current position in linked list
    logical :: found ! has the new entry been found in the existing list?
    ! ---------------------------------------------------------------------------

    Added = .false.
    found = .false.
    currPos => firstEntry ! The first entry of the list should be given here

    if(currPos%elem .le. 0) then
       ! If the element entry of the current entry is zero
       ! write into that position
       currPos%elem = entryPos
       nItems = 1
       Added = .true.
    else
      ! If element entry /= 0 then find the end of the list
      do while(associated(currPos%next))
        if( currPos%elem .eq. entryPos) found = .true.
        currPos => currPos%next
      enddo
      if( currPos%elem .eq. entryPos) found = .true.

      if ( .not. found) then
        ! at the end of the list, append new list item
        allocate(currPos%next)
        currPos => currPos%next
        ! And write to the item
        currPos%elem = entryPos
        nItems = nItems + 1
      endif
    endif

  end subroutine tem_appendLongListUnique
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry at the end of the integer list
  !! If the first entry is zero, write into that one
  !!
  subroutine tem_appendLongList(firstEntry, entryPos)
    ! ---------------------------------------------------------------------------
    !> linked list of resulting elements building the neighbor
    type(tem_longList),pointer :: firstEntry
    !> Add that element
    integer(kind=long_k),intent(in) :: entryPos
    ! ---------------------------------------------------------------------------
    type(tem_longList),pointer :: currPos    ! current position in linked list
    ! ---------------------------------------------------------------------------

    currPos => firstEntry

    if(currPos%elem .le. 0) then
       ! If the element entry of the current entry is zero
       ! write into that position
       currPos%elem = entryPos
    else
      ! If element entry /= 0 then find the end of the list
      do while(associated(currPos%next))
        currPos => currPos%next
      enddo
      ! at the end of the list, append new list item
      allocate(currPos%next)
      currPos => currPos%next
      ! And write to the item
      currPos%elem = entryPos
    endif

  end subroutine tem_appendLongList
! ****************************************************************************** !


! ****************************************************************************** !
  !> Destroy complete list
  !!
  pure subroutine tem_destroyLongList(ElemList)
    ! ---------------------------------------------------------------------------
    !> linked list of resulting elements building the neighbor
    type(tem_longList),pointer,intent(inout) :: ElemList
    ! ---------------------------------------------------------------------------
    type(tem_longList),pointer :: currPos    ! current position in linked list
    type(tem_longList),pointer :: lastPos    ! one before current position
    ! ---------------------------------------------------------------------------

    currPos  => ElemList

    ! Iterate as long as there are entries
    do while(associated(currPos%next))
      lastPos => currPos
      currPos => currPos%next
      deallocate(lastPos)
    enddo
    deallocate(currPos)

  end subroutine tem_destroyLongList
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry at the end of the integer list
  !! If the first entry is zero, write into that one
  !! Check, if the current entry already exists
  !! Count, how many elements there are in the list
  !!
  subroutine tem_appendIntListUnique( firstEntry, entryPos, nItems, added )
    ! ---------------------------------------------------------------------------
    !> linked list of resulting elements building the neighbor
    type(tem_intList),pointer :: firstEntry
    !> Add that element
    integer(kind=int_k),intent(in) :: entryPos
    !> how many items are in list
    integer,intent(inout) :: nItems
    !> has the current item been added?
    logical,intent(out),optional :: added
    ! ---------------------------------------------------------------------------
    type(tem_intList),pointer :: currPos    ! current position in linked list
    logical :: found ! has the new entry been found in the existing list?
    ! ---------------------------------------------------------------------------

    Added = .false.
    found = .false.
    currPos => firstEntry ! The first entry of the list should be given here

    if(currPos%elem .le. 0) then
       ! If the element entry of the current entry is zero
       ! write into that position
       currPos%elem = entryPos
       nItems = 1
       Added = .true.
    else
      ! If element entry /= 0 then find the end of the list
      do while(associated(currPos%next))
        if( currPos%elem .eq. entryPos) found = .true.
        currPos => currPos%next
      enddo
      if( currPos%elem .eq. entryPos) found = .true.

      if ( .not. found) then
        ! at the end of the list, append new list item
        allocate(currPos%next)
        currPos => currPos%next
        ! And write to the item
        currPos%elem = entryPos
        nItems = nItems + 1
      endif
    endif

  end subroutine tem_appendIntListUnique
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry at the end of the integer list
  !! If the first entry is zero, write into that one
  subroutine tem_appendIntList(firstEntry, entryPos)
    ! ---------------------------------------------------------------------------
    !> linked list of resulting elements building the neighbor
    type(tem_intList),pointer :: firstEntry
    !> Add that element
    integer(kind=int_k),intent(in) :: entryPos
    ! ---------------------------------------------------------------------------
    type(tem_intList),pointer :: currPos    ! current position in linked list
    ! ---------------------------------------------------------------------------

    currPos => firstEntry

    if(currPos%elem .le. 0) then
       ! If the element entry of the current entry is zero
       ! write into that position
       currPos%elem = entryPos
    else
      ! If element entry /= 0 then find the end of the list
      do while(associated(currPos%next))
        currPos => currPos%next
      enddo
      ! at the end of the list, append new list item
      allocate(currPos%next)
      currPos => currPos%next
      ! And write to the item
      currPos%elem = entryPos
    endif

  end subroutine tem_appendIntList
! ****************************************************************************** !


! ****************************************************************************** !
  !> Destroy complete list
  !!
  pure subroutine tem_destroyIntList(ElemList)
    ! ---------------------------------------------------------------------------
    !> linked list of resulting elements building the neighbor
    type(tem_intList),pointer, intent(inout) :: ElemList
    ! ---------------------------------------------------------------------------
    type(tem_intList),pointer :: firstPos   ! first entry in linked list
    type(tem_intList),pointer :: currPos    ! current position in linked list
    type(tem_intList),pointer :: lastPos    ! one before current position
    ! ---------------------------------------------------------------------------

    firstPos => ElemList
    currPos  => ElemList

    ! Iterate as long as there are entries
    do while(associated(currPos%next))
      lastPos => currPos
      currPos => currPos%next
      deallocate(lastPos)
    enddo

  end subroutine tem_destroyIntList
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry to an allocatable array 1d with long integer
  !! If the array is too small, reallocate with double size
  !!
  subroutine tem_appendIntLongArrayTo1dArray(Array, ArrayToAppend )
    ! ---------------------------------------------------------------------------
    !> array to append to
    integer(kind=long_k),intent(inout), allocatable :: Array(:)
    !> array to append
    integer(kind=long_k),intent(in) :: ArrayToAppend(:)
    ! ---------------------------------------------------------------------------
    integer(kind=long_k),allocatable :: tempArray(:)
    integer(kind=long_k) :: ArraySize,ArraySize2,ierr,NewSize
    ! ---------------------------------------------------------------------------

    ! Get size of array
    ArraySize  = size(Array,1)
    ArraySize2 = size(ArrayToAppend,1)

    NewSize = ArraySize + ArraySize2
    ! allocate temporary array with new size
    allocate(tempArray(NewSize),stat=ierr)
    ! Copy both arrays to temp array
    tempArray(1:ArraySize) = Array(1:ArraySize)
    tempArray(ArraySize+1:ArraySize+ArraySize2) = ArrayToAppend(1:ArraySize2)
    ! Deallocate Array
    deallocate(Array)
    ! Reallocate Array
    allocate(Array(NewSize),stat=ierr)
    Array(1:NewSize) = tempArray(1:NewSize)
    ! Deallocate temp array
    deallocate(tempArray)
    if(ierr .ne. 0) Write(*,*) 'Error in reallocating array'

  end subroutine tem_appendIntLongArrayTo1dArray
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry to an allocatable array 1d with long integer
  !! If the array is too small, reallocate with double size
  !!
  subroutine tem_appendIntLong1dArray(Array, Position, Value )
    ! ---------------------------------------------------------------------------
    !> array to append value to
    integer(kind=long_k),intent(inout), allocatable :: Array(:)
    !> position the value is appended to
    integer,intent(in) :: Position
    !> value to append
    integer(kind=long_k),intent(in) :: Value
    ! ---------------------------------------------------------------------------
    integer(kind=long_k),allocatable :: tempArray(:)
    integer :: ArraySize,ierr,NewSize
    logical :: sizeZero
    ! ---------------------------------------------------------------------------

    ! Get size of array
    ArraySize = size(Array,1)
    ! Compare position, where to store with size
    if(Position .gt. ArraySize) then ! If position is larger than size
      sizeZero = .false.
      if( ArraySize .eq. 0) then
         ArraySize = 1
         sizeZero  = .true.
      endif
      NewSize = max(ArraySize*2, Position)
      ! allocate temporary array with double size
      allocate(tempArray(NewSize), stat=ierr)
      ! Copy to temp array
      if( .not. sizeZero ) tempArray(1:ArraySize) = Array(1:ArraySize)
      ! Deallocate Array
      deallocate(Array)
      ! Reallocate Array
      allocate(Array(NewSize), stat=ierr)
      Array(1:ArraySize) = tempArray(1:ArraySize)
      ! Deallocate temp array
      deallocate(tempArray)
      if(ierr .ne. 0) Write(*,*) 'Error in reallocating array'
    endif

    Array(Position) = Value


  end subroutine tem_appendIntLong1dArray
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry to an allocatable array 1d with long integer
  !! If the array is too small, reallocate with double size
  !!
  subroutine tem_appendIntLong2dArray(Array, Position1, Position2, Value )
    ! ---------------------------------------------------------------------------
    !> array to append value to
    integer(kind=long_k),intent(inout), allocatable :: Array(:,:)
    !>
    integer(kind=long_k),intent(in) :: Position1
    !>
    integer(kind=long_k),intent(in) :: Position2
    !> value to append
    integer(kind=long_k),intent(in) :: Value
    ! ---------------------------------------------------------------------------
    integer(kind=long_k),allocatable :: tempArray(:,:)
    integer(kind=long_k) :: ArraySize1,ArraySize2,ierr
    integer(kind=long_k) :: NewSize1,NewSize2
    logical :: changeSize
    logical :: sizeZero
    ! ---------------------------------------------------------------------------

    ! Get size of array
    ArraySize1 = size(Array,1)
    ArraySize2 = size(Array,2)
    changeSize = .false.
    sizeZero   = .false.

    ! Compare position, where to store with size
    if(Position1 .gt. ArraySize1) then
      if( ArraySize1 .eq. 0) then
         ArraySize1 = 1
         sizeZero  = .true.
      endif
      NewSize1 = max(Position1, ArraySize1*2)
      changeSize = .true.
    else
      NewSize1 = ArraySize1
    endif
    if(Position2 .gt. ArraySize2) then
      if( ArraySize2 .eq. 0) then
         ArraySize2 = 1
         sizeZero  = .true.
      endif
      NewSize2 = max( Position2, ArraySize2*2 )
      changeSize = .true.
    else
      NewSize2 = ArraySize2
    endif

    if(changeSize) then
      ! allocate temporary array with double size
      allocate(tempArray(NewSize1, NewSize2),stat=ierr)
      ! Copy to temp array
      if(.not. sizeZero)                                                       &
        & tempArray( 1:ArraySize1,1:ArraySize2 ) =                             &
        &                                    Array( 1:ArraySize1,1:ArraySize2 )
      ! Deallocate Array
      deallocate(Array)

      ! Reallocate Array
      allocate(Array(NewSize1, NewSize2),stat=ierr)
      Array(1:ArraySize1,1:ArraySize2) = tempArray(1:ArraySize1,1:ArraySize2)
      ! Deallocate temp array
      deallocate(tempArray)
      if(ierr .ne. 0) Write(*,*) 'Error in reallocating array'
    endif

    Array(Position1,Position2) = Value


  end subroutine tem_appendIntLong2dArray
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry to an allocatable array 1d with single integer at the end
  !! If the array is too small, reallocate with double size
  !!
  subroutine tem_appendInt1dArrayUnique( Array, Value, nElems, Pos, Added )
    ! ---------------------------------------------------------------------------
    !> array to append value to
    integer(kind=int_k),intent(inout), allocatable :: Array(:)
    !> The value to add as an entry in Array
    integer(kind=int_k),intent(in) :: Value
    !> Number of entries in the array (changes, if added = .true.)
    integer,intent(inout) :: nElems
    !> position the value was appended
    integer,intent(out) :: Pos
    !> new entry in array added?
    logical,intent(out) :: Added
    ! ---------------------------------------------------------------------------
    integer,allocatable :: tempArray(:)
    integer :: ArraySize,ierr,NewSize, iPos
    logical :: found !< has the new entry been found in the existing list?
    ! ---------------------------------------------------------------------------

    added  = .false.
    found = .false.

    ! Get size of array
    ArraySize = size(Array,1)
    if( ArraySize .gt. 0) then
      ! check if value already exists
      do iPos = 1, nElems
        if(Array(iPos) .eq. Value ) then
          found = .true.
          Pos   = iPos
        endif
      enddo
      if(.not. found) then  ! Add to array
        nElems = nElems + 1
        if ( nElems .gt. ArraySize ) then ! Resize array
          NewSize = ArraySize*2
          allocate(tempArray(NewSize),stat=ierr)
          ! Copy to temp array
          tempArray(1:ArraySize) = Array(1:ArraySize)
          ! Deallocate Array
          deallocate(Array)
          ! Reallocate Array
          allocate(Array(NewSize),stat=ierr)
          Array(1:ArraySize) = tempArray(1:ArraySize)
          ! Deallocate temp array
          deallocate(tempArray)
          if(ierr .ne. 0) Write(*,*) 'Error in reallocating array'
        endif
        Array( nElems ) = Value
        added = .true.
        Pos = nElems
      endif ! not found yet
    else ! size zero. reallocate and add value.
      deallocate( Array )
      allocate(Array( InitialSize ))
      nElems = 1
      Array( nElems ) = Value
      added = .true.
      Pos = nElems
    endif

  end subroutine tem_appendInt1dArrayUnique
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry to an allocatable array 1d with single integer at the end
  !! If the array is too small, reallocate with double size
  !!
  subroutine tem_appendIntLong1dArrayUnique(Array, Value, nElems, Pos, Added )
    ! ---------------------------------------------------------------------------
    !> array to append value to
    integer(kind=long_k),intent(inout), allocatable :: Array(:)
    !> The value to add as an entry in Array
    integer(kind=long_k),intent(in) :: Value
    !> Number of entries in the array (changes, if added = .true.)
    integer,intent(inout) :: nElems
    !> position the value was appended
    integer,intent(out) :: Pos
    !> new entry in array added?
    logical,intent(out) :: Added
    ! ---------------------------------------------------------------------------
    integer(kind=long_k),allocatable :: tempArray(:)
    integer :: ArraySize, ierr, NewSize, iPos
    logical :: found ! has the new entry been found in the existing list?
    ! ---------------------------------------------------------------------------

    added  = .false.
    found = .false.

    ! Get size of array
    ArraySize = size(Array,1)
    if( ArraySize .gt. 0) then
      ! check if value already exists
      do iPos = 1, nElems
        if(Array(iPos) .eq. Value ) then
          found = .true.
          Pos   = iPos
        endif
      enddo
      if(.not. found) then  ! Add to array
        nElems = nElems + 1
        if ( nElems .gt. ArraySize ) then ! Resize array
          NewSize = ArraySize*2
          allocate(tempArray(NewSize), stat=ierr)
          ! Copy to temp array
          tempArray(1:ArraySize) = Array(1:ArraySize)
          ! Deallocate Array
          deallocate(Array)
          ! Reallocate Array
          allocate(Array(NewSize), stat=ierr)
          Array(1:ArraySize) = tempArray(1:ArraySize)
          ! Deallocate temp array
          deallocate(tempArray)
          if(ierr .ne. 0) Write(*,*) 'Error in reallocating array'
        endif
        Array( nElems ) = Value
        added = .true.
        Pos = nElems
      endif ! not found yet
    else ! size zero. reallocate and add value.
      deallocate( Array )
      allocate(Array( InitialSize ))
      nElems = 1
      Array( nElems ) = Value
      added = .true.
      Pos = nElems
    endif

  end subroutine tem_appendIntLong1dArrayUnique
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry to an allocatable array 1d with single integer
  !! If the array is too small, reallocate with double size
  !!
  subroutine tem_appendInt1dArray(Array, Position, Value )
    ! ---------------------------------------------------------------------------
    !> array to append value to
    integer(kind=int_k),intent(inout), allocatable :: Array(:)
    !> position the value is appended to
    integer,intent(in) :: Position
    !> value to append
    integer(kind=int_k),intent(in) :: Value
    ! ---------------------------------------------------------------------------
    integer,allocatable :: tempArray(:)
    integer :: ArraySize,ierr, NewSize
    logical :: sizeZero
    ! ---------------------------------------------------------------------------

    ! Get size of array
    ArraySize = size(Array,1)
    ! Compare position, where to store with size
    if(Position .gt. ArraySize) then ! If position is larger than size
      sizeZero = .false.
      if( ArraySize .eq. 0) then
         ArraySize = 1
         sizeZero  = .true.
      endif
      ! allocate temporary array with double size
      NewSize = max(ArraySize*2, Position)
      allocate(tempArray(NewSize),stat=ierr)
      ! Copy to temp array
      if( .not. sizeZero ) tempArray(1:ArraySize) = Array(1:ArraySize)
      ! Deallocate Array
      deallocate(Array)
      ! Reallocate Array
      allocate(Array(NewSize),stat=ierr)
      Array(1:ArraySize) = tempArray(1:ArraySize)
      ! Deallocate temp array
      deallocate(tempArray)
      if(ierr .ne. 0) Write(*,*) 'Error in reallocating array'
    endif

    Array(Position) = Value

  end subroutine tem_appendInt1dArray
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry to an allocatable array 1d with single integer
  !! If the array is too small, reallocate with double size
  !!
  subroutine tem_appendInt2dArray(Array, Position1, Position2, Value )
    ! ---------------------------------------------------------------------------
    !> array to append value to
    integer(kind=int_k),intent(inout), allocatable :: Array(:,:)
    !>
    integer,intent(in) :: Position1
    !>
    integer,intent(in) :: Position2
    !> value to append
    integer(kind=int_k),intent(in) :: Value
    ! ---------------------------------------------------------------------------
    integer,allocatable :: tempArray(:,:)
    integer :: ArraySize1,ArraySize2,ierr
    integer :: NewSize1,NewSize2
    logical :: changeSize
    logical :: sizeZero
    ! ---------------------------------------------------------------------------

    ! Get size of array
    ArraySize1 = size(Array,1)
    ArraySize2 = size(Array,2)
    changeSize = .false.
    sizeZero   = .false.

    ! Compare position, where to store with size
    if(Position1 .gt. ArraySize1) then
      if( ArraySize1 .eq. 0) then
         ArraySize1 = 1
         sizeZero  = .true.
      endif
      NewSize1 = max( Position1, ArraySize1*2 )
      changeSize = .true.
    else
      NewSize1 = ArraySize1
    endif
    if(Position2 .gt. ArraySize2) then
      if( ArraySize2 .eq. 0) then
         ArraySize2 = 1
         sizeZero  = .true.
      endif
      NewSize2 = max( ArraySize2*2, Position2 )
      changeSize = .true.
    else
      NewSize2 = ArraySize2
    endif

    if(changeSize) then
      ! allocate temporary array with double size
      allocate(tempArray(NewSize1, NewSize2),stat=ierr)
      ! Copy to temp array
      if(.not. sizeZero)                                                       &
        & tempArray( 1:ArraySize1,1:ArraySize2 ) =                             &
        &                                    Array( 1:ArraySize1,1:ArraySize2 )
      ! Deallocate Array
      deallocate(Array)
      ! Reallocate Array
      allocate(Array(NewSize1, NewSize2),stat=ierr)
      Array(1:ArraySize1,1:ArraySize2) = tempArray(1:ArraySize1,1:ArraySize2)
      ! Deallocate temp array
      deallocate(tempArray)
      if(ierr .ne. 0) Write(*,*) 'Error in reallocating array'
    endif

    Array(Position1,Position2) = Value

  end subroutine tem_appendInt2dArray
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry to an allocatable array 1d with single precision
  !! If the array is too small, reallocate with double size
  subroutine tem_appendSp1dArray(Array, Position, Value )
    ! ---------------------------------------------------------------------------
    !> array to append value to
    real(kind=single_k),intent(inout), allocatable :: Array(:)
    !> position the value is appended to
    integer,intent(in) :: Position
    !> value to append
    real(kind=single_k),intent(in) :: Value
    ! ---------------------------------------------------------------------------
    real(kind=single_k),allocatable :: tempArray(:)
    integer :: ArraySize,ierr, NewSize
    logical :: sizeZero
    ! ---------------------------------------------------------------------------

    ! Get size of array
    ArraySize = size(Array,1)
    ! Compare position, where to store with size
    if(Position .gt. ArraySize) then ! If position is larger than size
      sizeZero = .false.
      if( ArraySize .eq. 0) then
         ArraySize = 1
         sizeZero  = .true.
      endif
      NewSize = max(ArraySize*2, Position)
      ! allocate temporary array with double size
      allocate(tempArray(NewSize),stat=ierr)
      ! Copy to temp array
      if( .not. sizeZero ) tempArray(1:ArraySize) = Array(1:ArraySize)
      ! Deallocate Array
      deallocate(Array)
      ! Reallocate Array
      allocate(Array(NewSize),stat=ierr)
      Array(1:ArraySize) = tempArray(1:ArraySize)
      ! Deallocate temp array
      deallocate(tempArray)
      if(ierr .ne. 0) Write(*,*) 'Error in reallocating array'
    endif

    Array(Position) = Value

  end subroutine tem_appendSp1dArray
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry to an allocatable array 1d with single precision
  !! If the array is too small, reallocate with double size
  !!
  subroutine tem_appendSp2dArray(Array, Position1, Position2, Value )
    ! ---------------------------------------------------------------------------
    !> array to append value to
    real(kind=single_k),intent(inout), allocatable :: Array(:,:)
    !>
    integer,intent(in) :: Position1
    !>
    integer,intent(in) :: Position2
    !> value to append
    real(kind=single_k),intent(in) :: Value
    ! ---------------------------------------------------------------------------
    real(kind=single_k),allocatable :: tempArray(:,:)
    integer :: ArraySize1,ArraySize2,ierr
    integer :: NewSize1,NewSize2
    logical :: changeSize
    logical :: sizeZero
    ! ---------------------------------------------------------------------------

    ! Get size of array
    ArraySize1 = size(Array,1)
    ArraySize2 = size(Array,2)
    changeSize = .false.
    sizeZero   = .false.

    ! Compare position, where to store with size
    if(Position1 .gt. ArraySize1) then
      if( ArraySize1 .eq. 0) then
         ArraySize1 = 1
         sizeZero  = .true.
      endif
      NewSize1 = max( Position1, ArraySize1*2 )
      changeSize = .true.
    else
      NewSize1 = ArraySize1
    endif
    if(Position2 .gt. ArraySize2) then
      if( ArraySize2 .eq. 0) then
         ArraySize2 = 1
         sizeZero  = .true.
      endif
      NewSize2 = max( Position2, ArraySize2*2 )
      changeSize = .true.
    else
      NewSize2 = ArraySize2
    endif

    if(changeSize) then
      ! allocate temporary array with double size
      allocate(tempArray(NewSize1, NewSize2),stat=ierr)
      ! Copy to temp array
      if(.not. sizeZero)                                                       &
        & tempArray( 1:ArraySize1,1:ArraySize2 ) =                             &
        &                                    Array( 1:ArraySize1,1:ArraySize2 )
      ! Deallocate Array
      deallocate(Array)
      ! Reallocate Array
      allocate(Array(NewSize1, NewSize2),stat=ierr)
      Array(1:ArraySize1,1:ArraySize2) = tempArray(1:ArraySize1,1:ArraySize2)
      ! Deallocate temp array
      deallocate(tempArray)
      if(ierr .ne. 0) Write(*,*) 'Error in reallocating array'
    endif

    Array(Position1,Position2) = Value

  end subroutine tem_appendSp2dArray
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry to an allocatable array 1d with double precision
  !! If the array is too small, reallocate with double size
  !!
  subroutine tem_appendDp1dArray(Array, Position, Value )
    ! ---------------------------------------------------------------------------
    !> array to append value to
    real(kind=rk),intent(inout), allocatable :: Array(:)
    !> position the value is appended to
    integer,intent(in) :: Position
    !> value to append
    real(kind=rk),intent(in) :: Value
    ! ---------------------------------------------------------------------------
    real(kind=rk),allocatable :: tempArray(:)
    integer :: ArraySize,ierr, NewSize
    logical :: sizeZero
    ! ---------------------------------------------------------------------------

    ! Get size of array
    ArraySize = size(Array,1)
    ! Compare position, where to store with size
    if(Position .gt. ArraySize) then ! If position is larger than size
      sizeZero = .false.
      if( ArraySize .eq. 0) then
         ArraySize = 1
         sizeZero  = .true.
      endif
      NewSize = max(ArraySize*2, Position)
      ! allocate temporary array with double size
      allocate(tempArray(Newsize),stat=ierr)
      ! Copy to temp array
      if( .not. sizeZero ) tempArray(1:ArraySize) = Array(1:ArraySize)

      ! Deallocate Array
      deallocate(Array)
      ! Reallocate Array
      allocate(Array(NewSize),stat=ierr)
      Array(1:ArraySize) = tempArray(1:ArraySize)
      ! Deallocate temp array
      deallocate(tempArray)
      if(ierr .ne. 0) Write(*,*) 'Error in reallocating array'
    endif

    Array(Position) = Value

  end subroutine tem_appendDp1dArray
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry to an allocatable array 1d with double precision
  !! If the array is too small, reallocate with double size
  !!
  subroutine tem_appendDp2dArray(Array, Position1, Position2, Value )
    ! ---------------------------------------------------------------------------
    !> array to append value to
    real(kind=rk),intent(inout), allocatable :: Array(:,:)
    !>
    integer,intent(in) :: Position1
    !>
    integer,intent(in) :: Position2
    !> value to append
    real(kind=rk),intent(in) :: Value
    ! ---------------------------------------------------------------------------
    real(kind=rk),allocatable :: tempArray(:,:)
    integer :: ArraySize1,ArraySize2,ierr
    integer :: NewSize1,NewSize2
    logical :: changeSize
    logical :: sizeZero
    ! ---------------------------------------------------------------------------

    ! Get size of array
    ArraySize1 = size(Array,1)
    ArraySize2 = size(Array,2)
    changeSize = .false.
    sizeZero   = .false.

    ! Compare position, where to store with size
    if(Position1 .gt. ArraySize1) then
      if( ArraySize1 .eq. 0) then
         ArraySize1 = 1
         sizeZero  = .true.
      endif
      NewSize1 = max( Position1, ArraySize1*2)
      changeSize = .true.
    else
      NewSize1 = ArraySize1
    endif
    if(Position2 .gt. ArraySize2) then
      if( ArraySize2 .eq. 0) then
         ArraySize2 = 1
         sizeZero  = .true.
      endif
      NewSize2 = max( Position2, ArraySize2*2 )
      changeSize = .true.
    else
      NewSize2 = ArraySize2
    endif

    if(changeSize) then
      ! allocate temporary array with double size
      allocate(tempArray(NewSize1, NewSize2),stat=ierr)
      ! Copy to temp array
      if(.not. sizeZero)                                                       &
        & tempArray( 1:ArraySize1,1:ArraySize2 ) =                             &
        &                                    Array( 1:ArraySize1,1:ArraySize2 )
      ! Deallocate Array
      deallocate(Array)
      ! Reallocate Array
      allocate(Array(NewSize1, NewSize2),stat=ierr)
      Array(1:ArraySize1,1:ArraySize2) = tempArray(1:ArraySize1,1:ArraySize2)
      ! Deallocate temp array
      deallocate(tempArray)
      if(ierr .ne. 0) Write(*,*) 'Error in reallocating array'
    endif

    Array(Position1,Position2) = Value

  end subroutine tem_appendDp2dArray
! ****************************************************************************** !


! ****************************************************************************** !
  !> write out a message with the defined string+real content
  subroutine tem_horizontalSpacer( fUnit, before, after, toFile, str )
    ! ---------------------------------------------------------------------------
    !> output unit
    integer,intent(in)          :: fUnit
    !> optional white lines before print
    integer,intent(in),optional :: before
    !> optional white lines after print
    integer,intent(in),optional :: after
    !> output goes to debug file
    logical,intent(in),optional :: toFile
    !> String to be printed at the beginning of the spacer
    character(*),intent(in),optional :: str
    ! ---------------------------------------------------------------------------
    ! loop variable
    integer :: i
    ! local temporary variable for the output unit
    integer :: nUnit
    ! ---------------------------------------------------------------------------

    nUnit = fUnit

    if( present( toFile)) then
      if( toFile ) nUnit = 66
    endif

    if( present(before) ) then
      do i=1,before
        write(nUnit,*)
      enddo
    endif

    if( present(str) ) then
      write(nUnit, '(5a)') repeat('=', (80-len_trim(str))/2 -1), &
        &                  ' ',                                  &
        &                  trim(str),                            &
        &                  repeat(' ', 1+mod(len_trim(str), 2)), &
        &                  repeat('=',(80-len_trim(str))/2 -1)
    else
      write(nUnit, '(a)') repeat('-', 80)
    end if

    if( present(after ) ) then
      do i=1,after
        write(nUnit,*)
      enddo
    endif
  end subroutine tem_horizontalSpacer
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry to an allocatable array 1d with integer
  !! If the array is too small, reallocate with double size
  !!
  subroutine tem_resizeInt1dArray(Array, Newsize )
    ! ---------------------------------------------------------------------------
    !> array to resize
    integer(kind=int_k),intent(inout), allocatable :: Array(:)
    !> new size of the array
    integer,intent(in) :: Newsize
    ! ---------------------------------------------------------------------------
    integer,allocatable :: tempArray(:)
    integer :: ierr
    ! ---------------------------------------------------------------------------
    ! allocate temporary array with double size
    allocate(tempArray(NewSize), stat=ierr)
    ! Copy to temp array
    tempArray(1:NewSize) = Array(1:NewSize)
    ! Deallocate Array
    deallocate(Array)
    ! Reallocate Array
    allocate(Array(NewSize), stat=ierr)
    Array(1:NewSize) = tempArray(1:NewSize)
    ! Deallocate temp array
    deallocate(tempArray)
    if(ierr .ne. 0) Write(*,*) 'Error in reallocating array'

  end subroutine tem_resizeInt1dArray
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry to an allocatable array 2d with integer
  !! If the array is too small, reallocate with double size
  !!
  subroutine tem_resizeInt2dArray( Array, Newsize1, Newsize2 )
    ! ---------------------------------------------------------------------------
    !> array to resize
    integer(kind=int_k),intent(inout), allocatable :: Array(:,:)
    !> first new size
    integer,intent(in) :: Newsize1
    !> second new size
    integer,intent(in) :: Newsize2
    ! ---------------------------------------------------------------------------
    integer,allocatable :: tempArray(:,:)
    integer :: ierr
    ! ---------------------------------------------------------------------------

    ! allocate temporary array with new size
    allocate(tempArray(NewSize1,NewSize2),stat=ierr)
    if(ierr .ne. 0) Write(*,*) 'Error in allocating temp_array'
    ! Copy to temp array
    tempArray(1:NewSize1, 1:NewSize2) = Array(1:NewSize1, 1:NewSize2)
    ! Deallocate Array
    deallocate(Array)
    ! Reallocate Array
    allocate(Array(NewSize1, Newsize2), stat=ierr)
    ! Deallocate temp array
    Array(1:NewSize1,1:NewSize2) = tempArray(1:NewSize1, 1:NewSize2)
    deallocate(tempArray)
    if(ierr .ne. 0) Write(*,*) 'Error in reallocating array'

  end subroutine tem_resizeInt2dArray
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry to an allocatable array 1d with long integer
  !! If the array is too small, reallocate with double size
  !!
  subroutine tem_resizeIntLong1dArray(Array, Newsize )
    ! ---------------------------------------------------------------------------
    !> array to resize
    integer(kind=long_k),intent(inout), allocatable :: Array(:)
    !> new size of the array
    integer,intent(in) :: Newsize
    ! ---------------------------------------------------------------------------
    integer(kind=long_k),allocatable :: tempArray(:)
    integer(kind=long_k) :: ierr
    ! ---------------------------------------------------------------------------

    ! allocate temporary array with double size
    allocate(tempArray(NewSize), stat=ierr)
    ! Copy to temp array
    tempArray(1:NewSize) = Array(1:NewSize)
    ! Deallocate Array
    deallocate(Array)
    ! Reallocate Array
    allocate(Array(NewSize), stat=ierr)
    Array(1:NewSize) = tempArray(1:NewSize)
    ! Deallocate temp array
    deallocate(tempArray)
    if(ierr .ne. 0) Write(*,*) 'Error in reallocating array'

  end subroutine tem_resizeIntLong1dArray
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry to an allocatable array 2d with long integer
  !! If the array is too small, reallocate with double size
  !!
  subroutine tem_resizeIntLong2dArray( Array, Newsize1, Newsize2 )
    ! ---------------------------------------------------------------------------
    !> array to resize
    integer(kind=long_k),intent(inout), allocatable :: Array(:,:)
    !> first new size
    integer,intent(in) :: Newsize1
    !> second new size
    integer,intent(in) :: Newsize2
    ! ---------------------------------------------------------------------------
    integer(kind=long_k),allocatable :: tempArray(:,:)
    integer :: ierr
    ! ---------------------------------------------------------------------------

    ! allocate temporary array with new size
    allocate(tempArray(NewSize1,NewSize2),stat=ierr)
    ! Copy to temp array
    tempArray(1:NewSize1,1:NewSize2) = Array(1:NewSize1,1:NewSize2)
    ! Deallocate Array
    deallocate(Array)
    ! Reallocate Array
    allocate(Array(NewSize1,Newsize2),stat=ierr)
    ! Deallocate temp array
    Array(1:NewSize1,1:NewSize2) = tempArray(1:NewSize1,1:NewSize2)
    deallocate(tempArray)
    if(ierr .ne. 0) Write(*,*) 'Error in reallocating array'

  end subroutine tem_resizeIntLong2dArray
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry to an allocatable array 1d with integer
  !! If the array is too small, reallocate with double size
  !!
  subroutine tem_resizeDp1dArray( Array, Newsize )
    ! ---------------------------------------------------------------------------
    !> array to resize
    real(kind=rk),intent(inout), allocatable :: Array(:)
    !> new size of the array
    integer,intent(in) :: Newsize
    ! ---------------------------------------------------------------------------
    real(kind=rk),allocatable :: tempArray(:)
    integer :: ierr
    ! ---------------------------------------------------------------------------

    ! allocate temporary array with double size
    allocate(tempArray(NewSize), stat=ierr)
    ! Copy to temp array
    tempArray(1:NewSize) = Array(1:NewSize)
    ! Deallocate Array
    deallocate(Array)
    ! Reallocate Array
    allocate(Array(NewSize), stat=ierr)
    Array(1:NewSize) = tempArray(1:NewSize)
    ! Deallocate temp array
    deallocate(tempArray)
    if(ierr .ne. 0) Write(*,*) 'Error in reallocating array'

  end subroutine tem_resizeDp1dArray
! ****************************************************************************** !


! ****************************************************************************** !
  !> append an entry to an allocatable array 2d with integer
  !! If the array is too small, reallocate with double size
  !!
  subroutine tem_resizeDp2dArray( Array, Newsize1, Newsize2 )
    ! ---------------------------------------------------------------------------
    !> array to resize
    real(kind=rk),intent(inout), allocatable :: Array(:,:)
    !> first new size
    integer,intent(in) :: Newsize1
    !> second new size
    integer,intent(in) :: Newsize2
    ! ---------------------------------------------------------------------------
    real(kind=rk),allocatable :: tempArray(:,:)
    integer :: ierr
    ! ---------------------------------------------------------------------------

    ! allocate temporary array with new size
    allocate(tempArray(NewSize1,NewSize2),stat=ierr)
    ! Copy to temp array
    tempArray(1:NewSize1,1:NewSize2) = Array(1:NewSize1,1:NewSize2)
    ! Deallocate Array
    deallocate(Array)
    ! Reallocate Array
    allocate(Array(NewSize1,Newsize2),stat=ierr)
    ! Deallocate temp array
    Array(1:NewSize1,1:NewSize2) = tempArray(1:NewSize1,1:NewSize2)
    deallocate(tempArray)
    if(ierr .ne. 0) Write(*,*) 'Error in reallocating array'

  end subroutine tem_resizeDp2dArray
! ****************************************************************************** !


! ****************************************************************************** !
! ****************************************************************************** !


! ****************************************************************************** !
  !> return the position of a value in 'me', which is an array with sorted
  !! entries.
  !! if the value was not found,
  !!  - return 0 if nextifnotfound = .false.
  !!  - return position at the end if nextifnotfound = .true.
  !!
  function tem_positioninsorted_long( me, val, lower, upper ) result(pos)
    ! ---------------------------------------------------------------------------
    !> array to search in
    integer(kind=long_k), intent(in) :: me(:)
    !> value to look for
    integer(kind=long_k), intent(in) :: val
    !> lower search limit
    integer, intent(in), optional :: lower
    !> upper search limit
    integer, intent(in), optional :: upper
    !> position of val in the sorted list, 0 if not found
    integer :: pos
    ! ---------------------------------------------------------------------------
    logical :: retnext
    integer :: lb, ub
    integer :: mid
    integer(kind=long_k) :: lb_val, ub_val
    integer(kind=long_k) :: mid_val
    ! ---------------------------------------------------------------------------

    retnext = .false.

    lb = 1
    ub = size( me )

    if( present( lower ) ) lb = lower
    if( present( upper ) ) ub = upper

    pos = 0
    if (retnext) pos = lb

    ! binary search on sorted list
    do while(ub >= lb)

      lb_val = me(lb)

      if (val < lb_val) then
        if (retnext) pos = lb
        exit
      end if

      ub_val = me(ub)

      if (val > ub_val) then
        if (retnext) pos = ub+1
        exit
      end if

      mid = (lb+ub) / 2
      mid_val = me(mid)
      if (val == mid_val) then
        pos = mid
        exit
      end if
      if (val > mid_val) then
        lb = mid + 1
      else
        ub = mid - 1
      end if
    end do
  end function tem_positioninsorted_long
  !> return the position of a value in 'me', which is an array with sorted
  !! entries.
  !! if the value was not found,
  !!  - return 0 if nextifnotfound = .false.
  !!  - return position at the end if nextifnotfound = .true.
  !!
  function tem_positioninsorted_int( me, val, lower, upper ) result(pos)
    ! ---------------------------------------------------------------------------
    !> array to search in
    integer(kind=int_k), intent(in) :: me(:)
    !> value to look for
    integer(kind=int_k), intent(in) :: val
    !> lower search limit
    integer, intent(in), optional :: lower
    !> upper search limit
    integer, intent(in), optional :: upper
    !> position of val in the sorted list, 0 if not found
    integer :: pos
    ! ---------------------------------------------------------------------------
    logical :: retnext
    integer :: lb, ub
    integer :: mid
    integer(kind=int_k) :: lb_val, ub_val
    integer(kind=int_k) :: mid_val
    ! ---------------------------------------------------------------------------

    retnext = .false.

    lb = 1
    ub = size( me )

    if( present( lower ) ) lb = lower
    if( present( upper ) ) ub = upper

    pos = 0
    if (retnext) pos = lb

    ! binary search on sorted list
    do while(ub >= lb)

      lb_val = me(lb)

      if (val < lb_val) then
        if (retnext) pos = lb
        exit
      end if

      ub_val = me(ub)

      if (val > ub_val) then
        if (retnext) pos = ub+1
        exit
      end if

      mid = (lb+ub) / 2
      mid_val = me(mid)
      if (val == mid_val) then
        pos = mid
        exit
      end if
      if (val > mid_val) then
        lb = mid + 1
      else
        ub = mid - 1
      end if
    end do
  end function tem_positioninsorted_int
! ****************************************************************************** !


! ****************************************************************************** !
  !> print an array to the debugunit
  !!
  subroutine tem_printArray( me, itemLength, title, lineLength, nUnit )
    ! ---------------------------------------------------------------------------
    !> Array title in debug output for easy identification in the file
    character( len=* ),optional :: title
    !> long array to write to debug file
    integer(kind=long_k), intent(in) :: me(:)
    !> how many characters needs each item of the array to output
    integer,optional :: itemLength
    !> how long should the line be
    integer,optional :: lineLength
    !> to which unit to output
    integer :: nUnit
    ! ---------------------------------------------------------------------------
    integer :: iElem
    integer :: elemLength, meLength
    integer :: itemLength_loc, lineLength_loc
    character( len=128 ) :: buffer
    character( len=128 ) :: spacer
    ! ---------------------------------------------------------------------------

    if( present(linelength)) then
      linelength_loc = min( abs(linelength), 128 )
    else
      linelength_loc = 120
    endif
     if( present(itemlength)) then
      itemlength_loc = max( itemlength, 8 )
    else
      itemlength_loc =  8
    endif

    meLength = size( me(:) )
    if ( meLength > 0 ) then

      ! build spacer
      spacer = ''
      do iElem = 1, lineLength_loc -2
        write(spacer,'(2a)') trim(spacer),'-'
      enddo

      if( present(title)) then
        write(nUnit,*) trim(spacer)
        write(nUnit,'(2a,i0)') trim(title), ', size: ', size(me)
        ! write(nUnit,*) trim(spacer)
      endif
      buffer = ''
      elemLength = itemLength_loc + 3

      do iElem = 1, meLength
        write(buffer,'(2a,i8)') trim(buffer),' ',me( iElem )
        if( iElem == meLength .or. mod( iElem, elemLength ) == 0) then
          write(nUnit,*) trim(buffer)
          buffer = ''
        endif
      enddo
      write(nUnit,*) trim(spacer)
    end if ! meLength > 0

  end subroutine tem_printArray
! ****************************************************************************** !


  !> Read a file from a connected unit into a string.
  !!
  !! The connected file has to be opened for sequential formatted access.
  !! A string will be returned containing the characters read from the file.
  !! If there are potential problems arising, they are returned in the
  !! error code iError
  subroutine tem_file_to_string(funit, string, iError)
    ! ---------------------------------------------------------------------------
    !> File unit to read, has to be opened sequential and formatted.
    integer, intent(in) :: funit
    !> String to fill with the content of the file.
    character(len=*), intent(out) :: string
    !> Error code:
    !!
    !! 0 = no error
    !! 1 = end of string reached before end of file
    !! 2 = Unit not connected
    integer, intent(out) :: iError
    ! ---------------------------------------------------------------------------
    integer :: old_pos, length
    integer :: stringlen
    integer :: io
    logical :: nUnitOpened
    character(len=PathLen) :: loc_string
    ! ---------------------------------------------------------------------------
    string = ''
    stringlen = len(string)
    iError = 2
    inquire(unit=funit, opened=nUnitOpened)
    if (nUnitOpened) then
      ! Go to the start of the file.
      rewind(funit)
      old_pos = 0
      do
        ! Read all contents of the scratch file until eof is encountered
        read(funit, '(a)', iostat=io) loc_string
        if (io == 0) then
          length = len_trim(loc_string)+1
          ! Only proceed, if the line still fits into the string
          if (old_pos+length > stringlen) then
            string(old_pos+1:) = loc_string(:stringlen-old_pos-length+1)
            iError = 1
            EXIT
          end if
          string(old_pos+1:old_pos+length) = trim(loc_string)//new_line('A')
          old_pos = old_pos+length
        else
          ! Reached end of the file, exit the loop.
          iError = 0
          EXIT
        end if
      end do
    end if

  end subroutine tem_file_to_string
! ****************************************************************************** !

end module tem_tools_module
! ****************************************************************************** !
