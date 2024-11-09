! Copyright (c) 2018 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2018-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> Contains data_types and function for matrix operations
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
module tem_matrix_module
  use env_module,           only: rk, minLength, zeroLength
  use tem_aux_module,       only: tem_abort
  use tem_logging_module,   only: logUnit
  use tem_debug_module,     only: dbgUnit
  use tem_dyn_array_module, only: dyn_intArray_type, init, append, truncate, &
    &                             destroy
  use tem_grow_array_module, only: grw_logicalArray_type, init, append, &
    &                              truncate, destroy
  use tem_param_module,     only: c_x, c_y, c_z
  use tem_float_module,     only: operator(.feq.)

  implicit none

  private

  public :: tem_matrix_type
  public :: invert_matrix
  public :: tem_intpMatrixLSF_type
  public :: tem_matrix_dump
  public :: init, append, truncate, destroy, empty, placeAt

  integer, parameter :: maxIntp_order = 2
  ! Number of coefffs for quadratic interpolation
  ! (/ linear, quadratic /)
  !> For 1D stencil,
  !! 2 unknown coeffs: p(x)=a0+a1 x for linear 1st order interpolation
  !! 3 unknown coeffs: p(x)=a0+a1 x+a2 x^2 for quadratic 2nd order interpolation
  integer, dimension(maxIntp_order), parameter :: nCoeffs_1D = (/ 2,  3 /)
  !> For 2D stencil,
  !! 3 unknown coeffs for linear 1st order interpolation:
  !! p(x,y)=a0+a1 x+a2 y
  !! 6 unknown coeffs for quadratic 2nd order interpolation:
  !! p(x,y)=a0+a1 x+a2 y+a3 x^2+a4 y^2+ a5 xy
  integer, dimension(maxIntp_order), parameter :: nCoeffs_2D = (/ 3,  6 /)
  !> For 3D stencil,
  !! 4 unknown coeffs for linear 1st order interpolation:
  !! p(x,y,z)=a0+a1 x+a2 y+a3 z
  !! 10 unknown coeffs for quadratic 2nd order interpolation:
  !! p(x,y,z)=a0+a1 x+a2 y+a3 z+a4 x^2+a5 y^2+ a6 z^2+ a7 xy + a8 yz + a9 zx
  integer, dimension(maxIntp_order), parameter :: nCoeffs_3D = (/ 4, 10 /)

  !> This derived type encapsulates the definition of the matrix
  type tem_matrix_type
    !> inverted matrix to solve linear system of equation
    real(kind=rk), allocatable :: A(:,:)
    !> how many entries are in the 2d matrix?
    integer :: nEntries(2)
  end type tem_matrix_type

  !> growing array type for type(tem_matrix_type)
  type grw_matrixarray_type
    integer :: nvals = 0
    integer :: containersize = 0
    type(tem_matrix_type), allocatable :: val(:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga_matrix
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_ga_matrix
  end interface

  !> empty the entries  without changing arrays
  interface empty
    module procedure empty_ga_matrix
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga_matrix
  end interface

  !> insert an element at a given position
  interface placeat
    module procedure placeat_ga_matrix
    module procedure placeat_ga_matrix_vec
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_ga_matrix
    module procedure append_ga_matrix_vec
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga_matrix
  end interface


  !> This derived type encapsulates the definition of least square fit matrix
  !! for interpolation method which is required for every combination of
  !! available nSourceFromCoarser
  type tem_intpMatrixLSF_type
    type(grw_matrixArray_type) :: matArray
    !> Unique hash ID to identify different combination of available
    !! nSourceFromCoarser
    type(dyn_intArray_type) :: ID
    !> nCoeffs required for least square fit.
    !! Depends on nDims and order of interpolation
    integer:: nCoeffs
    !> For every matrix in matArray, store if its invertible or not
    !! to avoid rebuilding singular matrix
    type(grw_logicalArray_type) :: isInvertible
  end type tem_intpMatrixLSF_type

  interface init
    module procedure init_intpMatrixLSF
  end interface init

  interface append
    module procedure append_intpMatrixLSF
  end interface append

  interface truncate
    module procedure truncate_intpMatrixLSF
  end interface truncate

  interface destroy
    module procedure destroy_intpMatrixLSF
  end interface destroy

  interface assignment(=)
    module procedure copy_matrix
  end interface

contains

  subroutine init_ga_matrix(me, length)
    type(grw_matrixarray_type), intent(out) :: me !< dynamic array to init
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

  end subroutine init_ga_matrix

  subroutine destroy_ga_matrix(me)
    type(grw_matrixarray_type), intent(inout) :: me !< dynamic array to destroy

    me%containersize = 0
    me%nvals = 0
    if( allocated( me%val ) ) deallocate(me%val)

  end subroutine destroy_ga_matrix


  subroutine truncate_ga_matrix(me)
    !------------------------------------------------------------------------
    type(grw_matrixarray_type) :: me !< array to truncate
    !------------------------------------------------------------------------
    type(tem_matrix_type), allocatable :: tarray(:)
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

  end subroutine truncate_ga_matrix


  subroutine empty_ga_matrix(me)
    !------------------------------------------------------------------------
    type(grw_matrixarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------

    me%nvals = 0

  end subroutine empty_ga_matrix

  !> adds the value to a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! element at the requested position will be replaced.
  subroutine placeat_ga_matrix(me, val, pos, length)
    type(grw_matrixarray_type) :: me !< array to place the value into
    type(tem_matrix_type), intent(in) :: val !< value to place at the given position
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

  end subroutine placeat_ga_matrix


  !> adds the values starting from a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! elements starting from the requested position will be replaced up to
  !! the element at position `pos + size(val) - 1`.
  subroutine placeat_ga_matrix_vec(me, val, pos, length)
    type(grw_matrixarray_type) :: me !< array to append the value to
    type(tem_matrix_type), intent(in) :: val(:) !< values to append
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

  end subroutine placeat_ga_matrix_vec


  subroutine append_ga_matrix(me, val, length)
    type(grw_matrixarray_type) :: me !< array to append the value to
    type(tem_matrix_type), intent(in) :: val !< value to append
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

  end subroutine append_ga_matrix

  subroutine append_ga_matrix_vec(me, val, length)
    type(grw_matrixarray_type) :: me !< array to append the value to
    type(tem_matrix_type), intent(in) :: val(:) !< values to append
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

  end subroutine append_ga_matrix_vec


  subroutine expand_ga_matrix(me, pos, length)
    type(grw_matrixarray_type) :: me !< array to resize
    integer, intent(in), optional :: pos !< optional predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    type(tem_matrix_type), allocatable :: swpval(:)
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

  end subroutine expand_ga_matrix


! **************************************************************************** !
  !> This routine initialize interpolation matrix for least square fit
  subroutine init_intpMatrixLSF(me, length, nDims, order)
    ! --------------------------------------------------------------------------
    type(tem_intpMatrixLSF_type), intent(out) :: me
    integer, intent(in) :: length
    integer, intent(in) :: nDims
    integer, intent(in) :: order
    ! --------------------------------------------------------------------------
    call init(me = me%matArray, length = length)
    call init(me = me%ID, length = length)
    call init(me = me%isInvertible, length = length)

    ! Set coeffs required for each order
    if (order > 0 .and. order <= maxIntp_order) then
      select case(nDims)
      case(1)
        me%nCoeffs = nCoeffs_1D(order)
      case(2)
        me%nCoeffs = nCoeffs_2D(order)
      case(3)
        me%nCoeffs = nCoeffs_3D(order)
      end select
    else
      call tem_abort('Unsupported interpolation order')
    end if

  end subroutine init_intpMatrixLSF
! **************************************************************************** !

  ! ************************************************************************** !
  !> This routine builds up the matrix for least square fit used in
  !! linear and quadratic interpolation.
  !!
  !! Compute interpolation matrix for least square fit using stencil
  !! direction of available sources
  !! The parent of target childs coord is 0,0,0 so we could
  !! just use of stencil%cxDir to build up this matrix entries
  !! Every row in matrix is evaluated with coord of source element
  subroutine append_intpMatrixLSF(me, order, QQ, nDims, nSources, cxDirRK, &
    &                             neighDir, pos, success)
    ! --------------------------------------------------------------------------
    !> intpMatrix for LSF fill
    type(tem_intpMatrixLSF_type), intent(inout) :: me
    !> interpolation order calculated for current element depending on nSources
    !! if quadratic LSF matrix is singular fall back to linear
    integer, intent(inout) :: order
    !> Number of stencil directions
    integer, intent(in) :: QQ
    !> Number of dimensions
    integer, intent(in) :: nDims
    !> Number of sources from coarser found
    integer, intent(in) :: nSources
    !> Stencil directions
    real(kind=rk), intent(in) :: cxDirRK(3,QQ)
    !> direction in which sources are found
    integer, intent(in) :: neighDir(nSources)
    !> Pointer to position of interpolation matrix in growing array of matrix
    integer, intent(out) :: pos
    !> success if false if matrix is singular reduce interpolation order
    logical, intent(out) :: success
    ! --------------------------------------------------------------------------
    integer :: hashID
    type(tem_matrix_type) :: matLSF_tmp
    logical :: wasAdded
    integer :: iNeigh, iSrc
    ! --------------------------------------------------------------------------
write(dbgUnit(4),"(A,I0)") 'Inside append least square fit matrix for intp ' &
  &                     //'order: ', order
    ! hashID to identify unique combination of available nSources
    hashID = 0
    do iSrc = 1, nSources
      iNeigh = neighDir(iSrc)
      hashID = ibset(hashID, iNeigh)
    end do
write(dbgUnit(4),"(A,I0)") '   hashID: ', hashID
    call append( me       = me%ID,    &
      &          val      = hashID,   &
      &          pos      = pos,      &
      &          wasAdded = wasAdded  )
write(dbgUnit(4),"(A,L5)") '   wasAdded: ', wasAdded
write(dbgUnit(4),"(A,I5)") '        pos: ', pos

    if (wasAdded) then

      select case(order)
      case(1) ! linear
        call build_matrixLSF_linearIntp( me       = matLSF_tmp, &
          &                              QQ       = QQ,         &
          &                              nDims    = nDims,      &
          &                              nSources = nSources,   &
          &                              cxDirRK  = cxDirRK,    &
          &                              neighDir = neighDir,   &
          &                              nCoeffs  = me%nCoeffs, &
          &                              success  = success     )
      case(2) ! quadratic
        call build_matrixLSF_quadIntp( me       = matLSF_tmp, &
          &                            QQ       = QQ,         &
          &                            nDims    = nDims,      &
          &                            nSources = nSources,   &
          &                            cxDirRK  = cxDirRK,    &
          &                            neighDir = neighDir,   &
          &                            nCoeffs  = me%nCoeffs, &
          &                            success  = success     )
      case default
        write(logUnit(1),*) 'Unsupported interpolation order to build matrix ' &
          &                 //'for LSF: ', order
      end select

      call append( me = me%isInvertible, val = success )
      call append( me = me%matArray, val = matLSF_tmp )
    else
      ! hashID already exist then return success depends on matrix isInvertible
      if (me%isInvertible%val(pos)) then
        success = .true.
      else
        success = .false.
      end if
    end if
write(dbgUnit(4),"(A,L5)") '     success  : ', success

  end subroutine append_intpMatrixLSF
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine builds up the matrix for least square fit used in
  !! quadratic interpolation.
  !! We extract momentum information completely on the view of the source
  !! coordinate system
  !! Set the right hand side of the equation system
  !! Solve the problem, where b = rhs, x = coefficients
  !! A*x = b
  !! overdetermined, solve the least Square fit problem
  !! (A^T)A*x = (A^T)b
  !! x = ((A^T)A)^-1*(A^T)b
  !! Solve linear system of equation with inverted matrix.
  !! Size of matrix: (nCoeffs, QQ)
  !! matrix_LSF = ((A^T)A)^-1*(A^T)
  subroutine build_matrixLSF_quadIntp(me, QQ, nDims, nSources, cxDirRK, &
    &                                 neighDir, nCoeffs, success)
    ! --------------------------------------------------------------------------
    !> Matrix to fill
    type(tem_matrix_type), intent(out) :: me
    !> Number of stencil directions
    integer, intent(in) :: QQ
    !> Number of dimensions
    integer, intent(in) :: nDims
    !> Number of sources from coarser found
    integer, intent(in) :: nSources
    !> Stencil directions
    real(kind=rk), intent(in) :: cxDirRK(3,QQ)
    !> direction in which sources are found
    integer, intent(in) :: neighDir(nSources)
    !> nUnknown coeffs
    integer, intent(in) :: nCoeffs
    !> success if false if matrix is singular reduce interpolation order
    logical, intent(out) :: success
    ! --------------------------------------------------------------------------
    integer :: iDir, iSrc
    !> Each row represents a polynomial evaluated at coord of elements in
    ! stencil directions
    type(tem_matrix_type) :: tmp_matrix
    real(kind=rk) :: inv_AtA(nCoeffs,nCoeffs), AtA(nCoeffs,nCoeffs)
    integer :: errCode
    ! --------------------------------------------------------------------------
write(dbgUnit(4),"(A)") 'Inside build least square fit matrix for quadratic'
    ! me%A = ((A^T)A)^-1*(A^T)
    ! inv_AtA = ((A^T)A)^-1
    call alloc_matrix(tmp_matrix, nSources, nCoeffs)
    select case(nDims)
    case (1)
      do iSrc = 1, nSources
        iDir = neighDir(iSrc)
        tmp_matrix%A(iSrc,:) = polyQuadratic_1D( cxDirRK(:,iDir) )
      end do

    case (2)
      do iSrc = 1, nSources
        iDir = neighDir(iSrc)
        tmp_matrix%A(iSrc,:) = polyQuadratic_2D( cxDirRK(:,iDir) )
      end do

    case (3)
      do iSrc = 1, nSources
        iDir = neighDir(iSrc)
        tmp_matrix%A(iSrc,:) = polyQuadratic_3D( cxDirRK(:,iDir) )
      end do

    case default
      write(logUnit(1),*) 'Unknown nDims for quadratic interpolation'
      call tem_abort()
    end select

write(dbgUnit(4),"(A)") ' tmp matrix:'
call tem_matrix_dump(tmp_matrix, dbgUnit(4))
flush(dbgUnit(4))


    AtA =  matmul( transpose(tmp_matrix%A), tmp_matrix%A)
    inv_AtA = invert_matrix(AtA, errCode)
    if (errCode == 0) then
      ! matrix_LSF size is transpose of tmp_matrix
      call alloc_matrix(me, nCoeffs, nSources)
      me%A = matmul( inv_AtA, transpose(tmp_matrix%A) )
      success = .true.
    else
      ! singular matrix, reduce interpolation order
      call alloc_matrix(me, 1, 1)
      success = .false.
      return
    end if

write(dbgUnit(4),"(A)") ' matrix LSF:'
call tem_matrix_dump(me, dbgUnit(4))
flush(dbgUnit(4))
    !write(*,*) 'Matrix_LSF '
    !do iDir = 1, nCoeffs
    !  write(*,*) 'iDir ', iDir, me%A(iDir, :)
    !end do

  end subroutine build_matrixLSF_quadIntp
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine builds up the matrix for least square fit used in
  !! linear interpolation.
  subroutine build_matrixLSF_linearIntp(me, QQ, nDims, nSources, cxDirRK, &
    &                                 neighDir, nCoeffs, success)
    ! --------------------------------------------------------------------------
    !> Matrix to fill
    type(tem_matrix_type), intent(out) :: me
    !> Number of stencil directions
    integer, intent(in) :: QQ
    !> Number of dimensions
    integer, intent(in) :: nDims
    !> Number of sources from coarser found
    integer, intent(in) :: nSources
    !> Stencil directions
    real(kind=rk), intent(in) :: cxDirRK(3,QQ)
    !> direction in which sources are found
    integer, intent(in) :: neighDir(nSources)
    !> nUnknown coeffs
    integer, intent(in) :: nCoeffs
    !> success if false if matrix is singular reduce interpolation order
    logical, intent(out) :: success
    ! --------------------------------------------------------------------------
    integer :: iDir, iSrc
    !> Each row represents a polynomial evaluated at coord of elements in
    ! stencil directions
    type(tem_matrix_type) :: tmp_matrix
    real(kind=rk) :: inv_AtA(nCoeffs,nCoeffs), AtA(nCoeffs,nCoeffs)
    integer :: errCode
    ! --------------------------------------------------------------------------
write(dbgUnit(4),"(A)") 'Inside build least square fit matrix for linear'
    ! me%A = ((A^T)A)^-1*(A^T)
    ! inv_AtA = ((A^T)A)^-1
    call alloc_matrix(tmp_matrix, nSources, nCoeffs)
    select case(nDims)
    case (1)
      do iSrc = 1, nSources
        iDir = neighDir(iSrc)
        tmp_matrix%A(iSrc,:) = polyLinear_1D( cxDirRK(:,iDir) )
      end do

    case (2)
      do iSrc = 1, nSources
        iDir = neighDir(iSrc)
        tmp_matrix%A(iSrc,:) = polyLinear_2D( cxDirRK(:,iDir) )
      end do

    case (3)
      do iSrc = 1, nSources
        iDir = neighDir(iSrc)
        tmp_matrix%A(iSrc,:) = polyLinear_3D( cxDirRK(:,iDir) )
      end do

    case default
      write(logUnit(1),*) 'Unknown nDims for quadratic interpolation'
      call tem_abort()
    end select

write(dbgUnit(4),"(A)") ' tmp matrix:'
call tem_matrix_dump(tmp_matrix, dbgUnit(4))
flush(dbgUnit(4))

    AtA =  matmul( transpose(tmp_matrix%A), tmp_matrix%A)
    inv_AtA = invert_matrix(AtA, errCode)
    if (errCode == 0) then
      ! matrix_LSF size is transpose of tmp_matrix
      call alloc_matrix(me, nCoeffs, nSources)
      me%A = matmul( inv_AtA, transpose(tmp_matrix%A) )
      success = .true.
    else
      ! singular matrix, reduce interpolation order
      call alloc_matrix(me, 1, 1)
      success = .false.
      return
    end if

write(dbgUnit(4),"(A)") ' Matrix LSF:'
call tem_matrix_dump(me, dbgUnit(4))
flush(dbgUnit(4))

  end subroutine build_matrixLSF_linearIntp
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This function returns matrix entries for quadratic polynomial for 1D
  !! stencil
  pure function polyQuadratic_1D(cxDir) result (phi)
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: cxDir(3)
    real(kind=rk) :: phi(3)
    ! --------------------------------------------------------------------------
    phi(1) = 1.0_rk
    phi(2) = cxDir(c_x)
    phi(3) = cxDir(c_x)**2
  end function polyQuadratic_1D
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This function returns matrix entries for quadratic polynomial for 2D
  !! stencil
  pure function polyQuadratic_2D(cxDir) result (phi)
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: cxDir(3)
    real(kind=rk) :: phi(6)
    ! --------------------------------------------------------------------------
    phi(1) = 1.0_rk
    phi(2) = cxDir(c_x)
    phi(3) = cxDir(c_y)
    phi(4) = cxDir(c_x)**2
    phi(5) = cxDir(c_y)**2
    phi(6) = cxDir(c_x)*cxDir(c_y)
  end function polyQuadratic_2D
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This function returns matrix entries for quadratic polynomial for 3D
  !! stencil
  pure function polyQuadratic_3D(cxDir) result (phi)
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: cxDir(3)
    real(kind=rk) :: phi(10)
    ! --------------------------------------------------------------------------
    phi( 1) = 1.0_rk
    phi( 2) = cxDir(c_x)
    phi( 3) = cxDir(c_y)
    phi( 4) = cxDir(c_z)
    phi( 5) = cxDir(c_x)**2
    phi( 6) = cxDir(c_y)**2
    phi( 7) = cxDir(c_z)**2
    phi( 8) = cxDir(c_x)*cxDir(c_y)
    phi( 9) = cxDir(c_y)*cxDir(c_z)
    phi(10) = cxDir(c_z)*cxDir(c_x)

  end function polyQuadratic_3D
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This function returns matrix entries for Linear polynomial for 1D
  !! stencil
  pure function polyLinear_1D(cxDir) result (phi)
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: cxDir(3)
    real(kind=rk) :: phi(2)
    ! --------------------------------------------------------------------------
    phi(1) = 1.0_rk
    phi(2) = cxDir(c_x)
  end function polyLinear_1D
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This function returns matrix entries for Linear polynomial for 2D
  !! stencil
  pure function polyLinear_2D(cxDir) result (phi)
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: cxDir(3)
    real(kind=rk) :: phi(3)
    ! --------------------------------------------------------------------------
    phi(1) = 1.0_rk
    phi(2) = cxDir(c_x)
    phi(3) = cxDir(c_y)
  end function polyLinear_2D
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This function returns matrix entries for Linear polynomial for 3D
  !! stencil
  pure function polyLinear_3D(cxDir) result (phi)
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: cxDir(3)
    real(kind=rk) :: phi(4)
    ! --------------------------------------------------------------------------
    phi( 1) = 1.0_rk
    phi( 2) = cxDir(c_x)
    phi( 3) = cxDir(c_y)
    phi( 4) = cxDir(c_z)

  end function polyLinear_3D
  ! ************************************************************************** !

  ! ************************************************************************** !
  subroutine truncate_intpMatrixLSF(me)
    ! --------------------------------------------------------------------------
    type(tem_intpMatrixLSF_type), intent(inout) :: me
    ! --------------------------------------------------------------------------
    call truncate(me%matArray)
    call truncate(me%ID)
    call truncate(me%isInvertible)
  end subroutine truncate_intpMatrixLSF
  ! ************************************************************************** !

  ! ************************************************************************** !
  subroutine destroy_intpMatrixLSF(me)
    ! --------------------------------------------------------------------------
    type(tem_intpMatrixLSF_type), intent(inout) :: me
    ! --------------------------------------------------------------------------
    call destroy(me%matArray)
    call destroy(me%ID)
    call truncate(me%isInvertible)
  end subroutine destroy_intpMatrixLSF
  ! ************************************************************************** !

  ! ************************************************************************** !
  subroutine tem_matrix_dump(me, outUnit)
    ! --------------------------------------------------------------------------
    type(tem_matrix_type), intent(in) :: me
    integer, intent(in) :: outUnit
    ! --------------------------------------------------------------------------
    integer :: iRow
    ! --------------------------------------------------------------------------
    write(outUnit, "(A,i2,A,i2)") 'Matrix dimension: ', &
      &               me%nEntries(1), ' x ', me%nEntries(2)
    do iRow = 1, me%nEntries(1)
      write(outUnit, *) 'iRow ', iRow ,'Val ', me%A(iRow, :)
    end do

  end subroutine tem_matrix_dump
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine provides assignment operator of tem_matrix_type
  subroutine copy_matrix(left, right)
    ! --------------------------------------------------------------------------
    type(tem_matrix_type), intent(out) :: left
    type(tem_matrix_type), intent(in) :: right
    ! --------------------------------------------------------------------------
    left%nEntries = right%nEntries
    if (allocated(right%A)) then
      allocate(left%A(right%nEntries(1), right%nEntries(2)))
      left%A = right%A
    end if

  end subroutine copy_matrix
! **************************************************************************** !

! **************************************************************************** !
  !> This routine allocates matrix with given dimentions
  subroutine alloc_matrix( me, dim1, dim2 )
    ! --------------------------------------------------------------------------
    type( tem_matrix_type ) :: me
    integer, intent(in) :: dim1
    integer, intent(in) :: dim2
    ! --------------------------------------------------------------------------

    if ( dim1 > 0 .and. dim2 > 0 ) then
      me%nEntries(1) = dim1
      me%nEntries(2) = dim2
      allocate( me%A( dim1, dim2 ) )
      me%A = 0.0_rk
    else
      write( logUnit(1), "(A)" ) &
        &  'Failed to allocate matrix. Dimension is negative number.'
      call tem_abort()
    end if

  end subroutine alloc_matrix
! **************************************************************************** !


! **************************************************************************** !
  !> Returns the inverse of a matrix calculated by finding the LU
  !! decomposition.  Depends on LAPACK.
  !!
  function invert_matrix(A, errCode) result(Ainv)
    ! ---------------------------------------------------------------------------
    !> Matrix to invert
    real(kind=rk), dimension(:,:), intent(in) :: A
    !> If error code is present return error code and do not abort
    integer, optional, intent(out) :: errCode
    !> inverse of A
    real(kind=rk), dimension(size(A,1),size(A,2)) :: Ainv
    ! ---------------------------------------------------------------------------
    real(kind=rk), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info
    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI
    ! ---------------------------------------------------------------------------
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
      if (present(errCode)) then
        errCode = info
        write(logUnit(5),*) 'WARNING: Matrix is numerically singular!'
        return
      else
        call tem_abort('Matrix is numerically singular!')
      end if
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
      if (present(errCode)) then
        errCode = info
        write(logUnit(5),*) 'WARNING: Matrix inversion failed!'
        return
      else
        call tem_abort('Matrix inversion failed!')
      end if
    end if

    ! successfull
    if (present(errCode)) errCode = 0
  end function invert_matrix
! **************************************************************************** !



end module tem_matrix_module

