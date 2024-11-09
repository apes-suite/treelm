! Copyright (c) 2012, 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2014 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013-2014 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014-2015 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
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
! ****************************************************************************** !
! summary: module that holds datatypes and routines for coordinate transformations.
!!
module tem_coordinate_module

  ! include treelm modules
!  use atl_equation_module, only: scalar_variable, vector3_variable
  use tem_varSys_module,  only: tem_varSys_type
  use tem_logging_module, only: logUnit

  implicit none

  private

  public :: coordRotation_type
  public :: initcoordinateRotation
  public :: xToXAxes, yToXAxes, zToXAxes
  public :: dirToString

  !> constant representing the x direction
  integer, parameter :: xDir = 1
  !> constant representing the y direction
  integer, parameter :: yDir = 2
  !> constant representing the z direction
  integer, parameter :: zDir = 3

  !> Identity transformation
  integer, parameter :: xToXAxes = 1
  !> Transformation to transform a y axis into an x axis
  integer, parameter :: yToXAxes = 2
  !> Transformation to transform a z axis into an x axis
  integer, parameter :: zToXAxes = 3

  !>  summary: datatype to transform varibales given in one coordinate
  !! system to another coordinate system.
  !! This can be usefule e.g. for calculating the fluxes which is always
  !! done across the surface in positive x direction.
  type coordRotation_type
    !> The type of the rotation, see parameters above.
    integer :: rotationType
    !> Array of integers defining how to transform the variables of your
    !! state vector. Therefore the size of this array is the number of scalar
    !! variables of your equation system.
    integer, allocatable :: varTransformIndices(:)
    !> Array of integers defining how to transform the derivatives of your
    !! state vector. Therefore the size of this array is the total number
    !! of derivatives you store in your state vector.
    integer, allocatable :: derTransformIndices(:)
  end type coordRotation_type

  contains

! ****************************************************************************** !
  !>  summary: routine to specify a coordinate transformation for the state
  !! variables of your equation system.
  !!
  subroutine initCoordinateRotation( varSys, coordTrans, derivatives,          &
    &                                rotation, dimen )
    ! ---------------------------------------------------------------------------
    !> The variables to build the permutations for.
    type(tem_varSys_type), intent(in) :: varSys
    !> The rotation you want to obtain. Please have a look at the parameters
    !! of this module to find a valid input argument.
    integer, intent(in) :: coordTrans
    !> The number of derivatives (already multidimensional) you need for your
    !! equation. Zero means that we calculate cell values only, one means all
    !! first derivatives and so on.
    integer, intent(in) :: derivatives
    !> The coordinate rotation you want to initialize.
    type(coordRotation_type), intent(out) :: rotation
    !> The spatial dimension of the system
    integer, intent(in) :: dimen
    ! ---------------------------------------------------------------------------
    integer :: prevScalarVars, iVar
    ! ---------------------------------------------------------------------------

    rotation%rotationType = coordTrans

    ! get the number of scalar variables

    allocate( rotation%varTransformIndices(varSys%nScalars) )

    ! lets iterate over the variables and specify the rotation indices first.
    prevScalarVars = 0
    do iVar = 1, varSys%nStateVars
      call appendRotatedVariable( varSys%method%val(iVar)%nComponents, &
        &                         prevScalarVars, coordTrans, rotation )
      prevScalarVars = prevScalarVars + varSys%method%val(iVar)%nComponents
    end do

    ! iterate over all the derivatives and specifiy rotations for the
    ! derivatives.
    allocate( rotation%derTransformIndices(3*derivatives+1) )
    call appendDerivative( derivatives, coordTrans, rotation, dimen )

  end subroutine
! ****************************************************************************** !


! ****************************************************************************** !
  !>  summary: subroutine to append a rotated derivative to coordRotation type.
  !!
  subroutine appendDerivative( derivatives , coordTrans, rotation, dimen )
    ! ---------------------------------------------------------------------------
    !> The number of derivatives of your equation (inclunding the zeroth order
    !! derivative).
    integer, intent(in) :: derivatives
    !> The coordinate transformation you apply.
    integer, intent(in) :: coordTrans
    !> The coordinate rotation you want to initialize.
    type(coordRotation_type), intent(inout) :: rotation
    !> The spatial dimension of the system
    integer, intent(in) :: dimen
    ! ---------------------------------------------------------------------------

    select case(dimen)

    ! 3D
    case(3)
      select case(derivatives)
      case(1)
        select case(coordTrans)
        case(xToXAxes)
          rotation%derTransformIndices(1:4) = [1,2,3,4]
        case(yToXAxes)
          rotation%derTransformIndices(1:4) = [1,3,4,2]
        case(zToXAxes)
          rotation%derTransformIndices(1:4) = [1,4,2,3]
        case default
          write(logUnit(1),*)'unknown coordinate rotation for derivatives, '//   &
            &            'stopping....'
          stop
        end select
      case(0)
        select case(coordTrans)
        case(xToXAxes)
          rotation%derTransformIndices(1) = 1
        case(yToXAxes)
          rotation%derTransformIndices(1) = 1
        case(zToXAxes)
          rotation%derTransformIndices(1) = 1
        case default
          write(logUnit(1),*)'unknown coordinate rotation for derivatives, '//   &
            &            'stopping....'
          stop
        end select
      case default
        write(logUnit(1),*)'not able to create coordinate rotations for this '// &
          &            'order of derivatives, stopping...'
        stop
      end select

   ! 2D
   case(2)
      select case(derivatives)
      case(1)
        select case(coordTrans)
        case(xToXAxes)
          rotation%derTransformIndices(1:3) = [1,2,3]
        case(yToXAxes)
          rotation%derTransformIndices(1:3) = [1,3,2]
        case default
          write(logUnit(1),*)'unknown coordinate rotation for derivatives, '//   &
            &            'stopping....'
          stop
        end select
      case(0)
        select case(coordTrans)
        case(xToXAxes)
          rotation%derTransformIndices(1) = 1
        case(yToXAxes)
          rotation%derTransformIndices(1) = 1
        case default
          write(logUnit(1),*)'unknown coordinate rotation for derivatives, '//   &
            &            'stopping....'
          stop
        end select
      case default
        write(logUnit(1),*)'not able to create coordinate rotations for this '// &
          &            'order of derivatives, stopping...'
        stop
      end select
   ! 1D
   case(1)
     rotation%derTransformIndices(:) = 1
    case default
      write(logUnit(1),*)'not able to create coordinate rotations for this '// &
          &            'spatial dimension, stopping...'
      stop
    end select

  end subroutine
! ****************************************************************************** !


! ****************************************************************************** !
  !>  summary: routine to append a rotated varibale to corrdRoation type.
  !!
  subroutine appendRotatedVariable( nComponents, prevScalarVars, coordTrans,   &
    &                               rotation )
    ! ---------------------------------------------------------------------------
    !> nComponents of variable you want to append.
    integer, intent(in) :: nComponents
    !> The number of scalar variable you append before you append this variable.
    integer, intent(in) :: prevScalarVars
    !> The coordinate transformation you apply.
    integer, intent(in) :: coordTrans
    !> The coordinate rotation you want to initialize.
    type(coordRotation_type), intent(inout) :: rotation
    ! ---------------------------------------------------------------------------
    integer, allocatable :: varIndices(:)
    ! ---------------------------------------------------------------------------

    allocate( varIndices(nComponents) )

    select case(nComponents)
    case(1)
      varIndices = rotateScalar()
      rotation%varTransformIndices(prevScalarVars+1:prevScalarVars+nComponents)&
        &             = varIndices + prevScalarVars
    case(2)
      varIndices = rotateVector2(coordTrans)
      rotation%varTransformIndices(prevScalarVars+1:prevScalarVars+nComponents)&
        &             = varIndices + prevScalarVars
    case(3)
      varIndices = rotateVector3(coordTrans)
      rotation%varTransformIndices(prevScalarVars+1:prevScalarVars+nComponents)&
        &             = varIndices + prevScalarVars
    case(9)
      varIndices = rotateTensor3(coordTrans)
      rotation%varTransformIndices(prevScalarVars+1:prevScalarVars+nComponents)&
        &             = varIndices + prevScalarVars

    case default
      write(logUnit(1),*)'not able to transform this variable to a new '//     &
        &            'coordinate system, stopping....'
      stop
    end select

  end subroutine
! ****************************************************************************** !


! ****************************************************************************** !
  !>  summary: rotate a vector in 3D by a given rotation.
  !!
  function rotateVector3(coordTrans) result(rotationIndices)
    ! ---------------------------------------------------------------------------
    !> The coordinate transformation you apply.
    integer, intent(in) :: coordTrans
    !> Rotation indices for the given transformation
    integer :: rotationIndices(3)
    ! ---------------------------------------------------------------------------

    select case(coordTrans)
    case(xToXAxes)
      ! this is just the identity transformation
      rotationIndices = [1,2,3]
    case(yToXAxes)
      rotationIndices = [2,3,1]
    case(zToXAxes)
      rotationIndices = [3,1,2]
    case default
      write(logUnit(1),*)'unknown coordinate rotation, stopping....'
      stop
    end select

  end function
! ****************************************************************************** !

! ****************************************************************************** !
  !>  summary: rotate a tensor in 3D by a given rotation.
  !!
  function rotateTensor3(coordTrans) result(rotationIndices)
    ! ---------------------------------------------------------------------------
    !> The coordinate transformation you apply.
    integer, intent(in) :: coordTrans
    !> Rotation indices for the given transformation
    integer :: rotationIndices(9)
    ! ---------------------------------------------------------------------------

    select case(coordTrans)
    case(xToXAxes)
      ! this is just the identity transformation
      rotationIndices = [1,2,3,4,5,6,7,8,9]
    case(yToXAxes)
      rotationIndices = [5,6,4,8,9,7,2,3,1]
    case(zToXAxes)
      rotationIndices = [9,7,8,3,1,2,6,4,5]
    case default
      write(logUnit(1),*)'unknown coordinate rotation, stopping....'
      stop
    end select

  end function
! ****************************************************************************** !

! ****************************************************************************** !
  !>  summary: rotate a vector in 2D by a given rotation.
  !!
  function rotateVector2(coordTrans) result(rotationIndices)
    ! ---------------------------------------------------------------------------
    !> The coordinate transformation you apply.
    integer, intent(in) :: coordTrans
    !> Rotation indices for the given transformation
    integer :: rotationIndices(2)
    ! ---------------------------------------------------------------------------

    select case(coordTrans)
    case(xToXAxes)
      ! this is just the identity transformation
      rotationIndices = [1,2]
    case(yToXAxes)
      rotationIndices = [2,1]
    case(zToXAxes)
      ! rotation to z axes does not make sense.
      ! Nevertheless we define it for convenience.
      rotationIndices = [1,2]
    case default
      write(logUnit(1),*)'unknown coordinate rotation, stopping....'
      stop
    end select

  end function
! ****************************************************************************** !


! ****************************************************************************** !
  !>  summary: rotate a scalar?
  !!
  function rotateScalar() result(rotationIndices)
    ! ---------------------------------------------------------------------------
    !> The coordinate transformation you apply.
    ! integer, intent(in) :: coordTrans
    !> Rotation indices for the given transformation
    integer :: rotationIndices(1)
    ! ---------------------------------------------------------------------------
    rotationIndices = [1]

  end function rotateScalar
! ****************************************************************************** !


! ****************************************************************************** !
  !>  summary: function to convert a direction to a string.
  !!
  function dirToString(direction) result(dirAsChar)
    ! ---------------------------------------------------------------------------
    !> direction to convert
    integer, intent(in) :: direction
    !> direction as string
    character(len=1) :: dirAsChar
    ! ---------------------------------------------------------------------------
    select case(direction)
    case(xDir)
      dirAsChar = 'X'
    case(yDir)
      dirAsChar = 'Y'
    case(zDir)
      dirAsChar = 'Z'
    case default
      write(logUnit(1),*)'unknown direction for conversion to string.'
      write(logUnit(1),*)'stopping...'
      stop
    end select

  end function
! ****************************************************************************** !


end module tem_coordinate_module
! ****************************************************************************** !
