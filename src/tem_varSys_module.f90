! Copyright (c) 2014-2016, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014-2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014-2016, 2018-2019 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016-2017 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Daniel Fleischer <daniel.fleischer@student.uni-siegen.de>
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
!> The variable system allows the description of arbitrary quantities and their
!! relations.
!!
!! Its meant to enable the extraction of most of the data in the solver, but
!! also to describe additional quantities to be used by the solvers.
module tem_varSys_module
  use iso_c_binding,            only: c_ptr
  use env_module,               only: labelLen, rk, minLength, zeroLength
  use tem_time_module,          only: tem_time_type
  use tem_dyn_array_module,     only: dyn_labelArray_type, init, append, &
    &                                 empty, PositionOfVal
  use treelmesh_module,         only: treelmesh_type
  use tem_aux_module,           only: tem_abort
  use tem_logging_module,       only: logUnit, llerror, lldebug

  use aotus_module,             only: flu_state, aot_get_val
  use aot_table_module,         only: aot_table_open, aot_table_close, &
    &                                 aot_table_length, aot_get_val
  use aot_out_module,           only: aot_out_type, aot_out_open, &
    &                                 aot_out_close, aot_out_val, &
    &                                 aot_out_open_table, aot_out_close_table

  implicit none

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

  ! **************************************************************************** !
  !> Description of the method how to obtain a variable
  type tem_varSys_op_type

    !> Position of this variable in the variable system.
    integer :: mypos

    !> Position of state variable in the state array
    !! Currently used only in MUSUBI to access one dimensional state array
    integer, allocatable :: state_varPos(:)

    !> Position of auxiliary variable in the auxilied field array
    !! Currently used only in MUSUBI to access one dimensional auxiliary array
    !! In Musubi, auxField vars are conserved macroscopic quantities computed
    !! from PDF state
    integer, allocatable :: auxField_varPos(:)

    !> Number of components for this variable.
    integer :: nComponents

    !> Number of variables, that are needed as input for the operation
    !! to obtain the variable.
    integer :: nInputs

    !> Position of the input variables in the variable system.
    !!
    !! There are as many entries as nInputs.
    integer, allocatable :: input_varPos(:)

    !> Component index of the input variable in the variable system.
    !! It is used only when there is only one input variable.
    !! Index values must not be zero and > nComponents of input variable
    integer, allocatable :: input_varIndex(:)

    !> Data that is required by the get method.
    type(c_ptr) :: method_data

    !> Operation type
    character(len=labelLen) :: operType

    !> Function to actually obtain the variable at a given point.
    !!
    !! This is either a function accessing a state variable directly,
    !! a function returning a space time function evaluation or a
    !! derived quantity, that computes a new variable out of others.
    procedure(tem_varSys_proc_point), pointer :: get_point => null()

    !> Function to actually obtain the variable in a given element.
    procedure(tem_varSys_proc_element), pointer :: get_element => null()

    !> Function to set parameter in the data_type stored in method_data.
    procedure(tem_varSys_proc_setParams), pointer :: set_params => null()

    !> Function to get parameter in the data_type stored in method_data
    procedure(tem_varSys_proc_getParams), pointer :: get_params => null()

    !> Function to setup points set for boundaries and sources.
    !! Pointe set are stored in method_data level wise 1D growing array for
    !! dimension X,Y and Z.
    !!  * For solver variables, points are stored in solver container.
    !!  * For spacetime variables, points are stored in spacetime function.
    !!  * For operation variables, points are passed down to its
    !!    input_variable.
    procedure(tem_varSys_proc_setupIndices), pointer :: &
      &                                      setup_indices => null()

    !> Function to get value for point set stored in method_data
    !! for requested index in point set.
    !! This function either returns a pre-stored value or compute value
    !! depends on variable type and spacetime function.
    !! For time-independent spacetime function, values are computed
    !! in setupIndices and growing array of points are deleted
    procedure(tem_varSys_proc_getValOfIndex), pointer :: &
      &                                       get_valOfIndex => null()

  end type tem_varSys_op_type
  ! **************************************************************************** !


  !> growing array type for type(tem_varsys_op_type)
  type grw_varoparray_type
    integer :: nvals = 0
    integer :: containersize = 0
    type(tem_varsys_op_type), allocatable :: val(:)
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_ga_varop
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_ga_varop
  end interface

  !> empty the entries  without changing arrays
  interface empty
    module procedure empty_ga_varop
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_ga_varop
  end interface

  !> insert an element at a given position
  interface placeat
    module procedure placeat_ga_varop
    module procedure placeat_ga_varop_vec
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_ga_varop
    module procedure append_ga_varop_vec
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_ga_varop
  end interface



  ! **************************************************************************** !
  !> Description of the variable system.
  type tem_varSys_type
    !> A descriptive name for this system of variables.
    character(len=LabelLen) :: SystemName

    !> Number of variables in the state.
    integer :: nStateVars = 0

    !> Number of scalars in the state.
    !!
    !! This keeps track of the length of the state array.
    integer :: nScalars = 0

    !> Number of auxField variables
    integer :: nAuxVars = 0

    !> Number of scalars in the auxField
    !! This keeps track of the length of the auxField array
    integer :: nAuxScalars

    !> Definition of how to obtain a variable.
    type(grw_varOpArray_type) :: method

    !> List of variables in the system.
    type(dyn_labelArray_type) :: varname

  end type tem_varSys_type
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> A supporting data type to define a solver specific element evaluation for
  !! stfuns.
  type tem_varSys_solverData_evalElem_type
    !> Data from the solver for the evaluation.
    type(c_ptr) :: solver_bundle
    !> Callback function to set the appropriate function for stFun vars.
    procedure(tem_varsys_set_evalelem), pointer :: stFun_setter => null()
    !> Callback function to set the appropriate function for solver specific
    !! operation vars.
    procedure(tem_varsys_set_evalelem), pointer :: opVar_setter => null()
  end type tem_varSys_solverData_evalElem_type
  ! **************************************************************************** !


  ! **************************************************************************** !
  abstract interface
    !> Interface description for a variable access method (single point).
    !!
    !! To obtain values of a given variable, it is necessary to state the
    !! space and time at which the variable should be evaluated.
    !! The interface is vectorized and provides n values at n different space
    !! locations at the same time.
    !! Of course the variable system itself also needs to be passed in, to
    !! allow the computation of other derived quantities as needed.
    !! The method description itself is passed in automatically, and has not
    !! to be provided explicitly.
    subroutine tem_varSys_proc_point(fun, varsys, point, time, tree, nPnts, &
      &                              res)
      import :: tem_varSys_op_type, tem_varSys_type, tem_time_type, rk,        &
        &       treelmesh_type

      !> Description of the method to obtain the variables, here some preset
      !! values might be stored, like the space time function to use or the
      !! required variables.
      class(tem_varSys_op_type), intent(in) :: fun

      !> The variable system to obtain the variable from.
      type(tem_varSys_type), intent(in) :: varSys

      !> Three-dimensional coordinates at which the variable should be
      !! evaluated. Only useful for variables provided as space-time functions.
      real(kind=rk), intent(in) :: point(:,:)

      !> Point in time at which to evaluate the variable.
      type(tem_time_type), intent(in)  :: time

      !> global treelm mesh info
      type(treelmesh_type), intent(in) :: tree

      !> Number of values to obtain for this variable (vectorized access).
      integer, intent(in) :: nPnts

      !> Resulting values for the requested variable.
      !!
      !! Dimension: n requested entries x nComponents of this variable
      !! Access: (iElem-1)*fun%nComponents + iComp
      real(kind=rk), intent(out) :: res(:)
    end subroutine tem_varSys_proc_point


    !> Interface description for a variable access method (complete element).
    !!
    !! To obtain values of a given variable, it is necessary to state the
    !! treeID and time at which the variable should be evaluated.
    !! The interface is nDofs values to cover the all degrees of freedoms
    !! in the element.
    !! Of course the variable system itself also needs to be passed in, to
    !! allow the computation of other derived quantities as needed.
    !! The method description itself is passed in automatically, and has not
    !! to be provided explicitly.
    subroutine tem_varSys_proc_element(fun, varsys, elempos, time, tree, &
      &                                nElems, nDofs, res)
      import :: tem_varSys_op_type, tem_varSys_type, tem_time_type, rk, &
        &       treelmesh_type

      !> Description of the method to obtain the variables, here some preset
      !! values might be stored, like the space time function to use or the
      !! required variables.
      class(tem_varSys_op_type), intent(in) :: fun

      !> The variable system to obtain the variable from.
      type(tem_varSys_type), intent(in) :: varSys

      !> Position of element in tree%treeID to get the variable for.
      integer, intent(in) :: elempos(:)

      !> Point in time at which to evaluate the variable.
      type(tem_time_type), intent(in)  :: time

      !> global treelm mesh info
      type(treelmesh_type), intent(in) :: tree

      !> Number of elements to obtain for this variable (vectorized access).
      integer, intent(in) :: nElems

      !> Number of degrees of freedom within an element.
      integer, intent(in) :: nDofs

      !> Resulting values for the requested variable.
      !!
      !! Linearized array dimension:
      !! (nComponents of resulting variable) x (nDegrees of freedom) x (nElems)
      !! Access: (iElem-1)*fun%nComponents*nDofs +
      !!         (iDof-1)*fun%nComponents + iComp
      real(kind=rk), intent(out) :: res(:)
    end subroutine tem_varSys_proc_element


    !> Interface description for a variable to set parameter in data type
    !! stored in method_data.
    !!
    !! To set the parameter, provide a string (recommended in lua format).
    !! Depends on variable procedure, input string will be processed and
    !! information will be store in method_data.
    !! For operation variable, pass the information down to its
    !! input_variable.
    subroutine tem_varSys_proc_setParams(fun, varSys, instring)
      import :: tem_varSys_op_type, tem_varSys_type

      !> Description of the method to obtain the variables, here some preset
      !! values might be stored, like the space time function to use or the
      !! required variables.
      class(tem_varSys_op_type), intent(in) :: fun

      !> The variable system to obtain the variable from.
      type(tem_varSys_type), intent(in) :: varSys

      !> Input string with parameter to set in method_data
      character(len=*), intent(in) :: instring
    end subroutine tem_varSys_proc_setParams

    !> Interface description for a variable to get parameter in data type
    !! stored in method_data.
    !!
    !! To get the parameter, provide requested parameter via
    !! instring (recommended in lua format) and routine returns the parameter
    !! via outstring in lua format.
    !! For operation variable, pass the information down to its
    !! input_variable.
    subroutine tem_varSys_proc_getParams(fun, varSys, instring, outstring)
      import :: tem_varSys_op_type, tem_varSys_type

      !> Description of the method to obtain the variables, here some preset
      !! values might be stored, like the space time function to use or the
      !! required variables.
      class(tem_varSys_op_type), intent(in) :: fun

      !> The variable system to obtain the variable from.
      type(tem_varSys_type), intent(in) :: varSys

      !> Input string with parameter to set in method_data
      character(len=*), intent(in) :: instring

      !> Output string with requested parameter value from method_data
      character(len=*), intent(out) :: outstring
    end subroutine tem_varSys_proc_getParams

    !> Interface description for a variable to setup point sets and return
    !! indices of points in the 1D growing array.
    !!
    !! Before setupIndices, use setParams to set information about boundary
    !! or sources in the spacetime function type.
    !! Points are stored in the 1D growing array for each dimension in the
    !! level wise list in the method_data.
    !! offset_bit is required for each point on boundaries so the remote domain
    !! can use this offset to identify rank for the point.
    !! If offset_bit is not present then default is set to element center
    !! i.e offset_bit = achar(1+4+16).
    !!
    !! For state variable, this routine will also compute elemPos in the
    !! given level and local_coord for every point. Elempos and local_coord are
    !! stored in same place as growing array of points.
    !!
    !! For spacetime function variable, if function is time independent then
    !! values are computed and stored level wise after all points are
    !! added to spacetime function.
    subroutine tem_varSys_proc_setupIndices(fun, varSys, point, offset_bit,    &
      &                                     iLevel, tree, nPnts, idx)
      import :: tem_varSys_op_type, tem_varSys_type, rk, treelmesh_type

      !> Description of the method to obtain the variables, here some preset
      !! values might be stored, like the space time function to use or the
      !! required variables.
      class(tem_varSys_op_type), intent(in) :: fun

      !> The variable system to obtain the variable from.
      type(tem_varSys_type), intent(in) :: varSys

      !> List of space coordinate points to store as growing array in
      !! method_data
      real(kind=rk), intent(in) :: point(:,:)

      !> Offset bit encoded as character for every point.
      !!
      !! Offset integer coord(3) is converted into a character with
      !! offset_bit = achar( (coord(1)+1) + (coord(2)+1)*4 + (coord(3)+1)*16 )
      !! Backward transformation form character to 3 integer:
      !! coord(1) = mod(ichar(offset_bit),4) - 1
      !! coord(2) = mod(ichar(offset_bit),16)/4 - 1
      !! coord(3) = ichar(offset_bit)/16 - 1
      !!
      !! If not present default is to center i.e offset_bit = achar(1+4+16)
      character, optional, intent(in) :: offset_bit(:)

      !> Level to which input points belong to
      integer, intent(in) :: iLevel

      !> global treelm mesh info
      type(treelmesh_type), intent(in) :: tree

      !> Number of points to add in method_data of this variable
      integer, intent(in) :: nPnts

      !> Index of points in the growing array and variable val array.
      !! Size: n
      !!
      !! This must be stored in boundary or source depends on who
      !! calls this routine.
      !! This index is required to return a value using getValOfIndex.
      integer, intent(out) :: idx(:)
    end subroutine tem_varSys_proc_setupIndices


    !> Interface description for a variable to return a value at the given
    !! index position in the growing array points set stored in method_data.
    !!
    !! For spacetime function, if value is pre-stored, it will return a value
    !! at given index else if will evaluate a variable at index of a point
    !! for a given time and return a value.
    !!
    !! If index is not present and first is present then index is computed
    !! from first and number of return value (n).
    subroutine tem_varSys_proc_getValOfIndex(fun, varSys, time, iLevel, &
      &                                      idx, idxLen, nVals, res)
      import :: tem_varSys_op_type, tem_varSys_type, tem_time_type, rk

      !> Description of the method to obtain the variables, here some preset
      !! values might be stored, like the space time function to use or the
      !! required variables.
      class(tem_varSys_op_type), intent(in) :: fun

      !> The variable system to obtain the variable from.
      type(tem_varSys_type), intent(in) :: varSys

      !> Point in time at which to evaluate the variable.
      type(tem_time_type), intent(in) :: time

      !> Level on which values are requested
      integer, intent(in) :: iLevel

      !> Index of points in the growing array and variable val array to
      !! return.
      !! Size: most times nVals, if contiguous arrays are used it depends
      !! on the number of first indices
      integer, intent(in) :: idx(:)

      !> With idx as start index in contiguous memory,
      !! idxLength defines length of each contiguous memory
      !! Size: dependes on number of first index for contiguous array,
      !! but the sum of all idxLen is equal to nVals
      integer, optional, intent(in) :: idxLen(:)

      !> Number of values to obtain for this variable (vectorized access).
      integer, intent(in) :: nVals

      !> Resulting values for the requested variable.
      !!
      !! Dimension: n requested entries x nComponents of this variable
      !! Access: (iElem-1)*fun%nComponents + iComp
      real(kind=rk), intent(out) :: res(:)
    end subroutine tem_varSys_proc_getValOfIndex


    subroutine tem_varsys_set_evalelem(set_elem_eval, fun)
      import :: tem_varSys_op_type, tem_varSys_solverData_evalElem_type
      !> Description on how to set the element retrieval function for stfuns.
      !! and solver specific operation variables
      class(tem_varSys_solverData_evalElem_type), intent(in) :: set_elem_eval

      !> Description of the method to obtain the variables, here some preset
      !! values might be stored, like the space time function to use or the
      !! required variables.
      type(tem_varSys_op_type), intent(inout) :: fun
    end subroutine tem_varsys_set_evalelem

  end interface
  ! **************************************************************************** !


  ! **************************************************************************** !
  interface tem_varSys_load
    module procedure tem_varSys_load_vector
    module procedure tem_varSys_load_single
  end interface tem_varSys_load
  ! **************************************************************************** !

  ! **************************************************************************** !
  interface tem_varSys_dump
    module procedure tem_varSys_dump_vector
    module procedure tem_varSys_dump_single
  end interface tem_varSys_dump
  ! **************************************************************************** !

  ! **************************************************************************** !
  interface tem_varSys_out
    module procedure tem_varSys_out_vector
    module procedure tem_varSys_out_single
  end interface tem_varSys_out
  ! **************************************************************************** !

  interface assignment(=)
    module procedure copy_varOp
  end interface


contains


  subroutine init_ga_varop(me, length)
    type(grw_varoparray_type), intent(out) :: me !< dynamic array to init
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

  end subroutine init_ga_varop

  subroutine destroy_ga_varop(me)
    type(grw_varoparray_type), intent(inout) :: me !< dynamic array to destroy

    me%containersize = 0
    me%nvals = 0
    if( allocated( me%val ) ) deallocate(me%val)

  end subroutine destroy_ga_varop


  subroutine truncate_ga_varop(me)
    !------------------------------------------------------------------------
    type(grw_varoparray_type) :: me !< array to truncate
    !------------------------------------------------------------------------
    type(tem_varsys_op_type), allocatable :: tarray(:)
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

  end subroutine truncate_ga_varop


  subroutine empty_ga_varop(me)
    !------------------------------------------------------------------------
    type(grw_varoparray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------

    me%nvals = 0

  end subroutine empty_ga_varop

  !> adds the value to a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! element at the requested position will be replaced.
  subroutine placeat_ga_varop(me, val, pos, length)
    type(grw_varoparray_type) :: me !< array to place the value into
    type(tem_varsys_op_type), intent(in) :: val !< value to place at the given position
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

  end subroutine placeat_ga_varop


  !> adds the values starting from a given position inside the growing array.
  !!
  !! if the requested position is outside the current array bounds, the array
  !! will be resized accordingly. if it is inside the current array bounds, the
  !! elements starting from the requested position will be replaced up to
  !! the element at position `pos + size(val) - 1`.
  subroutine placeat_ga_varop_vec(me, val, pos, length)
    type(grw_varoparray_type) :: me !< array to append the value to
    type(tem_varsys_op_type), intent(in) :: val(:) !< values to append
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

  end subroutine placeat_ga_varop_vec


  subroutine append_ga_varop(me, val, length)
    type(grw_varoparray_type) :: me !< array to append the value to
    type(tem_varsys_op_type), intent(in) :: val !< value to append
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

  end subroutine append_ga_varop

  subroutine append_ga_varop_vec(me, val, length)
    type(grw_varoparray_type) :: me !< array to append the value to
    type(tem_varsys_op_type), intent(in) :: val(:) !< values to append
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

  end subroutine append_ga_varop_vec


  subroutine expand_ga_varop(me, pos, length)
    type(grw_varoparray_type) :: me !< array to resize
    integer, intent(in), optional :: pos !< optional predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    type(tem_varsys_op_type), allocatable :: swpval(:)
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

  end subroutine expand_ga_varop


  ! **************************************************************************** !
  !> Initialize a variable system.
  subroutine tem_varSys_init(me, systemName, length)
    ! --------------------------------------------------------------------------!
    !> Variable system to initialize.
    type(tem_varSys_type),           intent(out) :: me

    !> A descriptive name for this system of variables.
    character(len=*),                intent(in)  :: systemName

    !> Initial length for the arrays used in the variable system.
    integer,               optional, intent(in)  :: length
    ! --------------------------------------------------------------------------!

    me%systemName = trim(systemName)
    me%nStateVars = 0
    me%nScalars = 0
    me%nAuxVars = 0
    me%nAuxScalars = 0

    call init( me     = me%varname, &
      &        length = length      )

    call init( me     = me%method, &
      &        length = length     )

  end subroutine tem_varSys_init
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Append a new state variable to the variable system.
  !!
  !! Note that the actual data on how to access the state has to be encoded
  !! in the method_data pointer.
  !! The callers therefore have to ensure the consistency of those access
  !! definitions.
  !! nComponents has to be at least 1.
  subroutine tem_varSys_append_stateVar(me, varname, nComponents, method_data, &
    &                                   get_point, get_element, set_params,    &
    &                                   get_params, setup_indices,             &
    &                                   get_valOfIndex, pos, wasAdded)
    ! --------------------------------------------------------------------------!
    !> Variable system to append the state variable to.
    type(tem_varSys_type),             intent(inout) :: me

    !> Variable to append to the state.
    character(len=*),                  intent(in)    :: varname

    !> Number of components in this variable.
    integer, intent(in) :: nComponents

    !> Data that is required by the methods to obtain the variable.
    type(c_ptr), intent(in) :: method_data

    !> Procedure which allows the retrieval of the variable at given points.
    procedure(tem_varSys_proc_point), pointer :: get_point

    !> Procedure which allows the retrieval of the variable in an element.
    procedure(tem_varSys_proc_element), pointer :: get_element

    !> Procedure which allows to set parameter in method_data
    procedure(tem_varSys_proc_setParams), optional, pointer :: set_params

    !> Procedure which allows to get parameter in method_data
    procedure(tem_varSys_proc_getParams), optional, pointer :: get_params

    !> Procedure to setup growing array of points, variable value
    !! in method_data and return index of points set
    procedure(tem_varSys_proc_setupIndices), pointer :: setup_indices

    !> Procedure which allows to retrieval of the variable at point or val
    !! array index
    procedure(tem_varSys_proc_getValOfIndex), pointer :: get_valOfIndex

    !> Position of the variable in the system.
    integer,                 optional, intent(out)   :: pos

    !> Indicator, if the variable was actually added to the system.
    logical,                 optional, intent(out)   :: wasAdded
    ! --------------------------------------------------------------------------!
    logical :: newVar
    type(tem_varSys_op_type) :: state_access
    integer :: iComp, varPos
    ! --------------------------------------------------------------------------!

    call append( me       = me%varname, &
      &          val      = varname,    &
      &          pos      = varPos,     &
      &          wasAdded = newVar      )

    if (newVar) then
      state_access%mypos = varPos
      state_access%operType = 'state'
      state_access%nInputs = 0
      state_access%nComponents = nComponents
      allocate(state_access%state_varPos(nComponents))
      do iComp = 1, nComponents
        state_access%state_varPos(iComp) = me%nScalars + iComp
      end do
      allocate(state_access%input_varPos(0))

      state_access%method_data = method_data
      state_access%get_point => get_point
      state_access%get_element => get_element
      if (present(set_params)) then
        state_access%set_params => set_params
      else
        state_access%set_params => tem_varSys_setParams_dummy
      end if
      if (present(get_params)) then
        state_access%get_params => get_params
      else
        state_access%get_params => tem_varSys_getParams_dummy
      end if
      state_access%setup_indices => setup_indices
      state_access%get_valOfIndex => get_valOfIndex
      call append( me  = me%method,   &
        &          val = state_access )
      me%nStateVars = me%nStateVars+1
      me%nScalars = me%nScalars + nComponents
    end if

    if (present(pos)) pos = varPos
    if (present(wasAdded)) wasAdded = newVar

  end subroutine tem_varSys_append_stateVar
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Append a new auxiliaryField variable to the variable system.
  !!
  !! Note that the actual data on how to access the state has to be encoded
  !! in the method_data pointer.
  !! The callers therefore have to ensure the consistency of those access
  !! definitions.
  !! nComponents has to be at least 1.
  subroutine tem_varSys_append_auxFieldVar(me, varname, nComponents, &
    & method_data, get_point, get_element, set_params, get_params,   &
    & setup_indices, get_valOfIndex, pos, wasAdded)
    ! --------------------------------------------------------------------------!
    !> Variable system to append the state variable to.
    type(tem_varSys_type),             intent(inout) :: me

    !> Variable to append to the state.
    character(len=*),                  intent(in)    :: varname

    !> Number of components in this variable.
    integer, intent(in) :: nComponents

    !> Data that is required by the methods to obtain the variable.
    type(c_ptr), intent(in) :: method_data

    !> Procedure which allows the retrieval of the variable at given points.
    procedure(tem_varSys_proc_point), pointer :: get_point

    !> Procedure which allows the retrieval of the variable in an element.
    procedure(tem_varSys_proc_element), pointer :: get_element

    !> Procedure which allows to set parameter in method_data
    procedure(tem_varSys_proc_setParams), optional, pointer :: set_params

    !> Procedure which allows to get parameter in method_data
    procedure(tem_varSys_proc_getParams), optional, pointer :: get_params

    !> Procedure to setup growing array of points, variable value
    !! in method_data and return index of points set
    procedure(tem_varSys_proc_setupIndices), pointer :: setup_indices

    !> Procedure which allows to retrieval of the variable at point or val
    !! array index
    procedure(tem_varSys_proc_getValOfIndex), pointer :: get_valOfIndex

    !> Position of the variable in the system.
    integer,                 optional, intent(out)   :: pos

    !> Indicator, if the variable was actually added to the system.
    logical,                 optional, intent(out)   :: wasAdded
    ! --------------------------------------------------------------------------!
    logical :: newVar
    type(tem_varSys_op_type) :: auxField_access
    integer :: iComp, varPos
    ! --------------------------------------------------------------------------!

    call append( me       = me%varname, &
      &          val      = varname,    &
      &          pos      = varPos,     &
      &          wasAdded = newVar      )

    if (newVar) then
      auxField_access%mypos = varPos
      auxField_access%operType = 'auxfield'
      auxField_access%nInputs = 0
      auxField_access%nComponents = nComponents
      allocate(auxField_access%auxField_varPos(nComponents))
      do iComp = 1, nComponents
        auxField_access%auxField_varPos(iComp) = me%nAuxScalars + iComp
      end do
      allocate(auxField_access%input_varPos(0))

      auxField_access%method_data = method_data
      auxField_access%get_point => get_point
      auxField_access%get_element => get_element
      if (present(set_params)) then
        auxField_access%set_params => set_params
      else
        auxField_access%set_params => tem_varSys_setParams_dummy
      end if
      if (present(get_params)) then
        auxField_access%get_params => get_params
      else
        auxField_access%get_params => tem_varSys_getParams_dummy
      end if
      auxField_access%setup_indices => setup_indices
      auxField_access%get_valOfIndex => get_valOfIndex
      call append( me  = me%method,   &
        &          val = auxField_access )
      me%nAuxVars = me%nAuxVars+1
      me%nAuxScalars = me%nAuxScalars + nComponents
    end if

    if (present(pos)) pos = varPos
    if (present(wasAdded)) wasAdded = newVar

  end subroutine tem_varSys_append_auxFieldVar
  ! **************************************************************************** !

  ! **************************************************************************** !
  !> Append a new (non-state) variable to the variable system.
  !!
  !! Non-state variables might depend on other variables, and do not modify the
  !! internal state variable counter.
  subroutine tem_varSys_append_derVar(me, varname, operType, nComponents,    &
    &                                 input_varname, input_varindex,         &
    &                                 method_data, get_point, get_element,   &
    &                                 set_params, get_params, setup_indices, &
    &                                 get_valOfIndex, pos, wasAdded          )
    ! --------------------------------------------------------------------------!
    !> Variable system to append the state variable to.
    type(tem_varSys_type), intent(inout) :: me

    !> Variable to append to the state.
    character(len=*), intent(in) :: varname

    !> Operation type
    character(len=*), optional, intent(in) :: operType

    !> Number of components in this variable.
    integer, intent(in) :: nComponents

    !> List of variable names, this variable depends on.
    !!
    !! (Might be empty). The variable will only be appended, if all invars are
    !! already found in the variable system.
    character(len=*), optional, intent(in) :: input_varname(:)

    !> Component index to access from input variable nComponents
    !!
    !! (Might be empty). The variable will only be appended, if all index are
    !! available in input_varname
    integer, optional, intent(in) :: input_varIndex(:)

    !> Data that is required by the methods to obtain the variable.
    type(c_ptr), intent(in) :: method_data

    !> Procedure which allows the retrieval of the variable at given points.
    procedure(tem_varSys_proc_point), pointer :: get_point

    !> Procedure which allows the retrieval of the variable in an element.
    procedure(tem_varSys_proc_element), pointer :: get_element

    !> Procedure which allows to set parameter in method_data
    procedure(tem_varSys_proc_setParams), optional, pointer :: set_params

    !> Procedure which allows to get parameter in method_data
    procedure(tem_varSys_proc_getParams), optional, pointer :: get_params

    !> Procedure to setup growing array of points, variable value
    !! in method_data and return index of points set
    procedure(tem_varSys_proc_setupIndices), pointer :: setup_indices

    !> Procedure which allows to retrieval of the variable at point or val
    !! array index
    procedure(tem_varSys_proc_getValOfIndex), pointer :: get_valOfIndex

    !> Position of the variable in the system.
    integer, optional, intent(out)   :: pos

    !> Indicator, if the variable was actually added to the system.
    logical,  optional, intent(out)   :: wasAdded
    ! --------------------------------------------------------------------------!
    logical :: newVar, nonexisting_input
    type(tem_varSys_op_type) :: var_access
    integer :: iIn, nInputs, varPos, in_nComp, nVarIndex
    integer, allocatable :: inpos(:), input_varIndex_loc(:)
    ! --------------------------------------------------------------------------!

    if (present(input_varname)) then
      nInputs = size(input_varname)
    else
      nInputs = 0
    end if

    if (present(input_varindex)) then
      nVarIndex = size(input_varindex)
    else
      nVarIndex = 0
    end if

    allocate(inpos(nInputs))
    allocate(input_varIndex_loc(nVarIndex))

    nonexisting_input = .false.

    ! Check whether all input variables are present. If not, the variable can't
    ! be added, because variables it depends on are not present.
    do iIn=1,nInputs
      inpos(iIn) = PositionOfVal(me%varname, input_varname(iIn))
      nonexisting_input = (inpos(iIn) <= 0)
      if (nonexisting_input) then
        write(logUnit(llerror),*) 'Variable ' // trim(varname) &
          & // ' cant be added due to missing input variable '   &
          & // trim(input_varname(iIn))
        EXIT
      end if

      ! allocate and set input_varIndex only if nVarIndex > 0
      ! input_varIndex is used only for one nInputs
      if ( nVarIndex > 0 .and. nInputs == 1 ) then
        ! Check whether the component index requested from input_varname
        ! is available
        in_nComp = me%method%val(inpos(iIn))%nComponents
        if ( (maxval(input_varIndex) > in_nComp) .or. &
          &  (nVarIndex > in_nComp) ) then
          write(logUnit(llerror),*) 'WARNING: Component index requested ' &
            & // 'is not extractable from variable: '  &
            & // trim(input_varname(iIn))
          write(logUnit(llerror),*) ' Input variable name: ' &
            & // trim(input_varname(iIn))
          write(logUnit(llerror),*) ' nComponents of input_varname: ', &
            & in_nComp
          write(logUnit(llerror),*) ' Requested component index of ' &
            & // 'input var: ', input_varIndex
          EXIT
        end if
        ! Copy input_varIndex
        input_varIndex_loc = input_varIndex
      end if

    end do

    if (nonexisting_input) then

      varPos = 0
      newVar = .false.

    else

      ! All inputs are found in the variable system, thus this derived
      ! quantity can be added now.
      call append( me       = me%varname, &
        &          val      = varname,    &
        &          pos      = varPos,     &
        &          wasAdded = newVar      )

      if (newVar) then
        write(logUnit(lldebug),*) 'Variable ' // trim(varname) &
          & // ' was newly added'
        var_access%mypos = varPos
        if (present(operType)) then
          var_access%operType = trim(operType)
        else
          var_access%operType = 'derived'
        end if
        var_access%nInputs = nInputs
        var_access%nComponents = nComponents
        allocate(var_access%state_varPos(0))
        allocate(var_access%auxField_varPos(0))
        allocate(var_access%input_varPos(nInputs))
        var_access%input_varPos = inpos
        allocate(var_access%input_varIndex(nVarIndex))
        var_access%input_varIndex = input_varIndex_loc

        var_access%method_data = method_data
        var_access%get_point => get_point
        var_access%get_element => get_element
        if (present(set_params)) then
          var_access%set_params => set_params
        else
          var_access%set_params => tem_varSys_setParams_dummy
        end if
        if (present(get_params)) then
          var_access%get_params => get_params
        else
          var_access%get_params => tem_varSys_getParams_dummy
        end if
        var_access%setup_indices => setup_indices
        var_access%get_valOfIndex => get_valOfIndex
        call append( me  = me%method,   &
          &          val = var_access )
      else
        write(logUnit(lldebug),*) 'Variable ' // trim(varname) &
          & // ' is already present'
      end if

    end if

    if (present(pos)) pos = varPos
    if (present(wasAdded)) wasAdded = newVar

  end subroutine tem_varSys_append_derVar
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Load variable system(s) from a lua file
  !!
  subroutine tem_varSys_load_vector( me, conf, parent, key )
    ! ---------------------------------------------------------------------------
    !> The variable system to read in
    type(tem_varSys_type), intent(out), allocatable :: me(:)
    !> Lua handle connected to the script to read the table from
    type(flu_State) :: conf
    !> A parent table handle in which to look the current variables up
    integer, intent(in), optional :: parent
    !> load from this key
    character(len=*), optional, intent(in) :: key
    ! ---------------------------------------------------------------------------
    integer :: thandle, subhandle
    integer :: nSys, iSys
    character(len=LabelLen) :: local_key
    ! ---------------------------------------------------------------------------

    if( present( key )) then
      local_key = key
    else
      local_key = 'varsys'
    endif

    call aot_table_open( L       = conf,           &
      &                  parent  = parent,         &
      &                  thandle = thandle,        &
      &                  key     = trim(local_key) )
    !varsys table is present, load variable system
    if( thandle > 0 ) then
      ! Now the 'varsys' table is opened. Let's first check, if
      ! there is a single variable system directly inside, or if there are
      ! several, so we have to open a table first and check whether it has
      ! another table
      call aot_table_open( L       = conf,      &
        &                  parent  = thandle,   &
        &                  thandle = subhandle, &
        &                  pos     = 1          )
      ! more than one varSys
      if ( subhandle > 0 ) then
        call aot_table_close( L = conf, thandle = subhandle )
        nSys = aot_table_length( L = conf, thandle = thandle )
        allocate(me(nSys))
        ! load each varSys
        do iSys = 1, nSys
          call aot_table_open( L       = conf,      &
            &                  parent  = thandle,   &
            &                  thandle = subhandle, &
            &                  pos     = iSys       )
          call tem_varSys_load_single( me        = me(iSys),  &
            &                          conf      = conf,      &
            &                          parent    = subhandle, &
            &                          openTable = .false.    )
          call aot_table_close( L = conf, thandle = subhandle )
        end do
      else !single varSys
        nSys = 1
        allocate(me(nSys))
        call tem_varSys_load_single( me        = me(1),   &
          &                          conf      = conf,    &
          &                          parent    = thandle, &
          &                          openTable = .false.  )
      end if
    else
      write(logUnit(1),*) 'Error: Failed to load varsys table'
      call tem_abort()
    end if

    call aot_table_close( L = conf, thandle = thandle )

  end subroutine tem_varSys_load_vector
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> load varSys from lua file.
  !! Required for harvester to load varSys from tracking or restart header
  !! file.
  subroutine tem_varSys_load_single(me, conf, parent, key, openTable)
    ! --------------------------------------------------------------------------!
    !> varSys to read from the Lua script(conf) and fill
    type(tem_varSys_type), intent(out) :: me

    !> Lua handle connected to the script to read the table from
    type(flu_state) :: conf

    !> A parent table handle in which to look the current variable up
    integer, optional, intent(in) :: parent

    !> load varsys from this key
    character(len=*), optional, intent(in) :: key

    ! if varsys table is already opened then openTable = .false.
    logical, optional, intent(in) :: openTable
    ! --------------------------------------------------------------------------!
    integer :: syshandle, varhandle, varsubhandle, iVar, nVars, pos
    type(tem_varSys_op_type) :: method
    character(len=labelLen) :: varname
    integer :: iError
    integer, allocatable :: vError(:)
    character(len=LabelLen) :: local_key
    logical :: openTable_loc
    ! --------------------------------------------------------------------------!
    ! if varsys table is already opened then openTable = .false.
    if( present(openTable) ) then
      openTable_loc = openTable
    else
      openTable_loc = .true.
    end if

    if( present( key )) then
      local_key = key
    else
      local_key = 'varsys'
    endif

    if (openTable_loc) then
      write(logUnit(1),*) 'Loading varsys table '
      ! Try to open the varsys table
      call aot_table_open( L       = conf,           &
        &                  parent  = parent,         &
        &                  thandle = syshandle,      &
        &                  key     = trim(local_key) )
      write(logUnit(1),*) 'Syshandle ', syshandle
    else
      ! if table is open before just use the parent as syshandle to load
      ! variable info
      syshandle = parent
    end if

    !varsys table is present, load variable system
    if (syshandle > 0) then
      ! Get the variable system name
      call aot_get_val( L       = conf,          &
        &               thandle = syshandle,     &
        &               val     = me%systemname, &
        &               ErrCode = iError,        &
        &               key     = 'systemname'   )

      ! get nStateVars
      call aot_get_val( L       = conf,          &
        &               thandle = syshandle,     &
        &               val     = me%nStateVars, &
        &               ErrCode = iError,        &
        &               key     = 'nStateVars'   )

      call aot_get_val( L       = conf,        &
        &               thandle = syshandle,   &
        &               val     = me%nScalars, &
        &               ErrCode = iError,      &
        &               key     = 'nScalars'   )

      call aot_table_open( L       = conf,      &
        &                  parent  = syshandle, &
        &                  thandle = varhandle, &
        &                  key     = 'variable' )

      if (varhandle > 0) then
        nVars = aot_table_length( L = conf, thandle = varhandle )
        call init( me     = me%varname, &
          &        length = nVars       )

        call init( me     = me%method, &
          &        length = nVars      )

        do iVar = 1, nVars
          method%myPos = iVar
          call aot_table_open( L       = conf,         &
            &                  parent  = varhandle,    &
            &                  thandle = varsubhandle, &
            &                  pos     = iVar          )

          call aot_get_val( L       = conf,         &
            &               thandle = varsubhandle, &
            &               val     = varname,      &
            &               ErrCode = iError,       &
            &               key     = 'name'        )

          call append( me       = me%varname, &
            &          val      = varname,    &
            &          pos      = pos         )

          call aot_get_val( L       = conf,               &
            &               thandle = varsubhandle,       &
            &               val     = method%nComponents, &
            &               ErrCode = iError,             &
            &               key     = 'ncomponents'       )

          call aot_get_val( L         = conf,                &
            &               thandle   = varsubhandle,        &
            &               val       = method%state_varPos, &
            &               maxLength = method%nComponents,  &
            &               ErrCode   = vError,              &
            &               key       = 'state_varpos'       )

          call append( me = me%method, val = method )
          call aot_table_close( L = conf, thandle = varsubhandle)
        end do
        call aot_table_close( L = conf, thandle = varhandle)
      else
        write(logUnit(1),*) 'Error: Failed to load variable table within ' &
          &                 //'varsys table'
        call tem_abort()
      end if ! variable table
    else
      write(logUnit(1),*) 'Error: Failed to load varsys table single'
      call tem_abort()
    end if !varsys table
    call aot_table_close( L = conf, thandle = syshandle)

  end subroutine tem_varSys_load_single
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Dumps array of varSys to given unit
  subroutine tem_varSys_dump_vector(me, outUnit, dumpVarPos)
    ! ---------------------------------------------------------------------------
    !> variable to write into the lua file
    type(tem_varSys_type), intent(in) :: me(:)

    !> unit to write to
    integer, intent(inout) :: outUnit

    !> Position of variables to dump
    integer, optional, intent(in) :: dumpVarPos(:)
    ! ---------------------------------------------------------------------------
    ! aotus type handling the output to the file in lua format
    type(aot_out_type) :: conf
    ! ---------------------------------------------------------------------------
    call aot_out_open( put_conf = conf, outUnit = outUnit )
    call tem_varSys_out_vector( me, conf, dumpVarPos )
    call aot_out_close( put_conf = conf )

  end subroutine tem_varSys_dump_vector
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Dump single varSys to given unit
  subroutine tem_varSys_dump_single(me, outUnit, dumpVarPos)
    ! ---------------------------------------------------------------------------
    !> variable to write into the lua file
    type(tem_varSys_type), intent(in) :: me

    !> unit to write to
    integer, intent(inout) :: outUnit

    !> Position of variables to dump
    integer, optional, intent(in) :: dumpVarPos(:)
    ! ---------------------------------------------------------------------------
    ! aotus type handling the output to the file in lua format
    type(aot_out_type) :: conf
    ! ---------------------------------------------------------------------------
    call aot_out_open( put_conf = conf, outUnit = outUnit )
    call tem_varSys_out_single( me, conf, dumpVarPos )
    call aot_out_close( put_conf = conf )

  end subroutine tem_varSys_dump_single
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Allows the output of array of varSys to lua out
  subroutine tem_varSys_out_vector(me, conf, dumpVarPos)
    ! ---------------------------------------------------------------------------
    !> variable to write into the lua file
    type(tem_varSys_type), intent(in) :: me(:)

    !> aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf

    !> Position of variables to dump
    integer, optional, intent(in) :: dumpVarPos(:)
    ! ---------------------------------------------------------------------------
    integer :: iSys
    ! ---------------------------------------------------------------------------
    call aot_out_open_table( put_conf = conf, tname='varsys' )
    do iSys = 1, size(me)
      call tem_varSys_out_single( me         = me(iSys),   &
        &                         conf       = conf,       &
        &                         dumpVarPos = dumpVarPos, &
        &                         level      = 1           )
    end do
    call aot_out_close_table( put_conf = conf )

  end subroutine tem_varSys_out_vector
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Write the system of variables description into a Lua file.
  !! use the aotus_out functions for doing so, in order to obtain a neatly
  !! formatted lua file
  !!
  subroutine tem_varSys_out_single(me, conf, dumpVarPos, level)
    ! ---------------------------------------------------------------------------
    !> Variable system to write out
    type(tem_varSys_type), intent(in) :: me

    !> Aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf

    !> Position of variables to dump
    integer, optional, intent(in) :: dumpVarPos(:)

    !> to dump varSys with key or without key
    integer, optional, intent(in) :: level
    ! ---------------------------------------------------------------------------
    integer :: iVar
    ! number of variables to be dumped (nVals or nVars)
    integer :: nVarDump, pos, level_loc, nScalars, nStateVars
    integer, allocatable :: dumpVarPos_loc(:)
    ! --------------------------------------------------------------------------
    if (present(level)) then
      level_loc = level
    else
      level_loc = 0
    end if

    ! dumpVarPos is present, then dump only variable at this position in varSys
    if (present(dumpVarPos)) then
      nVarDump = size(dumpVarPos)
      allocate(dumpVarPos_loc(nVarDump))
      dumpVarPos_loc = dumpVarPos
      nStateVars = nVarDump
      nScalars = sum(me%method%val(dumpVarPos_loc(:))%nComponents)
    else
      nVarDump = me%varname%nVals
      allocate(dumpVarPos_loc(nVarDump))
      dumpVarPos_loc = me%method%val(1:nVarDump)%myPos
      nStateVars = me%nStateVars
      nScalars = me%nScalars
    end if


    if( level_loc == 0) then
      call aot_out_open_table( put_conf = conf, tname = 'varsys' )
    else
      call aot_out_open_table( put_conf = conf )
    end if

    call aot_out_val( put_conf = conf, val = trim(me%SystemName),              &
      &               vname = 'systemname')

    call aot_out_open_table( put_conf = conf, tname = 'variable' )
    do iVar = 1, nVarDump
      call aot_out_open_table( put_conf = conf )

      pos = dumpVarPos_loc(iVar)

      call aot_out_val( put_conf = conf,                                       &
        &               val      = trim(me%varname%val(pos)),                  &
        &               vname    = 'name' )
      call aot_out_val( put_conf = conf,                                       &
        &               val      = me%method%val(pos)%nComponents,             &
        &               vname    = 'ncomponents' )

      !call aot_out_val( put_conf = conf,                                       &
      !  &               val      = me%method%val(pos)%myPos,                   &
      !  &               vname    = 'mypos' )

      if (allocated(me%method%val(pos)%state_varPos)) then
        if (size(me%method%val(pos)%state_varPos)>0) then
          call aot_out_val( put_conf = conf,                                   &
            &               val      = me%method%val(pos)%state_varPos,        &
            &               vname    = 'state_varpos' )
        end if
      end if

      !if (me%method%val(pos)%nInputs > 0) then
      !  call aot_out_val( put_conf = conf,                                    &
      !    &               val      = me%method%val(pos)%input_varPos,         &
      !    &               vname    = 'input_varPos' )
      !end if

      call aot_out_close_table( put_conf = conf )
    end do

    call aot_out_close_table( put_conf = conf )
    call aot_out_val( put_conf = conf, val = nScalars, vname = 'nScalars' )
    call aot_out_val( put_conf = conf, val = nStateVars, vname = 'nStateVars' )
    call aot_out_val( put_conf = conf, val = me%nAuxScalars, &
      & vname = 'nAuxScalars' )
    call aot_out_val( put_conf = conf, val = me%nAuxVars, vname = 'nAuxVars' )
    call aot_out_close_table( put_conf = conf )

  end subroutine tem_varSys_out_single
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> A routine to evaluate chunk of elements for given list of variables
  !!
  !! If subTree is present, it will use map2Global from subTree else
  !! map2Global is created for current chunk of global mesh
  subroutine tem_get_element_chunk(varSys, varPos, elemPos, time, tree, nElems,&
    &                              nDofs, res)
    ! --------------------------------------------------------------------------!
    !> Variable system describing available data.
    type(tem_varsys_type), intent(in) :: varsys

    !> Position of variables to evaluate in varSys
    integer, intent(in) :: varPos(:)

    !> Mesh definition of the input data.
    type(treelmesh_type), intent(in) :: tree

    !> Time information for the current data.
    type(tem_time_type), intent(in) :: time

    !> Number of degrees of freedom.
    integer, intent(in) :: nDofs

    !> Position of treeID of the element to get the variable for.
    integer, intent(in) :: elemPos(:)

    !> number of elements to evaluate
    integer, intent(in) :: nElems

    !> Output data
    !! size: io_buffer_size
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------!
    integer :: maxComponents, nScalars, nComponents, elemsize, compOff
    integer :: iElem, iVar, iDof, nVars, res_size
    integer :: e_start, t_start, d_start
    real(kind=rk), allocatable :: tmpdat(:)
    ! --------------------------------------------------------------------------!
    ! number of variables to evaluate
    nVars = size(varPos)

    ! Number of scalars in current output
    nScalars = sum(varSys%method%val(varPos(:))%nComponents)

    ! Size of a single element
    elemsize = nScalars*nDofs

    ! Need to obtain the data variable for variable, and store it in an
    ! intermediate array, because all components should be put together in the
    ! res array.
    ! The temporary array therefore needs to be sufficiently large to store the
    ! maximal number of components.
    maxComponents = maxval(varSys%method%val(varPos(:))%nComponents)

    ! Using a temporary array to store the variables and transfer them to res
    ! in the correct ordering afterwards.
    allocate(tmpdat(nElems*maxComponents*nDofs))

    compOff = 0
    vars: do iVar=1, nVars

      ! get the number of components for variable iVar
      nComponents = varSys%method%val(varPos(iVar))%nComponents

      ! get the size of the needed part of the res array
      res_size = nElems * nDofs * nComponents

      ! derive the quantities for all the elements in the current chunk
      call varSys%method%val(varpos(iVar))%get_element(                     &
        &                                varSys  = varSys,                  &
        &                                elemPos = elemPos,                 &
        &                                time    = time,                    &
        &                                tree    = tree,                    &
        &                                nElems  = nElems,                  &
        &                                nDofs   = nDofs,                   &
        &                                res     = tmpdat(:res_size)        )

      ! copy the information to the right positions in the result array
      ! res contains results for all variables,
      ! tmpdat is only for one variable
      do iElem=1,nElems
        e_start = (iElem-1)*elemsize
        t_start = (iElem-1)*nComponents*nDofs
        do iDof=1,nDofs
          d_start = (iDof-1)*nScalars + compOff
          res( (e_start+d_start+1) : (e_start+d_start+nComponents) ) &
            &  = tmpdat( t_start + (iDof-1)*nComponents + 1 &
            &            :t_start + iDof*nComponents        )
        end do
      end do
      ! Increase the component offset for the next variables.
      compOff = compOff + nComponents
    end do vars

  end subroutine tem_get_element_chunk
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> A routine to evaluate chunk of points for given list of variables
  !!
  !! If subTree is present, it will use map2Global from subTree else
  !! map2Global is created for current chunk of global mesh
  subroutine tem_get_point_chunk(varSys, varPos, point, time, tree, nPnts, &
    &                            res)
    ! --------------------------------------------------------------------------!
    !> Variable system describing available data.
    type(tem_varsys_type), intent(in) :: varsys

    !> Position of variables to evaluate in varSys
    integer, intent(in) :: varPos(:)

    !> Mesh definition of the input data.
    type(treelmesh_type), intent(in) :: tree

    !> Time information for the current data.
    type(tem_time_type), intent(in) :: time

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Output data
    !! size: io_buffer_size
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------!
    integer :: maxComponents, nScalars, nComponents, compOff
    integer :: iPnt, iVar, nVars, res_size
    integer :: e_start, t_start
    real(kind=rk), allocatable :: tmpdat(:)
    ! --------------------------------------------------------------------------!
    ! number of variables to evaluate
    nVars = size(varPos)

    ! Number of scalars in current output
    nScalars = sum(varSys%method%val(varPos(:))%nComponents)

    ! Need to obtain the data variable for variable, and store it in an
    ! intermediate array, because all components should be put together in the
    ! res array.
    ! The temporary array therefore needs to be sufficiently large to store the
    ! maximal number of components.
    maxComponents = maxval(varSys%method%val(varPos(:))%nComponents)

    ! Using a temporary array to store the variables and transfer them to res
    ! in the correct ordering afterwards.
    allocate(tmpdat(nPnts*maxComponents))

    compOff = 0
    vars: do iVar=1, nVars

      ! get the number of components for variable iVar
      nComponents = varSys%method%val(varPos(iVar))%nComponents

      ! get the size of the needed part of the res array
      res_size = nPnts * nComponents

      ! derive the quantities for all the elements in the current chunk
      call varSys%method%val(varpos(iVar))%get_point(               &
        &                                varSys = varSys,           &
        &                                point  = point,            &
        &                                time   = time,             &
        &                                tree   = tree,             &
        &                                nPnts  = nPnts,            &
        &                                res    = tmpdat(:res_size) )

      ! copy the information to the right positions in the result array
      ! res contains results for all variables,
      ! tmpdat is only for one variable
      do iPnt=1,nPnts
        e_start = (iPnt-1)*nScalars + compOff
        t_start = (iPnt-1)*nComponents
        res( (e_start+1) : (e_start+nComponents) )         &
          &  = tmpdat( t_start + 1 : t_start + nComponents )
      end do
      ! Increase the component offset for the next variables.
      compOff = compOff + nComponents
    end do vars

  end subroutine tem_get_point_chunk
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> This function provides the assignment operator of tem_varSys_op_type
  subroutine copy_varOp(left, right)
    ! ---------------------------------------------------------------------------
    !> varSys op to copy to
    type(tem_varSys_op_type), intent(out) :: left
    !> varSys op to copy from
    type(tem_varSys_op_type), intent(in) :: right
    ! ---------------------------------------------------------------------------
    integer :: nVals
    ! ---------------------------------------------------------------------------
    left%myPos = right%myPos

    if (allocated(right%state_varPos)) then
      nVals = size(right%state_varPos)
      allocate(left%state_varPos(nVals))
      left%state_varPos = right%state_varPos
    end if
    if (allocated(right%auxField_varPos)) then
      nVals = size(right%auxField_varPos)
      allocate(left%auxField_varPos(nVals))
      left%auxField_varPos = right%auxField_varPos
    end if
    left%nComponents = right%nComponents
    left%nInputs = right%nInputs
    left%operType = right%operType

    if ( allocated(right%input_varPos) ) then
      allocate(left%input_varPos(size(right%input_varPos)))
      left%input_varPos = right%input_varPos
    end if

    if ( allocated(right%input_varIndex) ) then
      allocate(left%input_varIndex(size(right%input_varIndex)))
      left%input_varIndex = right%input_varIndex
    end if

    left%method_data = right%method_data
    left%get_element => right%get_element
    left%get_point => right%get_point
    left%set_params => right%set_params
    left%get_params => right%get_params
    left%setup_indices => right%setup_indices
    left%get_valOfIndex => right%get_valOfIndex

  end subroutine copy_varOp
  ! ************************************************************************ !
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Free a variable description.
  !!
  !! Deallocates arrays and nullifies pointers.
  !! Does not touch the method data. Thus you should free the method data
  !! before calling this routine.
  subroutine tem_free_varOp(fun)
    !> Variable method to free.
    type(tem_varSys_op_type), intent(inout) :: fun

    if (allocated(fun%state_varPos)) deallocate(fun%state_varPos)
    if (allocated(fun%auxField_varPos)) deallocate(fun%auxField_varPos)
    fun%nComponents = 0
    fun%nInputs = 0
    if (allocated(fun%input_varpos)) deallocate(fun%input_varpos)
    if (allocated(fun%input_varIndex)) deallocate(fun%input_varIndex)
    nullify(fun%get_point)
    nullify(fun%get_element)
    nullify(fun%set_params)
    nullify(fun%get_params)
    nullify(fun%setup_indices)
    nullify(fun%get_valOfIndex)

  end subroutine tem_free_varOp
  ! ************************************************************************ !
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Empty a variable system.
  !!
  !! This frees all variables and empties the method and varname arrays.
  !! Note, this does not free the memory used for method data.
  !! Any method data that is to be freed has to be dealt with beforehand.
  subroutine tem_empty_varsys(varsys)
    !> Variable system to empty.
    type(tem_varSys_type), intent(inout) :: varsys

    integer :: iMethod

    do iMethod=1,varsys%method%nVals
      call tem_free_varOp(varsys%method%val(iMethod))
    end do
    call empty(varsys%method)
    call empty(varsys%varname)
    varsys%nStateVars = 0
    varsys%nAuxVars = 0
    varsys%nScalars = 0
    varsys%nAuxScalars = 0
    varsys%systemName = ''

  end subroutine tem_empty_varsys
  ! ************************************************************************ !
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_varSys_getParams_dummy(fun, varSys, instring, outstring)
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Input string with parameter to set in method_data
    character(len=*), intent(in) :: instring

    !> Output string with requested parameter value from method_data
    character(len=*), intent(out) :: outstring

    ! Nothing to do here ...
    outstring = trim(instring)
    call tem_abort('Not implemented for '//trim(varSys%varName%val(fun%myPos)) &
      &         // ':is dummy getParams')
  end subroutine tem_varSys_getParams_dummy
  ! ************************************************************************ !

  ! ************************************************************************ !
  subroutine tem_varSys_setParams_dummy(fun, varSys, instring)
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Input string with parameter to set in method_data
    character(len=*), intent(in) :: instring

    ! Nothing to do here ...
    write(*,*) trim(instring)
    call tem_abort('Not implemented for '//trim(varSys%varName%val(fun%myPos)) &
      &         // ':is dummy setParams')
  end subroutine tem_varSys_setParams_dummy
  ! ************************************************************************ !

  ! ************************************************************************ !
  subroutine tem_varSys_setupIndices_dummy(fun, varSys, point, offset_bit,    &
    &                                      iLevel, tree, nPnts, idx)
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> List of space coordinate points to store as growing array in
    !! method_data
    real(kind=rk), intent(in) :: point(:,:)

    !> Offset bit encoded as character for every point.
    !!
    !! Offset integer coord(3) is converted into a character with
    !! offset_bit = achar( (coord(1)+1) + (coord(2)+1)*4 + (coord(3)+1)*16 )
    !! Backward transformation form character to 3 integer:
    !! coord(1) = mod(ichar(offset_bit),4) - 1
    !! coord(2) = mod(ichar(offset_bit),16)/4 - 1
    !! coord(3) = ichar(offset_bit)/16 - 1
    !!
    !! If not present default is to center i.e offset_bit = achar(1+4+16)
    character, optional, intent(in) :: offset_bit(:)

    !> Level to which input points belong to
    integer, intent(in) :: iLevel

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of points to add in method_data of this variable
    integer, intent(in) :: nPnts

    !> Index of points in the growing array and variable val array.
    !! Size: n
    !!
    !! This must be stored in boundary or source depends on who
    !! calls this routine.
    !! This index is required to return a value using getValOfIndex.
    integer, intent(out) :: idx(:)

    write(*,*) 'point(1,1):', point(1,1)
    if (present(offset_bit)) write(*,*) 'offset_bit(1):', offset_bit(1)
    write(*,*) 'iLevel:', iLevel
    write(*,*) 'tree%nElems:', tree%nElems
    write(*,*) 'nPnts:', nPnts
    idx = 0
    call tem_abort('Not implemented for '//trim(varSys%varName%val(fun%myPos)) &
      &           //':is dummy setupIndices')
  end subroutine tem_varSys_setupIndices_dummy
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_varSys_getValOfIndex_dummy(fun, varSys, time, iLevel, &
    &                                       idx, idxLen, nVals, res)

    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: most times nVals, if contiguous arrays are used it depends
    !! on the number of first indices
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: dependes on number of first index for contiguous array,
    !! but the sum of all idxLen is equal to nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)

    write(*,*) 'time%sim:', time%sim
    write(*,*) 'iLevel:', iLevel
    write(*,*) 'idx(1):', idx(1)
    if (present(idxLen)) write(*,*) 'offset_bit(1):', idxLen(1)
    write(*,*) 'nVals:', nVals
    res = 0.0_rk
    call tem_abort('Not implemented for '//trim(varSys%varName%val(fun%myPos)) &
      &           //':is dummy getValOfIndex')
  end subroutine tem_varSys_getValOfIndex_dummy
  ! ************************************************************************ !

  ! ************************************************************************ !
  subroutine tem_varSys_getPoint_dummy(fun, varsys, point, time, tree, nPnts, &
    &                                  res)
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)

    write(*,*) 'point(1,1): ', point(1,1)
    write(*,*) 'time%sim:', time%sim
    write(*,*) 'tree%nElems: ', tree%nElems
    write(*,*) 'nPnts:', nPnts
    res = 0.0_rk
    call tem_abort('Not implemented for '//trim(varSys%varName%val(fun%myPos)) &
      &           //':is dummy getPoint')
  end subroutine tem_varSys_getPoint_dummy
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_varSys_getElement_dummy(fun, varsys, elempos, time, tree, &
    &                                    nElems, nDofs, res)
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Position of element in tree%treeID to get the variable for.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of elements to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (nComponents of resulting variable) x (nDegrees of freedom) x (nElems)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)

    write(*,*) 'elempos(1): ', elempos(1)
    write(*,*) 'time%sim:', time%sim
    write(*,*) 'tree%nElems: ', tree%nElems
    write(*,*) 'nElems: ', nElems
    write(*,*) 'nDofs: ', nDofs
    res = 0.0_rk
    call tem_abort('Not implemented for '//trim(varSys%varName%val(fun%myPos)) &
      &           //':is dummy get_element')
  end subroutine tem_varSys_getElement_dummy
  ! ************************************************************************ !


  ! ************************************************************************** !
  subroutine tem_varSys_check_inArgs( fun, varSys, time, iLevel, idx, idxLen, &
    &                                 nVals, label                            )
    !------------------------------------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    integer, intent(in) :: nVals
    character(len=*), intent(in) :: label
    !------------------------------------------------------------------------!

    if (size(idx) /= nVals) then
      write(logunit(10),*) 'in ', trim(label), ' idx length /= nVals !'
      write(logunit(10),*) '  varname:', trim(varsys%varname%val(fun%mypos))
      write(logunit(10),*) '  time%sim:', time%sim
      write(logunit(10),*) '  iLevel:', iLevel
      write(logunit(10),*) '  nVals:', nVals
      write(logunit(10),*) '  size(idx):', size(idx)
      if (present(idxLen)) write(logunit(10),*) 'idxLen provided:', idxLen(1)
    end if
  end subroutine tem_varsys_check_inArgs
  ! ************************************************************************** !

end module tem_varSys_module
