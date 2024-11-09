! Copyright (c) 2014-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014, 2016, 2018-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016-2017 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016, 2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2018 Jana Gericke <jana.gericke@student.uni-siegen.de>
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
!> author: Kannan Masilamani
!! Contains varMap type definition to store requested variable position
!! in the global varSys and routine to create varMap from given variable
!! labels
!!
!! This is especially useful wherever a set of variables needs to be
!! provided, for example in tracking objects.
!! There is a routine to read a list of labels, and one routine to
!! look these labels up in a variable system, and store the position
!! of each variable in that system.
!!
module tem_varMap_module
  use env_module,                     only: rk, labelLen

  use aotus_module,                   only: flu_state, aoterr_WrongType
  use aot_table_module,               only: aot_table_open, aot_table_close, &
    &                                       aot_get_val, aot_exists

  use tem_grow_array_module,          only: grw_labelArray_type, &
    &                                       grw_intArray_type,   &
    &                                       init, append, truncate
  use tem_dyn_array_module,           only: dyn_labelArray_type, init, append, &
    &                                       truncate, PositionOfVal
  use tem_logging_module,             only: logUnit
  use tem_varSys_module,              only: tem_varSys_type, &
    &                                       tem_varSys_solverData_evalElem_type
  use tem_stringKeyValuePair_module,  only: tem_stringKeyValuePair_type,      &
    &                                       grw_stringKeyValuePairArray_type, &
    &                                       append
  use tem_spacetime_fun_module,       only: tem_spacetime_fun_type, &
    &                                       tem_load_spacetime,     &
    &                                       tem_spacetime_hash_id
  use tem_spacetime_var_module,       only: tem_varSys_append_stfun

  implicit none

  ! ****************************************************************************
  !> Contains variable labels and their positions in a varSys
  type tem_varMap_type
    !> Array of active variable labels
    type(grw_labelArray_type) :: varName
    !character(len=labelLen), allocatable :: varName(:)

    !> Position of appended variable in the global variable system
    type(grw_intArray_type) :: varPos

    !> number of scalars of variables in varPos
    integer :: nScalars
  end type tem_varMap_type
  ! ****************************************************************************


  ! ****************************************************************************
  !> Contains possible variable labels and their nComponents
  type tem_possible_variable_type
    !> Array of possible variable labels
    type(dyn_labelArray_type) :: varName

    !> Number of components of variable in dynamic array of varName
    type(grw_intArray_type) :: nComponents
  end type tem_possible_variable_type
  ! ****************************************************************************

  ! ****************************************************************************
  interface init
    module procedure init_possible_variable
  end interface init
  ! ****************************************************************************

  ! ****************************************************************************
  interface append
    module procedure append_possible_variable
  end interface append
  ! ****************************************************************************

  ! ****************************************************************************
  interface truncate
    module procedure truncate_possible_variable
  end interface truncate
  ! ****************************************************************************


  ! ****************************************************************************
  interface  tem_variable_loadMapping
    module procedure tem_variable_loadMapping_vector
    module procedure tem_variable_loadMapping_single
  end interface tem_variable_loadMapping
  ! ****************************************************************************


contains


  ! ****************************************************************************
  subroutine init_possible_variable( me, length )
    !---------------------------------------------------------------------------
    !> Possible variables
    type(tem_possible_variable_type), intent(out) :: me
    !> Minimum length to expand the array
    integer, intent(in) :: length
    !---------------------------------------------------------------------------
    ! Unique array is used to extract position of variable to get nComponents
    call init( me = me%varName, length = length )
    call init( me = me%nComponents, length = length )

  end subroutine init_possible_variable
  ! ****************************************************************************


  ! ****************************************************************************
  subroutine append_possible_variable( me, varName, nComponents, pos )
    !---------------------------------------------------------------------------
    !> Possible variables
    type(tem_possible_variable_type), intent(inout) :: me
    !> variable name
    character(len=*), intent(in) :: varName
    !> number of components
    integer, intent(in) :: nComponents
    integer, intent(out), optional :: pos
    !---------------------------------------------------------------------------
    integer :: internal_pos
    !---------------------------------------------------------------------------
    call append( me  = me%varName,    &
      &          val = trim(varName), &
      &          pos = internal_pos   )

    call append( me  = me%nComponents, &
      &          val = nComponents     )

    if( present( pos ) ) pos = internal_pos

  end subroutine append_possible_variable
  ! ****************************************************************************


  ! ****************************************************************************
  subroutine truncate_possible_variable( me )
    !---------------------------------------------------------------------------
    !> Possible variables
    type(tem_possible_variable_type), intent(inout) :: me
    !---------------------------------------------------------------------------
    call truncate( me  = me%varName )
    call truncate( me  = me%nComponents )

  end subroutine truncate_possible_variable
  ! ****************************************************************************



  ! ****************************************************************************
  !> Loads the variable mapping from a table defined by given key
  !! for the variable names defined in possVars list.
  !! A variable mapping is used to
  !! link a user defined variable to a variable expected from, e.g., an equation
  !! system. These mappings are stored in varDict, which basically is a
  !! dictionary, whereas the key contains the name of the expected variable and
  !! the value contains the name of the user defined variable in the variable
  !! table.
  subroutine tem_variable_loadMapping_vector( possVars, conf, parent, key,     &
    &                                         varDict, varSys                  )
    !---------------------------------------------------------------------------
    !> Possible variable names expected by the solver
    type(tem_possible_variable_type), intent(in) :: possVars
    !> lua config file
    type(flu_state), intent(in) :: conf
    !> optional parent handle
    integer, optional, intent(in) :: parent
    !> optional key for the table
    character(len=*), optional, intent(in) :: key
    !> The dictionary that contains the mapping between expected variables
    !! and the actual variables defined by the user.
    type(grw_stringKeyValuePairArray_type), intent(inout) :: varDict
    !> Variable system to append anonymous variables to.
    type(tem_varSys_type), intent(inout) :: varSys

    !---------------------------------------------------------------------------
    character(len=32) :: localKey
    integer :: iVar, srchandle
    integer :: ErrorCode(possVars%varname%nVals)
    !---------------------------------------------------------------------------

    if( present( key )) then
      localKey = key
    else
      localKey = 'source'
    endif
    write(logUnit(1),*) 'Loading the mappings from table: ' // trim(localkey)

    if (possVars%varname%nVals == 0) then
      write(logUnit(1),*) 'WARNING: No possible variables are defined!' &
        & // ' Definitions from ' // trim(localKey) // ' are not loaded.'
      return
    end if

    ! Try to open the mapping table
    call aot_table_open( L       = conf,      &
      &                  parent  = parent,    &
      &                  thandle = srchandle, &
      &                  key     = localkey   )

    ! if mapping table is defined
    if (srchandle > 0) then
      ! loop over all possible variables
      do iVar = 1, possVars%varname%nVals
        ! check if expected variable is defined in lua file
        ! and if variable exist then store in stringKeyValuePairArray
        call tem_variable_loadMapping_single(                  &
          &     expectedName = possVars%varName%val(iVar),     &
          &     conf         = conf,                           &
          &     thandle      = srchandle,                      &
          &     varDict      = varDict,                        &
          &     varSys       = varSys,                         &
          &     nComp        = possVars%nComponents%val(iVar), &
          &     ErrCode      = ErrorCode(iVar)                 )
      end do

      call aot_table_close( L = conf, thandle = srchandle)
    else
      write(logUnit(1),*) 'Table ' // trim(localKey) &
        & // ' could not be opened for reading variable mappings!'
    end if

  end subroutine tem_variable_loadMapping_vector
  ! ****************************************************************************


  ! ****************************************************************************
  !> Loads the variable mapping from a table for single expected name.
  !! A variable mapping is used to
  !! link a user defined variable to a variable expected from, e.g., an equation
  !! system. These mappings are stored in varDict, which basically is a
  !! dictionary, whereas the key contains the name of the expected variable and
  !! the value contains the name of the user defined variable in the variable
  !! table.
  subroutine tem_variable_loadMapping_single( expectedName, conf, thandle, &
    &                                         varDict, varSys, ncomp,      &
    &                                         solverData_evalElem, ErrCode )
    ! --------------------------------------------------------------------------
    !> Expected variable name from config file
    character(len=*), intent(in) :: expectedName
    !> lua config file
    type(flu_state), intent(in) :: conf
    !> table handle
    integer, intent(in) :: thandle
    !> The dictionary that contains the mapping between expected variables
    !! and the actual variables defined by the user.
    type(grw_stringKeyValuePairArray_type), intent(inout) :: varDict
    !> Variable system to append anonymous variables to.
    type(tem_varSys_type), intent(inout) :: varSys
    !> Number of components we expect for this variable.
    integer, optional, intent(in) :: nComp

    !> Routine to convert point information into an element representation.
    type(tem_varSys_solverData_evalElem_type), &
      &  optional, intent(in) :: solverData_evalElem

    !> Error code
    integer, optional, intent(out) :: ErrCode
    ! --------------------------------------------------------------------------
    character(len=labelLen) :: varname
    character(len=labelLen) :: anony_name
    logical :: varExist
    integer :: errorCode
    type(tem_stringKeyValuePair_type) :: kvp
    type(tem_spacetime_fun_type), pointer :: stfun(:)
    real(kind=rk) :: numtest
    logical :: check_stfun
    ! --------------------------------------------------------------------------

    varExist = aot_exists( L       = conf,        &
      &                    thandle = thandle,     &
      &                    key     = expectedName )

    check_stfun = varExist

    ! if variable exist load space time function for that variable
    if (varExist) then
      write(logUnit(3),* ) 'Expected variable '  &
        & // trim(expectedName) // ' is defined'
      ! First check, whether this variable definition is a number.
      ! (They also satisfy strings).
      ! We do not accept numbers as variable names, instead this
      ! will be read as constant stfun.
      call aot_get_val( L       = conf,         &
        &               thandle = thandle,      &
        &               key     = expectedName, &
        &               val     = numtest,      &
        &               ErrCode = errorCode     )
      if (btest(errorCode, aoterr_WrongType)) then
        ! Not a number, try to interpret it as a string.
        call aot_get_val( L       = conf,         &
          &               thandle = thandle,      &
          &               key     = expectedName, &
          &               val     = varname,      &
          &               ErrCode = errorCode     )
        if (errorCode == 0) then
          ! Found a string, use it to refer to a variable.
          write(logUnit(1),*) 'Corresponding variable for ' &
            & // trim(expectedName) // ' found: ' // trim(varname)
          check_stfun = .false.
          kvp%key = expectedName
          kvp%value = varname
          call append(me = varDict, val = kvp)
        end if
      end if
      ! If the variable is not defined yet, try to get it as an anonymous
      ! variable in a space-time function.
      if (check_stfun) then
        write(logUnit(3),*) 'Did not find a variable name for ' &
          &                 // trim(expectedName)
        write(logUnit(3),*) 'Trying to load it as an anonymous space' &
          &                 // ' time function'
        allocate(stfun(1))
        call tem_load_spacetime(me      = stfun(1),           &
          &                     conf    = conf,               &
          &                     parent  = thandle,            &
          &                     key     = trim(expectedName), &
          &                     nComp   = nComp,              &
          &                     errCode = errorCode           )
        ! use global mesh for anonymous variable
        stfun(1)%subTree%useGlobalMesh = .true.
        if (errorCode == 0) then
          kvp%key = expectedName
          anony_name = tem_spacetime_hash_id(me=stfun(1), conf=conf)
          kvp%value = 'anon_' // trim(anony_name)
          call tem_varSys_append_stFun(                   &
            & varSys              = varSys,               &
            & stfun               = stfun,                &
            & varname             = kvp%value,            &
            & nComp               = nComp,                &
            & evalType            = 'firstonly_asglobal', &
            & solverData_evalElem = solverData_evalElem   )
          call append(me = varDict, val = kvp)
          write(logUnit(1),*) 'Found a spacetime function for ' &
            &                 // trim(expectedName)
        else
          write(logUnit(1),* ) 'Expected variable ' // trim(expectedName) &
            & // ' is NOT defined'
          deallocate(stfun)
        end if
      end if
    else
      write(logUnit(1),* ) 'Expected variable ' // trim(expectedName) &
        &                  // ' is NOT defined'
      errorCode = -1
    end if

    if (present(ErrCode)) ErrCode = errorCode

  end subroutine tem_variable_loadMapping_single
  ! ***************************************************************************


  ! ****************************************************************************
  !> Creates a variable map. Therefore it looks for the variables stored in
  !! varname in the varSys. If found, the name and the corresponding position
  !! in the varSys is added to the resultung varMap.
  subroutine tem_create_varMap( varName, varSys, varMap )
    ! --------------------------------------------------------------------------
    !> array of variable labels
    character(len=*), intent(in) :: varName(:)

    !> variable system to look for varName
    type(tem_varSys_type), intent(in) :: varSys

    !> Contains position of varname in given varSys
    type(tem_varMap_type), intent(out) :: varMap
    ! --------------------------------------------------------------------------
    integer :: varPos, iVar
    ! --------------------------------------------------------------------------
    write(logUnit(10),*) 'Create varMap:'
    call init(varMap%varName)
    call init(varMap%varPos)
    do iVar = 1, size(varName)
      varPos = PositionOfVal( me  = varSys%varName,     &
        &                     val = trim(varName(iVar)) )
      if (varPos > 0) then
        write(logUnit(10),"(A,I0)") 'Position of ' // trim(varName(iVar)) &
          & // ': ', varPos
        call append( me = varMap%varPos, val = varPos )
        call append( me = varMap%varname, val = varName(iVar) )
      else
        write(logUnit(5),*) 'Variable ' // trim(varName(iVar)) // &
          &                 ' not found in varSys'
      end if
    end do
    call truncate(varMap%varName)
    call truncate(varMap%varPos)
    varMap%nScalars = sum( varSys%method%val(             &
      &                    varMap%varPos%val(:) )%nComponents )
    write(logUnit(10),"(A,I0)") 'nScalars in varMap: ', varMap%nScalars

  end subroutine tem_create_varMap
  ! ****************************************************************************

end module tem_varMap_module
