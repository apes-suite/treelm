! Copyright (c) 2014-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014-2016, 2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2016, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> author: Kannan Masilamani
!! A generic variable definition which can be loaded from a Lua file.
!!
!! This module provides the means to define new variables via the
!! configuration.
!!
!! Example: Variable to evaluate space time functions
!!
!!```lua
!! variable = {
!!   name = 'ref_density',
!!   ncomponents = 1,
!!   vartype = 'st_fun',
!!   st_fun = lua_fun_refDens
!! }
!!```
!!
!! Example: variable to perform certain operations
!!
!!```lua
!! variable = {
!!   name = 'density_diff',
!!   ncomponents = 3,
!!   vartype = 'operation',
!!   operation = {
!!     kind='difference',
!!     input_varname = { 'density', 'ref_density'}
!!   }
!! }
!!```
!!
!! Example: variable to perform operation extract
!!
!!```lua
!! variable = {
!!   name = 'vel_y',
!!   ncomponents = 1,
!!   vartype = 'operation',
!!   operation = {
!!     kind='extract',
!!     input_varname = {'velocity'},
!!     input_varindex = {2}
!!   }
!! }
!!```

module tem_variable_module
  use iso_c_binding, only: c_ptr, c_null_ptr
  use env_module,                    only: labelLen

  use aotus_module,                  only: flu_State, aot_get_val,           &
    &                                      aoterr_NonExistent, aoterr_Fatal, &
    &                                      aoterr_WrongType
  use aot_table_module,              only: aot_table_open,   &
    &                                      aot_table_close,  &
    &                                      aot_table_length, &
    &                                      aot_get_val
  use aot_out_module,                only: aot_out_type,       &
    &                                      aot_out_open,       &
    &                                      aot_out_close,      &
    &                                      aot_out_val,        &
    &                                      aot_out_open_table, &
    &                                      aot_out_close_table

  use tem_spacetime_fun_module,      only: tem_spacetime_fun_type, &
    &                                      tem_load_spacetime
  use tem_tools_module,              only: upper_to_lower,      &
    &                                      tem_horizontalSpacer
  use tem_logging_module,            only: logUnit
  use tem_aux_module,                only: tem_abort
  use tem_reduction_transient_module, only:                                &
    &                                 tem_reduction_transient_config_type, &
    &                                 tem_reduction_transient_load
  use tem_varSys_module,             only: tem_varSys_type, &
    &                                      tem_varSys_solverData_evalElem_type

  implicit none

  ! **************************************************************************** !
  !> Description of user defined variables.
  type tem_variable_type
    !> A name for this variable.
    character(len=LabelLen) :: label

    !> Number of components of this variable
    integer :: nComponents

    !> How to treat this variable: "operation" or "st_fun"
    character(len=labelLen) :: varType

    !> Type of Operation to perfom on this variable
    character(len=labelLen) :: operType

    !> additional info to load for reduction transient operation
    type(tem_reduction_transient_config_type) :: redTransConfig

    !> input variables names this variable depends on
    character(len=LabelLen), allocatable :: input_varName(:)

    !> Component index to extract from input_varName for operType = 'extract'
    !! NOTE: It is possible to extract of component index from only one
    !! variable so input_varName for operType="extract" must be single
    !! variable name
    integer, allocatable :: input_varIndex(:)

    !> space time functions
    type(tem_spacetime_fun_type), allocatable :: st_fun(:)

    !> The evaluation type to use. The evaluation type defines how the variable
    !! should evaluate a get_point or get_element-request when there are more
    !! than one space time functions that could fulfill this request.
    !!
    !! The standard evaluation type is 'add'. So every space time function is
    !! asked to fulfill the request and all results are added up to the final
    !! result returned back to the caller.
    !!
    !! Another evaluation type is 'last' which takes the sole value of the last
    !! space time function that was able to fulfill this request. For this
    !! evaluation type the ordering of the space time functions in the variable
    !! definition is important.
    character(len=labelLen) :: evalType

    !> A method to append the read solver specific variable to the varSys
    procedure(tem_append_solverVar_method), pointer :: append_solverVar &
      &                                                => NULL()
    type(c_ptr) :: solver_specifics = C_NULL_PTR
  end type tem_variable_type
  ! **************************************************************************** !


  interface
    subroutine tem_append_solverVar_method(var, varSys, pos, &
      &                                    solverData_evalElem)
      import :: tem_variable_type, &
        &       tem_varSys_type,   &
        &       tem_varSys_solverData_evalElem_type
      !> Data to describe the solver specific variable
      class(tem_variable_type), intent(in) :: var

      !> Variable system to append the variable to.
      type(tem_varSys_type), intent(inout) :: varSys

      !> Position of the appended variable in the varSys.
      integer, optional, intent(out) :: pos
      !> A setter routine that allows the caller to define a routine for the
      !! construction of an element representation.
      type(tem_varSys_solverData_evalElem_type), &
        &  optional, intent(in) :: solverData_evalElem
    end subroutine tem_append_solverVar_method

    subroutine tem_load_solverVar_method(L, parent, specifics, appender, iError)
      import :: flu_State, &
        &       c_ptr,     &
        &       tem_append_solverVar_method
      !> Lua script to load the variable data from.
      type(flu_State) :: L

      !> Parent table in the lua script to read the vairable from.
      integer, intent(in) :: parent

      !> Data to read for the solver specific variable
      type(c_ptr), intent(out) :: specifics

      !> Function pointer to use for appending the solver variable.
      procedure(tem_append_solverVar_method), pointer :: appender

      !> Indication whether the attempted reading was successful.
      integer, intent(out) :: iError

    end subroutine tem_load_solverVar_method
  end interface


  ! **************************************************************************** !
  interface tem_variable_load
    module procedure tem_variable_load_vector
    module procedure tem_variable_load_single
  end interface tem_variable_load
  ! **************************************************************************** !


  ! **************************************************************************** !
  interface tem_variable_dump
    module procedure tem_variable_dump_vector
    module procedure tem_variable_dump_single
  end interface tem_variable_dump
  ! **************************************************************************** !


  ! **************************************************************************** !
  interface tem_variable_out
    module procedure tem_variable_out_vector
    module procedure tem_variable_out_single
  end interface tem_variable_out
  ! **************************************************************************** !


  ! **************************************************************************** !
  interface assignment(=)
    module procedure copy_Var
  end interface
  ! **************************************************************************** !


contains


  ! **************************************************************************** !
  !> Load an array of variables from the configuration.
  subroutine tem_variable_load_vector( me, conf, parent, key, vError, nComp, &
    &                                  load_solvervar                        )
    ! --------------------------------------------------------------------------!
    !> The variable to read from the Lua script(conf) and fill
    type(tem_variable_type), allocatable, intent(out) :: me(:)

    !> Lua handle connected to the script to read the table from
    type(flu_state) :: conf

    !> A parent table handle in which to look the current variable up
    integer, optional, intent(in) :: parent

    !> key for array of variables
    character(len=*), optional, intent(in) :: key

    !> if Error .ne. 0 is variable is not loaded successfully.
    integer, allocatable, intent(out) :: vError(:)

    !> If the variable is expected to have a certain number of components,
    !! this can be provided with this argument.
    !!
    !! If the definition of the variable does not match this, we will fail
    !! loading the variable.
    integer, optional :: nComp

    !> A method to load solver specific variables.
    procedure(tem_load_solvervar_method), optional :: load_solvervar
    ! --------------------------------------------------------------------------!
    integer :: varhandle, nVars, varsubhandle, iVar, iError
    character(len=LabelLen) :: local_key
    ! --------------------------------------------------------------------------!
    call tem_horizontalSpacer(fUnit = logUnit(1))

    if( present( key )) then
      local_key = key
    else
      local_key = 'variable'
    endif

    ! Try to open the variable table
    call aot_table_open( L       = conf,            &
      &                  parent  = parent,          &
      &                  thandle = varhandle,       &
      &                  key     = trim(local_key ) )

    nVars = 0
    if (varhandle > 0) then
      ! Test whether the next thing is a table or not
      call aot_table_open( L       = conf,         &
        &                  parent  = varhandle,    &
        &                  thandle = varsubhandle, &
        &                  pos     = 1             )
      ! It is a table, so more than one variable is expected
      if (varsubhandle > 0) then
        call aot_table_close( L = conf, thandle = varsubhandle )
        nVars = aot_table_length( L = conf, thandle = varhandle )
        allocate(me(nVars))
        allocate(vError(nVars))

        do iVar = 1, nVars
          call aot_table_open( L       = conf,         &
            &                  parent  = varhandle,    &
            &                  thandle = varsubhandle, &
            &                  pos     = iVar          )

          call tem_variable_load_single( me             = me(iVar),      &
            &                            conf           = conf,          &
            &                            parent         = varsubhandle,  &
            &                            iError         = iError,        &
            &                            nComp          = nComp,         &
            &                            openTable      = .false.,       &
            &                            load_solvervar = load_solvervar )
          vError(iVar) = iError
          if (iError /= 0) then
            write(logUnit(1),*) 'Variable:'//trim(me(iVar)%label) &
              &        //' cannot be added to varSys'
          endif
          call aot_table_close( L = conf, thandle = varsubhandle )
        end do
      else ! it's not a table but a single variable
        nVars = 1
        allocate(me(nVars))
        allocate(vError(nVars))
        call tem_variable_load_single( me             = me(1),         &
          &                            conf           = conf,          &
          &                            parent         = varhandle,     &
          &                            iError         = iError,        &
          &                            nComp          = nComp,         &
          &                            openTable      = .false.,       &
          &                            load_solvervar = load_solvervar )
        vError(1) = iError
        if (iError /= 0) then
          write(logUnit(1),*) 'Variable:'//trim(me(iVar)%label) &
            &        //' cannot be added to varSys'
        endif
      end if
    else
      write(logUnit(1),*) 'Variable table not defined with key: ' &
        &                 //trim(local_key)
      allocate(me(nVars))
      allocate(vError(nVars))
    endif

    call aot_table_close( L = conf, thandle = varhandle )
    call tem_horizontalSpacer(fUnit = logUnit(1))

  end subroutine tem_variable_load_vector
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Reading a single variable from the Lua configuration.
  subroutine tem_variable_load_single( me, conf, parent, iError, key, nComp, &
    &                                  openTable, load_solvervar )
    ! --------------------------------------------------------------------------!
    !> The variable to read from the Lua script(conf) and fill
    type(tem_variable_type), intent(out) :: me

    !> Lua handle connected to the script to read the table from
    type(flu_state) :: conf

    !> A parent table handle in which to look the current variable up
    integer, intent(in) :: parent

    !> key for a single variable
    character(len=*), optional, intent(in) :: key

    !> if Error .ne. 0 is variable is not loaded successfully.
    integer, intent(out) :: iError

    !> If the variable is expected to have a certain number of components,
    !! this can be provided with this argument.
    !!
    !! If the definition of the variable does not match this, we will fail
    !! loading the variable.
    integer, optional :: nComp

    !> if variable table is already opened, set openTable = .false.
    logical, optional, intent(in) :: openTable

    !> A method to load solver specific variables.
    procedure(tem_load_solvervar_method), optional :: load_solvervar
    ! --------------------------------------------------------------------------!
    integer :: local_thandle, local_error
    logical :: openTable_loc
    character(len=LabelLen) :: local_key
    character(len=LabelLen) :: varname
    ! --------------------------------------------------------------------------!
    call tem_horizontalSpacer(fUnit = logUnit(1))
    ! if variable table is already opened then openTable = .false.
    if( present(openTable) ) then
      openTable_loc = openTable
    else
      openTable_loc = .true.
    end if

    local_thandle = parent
    if (openTable_loc) then
      if ( present(key) ) then

        local_key = trim(key)

        ! Attempt to read the variable as a single identifier.
        ! (Reference to another variable)
        call aot_get_val( L       = conf,    &
          &               thandle = parent,  &
          &               val     = varname, &
          &               ErrCode = iError,  &
          &               key     = key      )

        ! If this succeeds, take it and set the variable with the
        ! name of the key to be a combine operation (providing an alias).
        if (iError == 0) then
          me%label = trim(key)
          ! Negative number of components will inherit the components from
          ! the referred variable.
          me%nComponents = -1
          me%vartype = 'operation'
          me%opertype = 'combine'
          allocate(me%input_varName(1))
          me%input_varName(1) = varname
        end if

      else

        local_key = 'variable'

      end if

      call aot_table_open( L       = conf,          &
        &                  thandle = local_thandle, &
        &                  parent  = parent,        &
        &                  key     = local_key      )
    end if

    ! Get the name of the variable
    call aot_get_val( L       = conf,          &
      &               thandle = local_thandle, &
      &               val     = me%label,      &
      &               ErrCode = iError,        &
      &               key     = 'name',        &
      &               pos     = 1,             &
      &               default = key            )

    if (iError /= 0) then
      write(logUnit(1),*) 'Unable to load "name" with pos and label.'
      return
    end if

    write(logUnit(1),*) 'Loading variable ', trim(me%label)

    ! Get the number of components for this variable
    call aot_get_val( L       = conf,           &
      &               thandle = local_thandle,  &
      &               val     = me%nComponents, &
      &               ErrCode = iError,         &
      &               key     = 'ncomponents',  &
      &               pos     = 2,              &
      &               default = nComp           )

    if( iError /= 0 ) then
      write(logUnit(1),*) 'No ncomponents specified for variable ' &
        &                 //trim( me%label )
      return
    end if

    write(logUnit(5),*) '  nComponents ', me%nComponents

    ! Do not proceed, if the number of provided components
    ! does not match the number of expect components.
    if (present(nComp)) then
      if (nComp /= me%nComponents) RETURN
    end if

    ! load variable type
    call aot_get_val( L       = conf,           &
      &               thandle = local_thandle,  &
      &               val     = me%varType,     &
      &               ErrCode = iError,         &
      &               default = 'none',         &
      &               key     = 'vartype'       )

    select case(trim(me%varType))
    case('operation')
      call load_variable_operation( me     = me,            &
        &                           conf   = conf,          &
        &                           parent = local_thandle, &
        &                           iError = iError         )
      ! if operation table is not loaded successfully. this variable
      ! cannot be added to variable system
      if (iError /= 0) return

    case('st_fun')
      ! In case we have a space time function, we try to read the evaluation
      ! type, This is only needed for variables of type space time function as
      ! they can have several space time functions providing values and thus
      ! there are use cases where these space time functions have to be merged
      ! differently.
      ! The default is to add all values, if several space time functions
      ! provide values for a given request.
      call aot_get_val( L       = conf,          &
        &               thandle = local_thandle, &
        &               val     = me%evalType,   &
        &               ErrCode = local_Error,   &
        &               key     = 'evaltype',    &
        &               default = 'add'          )
      me%evalType = upper_to_lower(me%evalType)

      write(logUnit(1),*) 'loading the spacetime functions for variable '//    &
        &             trim( me%label )
      call tem_load_spacetime( me     = me%st_fun,       &
        &                      conf   = conf,            &
        &                      parent = local_thandle,   &
        &                      nComp  = me%nComponents,  &
        &                      key    = 'st_fun' )

      if (.not. allocated(me%st_fun)) then
        write(logUnit(1),*) 'Error: no stfun found for '//trim(me%label)
        write(logUnit(1),*) 'If you define a variable with a stfun, this has'
        write(logUnit(1),*) 'to be properly defined.'
        write(logUnit(1),*) 'Check your variable definition!'
        call tem_abort()
      else
        if (size(me%st_fun) == 0) then
          write(logUnit(1),*) 'Error: no stfun found for '//trim(me%label)
          write(logUnit(1),*) 'If you define a variable with a stfun, this has'
          write(logUnit(1),*) 'to be properly defined.'
          write(logUnit(1),*) 'Check your variable definition!'
          call tem_abort()
        end if
      end if

    case default
      iError = 1
      if (present(load_solverVar)) then
        call load_solverVar(L         = conf,                &
          &                 parent    = local_thandle,       &
          &                 specifics = me%solver_specifics, &
          &                 appender  = me%append_solvervar, &
          &                 iError    = iError               )
      end if
      if (iError /= 0) then
        call tem_abort( 'Error: varType '                                       &
          & // trim(me%varType)                                                 &
          & // ' not supported! Supported varType are "st_fun" and "operation"' )
      end if
    end select

    if (present(key)) call aot_table_close( L = conf, thandle = local_thandle)

  end subroutine tem_variable_load_single
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Loading of a variable, that is defined in terms of an operation on other
  !! variables.
  subroutine load_variable_operation( me, conf, parent, iError)
    ! --------------------------------------------------------------------------!
    !> The variable to read from the Lua script(conf) and fill
    type(tem_variable_type), intent(inout) :: me

    !> Lua handle connected to the script to read the table from
    type(flu_state) :: conf

    !> A parent table handle in which to look the current operation table
    integer, intent(in) :: parent

    !> iError .ne. 0 if operation is not loaded successfully.
    integer, intent(out) :: iError
    ! --------------------------------------------------------------------------!
    integer :: oper_handle, iIn, nInvar
    integer, allocatable :: vError(:)
    ! --------------------------------------------------------------------------!

    call aot_table_open( L       = conf,          &
      &                  parent  = parent,        &
      &                  thandle = oper_handle,   &
      &                  key     = 'operation'    )
    if (oper_handle == 0) then
      iError = -1
      write(logUnit(1),*) 'operation table not defined.'
      return
    end if
    write(logUnit(5),*) '  Loading operation table '

    ! get operation kind
    call aot_get_val( L       = conf,        &
      &               thandle = oper_handle, &
      &               val     = me%operType, &
      &               ErrCode = iError,      &
      &               key     = 'kind'       )

    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*) 'FATAL Error occured, while retrieving "kind" '//    &
        &                 'from operation table :'
      if ( btest( iError, aotErr_NonExistent ))                                &
        & write(logUnit(1),*)'Variable not existent!'
      if (btest(iError, aoterr_WrongType))                                     &
        & write(logUnit(1),*)'Variable has wrong type!'
    end if
    write(logUnit(5),*) '    kind = '//trim(me%operType)

    ! Get the dependent variables
    call aot_get_val( val       = me%input_varName, &
      &               ErrCode   = vError,           &
      &               maxlength = 100,              &
      &               L         = conf,             &
      &               thandle   = oper_handle,      &
      &               key       = 'input_varname'   )

    if ( any(btest(vError, aoterr_Fatal)) ) then

      write(logUnit(1),*) 'FATAL Error occured, while retrieving '
      write(logUnit(1),*) '"input_varName" for operation' // trim(me%operType)
      call tem_abort()

    end if

    ! number of variables in input_varname
    nInVar = size(me%input_varName)

    write(logUnit(5),*) '    input_varname:'
    do iIn = 1, nInVar
      write(logUnit(5), *) '      ', iIN, trim(me%input_varName(iIn))
    end do

    deallocate(vError)

    ! Load operation specific infos
    select case(trim(me%operType))
    case ('reduction_transient')
      call tem_reduction_transient_load(me     = me%redTransConfig, &
        &                               conf   = conf,              &
        &                               parent = oper_handle        )

    case ('extract')
      allocate(me%input_varIndex(me%nComponents))
      allocate(vError(me%nComponents))
      ! Extraction of component index from input_varname is possible
      ! only when there is only one variable name in input_varname
      if ( nInVar == 1 ) then
        call aot_get_val( L         = conf,               &
          &               thandle   = oper_handle,        &
          &               val       = me%input_varIndex, &
          &               ErrCode   = vError,             &
          &               key       = 'input_varindex'   )

        if ( any(btest(vError, aoterr_Fatal)) .or. &
          &  any(me%input_varIndex <= 0) ) then

          write(logUnit(1),*) 'FATAL Error occured, while retrieving '
          write(logUnit(1),*) '"input_varindex" for extract operation'
          write(logUnit(1),*) 'input_varindex: ',me%input_varIndex
          call tem_abort()

        end if
      else
        write(logUnit(1),*) 'Error loading component index for operation'//&
          &                 ' kind extract'
        write(logUnit(1),*) 'Size of input_varname: ', nInVar
        write(logUnit(1),*) 'Extract component index from input_varname'//&
          &                 'works only when size(input_varname) == 1'
        call tem_abort()
      end if

      write(logUnit(5),*) '     Component index to extract:', me%input_varIndex
    end select
    call aot_table_close( L = conf, thandle = oper_handle)

  end subroutine load_variable_operation
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Dumps array of variable to given unit
  subroutine tem_variable_dump_vector(me, outUnit)
    ! ---------------------------------------------------------------------------
    !> variable to write into the lua file
    type(tem_variable_type), intent(in) :: me(:)
    !> unit to write to
    integer, intent(in) :: outUnit
    ! ---------------------------------------------------------------------------
    ! aotus type handling the output to the file in lua format
    type(aot_out_type) :: conf
    ! ---------------------------------------------------------------------------
    call aot_out_open( put_conf = conf, outUnit = outUnit )
    call tem_variable_out_vector( me, conf )
    call aot_out_close( put_conf = conf )

  end subroutine tem_variable_dump_vector
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Dump single variable to given unit
  subroutine tem_variable_dump_single(me, outUnit)
    ! ---------------------------------------------------------------------------
    !> variable to write into the lua file
    type(tem_variable_type), intent(in) :: me
    !> unit to write to
    integer, intent(in) :: outUnit
    ! ---------------------------------------------------------------------------
    ! aotus type handling the output to the file in lua format
    type(aot_out_type) :: conf
    ! ---------------------------------------------------------------------------
    call aot_out_open( put_conf = conf, outUnit = outUnit )
    call tem_variable_out_single( me, conf )
    call aot_out_close( put_conf = conf )

  end subroutine tem_variable_dump_single
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Allows the output of array of variable to lua out
  subroutine tem_variable_out_vector(me, conf)
    ! ---------------------------------------------------------------------------
    !> variable to write into the lua file
    type(tem_variable_type), intent(in) :: me(:)
    !> aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    ! ---------------------------------------------------------------------------
    integer :: iVar
    ! ---------------------------------------------------------------------------
    call aot_out_open_table( put_conf = conf, tname='variable' )
    do iVar = 1,size(me)
      call tem_variable_out_single( me(iVar), conf, level=1 )
    end do
    call aot_out_close_table( put_conf = conf )

  end subroutine tem_variable_out_vector
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Allows the output of the single variable to lua out.
  !!
  !! The data is written into the file, the lunit is connected to.
  !! It is formatted as a Lua table.
  !!
  subroutine tem_variable_out_single(me, conf, level)
    ! ---------------------------------------------------------------------------
    !> variable to write into the lua file
    type(tem_variable_type), intent(in) :: me
    !> aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    !> to dump variable with key or without key
    integer, optional, intent(in) :: level
    ! ---------------------------------------------------------------------------
    integer :: level_loc
    ! ---------------------------------------------------------------------------
    if (present(level)) then
      level_loc = level
    else
      level_loc = 0
    end if

    if( level_loc == 0) then
      call aot_out_open_table( put_conf = conf, tname = 'variable' )
    else
      call aot_out_open_table( put_conf = conf )
    end if

    call aot_out_val( put_conf = conf,                                         &
      &               val      = trim(me%label),                               &
      &               vname    = 'name' )
    call aot_out_val( put_conf = conf,                                         &
      &               val      = me%nComponents,                               &
      &               vname    = 'ncomponents' )

    if (trim(me%varType) /= '') then
      call aot_out_val( put_conf = conf,                                       &
        &               val      = me%varType,                                 &
        &               vname    = 'vartype' )

      if( trim(me%varType) == 'operation' ) then
        ! write the operation kind and dependent variables
        call aot_out_open_table( put_conf = conf, tname='operation' )
        call aot_out_val( put_conf = conf,        &
          &               val      = me%operType, &
          &               vname    = 'kind'       )
        call aot_out_val( put_conf = conf,             &
          &               val      = me%input_varName, &
          &               vname    = 'input_varname'   )
        if (allocated(me%input_varIndex)) then
          call aot_out_val( put_conf = conf,              &
            &               val      = me%input_varIndex, &
            &               vname    = 'input_varindex'   )
        end if
        call aot_out_close_table( put_conf = conf )
      end if
    end if
    call aot_out_close_table( put_conf = conf )

  end subroutine tem_variable_out_single
  ! **************************************************************************** !




  ! **************************************************************************** !
  !> This routine provide funtionality to copy variable type
  subroutine copy_Var(left, right)
    ! ---------------------------------------------------------------------------
    !> variable to copy to
    type(tem_variable_type), intent(out) :: left
    !> variable to copy from
    type(tem_variable_type), intent(in) :: right
    ! ---------------------------------------------------------------------------
    left%label = right%label
    left%nComponents = right%nComponents
    left%varType = right%varType
    left%operType = right%operType
    if (allocated(right%input_varName)) then
      allocate(left%input_varName(size(right%input_varName)))
      left%input_varName = right%input_varName
    end if
    if (allocated(right%input_varIndex)) then
      allocate(left%input_varIndex(size(right%input_varIndex)))
      left%input_varIndex = right%input_varIndex
    end if
    if (allocated(right%st_fun)) then
      allocate(left%st_fun(size(right%st_fun)))
      left%st_fun = right%st_fun
    end if
  end subroutine copy_Var
  ! **************************************************************************** !


end module tem_variable_module
