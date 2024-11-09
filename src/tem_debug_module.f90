! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2013, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2014-2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2013-2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Daniel Petró <daniel.petro@student.uni-siegen.de>
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
!> author: Manuel Hasert
!! Provide a general debug environment.
!! See [[tem_load_debug]] for configuration options.
!!
module tem_debug_module

  ! include treelm modules
  use env_module,              only: labelLen, PathLen, my_status_string, &
    &                                long_k, solSpecLen
  use tem_tools_module,        only: tem_horizontalSpacer
  use tem_logging_module,      only: tem_logging_type, tem_logging_load,       &
    &                                tem_logging_init, tem_logging_isActive,   &
    &                                tem_last_lu, logUnit, tem_toStr

  ! include aotus modules
  use flu_binding,      only: flu_State
  use aot_table_module, only: aot_get_val, aot_table_open, aot_table_close

  implicit none

  private

  !> The length of the buffer to create a string representation of arrays.
  integer, parameter :: buffer_length = 128

  integer, public :: dbgUnit(0:tem_last_lu)

  public :: tem_debug_type
  public :: tem_load_debug
  public :: tem_debug_load_main
  public :: tem_print_array
  public :: tem_reportStatus
  public :: tem_debug
  public :: main_debug

  !> Debug option definitions
  type tem_debug_type
    !> debug mode activated?
    logical :: active     = .false.
    !> open separate file for each process
    logical :: debugFiles = .false.
    !> folder to store the debug mesh
    character(len=PathLen) :: debugMesh
    !> unit to write in
    integer :: unit  = -1
    !> output debug output after each step in compute
    logical :: checkEachAlgorithmicStep  = .false.
    !> write element state information to the debugFiles
    logical :: dumpLevelwiseState        = .false.
    !> write halo state information to the debugFiles
    logical :: dumpHaloState             = .false.
    !> write all the required and generated treeIDs in a level-wise manner to
    !! the debug file this includes fluid, ghost and halo elements and can
    !! serve as a detailed debug output
    logical :: dumpTreeIDlists   = .false.
    !> write all the property bits to disk for all elements
    logical :: dumpPropBits      = .false.
    !> write all dependencies of ghost elements to disk
    logical :: dumpDependencies  = .false.
    !> write all dependencies of ghost elements to disk
    logical :: debugDependencies  = .false.
    !> check if the dependencies are correct by cross checking
    logical :: checkDependencies  = .false.
    !> write auxiliary lists to disk such as boundary element lists
    logical :: dumpAuxLists      = .false.
    logical :: unitTests         = .false.
    !> write out source debug statements to debug unit
    logical :: dumpSource        = .false.
    !> write out restart debug info
    logical :: debugRestart  = .false.
    !> trace memory consumption
    logical :: traceMemory  = .false.
    !> Check the state array for NaNs
    logical :: checkNaNs = .false.

    !> Dump boundary elements information
    logical :: dumpBoundaries = .false.

    !> A logger to describe the output capabilities of this debug object.
    type(tem_logging_type) :: logger
  end type tem_debug_type

  type(tem_debug_type) :: main_debug

  interface tem_debug
    module procedure tem_debug_for
    module procedure tem_debug_main
  end interface tem_debug

contains

! ****************************************************************************** !
  !> Read the debug configuration into the debug type 'me'
  !! The debug definition is placed in the main level of the musubi.lua file
  !! It can look like
  !!~~~~~~~~~~~~~~~~~~~~~{.lua}
  !! debug = {
  !!         debugMode = true,        -- default= false
  !!         -- Use a logging object to track output in individual files.
  !!         -- The output messages to this logging object are usually
  !!         -- generated by tem_debug calls.
  !!         logging = { level = 1, -- how detailed should the log be?
  !!                     filename = 'dbg', -- to which file to write the log
  !!                     root_only = .false. } -- should only root write msgs?
  !!
  !!         debugFiles = true,       -- default= false
  !!         -- What to dump into debugFiles? --
  !!          dumpTreeIDs = true,      -- default= false
  !!          dumpPropBits = true,     -- default= false
  !!          dumpAuxLists = true,     -- default= false
  !!          dumpDependencies = true, -- default= false
  !!          dumpState = true,        -- default= false
  !!          dumpHaloState = true,    -- default= false
  !!         --  end debugFiles  --
  !!         debugDependencies = true, -- default= false
  !!         checkDependencies = true, -- default= false
  !!         checkNans = true,         -- default= false
  !!         checkSteps = true,        -- default= false
  !!         debugMesh = 'dbg/mesh_',  -- default= ''
  !!         debugSource = true,       -- default= false
  !!         debugRestart = true,      -- default= false
  !!         traceMemory = true,       -- default= false
  !!        }
  !!~~~~~~~~~~~~~~~~~~~~~
  !! Possible Options are
  !! - active = {true, false}\n
  !!   activate or deactivate the complete debug mode. If deactivated, all
  !!   subsequent definitions are ignored
  !! @todo: Check the following reference
  !! - [[tem_debug_module:tem_debug_type.debugFiles]] `= {true, false}`
  !!   open debug files. They can be accessed by writing to the unit dbgUnit.
  !!   the dbgUnit is given in the tem_debug_module and simply needs to be included
  !!   into the use statements
  !!```
  !!     use tem_debug_module, only: dbgUnit
  !!```
  !! - [[tem_subtree_module:tem_write_debugmesh]] = {true, false}
  !! - unit = {integer}
  !!   debug unit to write to
  !! - checkEachAlgorithmicStep = {true, false}
  !!   output debug output after each step in compute
  !! - dumpLevelwiseState = {true, false}
  !!   write element state information to the debugFiles
  !! - dumpHaloState = {true, false}
  !!   write halo state information to the debugFiles
  !! - [[tem_construction_module:tem_dumptreeidlists]] = {true, false}
  !!   write all the treeIDs fluid, ghost, halo for each level in a nice way to
  !!   the debug files
  !! - dumpPropBits = {true, false}
  !!   write all the property bits to disk for all elements
  !! - dumpDependencies = {true, false}
  !!   write the dependencies of the ghost elements to disk
  !! - checkDependencies = {true, false}
  !!   check if the dependencies are correct by cross checking
  !! - dumpAuxLists = {true, false}
  !!   write the auxiliary lists to the debug file, such as boundary element
  !!   lists
  !! - dumpState = {true, false}
  !!   write the contents of the state vectors to the debug file for all the
  !!   elements of a process
  !! - dumpSource = {true, false}
  !!   write out source debug statements to debug unit
  !! - dumpRestart = {true, false}
  !!   write out restart debug info
  !! - traceMemory = {true, false}
  !!   trace memory consumption
  !! - checkNaNs = {true, false}
  !!   check the state array for NaNs
  !!
  subroutine tem_load_debug(me, conf, rank)
    ! ---------------------------------------------------------------------------
    !> debug type to store information
    type(tem_debug_type) :: me
    !> lua state
    type(flu_state) :: conf
    !> Rank of the calling process
    integer, intent(in) :: rank
    ! ---------------------------------------------------------------------------
    integer :: myrank
    integer :: tc_handle
    integer :: log_tab
    integer :: iError
    ! ---------------------------------------------------------------------------

    ! tables inside the boundary table
    me%active     = .false.

    myrank = rank

    ! Read the number of trackings in the lua file
    call aot_table_open( L = conf, thandle = tc_handle, key = 'debug' )

    call aot_table_open( L       = conf,                                       &
      &                  parent  = tc_handle,                                  &
      &                  key     = 'logging',                                  &
      &                  thandle = log_tab )

    if (log_tab /= 0) then
      call tem_logging_load( conf    = conf,                                   &
        &                    thandle = log_tab,                                &
        &                    rank    = myrank,                                 &
        &                    me      = me%logger )
    else
      ! Deactivate the logging, if there is no configuration for it.
      ! (Setting root_only to true and rank to -1 results in non participating
      !  processes for all processes.)
      call tem_logging_init( me        = me%logger,                            &
        &                    level     = 0,                                    &
        &                    root_only = .true.,                               &
        &                    rank      = -1 )
    end if

    call aot_table_close(L = conf, thandle = log_tab)

    call aot_get_val( L       = conf,                                          &
      &               thandle = tc_handle,                                     &
      &               val     = me%active,                                     &
      &               ErrCode = iError,                                        &
      &               key     = 'debugMode',                                   &
      &               default = .false. )
    if( me%active ) then
      write(logUnit(1),*)
      write(logUnit(1),*) '          ******** DEBUG MODE ACTIVATED **********'
      write(logUnit(1),*)
    endif

    call aot_get_val( L       = conf,                                          &
      &               thandle = tc_handle,                                     &
      &               val     = me%debugFiles,                                 &
      &               ErrCode = iError,                                        &
      &               key     = 'debugFiles',                                  &
      &               default = .false. )

    if( me%debugFiles ) then
      write(logUnit(1),*) 'Opening debug files for each process ...'
    endif

    if( me%debugFiles ) then
      call aot_get_val( L       = conf,                                        &
        &               thandle = tc_handle,                                   &
        &               val     = me%dumpTreeIDlists,                          &
        &               ErrCode = iError,                                      &
        &               key     = 'dumpTreeIDs',                               &
        &               default = .false. )
      if( me%dumpTreeIDlists ) then
        write(logUnit(1),*)' Writing a complete level-wise dump of the treeIDs'&
          &            //' to debug file '
      endif

      call aot_get_val( L       = conf,                                        &
        &               thandle = tc_handle,                                   &
        &               val     = me%dumpPropBits,                             &
        &               ErrCode = iError,                                      &
        &               key     = 'dumpPropBits',                              &
        &               default = .false. )
      if( me%dumpPropBits ) then
        write(logUnit(1),*)' Writing a complete level-wise PropertyBits'//     &
          &            ' to debug file '
      endif

       call aot_get_val( L       = conf,                                       &
         &               thandle = tc_handle,                                  &
         &               val     = me%dumpAuxLists,                            &
         &               ErrCode = iError,                                     &
         &               key     = 'dumpAuxLists',                             &
         &               default = .false. )
      if( me%dumpAuxLists ) then
        write(logUnit(1),*)' Writing auxiliary lists to disk'//' to debug file '
      endif

      call aot_get_val( L       = conf,                                        &
        &               thandle = tc_handle,                                   &
        &               val     = me%dumpDependencies,                         &
        &               ErrCode = iError,                                      &
        &               key     = 'dumpDependencies',                          &
        &               default = .false. )
      if( me%dumpDependencies ) then
        write(logUnit(1),*) ' Writing dependencies to debug file '
      endif

      call aot_get_val( L       = conf,                                        &
        &               thandle = tc_handle,                                   &
        &               val     = me%dumpLevelwiseState,                       &
        &               ErrCode = iError,                                      &
        &               key     = 'dumpState',                                 &
        &               default = .false. )
      if( me%dumpLevelwiseState ) then
        write(logUnit(1),*)' Writing state vector level-wise for each element' &
          &                //' to disk, also ghosts '
      endif
      call aot_get_val( L       = conf,                                        &
        &               thandle = tc_handle,                                   &
        &               val     = me%dumpHaloState,                            &
        &               ErrCode = iError,                                      &
        &               key     = 'dumpHaloState',                             &
        &               default = .true. )
    endif

    call aot_get_val( L       = conf,                                          &
      &               thandle = tc_handle,                                     &
      &               val     = me%debugDependencies,                          &
      &               ErrCode = iError,                                        &
      &               key     = 'debugDependencies',                           &
      &               default = .false. )
    if( me%dumpDependencies ) then
      write(logUnit(1),*) ' Debugging Dependencies:        on'
    endif

    call aot_get_val( L       = conf,                                          &
      &               thandle = tc_handle,                                     &
      &               val     = me%checkDependencies,                          &
      &               ErrCode = iError,                                        &
      &               key     = 'checkDependencies',                           &
      &               default = .false. )
    if( me%dumpDependencies ) then
      write(logUnit(1),*) ' Check horizontal Dependencies: on '
    endif

    call aot_get_val( L       = conf,                                          &
      &               thandle = tc_handle,                                     &
      &               val     = me%checkNaNs,                                  &
      &               ErrCode = iError,                                        &
      &               key     = 'checkNans',                                   &
      &               default = .false. )
    if( me%checkNaNs ) then
      write(logUnit(1),*) ' Checking for NaNs:             on '
    endif

    call aot_get_val( L       = conf,                                          &
      &               thandle = tc_handle,                                     &
      &               val     = me%checkEachAlgorithmicStep,                   &
      &               ErrCode = iError,                                        &
      &               key     = 'checkSteps',                                  &
      &               default = .false. )
    if( me%checkEachAlgorithmicStep ) then
      write(logUnit(1),*)' Checking each algorithmic step and advancing time'//&
        &            ' afterwards '
    endif

    call aot_get_val( L       = conf,                                          &
      &               thandle = tc_handle,                                     &
      &               val     = me%debugMesh,                                  &
      &               ErrCode = iError,                                        &
      &               key     = 'debugMesh',                                   &
      &               default = '' )

    call aot_get_val( L       = conf,                                          &
      &               thandle = tc_handle,                                     &
      &               val     = me%dumpSource,                                 &
      &               ErrCode = iError,                                        &
      &               key     = 'debugSource',                                 &
      &               default = .false. )

    call aot_get_val( L       = conf,                                          &
      &               thandle = tc_handle,                                     &
      &               val     = me%debugRestart,                               &
      &               ErrCode = iError,                                        &
      &               key     = 'debugRestart',                                &
      &               default = .false. )

    ! Trace memory consumption?
    ! Must be done explicitely in the solvers
    call aot_get_val( L       = conf,                                          &
      &               thandle = tc_handle,                                     &
      &               val     = me%traceMemory,                                &
      &               ErrCode = iError,                                        &
      &               key     = 'traceMemory',                                 &
      &               default = .false. )
    if( me%traceMemory ) then
      write(logUnit(1),*)' Tracing memory consumption (if implemented in the'//&
        &            ' solver) ... '
    end if

    call aot_get_val( L       = conf,              &
      &               thandle = tc_handle,         &
      &               val     = me%dumpBoundaries, &
      &               ErrCode = iError,            &
      &               key     = 'dumpBoundaries',  &
      &               default = .false. )

!    call tem_loadRestartConfig( controller = me%restart%controller, conf = conf, &
!                           key = 'debugStates', parent_table = tc_handle  )

    call aot_table_close(L=conf, thandle=tc_handle)
    call tem_horizontalSpacer( fUnit = logUnit(1) )

  end subroutine tem_load_debug
! ****************************************************************************** !


! ****************************************************************************** !
  !> Load the main debugger object
  !!
  subroutine tem_debug_load_main( conf, rank )
    ! ---------------------------------------------------------------------------
    !> Lua state, to get the configuration from.
    type(flu_state) :: conf
    !> Rank of the calling process.
    integer, optional, intent(in) :: rank
    ! ---------------------------------------------------------------------------

    call tem_load_debug( me = main_debug, conf = conf, rank = rank )

    ! Fill shorthand for the debug units.
    dbgUnit = main_debug%logger%funit

  end subroutine tem_debug_load_main
! ****************************************************************************** !


! ****************************************************************************** !
  !> Write a message to a dedicated debug logger.
  !!
  subroutine tem_debug_for( me, level, msg )
    ! ---------------------------------------------------------------------------
    !> Debug type including the logger to use.
    type(tem_debug_type), intent(in) :: me
    !> Log-Level of this message.
    integer, intent(in) :: level
    !> Message to write.
    character(len=*), intent(in) :: msg
    ! ---------------------------------------------------------------------------

    write(me%logger%funit(level),*) msg

  end subroutine tem_debug_for
! ****************************************************************************** !


! ****************************************************************************** !
  !> Write a message to the main debug logger.
  !!
  subroutine tem_debug_main(level, msg)
    ! ---------------------------------------------------------------------------
    !> Log-Level of this message.
    integer, intent(in) :: level
    !> Message to write.
    character(len=*), intent(in) :: msg
    ! ---------------------------------------------------------------------------

    write(dbgUnit(level),*) msg

  end subroutine tem_debug_main
! ****************************************************************************** !


! ****************************************************************************** !
  !> print an array to the debugunit
  !!
  subroutine tem_reportStatus( level, text, debug, string )
    ! ---------------------------------------------------------------------------
    !> level for debug output
    integer, intent(in) :: level
    !> Array title in debug output for easy identification in the file
    character( len=* ) :: text
    !> optional debug type
    type(tem_debug_type), optional, intent(in) :: debug
    !> optional additional string extending the title
    character( len=* ), optional :: string
    ! ---------------------------------------------------------------------------
    character( len=labelLen ) :: traceString
    character(len=labelLen) :: stat_str
    integer :: nUnit
    logical :: isActive
    ! ---------------------------------------------------------------------------
    if( present( debug ))then
      nUnit = debug%logger%funit(level)
      isActive = tem_logging_isActive( debug%logger, level )
    else
      nUnit = main_debug%logger%funit(level)
      isActive = tem_logging_isActive( main_debug%logger, level )
    end if

    if( present ( string )) then
      traceString = string
    else
      traceString = 'VmHWM'
    endif

    if( isActive ) then
      stat_str = my_status_string(trim(traceString))
      write(nUnit,*) trim(stat_str), '  ',  trim(text)
    end if

  end subroutine tem_reportStatus
! ****************************************************************************** !


! ****************************************************************************** !
  !> print an array to the debugunit
  !!
  subroutine tem_print_array( me, nVals, itemLength, title, outUnit )
    ! ---------------------------------------------------------------------------
    !> Array title in debug output for easy identification in the file
    character( len=* ),optional :: title
    !> number of values in array me
    integer, intent(in) :: nVals
    !> long array to write to debug file
    integer(kind=long_k), intent(in) :: me(:)
    !> how many characters needs each item of the array to output
    integer,optional :: itemLength
    !> output unit
    integer, intent(in) :: outUnit
    ! ---------------------------------------------------------------------------
    integer :: iVal, iLine, nLines
    integer :: itemLength_loc
    character( len=buffer_length ) :: buffer
    ! ---------------------------------------------------------------------------

    if( present(itemlength)) then
      itemlength_loc = itemlength
    else
      itemlength_loc = 8
    endif

    write(outUnit,*) ''
    write(outUnit,*) '------------------------------------------------'

    if( present(title)) then
      write(outUnit,"(A,I0)") '   '//trim(title)//', nVals: ', nVals
      write(outUnit,*) ''
    endif

    buffer = ''

    nLines = nVals / itemLength_loc + 1

    do iLine = 1, nLines - 1
      buffer = ''
      do iVal = 1, itemLength_loc
        call append_to_buffer( buffer, me(iVal) )
      end do
      write(outUnit,*) trim(buffer)
    end do

    buffer = ''
    do iVal = 1, mod( nVals, itemLength_loc )
      call append_to_buffer( buffer, me(iVal) )
    end do
    write(outUnit,*) trim(buffer)

    write(outUnit,*) '------------------------------------------------'
    write(outUnit,*) ''

    ! do iVal = 1, nVals
    !   buffer = trim(buffer)//'  '//tem_toStr( me(iVal), logger )
    !   if( iVal == nVals .or. mod( iElem, itemLength_loc ) == 0) then
    !     write(outUnit,*) trim(buffer)
    !     buffer = ''
    !   endif
    ! enddo

    contains

      subroutine append_to_buffer(buffer, val )
        character( len=buffer_length ), intent(inout) :: buffer
        integer(kind=long_k), intent(in) :: val
        character(len=SolSpecLen) :: str

        str = tem_toStr( val, main_debug%logger )
        if ( (2 + len_trim(buffer) + len_trim(str)) > buffer_length) then
          write(logunit(1),*) "Can't append value to output in tem_print_array."
          write(logunit(1),*) "  Internal buffer too short."
        else
          buffer = trim(buffer) // '  ' // str
        endif

      end subroutine append_to_buffer

  end subroutine tem_print_array
! ****************************************************************************** !


end module tem_debug_module
! ****************************************************************************** !
