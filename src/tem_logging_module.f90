! Copyright (c) 2013, 2016, 2018-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2015, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
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
!! Providing a logging infrastructure to easily control
!! the verbosity of output.
!!
!! The internal log_level can be set at runtime and output
!! will only happen, if the log_level is higher than that
!! of the messages to write.
!!
!! The numbering of levels is limited to 20, the higher the number the less
!! important the message is:
!!
!! -  1 A message that should always appear (except the logger is explicitly
!!     turned off)
!! -  3 Informational warnings of less importance
!! -  6 Detailed information
!! - 10 Very detailed debugging information
!!
!! Messages should be defined with positive integers between 1 and 10 and the
!! list above provides some idea on which levels to use.
!! Log levels beyond 10 should only be used for temporary outputs during
!! development.
!!
!! There is a primary logging data structure provided by the module, which
!! will provide output by the root rank usually to the standard output.
!!
!! To write a message to the log, use its funit.
!! The funit is an array of file units, that are either connected to the
!! configured output file or to the null device, depending on the configured
!! log level.
!! All units beyond the log level are connected to the null device, while
!! those up that level are connected to the file.
!! For the primary logger a shorthand is provided by the module variable
!! logUnit.
!! Thus, to write some message on log level 4 to the primary logger
!! you do something like this:
!!
!! write(logUnit(4),*) 'some real: ', x
!!
!! This message will only appear, if the log level in the configuration is
!! set to 4 or higher.
!!
!! Unfortunately debug outputs are costly, even if written to the null device
!! they might heavily affect performance.
!! Thus you should not put log messages into hot code parts, that is mostly
!! into long or deeply nested loops.
!! If there is the need for debug output in such parts, there is the possibility
!! to enclose the calls in CoCo preprocessing commands:
!!
!!```lua
!! ?? IF (loglvl >= 3) THEN
!!   write(logunit(3),*) 'We only write this if compilation includes this.'
!! ?? ENDIF
!!```
!!
!! As loglvl needs to be always defined you currently need to include the
!! logMacros.inc file everywhere you want to use that.
!! However, it might be possible to shift the default definition of loglvl=0
!! to the setfile (needs to be checked).
!! To set a certain loglvl at compile time it is then necessary to set the
!! COCOFLAGS environment variable:
!!
!!```lua
!! export COCOFLAGS='-Dloglvl=3'
!!```
!!
!! In addition there are some helping routines defined, that help to convert
!! data to strings.
!! To convert a number or logical to a string use tem_toStr.
!! For example you can write a message to the primary log with the real "ar"
!! like this:
!!```fortran
!! call tem_log(level=1, msg='some real: '//trim(tem_toStr(ar))//' !')
!!```
!!
module tem_logging_module

  ! include treelm modules
  use env_module, only: pathLen, SolSpecLen, single_k, double_k, int_k, long_k, &
    &                   newUnit, stdOutUnit, tem_connect_toNull

  ! include aotus modules
  use aotus_module,     only: aot_get_val, flu_State
  use aot_table_module, only: aot_table_open, aot_table_close

  implicit none

  private

  !> The last logging unit, defining the length of log unit arrays.
  integer, parameter, public :: tem_last_lu = 21

  integer, public :: logUnit(0:tem_last_lu)

  integer, parameter :: form_len = 10

  !> A message that should always appear (except the logger is explicitly
  !! turned off)
  integer, parameter, public :: llerror = 1
  !> Informational warnings of less importance
  integer, parameter, public :: llwarning = 3
  !> Detailed information
  integer, parameter, public :: llinfo = 6
  !> Very detailed debugging information
  integer, parameter, public :: lldebug = 10

  type tem_logging_type
    integer :: log_level
    integer :: funit(0:tem_last_lu)
    logical :: participating
    character(len=form_len) :: real_form
    character(len=form_len) :: int_form
  end type tem_logging_type

  interface tem_logging_isActive
    module procedure tem_logging_isActive_for
    module procedure tem_logging_isActive_primary
  end interface tem_logging_isActive

  interface tem_logging_unit
    module procedure tem_logging_unit_for
    module procedure tem_logging_unit_primary
  end interface tem_logging_unit

  interface tem_log_write
    module procedure tem_log_write_for
    module procedure tem_log_write_primary
  end interface tem_log_write

  interface tem_log
    module procedure tem_log_for
    module procedure tem_log_primary
  end interface


  interface tem_logging_init
    module procedure tem_logging_init_logger
    module procedure tem_logging_init_primary
  end interface

  interface tem_toStr
    module procedure tem_r2str
    module procedure tem_d2str
    module procedure tem_i2str
    module procedure tem_l2str
    module procedure tem_b2str
    module procedure tem_r2str_arr
    module procedure tem_d2str_arr
    module procedure tem_i2str_arr
    module procedure tem_l2str_arr
    module procedure tem_b2str_arr
  end interface

  public :: tem_logging_type
  public :: tem_logging_isActive
  public :: tem_logging_unit
  public :: tem_log_write
  public :: tem_log
  public :: tem_toStr
  public :: tem_logging_init
  public :: tem_logging_load
  public :: tem_logging_load_primary
  public :: tem_logging_init_primary

  type(tem_logging_type) :: primary

  contains

! ****************************************************************************** !
  !> Initialize a logging data type.
  !!
  subroutine tem_logging_init_logger( me, level, rank, filename, root_only,    &
    &                                 real_form, int_form )
    ! ---------------------------------------------------------------------------
    !> Logger to initialize
    type(tem_logging_type), intent(out) :: me
    !> Level of output to log with this logger
    integer, intent(in) :: level
    !> Rank of the process executing the initialization.
    integer, intent(in) :: rank

    !> File to write output to, default is null device.
    !!
    !! If this is empty or not provided, the output will be to the
    !! null device.
    !! To send messages to the stdout set this parameter to
    !! '/stdout:'.
    character(len=*), optional, intent(in) :: filename

    !> Indicate if root only should write messages
    !! Default is true.
    logical,          optional, intent(in) :: root_only

    !> Format to write real numbers.
    !! Default is 'EN12.3'
    character(len=*), optional, intent(in) :: real_form

    !> Format to write integer numbers.
    !! Default is 'I0'
    character(len=*), optional, intent(in) :: int_form
    ! ---------------------------------------------------------------------------
    logical :: root_out
    character(len=7) :: rankstamp
    character(len=pathLen) :: fname
    logical             :: nUnitOpened
    integer             :: UnitNumber
    logical             :: file_exists
    ! ---------------------------------------------------------------------------

    me%log_level = level

    if (present(root_only)) then
      root_out = root_only
    else
      root_out = .true.
    end if

    if (root_out) then
      me%participating = (rank == 0)
    else
      me%participating = .true.
    end if

    if (present(filename)) then
      if (trim(filename) == '/stdout:') then
        ! Only root should write to stdout.
        root_out = .true.
        if (rank == 0) then
          fname = filename
          me%participating = .true.
        else
          fname = ''
          me%participating = .false.
        end if
      else
        fname = trim(filename)
      end if
    else
      fname = ''
    end if

    if (me%participating .and. (trim(adjustl(fname)) /= '')) then
      ! If I am participating and the filename is not actually empty,
      ! proceed connecting to appropriate external files.
      if (trim(adjustl(fname)) == '/stdout:') then
        ! Messages should be written to the standard output.
        me%funit(0) = stdOutUnit
      else
        ! Messages should be written to some file.
        if (root_out) then
          rankstamp = ''
        else
          write(rankstamp, '(a1,I6.6)') '.', rank
        end if
        ! changes for dynamic load balancing
        ! check if the file exists
        inquire( file  = trim(fname)//trim(rankstamp),                         &
          &      exist = file_exists )
        if( file_exists )then
          ! in case the file exists, check wether it is already opened somewhere
          ! else (dyn load balancing)
          inquire( file   = trim(fname)//trim(rankstamp),                      &
            &      opened = nUnitOpened,                                       &
            &      number = UnitNumber )
          if (nUNitOpened) then
            me%funit(0) = UnitNumber
          else
            me%funit(0) = newunit()
            open( unit     = me%funit(0),                                      &
              &   file     = trim(fname)//trim(rankstamp),                     &
              &   position = 'APPEND',                                         &
              &   status   = 'OLD' )
          end if
        else
          me%funit(0) = newunit()
          open( unit     = me%funit(0),                                        &
            &   file     = trim(fname)//trim(rankstamp),                       &
            &   status   = 'REPLACE' )
        end if
      end if
    else
      ! Output should be deactivated for this logger.
      ! Connect its unit to the null device.
      call tem_connect_toNull(me%funit(0))
    end if

    ! Always connect the last logging unit to the null device.
    call tem_connect_toNull(me%funit(tem_last_lu))

    ! Set all units according to the configured logging level:
    me%funit(1:level) = me%funit(0)
    me%funit(level+1:tem_last_lu-1) = me%funit(tem_last_lu)

    if (present(real_form)) then
      me%real_form = real_form
    else
      me%real_form = 'EN12.3'
    end if

    if (present(int_form)) then
      me%int_form = int_form
    else
      me%int_form = 'I0'
    end if

    !! @todo: this should move to a proper place and have a proper format!
    !if (present(root_only)) then
    !  write(me%fUnit(0),*) 'rank= ', rank, 'filename= ', trim(filename), &
    !    &                                'root_only= ', root_only
    !end if
    !write(me%fUnit(0),*) rank, root_out, 'log unit ', me%fUnit

  end subroutine tem_logging_init_logger
! ****************************************************************************** !


! ****************************************************************************** !
  !> Initialize the primary logger (default to stdout
  !! instead of null device).
  !!
  subroutine tem_logging_init_primary( level, rank, filename, &
    &                                  root_only,             &
    &                                  real_form, int_form    )
    ! ---------------------------------------------------------------------------
    !> Level of output to log with this logger
    integer, intent(in) :: level
    !> Rank of the process executing the initialization.
    integer, intent(in) :: rank

    !> Indication whether only root should print log messages.
    logical, optional, intent(in) :: root_only

    !> File to write output to, default is standard out.
    !!
    !! To send messages to the stdout set this parameter to
    !! '/stdout:'.
    character(len=*), optional, intent(in) :: filename

    !> Format to write real numbers.
    !! Default is 'EN12.3'
    character(len=*), optional, intent(in) :: real_form

    !> Format to write integer numbers.
    !! Default is 'I0'
    character(len=*), optional, intent(in) :: int_form
    ! ---------------------------------------------------------------------------
    character(len=PathLen) :: fname
    ! ---------------------------------------------------------------------------

    if (present(filename)) then
      if (trim(filename) == '') then
        fname = '/stdout:'
      else
        fname = adjustl(filename)
      end if
    else
      fname = '/stdout:'
    end if

    call tem_logging_init_logger( me        = primary,     &
      &                           level     = level,       &
      &                           rank      = rank,        &
      &                           root_only = root_only,   &
      &                           filename  = trim(fname), &
      &                           real_form = real_form,   &
      &                           int_form  = int_form     )

    ! Set the shorthand module variable for the primary file units:
    logUnit = primary%funit

  end subroutine tem_logging_init_primary
! ****************************************************************************** !


! ****************************************************************************** !
  !> Load logging configuration from a Lua table and return the initialized
  !! logger.
  !!
  !! The table provided via the thandle has to contain all settings to describe
  !! the logging object.
  !! For all settings there are defaults, and the routine will silently assume
  !! those in case of failing to read a setting for any reason.
  !!
  !! For the primary logger only level is mandatory,
  !! filename, real_form and int_form are also used but optional.
  !! In contrast to other logging facilities, the primary logging defaults to
  !! stdout and not to the null device.
  !!
  subroutine tem_logging_load( conf, thandle, rank, me )
    ! ---------------------------------------------------------------------------
    !>
    type(flu_State) :: conf
    !>
    integer, intent(in) :: thandle
    !>
    integer, intent(in) :: rank
    !>
    type(tem_logging_type), optional, intent(out) :: me
    ! ---------------------------------------------------------------------------
    integer :: level
    character(len=pathLen) :: filename
    logical :: root_only
    character(len=form_len) :: real_form
    character(len=form_len) :: int_form
    integer :: iError
    ! ---------------------------------------------------------------------------

    call aot_get_val( level, ErrCode = iError,                                 &
      &               L = conf, thandle = thandle, key = 'level',              &
      &               default = 1 )
    call aot_get_val( filename, ErrCode = iError,                              &
      &               L = conf, thandle = thandle, key = 'filename',           &
      &               default = '' )
    call aot_get_val( root_only, ErrCode = iError,                             &
      &               L = conf, thandle = thandle, key = 'root_only',          &
      &               default = .true. )
    call aot_get_val( real_form, ErrCode = iError,                             &
      &               L = conf, thandle = thandle, key = 'real_form',          &
      &               default = 'EN12.3' )
    call aot_get_val( int_form, ErrCode = iError,                              &
      &               L = conf, thandle = thandle, key = 'int_form',           &
      &               default = 'I0' )


    if (present(me)) then
      call tem_logging_init_logger( me        = me,                            &
        &                           level     = level,                         &
        &                           rank      = rank,                          &
        &                           filename  = filename,                      &
        &                           root_only = root_only,                     &
        &                           real_form = real_form,                     &
        &                           int_form  = int_form   )
    else
      call tem_logging_init_primary( level     = level,                        &
        &                            rank      = rank,                         &
        &                            root_only = root_only,                    &
        &                            filename  = filename,                     &
        &                            real_form = real_form,                    &
        &                            int_form  = int_form )
    end if

  end subroutine tem_logging_load
! ****************************************************************************** !


! ****************************************************************************** !
  !> Load the primary logger from a Lua script under its default name of
  !! logging as global variable.
  !!
  !! If it is not set by the user, a default logging level of 1 is assumed.
  !!
  subroutine tem_logging_load_primary( conf, rank )
    ! ---------------------------------------------------------------------------
    !>
    type(flu_state) :: conf
    !>
    integer :: rank
    ! ---------------------------------------------------------------------------
    integer :: thandle
    ! ---------------------------------------------------------------------------

    call aot_table_open(L = conf, thandle = thandle, key = 'logging')
    if (thandle == 0) then
      call tem_logging_init_primary(level = 1, rank = rank)
    else
      call tem_logging_load(conf = conf, thandle = thandle, rank = rank)
    end if
    call aot_table_close(L = conf, thandle = thandle)

  end subroutine tem_logging_load_primary
! ****************************************************************************** !


! ****************************************************************************** !
  !> Check if the given logger is active for the given level.
  !!
  function tem_logging_isActive_for(me, level) result(isActive)
    ! ---------------------------------------------------------------------------
    !>
    type(tem_logging_type), intent(in) :: me
    !>
    integer, intent(in) :: level
    !>
    logical :: isActive
    ! ---------------------------------------------------------------------------

    isActive = (me%participating .and. (me%log_level >= level))

  end function tem_logging_isActive_for
! ****************************************************************************** !


! ****************************************************************************** !
  !> Check if the primary logger is active for the given level.
  !!
  function tem_logging_isActive_primary(level) result(isActive)
    ! ---------------------------------------------------------------------------
    !>
    integer, intent(in) :: level
    !>
    logical :: isActive
    ! ---------------------------------------------------------------------------

    isActive = (primary%participating .and. (primary%log_level >= level))

  end function tem_logging_isActive_primary
! ****************************************************************************** !


! ****************************************************************************** !
  !> Return the unit provided by a given log
  !!
  function tem_logging_unit_for(me) result(unit)
    ! ---------------------------------------------------------------------------
    !>
    type(tem_logging_type), intent(in) :: me
    !>
    integer :: unit
    ! ---------------------------------------------------------------------------

    unit = me%funit(0)

  end function tem_logging_unit_for
! ****************************************************************************** !


! ****************************************************************************** !
  !> Return the unit provided by the primary log
  !!
  function tem_logging_unit_primary() result(unit)
    ! ---------------------------------------------------------------------------
    !>
    integer :: unit
    ! ---------------------------------------------------------------------------

    unit = primary%funit(0)

  end function tem_logging_unit_primary
! ****************************************************************************** !


! ****************************************************************************** !
  !> Write msg unconditionally to the logger given in me.
  !!
  subroutine tem_log_write_for(me, msg)
    ! ---------------------------------------------------------------------------
    !>
    type(tem_logging_type), intent(in) :: me
    !>
    character(len=*), intent(in) :: msg
    ! ---------------------------------------------------------------------------

    write(me%funit(0), '(a)') msg

  end subroutine tem_log_write_for
! ****************************************************************************** !


! ****************************************************************************** !
  !> Write msg unconditionally to the primary logger.
  !!
  subroutine tem_log_write_primary(msg)
    ! ---------------------------------------------------------------------------
    !>
    character(len=*), intent(in) :: msg
    ! ---------------------------------------------------------------------------

    write(primary%funit(0), '(a)') msg

  end subroutine tem_log_write_primary
! ****************************************************************************** !


! ****************************************************************************** !
  !> Log a message in the given logger.
  !!
  subroutine tem_log_for( me, level, msg )
    ! ---------------------------------------------------------------------------
    !>
    type(tem_logging_type), intent(in) :: me
    !>
    integer, intent(in) :: level
    !>
    character(len=*), intent(in) :: msg
    ! ---------------------------------------------------------------------------

    write(me%funit(level),'(a)') msg

  end subroutine tem_log_for
! ****************************************************************************** !


! ****************************************************************************** !
  !> Log a message in the primary logger.
  !!
  subroutine tem_log_primary( level, msg )
    ! ---------------------------------------------------------------------------
    !>
    integer, intent(in) :: level
    !>
    character(len=*), intent(in) :: msg
    ! ---------------------------------------------------------------------------

    write(primary%funit(level),'(a)') msg

  end subroutine tem_log_primary
! ****************************************************************************** !


! ****************************************************************************** !
  !> Convert a real to a string according to the format provided in the logger.
  !!
  function tem_r2str( val, logger ) result(str)
    ! ---------------------------------------------------------------------------
    !>
    real(kind=single_k), intent(in) :: val
    !>
    type(tem_logging_type), optional, intent(in) :: logger
    !>
    character(len=SolSpecLen) :: str
    ! ---------------------------------------------------------------------------
    character(len=form_len) :: form
    ! ---------------------------------------------------------------------------

    if (present(logger)) then
      form = logger%real_form
    else
      form = primary%real_form
    end if
    write(str, '(' // trim(form) // ')') val

  end function tem_r2str
! ****************************************************************************** !


! ****************************************************************************** !
  !> Converts a real "array" to a single string according to the format provided
  !! in the logger.
  !!
  function tem_r2str_arr( val, sep, logger) result(str)
    ! ---------------------------------------------------------------------------
    !> array to convert
    real(kind=single_k), intent(in) :: val(:)
    !> seperator between array elements
    character(len=*), intent(in)    :: sep
    !> logger type which provides output format
    type(tem_logging_type), optional, intent(in) :: logger
    !> output string
    character(len=SolSpecLen) :: str
    ! ---------------------------------------------------------------------------
    character(len=SolSpecLen)     :: temp_str
    character(len=form_len)       :: form
    integer                       :: iter, idx, offset, off_sep
    integer                       :: instrlen, maxstrlen
    ! ---------------------------------------------------------------------------

    if (present(logger)) then
      form = logger%real_form
    else
      form = primary%real_form
    end if

    maxstrlen = 10
    instrlen = size(val)
    off_sep = len(sep)
    offset = 0

    str = ''
    idx = 1
    ! If the length of input string is within limits
    if( instrlen .le. maxstrlen ) then
      do iter=1, instrlen
        ! Convert the ith element of array into character
        write(temp_str, '(' // trim(form) // ')') val(iter)
        ! Append the above obtained character to string followed by separator
        offset = len_trim(temp_str)
        str(idx:idx+offset-1) = trim(temp_str(1:offset) )
        if( iter .ne. instrlen ) then
          str(idx+offset:idx+offset+off_sep) = sep(1:off_sep)
        end if
        idx = idx + offset + off_sep + 1
      end do
      ! Clip the final string to removed unwanted stuff
      str = str(1:instrlen*(offset+off_sep+1)-off_sep-1)
    ! If not then print 1,2,3,...,last
    else
      do iter=1, 3
        write(temp_str, '(' // trim(form) // ')') val(iter)
        offset = len_trim(temp_str)
        str(idx:idx+offset-1) = trim(temp_str(1:offset) )
        str(idx+offset:idx+offset+off_sep) = sep(1:off_sep)
        idx = idx + offset + off_sep + 1
      end do
      ! Now add ,..., and the last entry of array
      str(idx:idx+2)           = '...'
      str(idx+3:idx+3+off_sep) = sep
      write(temp_str, '(' // trim(form) // ')') val(instrlen)
      offset = len_trim(temp_str)
      str(idx+4+off_sep:idx+4+off_sep+offset) = trim(temp_str(1:offset) )
      ! Clip the final string to removed unwanted stuff
      str = str(1:4*(offset+off_sep)-off_sep+4)
    end if

  end function tem_r2str_arr
! ****************************************************************************** !


! ****************************************************************************** !
  !> Converts a double to a string according to the format provided in the
  !! logger.
  !!
  function tem_d2str( val, logger ) result(str)
    ! ---------------------------------------------------------------------------
    !>
    real(kind=double_k), intent(in) :: val
    !>
    type(tem_logging_type), optional, intent(in) :: logger
    !>
    character(len=SolSpecLen) :: str
    ! ---------------------------------------------------------------------------
    character(len=form_len) :: form
    ! ---------------------------------------------------------------------------

    if (present(logger)) then
      form = logger%real_form
    else
      form = primary%real_form
    end if
    write(str, '(' // trim(form) // ')') val

  end function tem_d2str
! ****************************************************************************** !


! ****************************************************************************** !
  !> Converts an array of doubles to a string according to the format provided
  !! in the logger.
  !!
  function tem_d2str_arr(val, sep, logger) result(str)
    ! ---------------------------------------------------------------------------
    !>
    real(kind=double_k), intent(in) :: val(:)
    !>
    character(len=*), intent(in)    :: sep
    !>
    type(tem_logging_type), optional, intent(in) :: logger
    !>
    character(len=SolSpecLen) :: str
    ! ---------------------------------------------------------------------------
    character(len=SolSpecLen)     :: temp_str
    character(len=form_len)       :: form
    integer                       :: iter, idx, offset, off_sep
    ! ---------------------------------------------------------------------------
    if (present(logger)) then
      form = logger%real_form
    else
      form = primary%real_form
    end if

    off_sep = len(sep)
    offset = 0

    str = ''
    idx = 1
    do iter=1, size(val)
      ! Convert the ith element of array into character
      write(temp_str, '(' // trim(form) // ')') val(iter)
      ! Append the above obtained character to string followed by separator
      offset = len_trim(temp_str)
      str(idx:idx+offset-1) = trim(temp_str(1:offset) )
      if( iter .ne. size(val) ) then
        str(idx+offset:idx+offset+off_sep) = sep(1:off_sep)
      end if
      idx = idx + offset + off_sep + 1
    end do
    ! Clip the final string to removed unwanted stuff
    str = str(1:size(val)*(offset+off_sep+1)-off_sep-1)

  end function tem_d2str_arr
! ****************************************************************************** !


! ****************************************************************************** !
  !> Converts an integer to a string according to the format provided in the
  !! logger.
  !!
  function tem_i2str( val, logger ) result(str)
    ! ---------------------------------------------------------------------------
    !>
    integer(kind=int_k), intent(in) :: val
    !>
    type(tem_logging_type), optional, intent(in) :: logger
    !>
    character(len=SolSpecLen) :: str
    ! ---------------------------------------------------------------------------
    character(len=form_len) :: form
    ! ---------------------------------------------------------------------------

    if (present(logger)) then
      form = logger%int_form
    else
      form = primary%int_form
    end if
    write(str, '(' // trim(form) // ')') val

  end function tem_i2str
! ****************************************************************************** !


! ****************************************************************************** !
  !> Converts an array of integers to a string according to the format provided
  !! in the logger.
  !!
  function tem_i2str_arr(val, sep, logger) result(str)
    ! ---------------------------------------------------------------------------
    !> array to convert
    integer(kind=int_k), intent(in)          :: val(:)
    !>
    character(len=*), intent(in) :: sep
    !>
    type(tem_logging_type), optional, intent(in) :: logger
    !>
    character(len=SolSpecLen) :: str
    ! ---------------------------------------------------------------------------
    character(len=SolSpecLen)     :: temp_str
    character(len=form_len)       :: form
    integer                       :: iter, idx, offset, off_sep
    ! ---------------------------------------------------------------------------
    if (present(logger)) then
      form = logger%int_form
    else
      form = primary%int_form
    end if

    off_sep = len(sep)
    offset = 0

    str = ''
    idx = 1
    do iter=1, size(val)
      ! Convert the ith element of array into character
      write(temp_str, '(' // trim(form) // ')') val(iter)
      ! Append the above obtained character to string followed by separator
      offset = len_trim(temp_str)
      str(idx:idx+offset-1) = trim(temp_str(1:offset) )
      if( iter .ne. size(val) ) then
        str(idx+offset:idx+offset+off_sep) = sep(1:off_sep)
      end if
      idx = idx + offset + off_sep + 1
    end do
    ! Clip the final string to removed unwanted stuff
    str = str(1:size(val)*(offset+off_sep+1)-off_sep-1)

  end function tem_i2str_arr
! ****************************************************************************** !


! ****************************************************************************** !
  !> Converts a long to a string according to the format provided in the
  !! logger.
  !!
  function tem_l2str( val, logger ) result(str)
    ! ---------------------------------------------------------------------------
    !>
    integer(kind=long_k), intent(in) :: val
    !>
    type(tem_logging_type), optional, intent(in) :: logger
    !>
    character(len=SolSpecLen) :: str
    ! ---------------------------------------------------------------------------
    character(len=form_len) :: form
    ! ---------------------------------------------------------------------------

    if (present(logger)) then
      form = logger%int_form
    else
      form = primary%int_form
    end if
    write(str, '(' // trim(form) // ')') val

  end function tem_l2str
! ****************************************************************************** !


! ****************************************************************************** !
  !> Converts an array of longs to a string according to the format provided
  !! in the logger.
  !!
  function tem_l2str_arr(val, sep, logger) result(str)
    ! ---------------------------------------------------------------------------
    !>
    integer(kind=long_k), intent(in) :: val(:)
    !>
    character(len=*), intent(in)     :: sep
    !>
    type(tem_logging_type), optional, intent(in) :: logger
    !>
    character(len=SolSpecLen) :: str
    ! ---------------------------------------------------------------------------
    character(len=SolSpecLen)     :: temp_str
    character(len=form_len)       :: form
    integer                       :: iter, idx, offset, off_sep
    ! ---------------------------------------------------------------------------
    if (present(logger)) then
      form = logger%int_form
    else
      form = primary%int_form
    end if

    off_sep = len(sep)
    offset = 0

    str = ''
    idx = 1
    do iter=1, size(val)
      ! Convert the ith element of array into character
      write(temp_str, '(' // trim(form) // ')') val(iter)
      ! Append the above obtained character to string followed by separator
      offset = len_trim(temp_str)
      str(idx:idx+offset-1) = trim(temp_str(1:offset) )
      if( iter .ne. size(val) ) then
        str(idx+offset:idx+offset+off_sep) = sep(1:off_sep)
      end if
      idx = idx + offset + off_sep + 1
    end do
    ! Clip the final string to removed unwanted stuff
    str = str(1:size(val)*(offset+off_sep+1)-off_sep-1)

  end function tem_l2str_arr
! ****************************************************************************** !


! ****************************************************************************** !
  !> Converts a bool to a string.
  !!
  function tem_b2str( val ) result(str)
    ! ---------------------------------------------------------------------------
    !>
    logical, intent(in) :: val
    !>
    character(len=SolSpecLen) :: str
    ! ---------------------------------------------------------------------------
    character(len=SolSpecLen) :: tmp
    ! ---------------------------------------------------------------------------

    write(tmp, *) val
    str = adjustl(tmp)

  end function tem_b2str
! ****************************************************************************** !


! ****************************************************************************** !
  !> Converts an array of booleans to a string.
  !!
  function tem_b2str_arr( val, sep) result(str)
    ! ---------------------------------------------------------------------------
    !>
    logical, intent(in) :: val(:)
    !>
    character(len=*), intent(in) :: sep
    !>
    character(len=SolSpecLen) :: str
    ! ---------------------------------------------------------------------------
    character(len=SolSpecLen)     :: tmp, temp_str
    integer                       :: iter, idx, offset, off_sep
    ! ---------------------------------------------------------------------------
    off_sep = len(sep)
    idx = 1
    offset = 0

    str = ''
    do iter=1, size(val)
      ! Convert the ith element of array into character
      write(tmp, *) val(iter)
      temp_str = adjustl(tmp)
      ! Append the above obtained character to string followed by separator
      offset = len_trim(temp_str)
      str(idx:idx+offset-1) = trim(temp_str(1:offset) )
      if( iter .ne. size(val) ) then
        str(idx+offset:idx+offset+off_sep) = sep(1:off_sep)
      end if
      idx = idx + offset + off_sep + 1
    end do
    ! Clip the final string to removed unwanted stuff
    str = str(1:size(val)*(offset+off_sep)-off_sep-1)

  end function tem_b2str_arr
! ****************************************************************************** !


end module tem_logging_module
! ****************************************************************************** !
