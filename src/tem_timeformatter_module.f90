! Copyright (c) 2025 Harald Klimach <harald.klimach@dlr.de>
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
! **************************************************************************** !
!
!> Formatting of time points.
!!
!! The formatters alllow the specification of how time-stamps for given time
!! objects are to be formatted.
!!
module tem_timeformatter_module
  use env_module, only: labelLen

  use aotus_module, only: flu_State, &
    &                     aot_get_val
  use aot_table_module, only: aot_table_open, &
    &                         aot_table_close

  use tem_time_module, only: tem_time_type

  implicit none

  private

  public :: tem_timeformatter_type
  public :: tem_timeformatter_load
  public :: tem_timeformatter_init

  character(len=8), parameter :: default_form = '(EN12.3)'

  type tem_timeformatter_type
    character(len=labelLen) :: timeform = default_form
    procedure(timestamp), pointer, pass(formatter) :: stamp => null()
  end type tem_timeformatter_type

  abstract interface
    function timestamp(formatter, time) result(stamp)
      import :: tem_time_type, tem_timeformatter_type, labelLen
      class(tem_timeformatter_type), intent(in) :: formatter
      type(tem_time_type), intent(in) :: time
      character(len=LabelLen) :: stamp
    end function timestamp
  end interface


contains


  ! ************************************************************************ !
  !> Reading a timeformatter description from a Lua script given by conf.
  !!
  !! Initializing a timeformatter
  !! * timeform defines the formatting string to be used for the timestamp
  !!   this defaults to `default_form`
  !! * stamp defines the routine to use for the timestamp generation and
  !!   defaults to `tem_timeformatter_sim_stamp`
  function tem_timeformatter_init(timeform, stamp) result(formatter)
    character(len=*), optional, intent(in) :: timeform
    procedure(timestamp), optional :: stamp
    type(tem_timeformatter_type) :: formatter

    if (present(timeform)) then
      formatter%timeform = timeform
    else
      formatter%timeform = default_form
    end if

    if (present(stamp)) then
      formatter%stamp => stamp
    else
      formatter%stamp => tem_timeformatter_sim_stamp
    end if

  end function tem_timeformatter_init
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Reading a timeformatter description from a Lua script given by conf.
  !!
  !! The timeformatter description has to be provided in the table given by
  !! thandle.
  !! There are two settings for the formatting:
  !! - use_iter: use the number of iterations rather than simulated time
  !!             Default: false, but may be overwritten by caller.
  !! - simform: Fortran formatting string for the real number of the
  !!            simulated time.
  !!            Default: '(EN12.3)'
  !! Thus, the configuration looks like:
  !!
  !! \verbatim
  !! timeformat = { use_iter = false, simform = '(EN12.3)' }
  !! \end verbatim
  !!
  !! If timeformat is not a table, the value is attempted to be read as a
  !! string for the simform, so it is also possible to write:
  !!
  !! \verbatim
  !! timeformat = '(EN15.6)'
  !! \end verbatim
  !!
  !! The key to be used for reading these settings could also be changed by
  !! the caller by passing the key to be used in the key argument.
  subroutine tem_timeformatter_load(me, conf, key, parent, use_iter_default)
    ! -------------------------------------------------------------------- !
    !> Time to be read from the Lua script
    type(tem_timeformatter_type), intent(out) :: me

    !> Handle to the Lua script.
    type(flu_state), intent(inout) :: conf

    !> Name of the table containing the time definition. Default: 'time'.
    character(len=*), intent(in), optional :: key

    !> Handle to the parent table.
    integer, intent(in), optional :: parent

    !> Default for the use_iter setting
    logical, intent(in), optional :: use_iter_default
    ! -------------------------------------------------------------------- !
    integer :: iErr
    character(len=labelLen) :: loc_key
    character(len=labelLen) :: timeform
    integer :: thandle
    logical :: use_iter
    logical :: loc_iter_default
    ! -------------------------------------------------------------------- !

    if ( present(key) ) then
      loc_key = trim(key)
    else
      loc_key = 'timeformat'
    end if

    if ( present(use_iter_default) ) then
      loc_iter_default = use_iter_default
    else
      loc_iter_default = .false.
    end if

    call aot_table_open(L       = conf,    &
      &                 parent  = parent,  &
      &                 thandle = thandle, &
      &                 key     = loc_key  )

    use_iter = loc_iter_default

    if (thandle /= 0) then
      ! The timeformat is given as a table load its components accordingly.
      call aot_get_val(L       = conf,             &
        &              thandle = thandle,          &
        &              val     = use_iter,         &
        &              key     = 'use_iter',       &
        &              default = loc_iter_default, &
        &              ErrCode = iErr              )

      call aot_get_val(L       = conf,         &
        &              thandle = thandle,      &
        &              val     = timeform,     &
        &              key     = 'simform',    &
        &              default = default_form, &
        &              ErrCode = iErr          )

    else
      ! The time is not given as a table, try to interpret it as a setting for
      ! the simform.
      call aot_get_val(L       = conf,         &
        &              thandle = parent,       &
        &              key     = loc_key,      &
        &              val     = timeform,     &
        &              default = default_form, &
        &              ErrCode = iErr          )
    end if

    call aot_table_close(conf, thandle)

    if (use_iter) then
      me = tem_timeformatter_init(timeform = '(I0)',                      &
        &                         stamp    = tem_timeformatter_iter_stamp )
    else
      me = tem_timeformatter_init(timeform = timeform,                   &
        &                         stamp    = tem_timeformatter_sim_stamp )
    end if

  end subroutine tem_timeformatter_load
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Generate a time stamp from the simulation time in the given time
  !! definition.
  !!
  !! This basically generates a string identifying the solution time and
  !! writing it in a meaningful format, so that it can be easily recognized.
  function tem_timeformatter_sim_stamp(formatter, time) result(timeStamp)
    ! -------------------------------------------------------------------- !
    !> Formatting object to generate the time stamp
    class(tem_timeformatter_type), intent(in) :: formatter
    !> Time definition to create the stamp off.
    type(tem_time_type), intent(in) :: time

    !> String representation of the given simulation time.
    character(len=labelLen) :: timeStamp
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    write(timeStamp, formatter%timeform) time%sim

    ! remove leading empty spaces in the timestamp
    timeStamp = adjustl(timeStamp)

  end function  tem_timeformatter_sim_stamp
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Generate a time stamp from the iteration in the given time
  !! definition.
  function tem_timeformatter_iter_stamp(formatter, time) result(timeStamp)
    ! -------------------------------------------------------------------- !
    !> Formatting object to generate the time stamp
    class(tem_timeformatter_type), intent(in) :: formatter
    !> Time definition to create the stamp off.
    type(tem_time_type), intent(in) :: time

    !> String representation of the given simulation time.
    character(len=labelLen) :: timeStamp
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    write(timeStamp, formatter%timeform) time%iter

    ! remove leading empty spaces in the timestamp
    timeStamp = adjustl(timeStamp)

  end function  tem_timeformatter_iter_stamp
  ! ************************************************************************ !

end module tem_timeformatter_module
