! Copyright (c) 2013-2014, 2019-2021 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2014 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2015, 2017, 2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
! **************************************************************************** !
!> This module provides the definition of various abort criteria upon which
!! a simulation should be stopped.
!! Note that solvers may extend this table and have their own set of
!! addititonal abort parameters to be set in this table.
!!
!! There are two primary options that may be set:
!!
!! * `stop_file` this denotes a file that will be checked for existence in the
!!   current working directory. If it exists the simulation will be stopped (if
!!   restarts are active, a restart will be written).
!!   If the file is empty the file will be deleted upon completion of the
!!   application. Otherwise, if the file is not empty it will be kept around.
!! * `steady_state` this is a boolean to indicate whether the simulation is to
!!   halt when a steady state is reached. If this is true a convergence table,
!!   where the condition for a steady state may be defined, is read. See
!!   [[tem_convergence_module]] for details.
!!
!! If the abort criteria table is not provided, the defaults for both settings
!! will be used, which is equivalent to the following definition:
!!
!!```lua
!!  abort_criteria = {stop_file = '', steady_state = false}
!!```
!!
!! That is, no stop files to look for, and no check whether a steady state is
!! reached.
!!
!! A more complete example with a check for steady state and a convergence
!! table (more details in [[tem_convergence_module]] could look like this:
!!
!!```lua
!!  abort_criteria = {
!!    stop_file = 'stop',
!!    steady_state = true,
!!    convergence = {
!!      variable = {'pressure', 'velocity'},
!!      shape = {kind = 'all'},
!!      time_control = {
!!        min = {iter = 0},
!!        max = {iter = tmax},
!!        interval = {iter = 10}
!!      },
!!      reduction = { 'average', 'average' },
!!      norm = 'average',
!!      nvals = 100,
!!      absolute = true,
!!      condition = {
!!         { threshold = 1.e-15, operator = '<=' },
!!         { threshold = 1.e-12, operator = '<=' }
!!      }
!!    }
!!  }
!!```
!!
!! This results in the application to look for a file named `stop` in the
!! working directory and will abort the run if it is found to exist.
!! The application will check for a steady state, where a steady state is
!! assumed when the average of pressure and velocity across the complete
!! domain does not deviate by more of 1.e-15 for pressure and not more than
!! 1.e-12 for the velocity.
!!
module tem_abortCriteria_module
  use env_module,             only: rk, labelLen, newunit
  use tem_convergence_module, only: tem_convergence_type, tem_convergence_load

  use aotus_module, only: flu_State, aot_get_val
  use aot_table_module, only: aot_table_open, &
    &                         aot_table_close

  use aot_out_general_module, only: aot_out_type,        &
    &                               aot_out_close_table, &
    &                               aot_out_open_table

  use aot_out_module, only: aot_out_val

  implicit none

  private

  public :: tem_abortCriteria_type
  public :: tem_abortCriteria_new
  public :: tem_abortCriteria_load
  public :: tem_abortCriteria_out
  public :: tem_abortCriteria_dump
  public :: tem_stop_file_exists

  !> Abstract type to describe solver specific abort criteria.
  !!
  !! Solvers may extend this type and pass it to load additional,
  !! solver specific criteria.
  type, abstract, public :: tem_solverAborts_type
  contains
    procedure(load_aborts), deferred :: load
  end type tem_solverAborts_type

  abstract interface

    !> Loading additional parameters for solver specific abort criteria
    !! from the configuration file.
    subroutine load_aborts(me, conf, abort_table)
      use aotus_module, only: flu_State
      import tem_solverAborts_type

      !> The solver specific type to hold additional abort parameters.
      class(tem_solverAborts_type), intent(inout) :: me

      !> Handle to Lua configuration file to load parameters from.
      type(flu_state), intent(in) :: conf

      !> Handle to the abort criteria table to read from.
      integer, intent(in) :: abort_table
    end subroutine load_aborts

  end interface

  !> Definition of the various abort criteria.
  !!
  !! Currently we only have two in addition to the time controlled and
  !! erroneous aborts.
  !! Solvers may pass an additional type to load extra parameters for
  !! aborts from the abortcriteria table.
  type tem_abortCriteria_type
    !> A file which should cause the simulation to stop.
    !! Default: ''.
    !!
    !! If this is a non empty string, the solver will stop at the next
    !! opportunity if it detects a file with the name provided here in the
    !! current working directory. Thus the simulation could be stopped by
    !! doing a "touch stop" in the working directory of the application.
    !! Such an empty file will be deleted, after it is detected.
    !! If you want the stop file to stay on the file system, there has to
    !! be something in it, which can be achieved by "echo keep > stop"
    !! for example.
    character(len=labelLen) :: stop_file

    !> Should the simulation be checked for a steady state convergence and
    !! stop if it is detected? Default: .false.
    logical :: steady_state

    !> Convergence conditions for steady state check.
    !! Filled only when steady_state is True
    type(tem_convergence_type), allocatable :: convergence(:)

  end type tem_abortCriteria_type


contains


  ! ************************************************************************ !
  !> Define new abortCriteria.
  !!
  !! A new abortCriteria object will be filled according to the parameters
  !! passed into the function.
  function tem_abortCriteria_new(stop_file, steady_state) result(ac)
    ! -------------------------------------------------------------------- !
    !> Name of the stop file to react on. Default=''.
    !!
    !! Any non-empty string activates this criterion.
    character(len=*), optional, intent(in) :: stop_file

    !> Flag to indicate if the simulation should stop upon reaching a steady
    !! state. What a steady state exactly is has to be defined in the solver.
    !! Default: .false.
    logical, optional, intent(in) :: steady_state

    !> A new variable of abortCriteria filled with the values provided as
    !! arguments.
    type(tem_abortCriteria_type) :: ac
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    if (present(stop_file)) then
      ac%stop_file = trim(stop_file)
    else
      ac%stop_file = ''
    end if

    if (present(steady_state)) then
      ac%steady_state = steady_state
    else
      ac%steady_state = .false.
    end if

    allocate(ac%convergence(0))

  end function tem_abortCriteria_new
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Load the abortCriteria from a given configuration.
  !!
  !! The abort_critera are defined in a table as follows:
  !!
  !!```lua
  !! abort_criteria = {stop_file='', steady_state=false}
  !!```
  !!
  !! stop_file indicates, which file should be checked for
  !! to stop the execution. A typical setting would for example
  !! be stop_file = 'stop'. If the string is empty, no checks
  !! are performed. The default is an empty string.
  !! Empty stop files will be deleted after they are encountered.
  !! Non-empty ones are kept.
  !!
  !! steady_state indicates if the simulation should stop when
  !! a steady_state solution is found. Default ist false.
  !!
  !! If steady_state is True then load convergence table for condition to
  !! check for steady state
  !!
  !! Solvers may pass an additional solverAborts for specific abort parameters
  !! to be filled from the abort_criteria table.
  subroutine tem_abortCriteria_load(me, conf, parent, key, solverAborts)
    ! -------------------------------------------------------------------- !
    !> Abort criteria to load from the Lua table.
    type(tem_abortCriteria_type), intent(out) :: me

    !> Handle for the Lua script.
    type(flu_state) :: conf

    !> Parent table to read from.
    integer, intent(in), optional :: parent

    !> Name of the time control table. Default: 'time_control'
    character(len=*), intent(in), optional :: key

    !> Solver specific abort criteria to load.
    class(tem_solverAborts_type), intent(inout), optional :: solverAborts
    ! -------------------------------------------------------------------- !
    character(len=labelLen) :: loc_key
    integer :: thandle
    integer :: iErr
    ! -------------------------------------------------------------------- !

    loc_key = 'abort_criteria'
    if (present(key)) loc_key = key

    call aot_table_open( L       = conf,    &
      &                  parent  = parent,  &
      &                  thandle = thandle, &
      &                  key     = loc_key  )

    if (thandle /= 0) then

      call aot_get_val( L       = conf,         &
        &               thandle = thandle,      &
        &               val     = me%stop_file, &
        &               key     = 'stop_file',  &
        &               default = '',           &
        &               ErrCode = iErr          )

      call aot_get_val( L       = conf,            &
        &               thandle = thandle,         &
        &               val     = me%steady_state, &
        &               key     = 'steady_state',  &
        &               default = .false.,         &
        &               ErrCode = iErr             )

      if (me%steady_state) then
        call tem_convergence_load( me           = me%convergence, &
          &                        conf         = conf,           &
          &                        parent       = thandle,        &
          &                        steady_state = me%steady_state )
      end if

      if (present(solverAborts)) then
        call solverAborts%load( conf        = conf,   &
          &                     abort_table = thandle )
      end if

    else

      me%stop_file = ''
      me%steady_state = .false.

    end if

    ! If no steady state is defined
    if (.not. me%steady_state) allocate(me%convergence(0))

    call aot_table_close( L       = conf,   &
      &                   thandle = thandle )

  end subroutine tem_abortCriteria_load
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Saves the abortCriteria to a given configuration.
  !!
  !! For further information, see TEM_abortCriteria_load()
  subroutine tem_abortCriteria_out(me, conf, key)
    ! -------------------------------------------------------------------- !
    !> The abortCriteria to write out as a Lua table.
    type(tem_abortCriteria_type), intent(in) :: me

    !> Handle for the Lua script to write to.
    type(aot_out_type), intent(inout) :: conf

    !> A name for the table to write the abortCriteria to.
    !! Default: 'abort_criteria'.
    character(len=*), optional :: key
    ! -------------------------------------------------------------------- !
    character(len=labelLen) :: loc_key
    ! -------------------------------------------------------------------- !

    loc_key = 'abort_criteria'
    if (present(key)) loc_key = key

    call aot_out_open_table( put_conf = conf,   &
      &                      tname    = loc_key )

    call aot_out_val( put_conf = conf,   &
      &              val = me%stop_file, &
      &              vname = 'stop_file' )

    call aot_out_val( put_conf = conf,       &
      &               val = me%steady_state, &
      &               vname = 'steady_state' )

    call aot_out_close_table(put_conf = conf)

  end subroutine tem_abortCriteria_out
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Dump abort criteria information to the specified outUnit.
  subroutine tem_abortCriteria_dump(me, outUnit)
    ! -------------------------------------------------------------------- !
    !> Abort criteria settings to write on outUnit.
    type(tem_abortCriteria_type), intent(in) :: me

    !> File unit to write the settings to.
    integer, intent(in) :: outUnit
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    write(outUnit,*) ' stop_file: '//trim(me%stop_file)
    write(outUnit,*) ' steady_state: ', me%steady_state

  end subroutine tem_abortCriteria_dump
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Check if the stop file exists.
  !!
  !! The check is only done by the root process and only if the stop_file
  !! setting is not empty.
  !! If the stop file exists, but is empty it is deleted after probing its
  !! existence. Non-empty files are kept.
  !! Thus, you can create a stop file that is to be deleted upon
  !! encountering with: touch stop
  !! While one, that should be kept can be created by:
  !! echo keep > stop
  function tem_stop_file_exists(abortCriteria, rank) result(sf_exists)
    ! -------------------------------------------------------------------- !
    !> Abort criteria settings to use in this check for a stop file.
    type(tem_abortCriteria_type), intent(in) :: abortCriteria

    !> Rank of the probing process, only rank==0 actually checks for the file.
    integer, intent(in) :: rank

    !> Result that indicates, if the stop files exists.
    logical :: sf_exists
    ! -------------------------------------------------------------------- !
    integer :: fu
    integer :: ios
    character(len=labelLen) :: probe
    ! -------------------------------------------------------------------- !

    sf_exists = .false.

    if (trim(abortCriteria%stop_file) /= '') then
      if (rank == 0) then
        inquire( file  = trim(abortCriteria%stop_file), &
          &      exist = sf_exists                      )

        if (sf_exists) then
          fu = newunit()
          open( unit   = fu,                            &
            &   file   = trim(abortCriteria%stop_File), &
            &   status = 'old'                          )

          read(fu,'(a)', iostat=ios) probe
          if (ios < 0) then
            close( unit   = fu,      &
              &    status = 'DELETE' )
          else
            close(fu)
          end if

        end if
      end if
    end if

  end function tem_stop_file_exists
  ! ************************************************************************ !

end module tem_abortCriteria_module
