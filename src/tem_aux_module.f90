! Copyright (c) 2011-2016,2021 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2012 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2014 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2015, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
! ---------------------------------------------------------------------------- !
!> Some auxilary routines, providing
!! frequently needed common tasks.
module tem_aux_module

  ! include treelm modules
  use mpi
  use env_module,              only: rk, labelLen, pathLen, my_status_int, &
    &                                env_nu => newunit
  use tem_comm_env_module,     only: tem_comm_env_type
  use tem_lua_requires_module, only: tem_require_track_rq,  &
    &                                tem_get_required_lua,  &
    &                                tem_pop_from_track_rq, &
    &                                tem_require_rq_store,  &
    &                                tem_push_to_rq_store
  use tem_tools_module,        only: tem_horizontalSpacer, upper_to_lower
  use tem_logging_module,      only: logUnit
  use soi_revision_module ! Providing parameters on the compilation environment

  ! include aotus modules
  use flu_binding,  only: flu_State, cbuf_type, flu_free_cbuf
  use aotus_module, only: open_config_buffer, open_config_file, aot_get_val, &
    &                     aoterr_Fatal, aoterr_NonExistent, aoterr_WrongType
  use aot_table_module, only: aot_table_length, aot_table_open, aot_table_close

  implicit none

  private

  public :: tem_open_distconf
  public :: tem_open_distconf_array
  public :: tem_open
  public :: tem_abort
  public :: tem_unit_close
  public :: tem_checkLabel
  public :: tem_print_execInfo
  public :: tem_global_vmhwm
  public :: utc_date_string
  public :: check_mpi_error
  public :: check_aot_error

  interface check_aot_error
    module procedure check_aot_error_scalar
    module procedure check_aot_error_vector
  end interface check_aot_error


contains


  ! ------------------------------------------------------------------------ !
  !> Read a Lua file on the first process and distribute it to all.
  !!
  !! @todo HK: Maybe deprecate and remove this routine in favor of
  !! TEM_open_distconf_array to avoid code duplication? Or keep it around and
  !! put a generic interface in place?
  !!
  !! This is a drop in replacement for open_config_file from Aotus and allows
  !! the scalable processing of Lua files, as they are read by a single process
  !! and then streamed to all in proc.
  !! There should be no restrictions on the Lua scripts themselves in this
  !! method, as it uses an overloading of the require mechanism in Lua itself to
  !! replace the file searches by lookups of buffered Lua code snippets.
  !! The execution of the Lua script itself is not changed.
  !!
  subroutine tem_open_distconf(L, fileName, proc)
    ! -------------------------------------------------------------------- !
    type(flu_State) :: L !< Opened Lua state with the loaded script.
    character(len=*), intent(in) :: fileName !< Name of the file to open.
    !> Process description to use.
    type(tem_comm_env_type), intent(in) :: proc
    ! -------------------------------------------------------------------- !
    type(cbuf_type) :: scriptbuf
    type(cbuf_type) :: modbuf
    integer :: iError
    integer :: comm
    integer :: nProcs
    integer :: nFiles
    integer :: bufsize
    integer :: iFile
    character(len=labelLen), allocatable :: req_label(:)
    character(len=pathLen), allocatable  :: req_file(:)
    ! -------------------------------------------------------------------- !

    comm = proc%comm
    nProcs = proc%comm_size

    if (nProcs > 1) then

      if ( proc%isRoot ) then
        ! Only rank 0 reads and executes the config file, while doing so, it
        ! keeps track of all the required files.
        call tem_require_track_rq(L)
        call open_config_file(L, trim(filename), buffer=scriptbuf)
        call tem_get_required_Lua( L, fileList = req_file, &
          &                        labelList = req_label   )
        nFiles = size(req_Label)
      else
        ! Open the configuration on all other processes with the rq_store
        ! loaded:
        call tem_require_rq_store(L)
        ! this is necessary to suppress valgrind debug output
        nFiles = 0
      end if

      ! Broadcast the number of required files
      ! MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
      call MPI_Bcast(nFiles, 1, MPI_INTEGER, proc%root, comm, iError)

      if ( .not. proc%isRoot ) allocate(req_label(nFiles))
      if ( .not. proc%isRoot ) allocate(req_file(nFiles))

      ! Broadcast the module names
      call MPI_Bcast( req_label, nFiles*labelLen, MPI_CHARACTER, proc%root, &
        &             comm, iError                                          )
      ! Broadcast the file names
      call MPI_Bcast( req_file, nFiles*pathLen, MPI_CHARACTER, proc%root, &
        &             comm, iError                                        )

      ! Now go on opening all required files
      do iFile=1,nFiles
        if ( proc%isRoot ) then
          call tem_pop_from_track_rq(L, trim(req_label(iFile)), modbuf)
          bufsize = size(modbuf%buffer)
        end if
        ! Broadcast the loaded script to all processes.
        ! MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
        call MPI_Bcast(bufsize, 1, MPI_INTEGER, proc%root, comm, iError)
        if ( .not. proc%isRoot ) allocate(modbuf%buffer(bufsize))
        call MPI_Bcast( modbuf%buffer, bufsize, MPI_CHARACTER, proc%root, &
          &             comm, iError                                      )
        if ( .not. proc%isRoot ) then
          call tem_push_to_rq_store( L,                                 &
            &                        modname  = trim(req_label(iFile)), &
            &                        filename = trim(req_file(iFile)),  &
            &                        buffer   = modbuf%buffer           )
          deallocate(modbuf%buffer)
        else
          call flu_free_cbuf(modbuf)
        end if
      end do

      ! Broadcast the loaded script to all processes.
      if ( proc%isRoot ) bufsize = size(scriptbuf%buffer)
      ! MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
      call MPI_Bcast(bufsize, 1, MPI_INTEGER, proc%root, comm, iError)
      if ( .not. proc%isRoot ) allocate(scriptbuf%buffer(bufsize))
      call MPI_Bcast( scriptbuf%buffer, bufsize, MPI_CHARACTER, proc%root, &
        &             comm, iError                                         )

      if ( .not. proc%isRoot ) then
        call open_config_buffer(L = L, buffer = scriptbuf%buffer)
        deallocate(scriptbuf%buffer)
      else
        call flu_free_cbuf(scriptbuf)
      end if

    else

      ! Only a single process, no need for broadcasting.
      call open_config_file(L = L, filename = trim(fileName))

    end if

  end subroutine tem_open_distconf
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Open an array of Lua handles.
  !!
  !! This is a drop in replacement for open_config_file from Aotus and allows
  !! the scalable processing of Lua files, as they are read by a single process
  !! and then streamed to all in proc.
  !! There should be no restrictions on the Lua scripts themselves in this
  !! method, as it uses an overloading of the require mechanism in Lua itself to
  !! replace the file searches by lookups of buffered Lua code snippets.
  !! The execution of the Lua script itself is not changed.
  !!
  !! This variant of the routine opens an array of handles with the same script.
  !! It is used to provide individual Lua states, that can be used independently
  !! later on, for example in the evaluation of Lua functions for boundary
  !! conditions.
  !!
  subroutine tem_open_distconf_array(L, fileName, proc)
    ! -------------------------------------------------------------------- !
    type(flu_State) :: L(:) !< Opened Lua state with the loaded script.
    character(len=*), intent(in) :: fileName !< Name of the file to open.
    type(tem_comm_env_type) :: proc !< Process description to use.
    ! -------------------------------------------------------------------- !
    type(cbuf_type) :: scriptbuf
    type(cbuf_type) :: modbuf
    integer :: iError
    integer :: comm
    integer :: nProcs
    integer :: nFiles
    integer :: nLuaStates
    integer :: bufsize
    integer :: iFile
    integer :: iState
    integer :: Lua_lb
    character(len=labelLen), allocatable :: req_label(:)
    character(len=pathLen), allocatable  :: req_file(:)
    ! -------------------------------------------------------------------- !

    nLuaStates = size(L)

    comm = proc%comm
    nProcs = proc%comm_size

    if ( proc%isRoot ) then
      Lua_lb = 2
    else
      Lua_lb = 1
    end if

    if (nProcs > 1) then
      nFiles = 0

      if ( proc%isRoot ) then
        ! Only rank 0 reads and executes the config file for the first Lua
        ! state, while doing so, it keeps track of all the required files.
        call tem_require_track_rq(L(1))
        call open_config_file(L(1), trim(filename), buffer=scriptbuf)
        call tem_get_required_Lua( L(1), fileList = req_file, &
          &                        labelList = req_label      )
        nFiles = size(req_Label)
      end if
      ! Open the configuration on all other processes and the remaining local
      ! Lua states with the rq_store loaded:
      do iState=Lua_lb,nLuaStates
        call tem_require_rq_store(L(iState))
      end do

      ! Broadcast the number of required files
      ! MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
      call MPI_Bcast(nFiles, 1, MPI_INTEGER, 0, comm, iError)

      if ( .not. proc%isRoot ) allocate(req_label(nFiles))
      if ( .not. proc%isRoot ) allocate(req_file(nFiles))

      ! Broadcast the module names
      call MPI_Bcast( req_label, nFiles*labelLen, MPI_CHARACTER, 0, &
        &             comm, iError                                  )
      ! Broadcast the file names
      call MPI_Bcast( req_file, nFiles*pathLen, MPI_CHARACTER, 0, &
        &             comm, iError                                )

      ! Now go on opening all required files
      do iFile=1,nFiles
        if (proc%isRoot) then
          call tem_pop_from_track_rq(L(1), trim(req_label(iFile)), modbuf)
          bufsize = size(modbuf%buffer)
        end if
        ! Broadcast the loaded script to all processes.
        ! MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
        call MPI_Bcast(bufsize, 1, MPI_INTEGER, 0, comm, iError)
        if (.not. proc%isRoot) allocate(modbuf%buffer(bufsize))
        call MPI_Bcast(modbuf%buffer, bufsize, MPI_CHARACTER, 0, comm, iError)
        do iState=Lua_lb,nLuaStates
          call tem_push_to_rq_store( L(iState),                         &
            &                        modname  = trim(req_label(iFile)), &
            &                        filename = trim(req_file(iFile)),  &
            &                        buffer   = modbuf%buffer           )
        end do
        if (proc%isRoot) then
          call flu_free_cbuf(modbuf)
        else
          deallocate(modbuf%buffer)
        end if
      end do

      ! Broadcast the loaded script to all processes.
      if ( proc%isRoot ) bufsize = size(scriptbuf%buffer)
      ! MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
      call MPI_Bcast(bufsize, 1, MPI_INTEGER, 0, comm, iError)
      if ( .not. proc%isRoot ) allocate(scriptbuf%buffer(bufsize))
      call MPI_Bcast(scriptbuf%buffer, bufsize, MPI_CHARACTER, 0, comm, iError)

      do iState=Lua_lb,nLuaStates
        call open_config_buffer(L = L(iState), buffer = scriptbuf%buffer)
      end do
      if ( proc%isRoot ) then
        call flu_free_cbuf(scriptbuf)
      else
        deallocate(scriptbuf%buffer)
      end if

    else

      ! Only a single process, no need for broadcasting.
      do iState=1,nLuaStates
        call open_config_file(L = L(iState), filename = trim(fileName))
      end do

    end if

  end subroutine tem_open_distconf_array
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Wrapper around Fortran open of files to take care of errors and improve
  !! the error message in case the opening goes wrong.
  !!
  !! Use newunit to let tem_open provide a new file unit for the opened file.
  subroutine tem_open(file, unit, newunit, status, position, action, form, &
    &                 access, recl)
    character(len=*), intent(in) :: file
    character(len=*), intent(in), optional :: status
    character(len=*), intent(in), optional :: position
    character(len=*), intent(in), optional :: action
    character(len=*), intent(in), optional :: form
    character(len=*), intent(in), optional :: access
    integer, intent(in), optional :: recl
    integer, intent(in), optional :: unit
    integer, intent(out), optional :: newunit
    ! -------------------------------------------------------------------- !
    character(len=labelLen) :: loc_status
    character(len=labelLen) :: loc_position
    character(len=labelLen) :: loc_action
    character(len=labelLen) :: loc_form
    character(len=labelLen) :: loc_access
    integer :: stat
    integer :: funit
    ! -------------------------------------------------------------------- !

    ! Defaults:
    loc_status   = 'unknown'
    loc_position = 'asis'
    loc_action   = 'readwrite'
    loc_form     = 'formatted'
    loc_access   = 'sequential'

    if (present(status)) loc_status = upper_to_lower(status)
    if (present(position)) loc_position = upper_to_lower(position)
    if (present(action)) loc_action = upper_to_lower(action)
    if (present(access)) loc_access = upper_to_lower(access)

    ! Stream IO is by default unformatted.
    if (loc_access == 'stream') loc_form = 'unformatted'

    if (present(form)) loc_form = upper_to_lower(form)

    if (present(unit)) then
      funit = unit
    else
      funit = env_nu()
      if (present(newunit)) then
        newunit = funit
      end if
    end if

    rl_provided: if (present(recl)) then

      pos_provided: if (present(position)) then
        open( unit     = funit,            &
          &   file     = file,             &
          &   action   = trim(loc_action), &
          &   access   = loc_access,       &
          &   status   = loc_status,       &
          &   position = loc_position,     &
          &   form     = loc_form,         &
          &   recl     = recl,             &
          &   iostat   = stat              )
      else pos_provided
        open( unit     = funit,            &
          &   file     = file,             &
          &   action   = trim(loc_action), &
          &   access   = loc_access,       &
          &   status   = loc_status,       &
          &   form     = loc_form,         &
          &   recl     = recl,             &
          &   iostat   = stat              )
      end if pos_provided

    else rl_provided

      seqpos: if ( (loc_access == 'sequential') .and. present(position)) then
        open( unit     = funit,            &
          &   file     = file,             &
          &   action   = trim(loc_action), &
          &   access   = loc_access,       &
          &   status   = loc_status,       &
          &   position = loc_position,     &
          &   form     = loc_form,         &
          &   iostat   = stat              )
      else seqpos
        open( unit     = funit,            &
          &   file     = file,             &
          &   action   = trim(loc_action), &
          &   access   = loc_access,       &
          &   status   = loc_status,       &
          &   form     = loc_form,         &
          &   iostat   = stat              )
      end if seqpos

    end if rl_provided

    if (stat /= 0) then
      write(logUnit(1), *) 'Could not open file!'
      write(logUnit(1), *) 'iostat=', stat
      write(logUnit(1), *) 'File:     ' // trim(file)
      if (present(action))   write(logUnit(1), *) 'Action:   ' // trim(action)
      if (present(form))     write(logUnit(1), *) 'Form:     ' // trim(form)
      if (present(access))   write(logUnit(1), *) 'Access:   ' // trim(access)
      if (present(status))   write(logUnit(1), *) 'Status:   ' // trim(status)
      if (present(recl))     write(logUnit(1), *) 'Recl:     ', recl
      if (present(position)) write(logUnit(1), *) 'Position: ' // trim(position)
      write(logUnit(1), *) 'Aborting...'
      call tem_abort()
    end if

  end subroutine tem_open
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Abort the program with finalization of the environment
  !!
  subroutine tem_abort( errorMsg )
    ! -------------------------------------------------------------------- !
    !> An optional error message to print a reason for the abort.
    character(len=*), intent(in), optional :: errorMsg
    ! -------------------------------------------------------------------- !
    integer :: iError
    ! -------------------------------------------------------------------- !

    ! @todo JZ: commented out the the tem_finalize here: In case that one rank
    ! while the other ranks are still waiting for communication the solver
    ! will not terminate.
    !call tem_finalize()
    if( present( errorMsg ) ) write(logUnit(1),*) errorMsg
    write(logUnit(1),*)
    write(logUnit(1),*) " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(logUnit(1),*) "               Aborting. "
    write(logUnit(1),*) " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(logUnit(1),*)
    flush(logUnit(1))
    call mpi_abort(MPI_COMM_WORLD, 1, iError)

  end subroutine tem_abort
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Check, if a unit is open, and close it
  !!
  subroutine tem_unit_close(me)
    ! -------------------------------------------------------------------- !
    !> the restart type to close
    integer :: me
    ! -------------------------------------------------------------------- !
    logical :: nUnitOpened
    ! -------------------------------------------------------------------- !
    ! Check, if any open units have to be closed
    if ( me >= 0 ) then
      ! unit has be to be >= 0
      inquire(unit=me, opened=nUnitOpened)
      if (nUnitOpened) then
        close( me )
      end if
    end if

  end subroutine tem_unit_close
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> check whether the current label is already given
  !!
  subroutine tem_checkLabel(label, nLabels)
    ! -------------------------------------------------------------------- !
    !> holding array of labels with label(n) contains current label
    character(len=*), intent(in) :: label(:)
    !> Number of schemes already existing in the scheme array
    !! (they are being added currently, so we only have to compare against
    !! the ones coming before the current one, up to nSchemes-1)
    integer, intent(in) :: nLabels
    ! -------------------------------------------------------------------- !
    integer :: iLabel
    ! -------------------------------------------------------------------- !

    write(logUnit(1),*)' Comparing the labels...'
    do ilabel = 1, nLabels - 1
      ! Check the labels here
      if ( trim(label(ilabel)) == trim(label(nLabels)) ) then
        write(logUnit(1),*) 'Error: identical label have been encountered.'
        write(logUnit(1),*) '       Please specify a unique label for ' &
          &                 // 'multiple tables in the config file.'
        write(logUnit(1),*) "Example: scheme = {{ label = 'scheme1' , ...  }, "
        write(logUnit(1),*) "                   { label = 'scheme2' , ...  }} "
        write(logUnit(1),*) "                                                 "
        write(logUnit(1),*) "Stopping...       "
        call tem_abort()
       end if
    end do

  end subroutine tem_checkLabel
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Print information of the executable
  !!
  subroutine tem_print_execInfo()
    ! -------------------------------------------------------------------- !
    integer :: flagline
    ! -------------------------------------------------------------------- !

    call tem_horizontalSpacer(fUnit = logUnit(1))
    write(logUnit(1),*) '| INFORMATION ON THE EXECUTABLE'
    write(logUnit(1),*) '| Revision of the code in this executable: ' &
      &                 // trim(soi_solver_revision)
    write(logUnit(1),*) '|'
    write(logUnit(1),*) '| Compiled with '//trim(soi_FC_name) &
      &                 // ' in version ' // trim(soi_FC_version)
    write(logUnit(1),*) '| Using the command ' // trim(soi_FC_command)
    write(logUnit(1),*) '| And the following flags:'
    do flagline = 1, soi_FC_nFlagLines
      write(logUnit(1),*) '| ' // trim(soi_FC_flags(flagline))
    end do
    write(logUnit(1),*) '|'
    write(logUnit(1),*) '| Build date: ' // trim(soi_build_date)
    call tem_horizontalSpacer(fUnit = logUnit(1))

  end subroutine tem_print_execInfo
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Obtain the memory status from all processes (min, max, avg).
  !!
  !! Find min, max and average high water mark of the virtual memory usage
  !! across all processes (MPI_COMM_WORLD) on rank 0.
  !! Results are in Megabytes, and the resulting array contains min, max, avg
  !! in this order.
  function tem_global_vmhwm() result(hwm)
    real(kind=rk) :: hwm(3)
    ! -------------------------------------------------------------------- !
    integer :: myhwm, minhwm, maxhwm
    integer :: nProcs
    integer :: iError
    real :: myMB
    real :: sumhwm
    ! -------------------------------------------------------------------- !

    call MPI_Comm_Size(MPI_COMM_WORLD, nProcs, iError)
    myhwm = my_status_int('VmHWM:')
    call MPI_Reduce( myhwm, minhwm, 1, MPI_INTEGER, MPI_MIN, 0, &
      &              MPI_COMM_WORLD, iError                     )
    call MPI_Reduce( myhwm, maxhwm, 1, MPI_INTEGER, MPI_MAX, 0, &
      &              MPI_COMM_WORLD, iError                     )
    myMB = real(myhwm)/1024.0
    call MPI_Reduce( myMB, sumhwm, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, &
      &              iError                                                 )

    hwm(1) = real(minhwm, kind=rk)/1024.0_rk
    hwm(2) = real(maxhwm, kind=rk)/1024.0_rk
    hwm(3) = sumhwm / real(nProcs, kind=rk)

  end function tem_global_vmhwm
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Function to create a human readable UTC date string.
  !!
  !! The resulting string has 26 characters.
  !!
  function utc_date_string() result(dat_string)
    ! -------------------------------------------------------------------- !
    character(len=26) :: dat_string
    ! -------------------------------------------------------------------- !
    integer, parameter :: year = 1
    integer, parameter :: month = 2
    integer, parameter :: day = 3
    integer, parameter :: hour = 5
    integer, parameter :: minute = 6
    integer, parameter :: utc_diff = 4
    character(len=9) :: u_off_string
    integer :: off_min, off_hour
    integer :: dat(8)
    ! -------------------------------------------------------------------- !

    call date_and_time(values=dat)
    off_min = mod(dat(utc_diff),60)
    off_hour = dat(utc_diff)/60
    if (dat(utc_diff) >= 0) then
      write(u_off_string,'(a4,i2.2,a1,i2.2)') 'UTC+', off_hour, ':', off_min
    else
      write(u_off_string,'(a3,i3.2,a1,i2.2)') 'UTC', off_hour, ':', off_min
    end if
    write(dat_string,'(i4,4(a1,i2.2),a10)') dat(year), '-', &
      &     dat(month), '-', dat(day), ' ', dat(hour), ':', &
      &     dat(minute), ' '//u_off_string

  end function utc_date_string
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine check_mpi_error( iError, event_string )
    integer, intent(in) :: iError
    character(len=*), intent(in) :: event_string

    character(len=100) :: IOError
    integer :: resultlen = 100
    integer :: ErrErr

    if (iError /= MPI_SUCCESS) then
      call MPI_ERROR_STRING(iError, IOError, resultlen, ErrErr)
      write(logUnit(0),*) 'MPI Error when '//trim(event_string),': ' &
        &                 //trim(IOError)
      call tem_abort()
    end if
  end subroutine check_mpi_error
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Auxiliary subroutine to check on errors from attempting to get values
  !! from the Lua script with aot_get_val.
  !!
  !! If a fatal error was encountered, the routine aborts the program!
  subroutine check_aot_error_scalar( iError, key, event_string )
    !> aoterr code to interpret (returned by aot_get_val)
    integer, intent(in) :: iError
    !> Lua key that was attempted to be read
    character(len=*), intent(in) :: key
    !> Optional event string to describe the circumstances
    character(len=*), intent(in), optional :: event_string

    if (btest(iError, aoterr_Fatal)) then
      if (present(event_string)) then
        write(logUnit(0),*) 'Lua Error while ' // trim(event_string) // '!'
      else
        write(logUnit(0),*) 'Encountered Lua Error for ' // trim(key) // '!'
      end if
      if (btest(iError, aoterr_NonExistent)) then
        write(logUnit(0),*) ' -> expected value "'//trim(key)//'" not found!'
      end if
      if (btest(iError, aoterr_WrongType)) then
        write(logUnit(0),*) ' -> value "'//trim(key)//'" has wrong type!'
      end if
      write(logUnit(0),*) 'Aborting... '
      call tem_abort()
    end if
  end subroutine check_aot_error_scalar
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Auxiliary subroutine to check on errors from attempting to get an array
  !! of values from the Lua script with aot_get_val.
  !!
  !! If a fatal error was encountered, the routine aborts the program!
  subroutine check_aot_error_vector( iError, key, event_string )
    !> aoterr code to interpret (returned by aot_get_val)
    integer, intent(in) :: iError(:)
    !> Lua key that was attempted to be read
    character(len=*), intent(in) :: key
    !> Optional event string to describe the circumstances
    character(len=*), intent(in), optional :: event_string

    integer :: iComp

    if (any(btest(iError, aoterr_Fatal))) then
      if (present(event_string)) then
        write(logUnit(0),*) 'Lua Error while ' // trim(event_string) // '!'
      else
        write(logUnit(0),*) 'Encountered Lua Error for ' // trim(key) // '!'
      end if
      do iComp=lbound(iError,1),ubound(iError,1)
        if (btest(iError(iComp), aoterr_NonExistent)) then
          write(logUnit(0),"(a,i0,a)") ' -> expected value "'//trim(key)//'[', &
            &                        iComp, ']" not found!'
        end if
        if (btest(iError(iComp), aoterr_WrongType)) then
          write(logUnit(0),"(a,i0,a)") ' -> value "'//trim(key)//'[', iComp, &
            &                        ']" has wrong type!'
        end if
      end do
      write(logUnit(0),*) 'Aborting... '
      call tem_abort()
    end if
  end subroutine check_aot_error_vector

end module tem_aux_module
! ---------------------------------------------------------------------------- !
