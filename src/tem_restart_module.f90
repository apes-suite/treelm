! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2012 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011-2012 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2011-2019, 2021 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2012, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012-2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2014 Julia Moos <julia.moos@student.uni-siegen.de>
! Copyright (c) 2014, 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2017 Raphael Haupt <Raphael.Haupt@student.uni-siegen.de>
! Copyright (c) 2019 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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
!> author: Simon Zimny, Harald Klimach
!! This module provides the main IO facilities to read and write large amounts
!! of data during the simulation.
!!
!! Especially, it enables the dumping of the simulation status
!! (e.g. the state vector) to disk and reading the information later on to
!! restart the simulation at a certain point of time.
!! Additionally, information like solver, version, number of elements, etc.
!! are stored in a seperate header file.
!!
!! The general procedure is as follows:
!!
!! - read the configuration from the lua file with tem_load_restart
!! - initialize the read/write restart with tem_init_restart
!! - read restart
!! - write restart
!!
!! Find more information about usage of restart in the section '[restart_usage]'
!!
!! The Restart / Harvester Format
!! ---
!!
!! @todo Add an extensive description of the restart format here
!!
!! The Restart Usage
!! ---
!!
!! To use the restart functionality the user has to define the io buffersize
!!
!!```lua
!! -- IO buffer size in MB (default = 8)
!! io_buffer_size = 1
!!```
!!
!! and the restart settings in the lua configuration file.
!!
!!```lua
!! -- Restart settings
!! restart = { read  = 'restart/gaussian_lastHeader.lua',
!!                     -- file to restart from
!!             write = 'restart/',
!!                     -- prefix to write the files to
!!             time_control = { min = 0, max = 10, interval = 10,
!!                              align_trigger = false,
!!             }
!!             -- timing definitions (either iterations or simulation time)
!!           }
!!```
!! Another option to read restart is from initial conditions.
!! The user then has to define restart data to be read into a variable in this
!! way:
!!```lua
!! initial_condition = { file = 'restart/channel_lastHeader.lua',
!!   depend = {
!!     variable = { {'pdf', 19} },
!!     usescheme = 'mini-channel'
!!   }
!! }
!!```
!!
!! If only the write option shall be used the identifier read has to be removed
!! and vice versa.
module tem_restart_module

  ! include treelm modules
  use mpi
  use env_module,              only: LabelLen, PathLen, rk_mpi, long_k, rk, &
    &                                tem_create_endianSuffix, pathSep,      &
    &                                io_buffer_size, globalMaxLevels,       &
    &                                long_k_mpi, newUnit
  use treelmesh_module,        only: treelmesh_type, load_tem, dump_treelmesh
  use tem_global_module,       only: tem_global_type, tem_mesh_out
  use tem_aux_module,          only: tem_abort, tem_open_distconf,          &
    &                                tem_open, check_mpi_error
  use tem_comm_env_module,     only: tem_comm_env_type
  use tem_solveHead_module,    only: tem_solveHead_type, tem_solverTag
  use tem_timeControl_module,  only: tem_timeControl_type, &
    &                                tem_timeControl_load, &
    &                                tem_timeControl_align_trigger, &
    &                                tem_timeControl_dump
  use tem_time_module,         only: tem_time_sim_stamp, tem_time_type, &
    &                                tem_time_out, tem_time_load,       &
    &                                tem_time_set_clock, tem_time_dump, &
    &                                tem_time_reset
  use tem_tools_module,        only: tem_horizontalSpacer
  use tem_varSys_module,       only: tem_varSys_type, tem_varSys_out,  &
    &                                tem_varSys_load, tem_varSys_dump, &
    &                                tem_get_element_chunk
  use tem_varMap_module,       only: tem_varMap_type
  use tem_subTree_type_module, only: tem_subTree_type, tem_dump_subTree
  use tem_logging_module,      only: logUnit
  use tem_debug_module,        only: dbgUnit
  ! use tem_transient_reduction_module,                               &
  !   &                          only: tem_transient_reduction_open,  &
  !   &                                tem_transient_reduction_apply, &
  !   &                                tem_transient_reduction_close, &
  !   &                                tem_transient_reduction_type

  ! include aotus modules
  use aotus_module,     only: flu_State, aot_get_val, close_config
  use aot_out_module,   only: aot_out_type, aot_out_val,            &
    &                         aot_out_open, aot_out_close,          &
    &                         aot_out_open_table, aot_out_close_table
  use aot_table_module, only: aot_table_open, aot_table_close

  implicit none

  private

  public :: tem_restart_type
  public :: tem_restartControl_type
  public :: tem_restartHeader_type
  public :: tem_load_restart
  public :: tem_init_restart
  public :: tem_init_restart_alloc
  public :: tem_init_restart_create_types
  public :: tem_restart_finalize
  public :: tem_restart_writeData
  public :: tem_restart_writeHeader
  public :: tem_restart_readData
  public :: tem_restart_readData_single
  public :: tem_restart_openRead
  public :: tem_restart_openRead_single
  public :: tem_restart_closeRead
  public :: tem_restart_closeRead_single
  public :: tem_restart_openWrite
  public :: tem_restart_closeWrite
  public :: tem_restart_getTotalChunks
  public :: tem_restart_dump_data

  !> Control the behavior of the restart, like at which point in time etc.
  type tem_restartControl_type
    !> Do a normal initialization if the read restart file is not found?
    logical :: init_on_missing
    !> is the restart read active?
    logical :: readRestart
    !> is the restart write active?
    logical :: writeRestart
    !> read the restart file from following file
    character(len=PathLen) :: readFileName
    !> write the restart file into a file with the following prefix
    character(len=PathLen) :: writePrefix
    !> control about when to do the restart
    type(tem_timeControl_type) :: timeControl
  end type tem_restartControl_type

  !> Define quantities like the prefix, the mesh and the timestamp
  type tem_restartHeader_type
    !> a unique tag for the solver including version
    character(len=LabelLen) :: solverTag = 'unknown'
    !> solver creating the input data (Ateles, Musubi)
    character(len=LabelLen) :: solverKind
    !> the simulation name
    character(len=LabelLen) :: simName = ''
    !> the solver config file name
    character(len=pathLen) :: solverConfigFile = ''
    !> name of the binary files (this variable is set for every timestamp
    !! differently, used for read and write)
    character(len=PathLen) :: binName
    !> prefix for the binary files (includes directory and file prefix)
    character(len=PathLen) :: binPrefix
    !> prefix for the header files (includes directory and file prefix)
    character(len=PathLen) :: headerPrefix
    !> the current time generated from tem_time_module
    character(len=16) :: timestamp
    !> mesh directory
    ! character(len=pathLen) :: meshDir
    !> Number of elements as read from the restart header
    integer :: nElems
    !> variable system dumped in restart header
    type(tem_varSys_type) :: varSys
  end type tem_restartHeader_type


  type tem_file_layout_type
    !> Number of degrees of freedom for each of the scalars entries in
    !! the varibale system in the file.
    integer :: nDofs = 1

    !> Number of local chunks required to fit complete data in.
    !! Set in routine: tem_restart_getTotalChunks
    integer :: nChunks

    !> Globally maximal number of chunks required.
    !! Set in routine: tem_restart_getTotalChunks
    integer :: maxnChunks

    !> Number of elements that fit into the buffer.
    !! Set in routine: tem_restart_getTotalChunks
    integer :: chunkSize

    !> Handle for the MPI type to describe the view for IO in binary IO
    integer :: ftype

    !> Handle for the MPI type describing the vector of data in each element
    integer :: vectype

    !> start of view for MPI_SET_VIEW
    integer(kind=MPI_OFFSET_KIND)     :: displacement
  end type

  !> The restart type defining everything related to the
  !! disk input/output
  type tem_restart_type
    !> communicator for the processes participating in this restart (might be
    !! only a subset of the global communicator)
    type(tem_comm_env_type) :: comm
    !> actual number of elements in the current chunk (= chunkSize or
    ! tree%nElems-(nChunks-1)*ChunkSize)
    integer :: nChunkElems

    !> Description of the data layout to use when reading a file.
    type(tem_file_layout_type) :: read_file

    !> Description of the data layout to use when writing a file.
    type(tem_file_layout_type) :: write_file

    !> Control the behavior of the restart, like at which point in time etc.
    type(tem_restartControl_type) :: controller
    !> Define quantities like the prefix, the mesh and the timestamp
    type(tem_restartHeader_type)  :: header
    !> unit integer to write binary data to
    integer :: binaryUnit
    !> name and position of variables in global variable system
    type(tem_varMap_type) :: varMap
    !!> number of scalars of variables in varPos
    integer :: nScalars
    !> scratch file unit contains solver specific info in dump in restart header
    !! This file should contain the information in form of a Lua script.
    integer :: solSpec_unit = -1
    !> The time when the last restart file was written.
    type(tem_time_type) :: lastWritten
  end type tem_restart_type


contains


  ! ************************************************************************ !
  !> Read all necessary information for the restart from the lua config file.
  !!
  !! Include this routine into your general configuration load routine.
  !! The configuration looks as follows
  !!```lua
  !! restart = { read = 'restart/lastHeader.lua', -- Which file to restart from,
  !!                                              -- if any
  !!             write = 'restart/', -- Where to write the restart files to,
  !!                                 -- if any
  !!             time = { min = 0, max = 10, interval = 10} -- when to output
  !!             }
  !!```
  !! Here, the restart is loaded from `restart/lastHeader.lua` and reads in the
  !! related data and configuration.
  !! Restart files are written out in `restart/` folder
  !!
  subroutine tem_load_restart( me, conf, tree, timing, globProc, parent_table, &
    &                          key )
    ! -------------------------------------------------------------------- !
    !> restart type to be filled
    type(tem_restart_type), intent(inout)  :: me
    !> lua configuration file
    type(flu_state)                        :: conf
    !> mesh, provided in treelm format
    type(treelmesh_type), intent(inout)    :: tree
    !> the timing for re-setting the times
    type(tem_time_type), intent(inout)     :: timing
    !> Global process communicator env
    type( tem_comm_env_type ), intent(in)  :: globProc
    !> optional parent handle
    integer, optional, intent(in)          :: parent_table
    !> optional key for table
    character(len=*), optional, intent(in) :: key
    ! -------------------------------------------------------------------- !
    character(len=32) :: localKey
    logical :: readexists
    integer :: restart_table
    integer :: iError
    ! -------------------------------------------------------------------- !

    if (present(key)) then
      ! The table to look for is not named restart, look for this different
      ! key.
      localKey = key
    else
      ! Use the default name restart for the table.
      localKey = 'restart'
      ! Set current folder as default prefix for writing
      me%controller%writePrefix = '.'//pathSep
    end if

    me%controller%readRestart  = .false.
    me%controller%writeRestart = .false.

    ! Attempt to open the restart table (within another table, if a parent is
    ! given).
    call aot_table_open( L       = conf,          &
      &                  thandle = restart_table, &
      &                  parent  = parent_table,  &
      &                  key     = trim(localKey) )

    ! Initialize the last written time to 0.
    call tem_time_reset(me%lastWritten)
    me%lastWritten%iter = -1

    ! If the restart table is present, the parameters are loaded.
    ! In case of dynamic load balancing, parameters are loaded
    ! in a different manner i.e. timings etc. are read from balance table
    if (restart_table .ne. 0 ) then
      call tem_horizontalSpacer(fUnit = logUnit(1))
      write(logUnit(1),*) 'Loading restart ...'
      ! Successfully opened the table.

      ! First we get all the informations of the read table.
      ! Reading the filename to read the restart data from.
      call aot_get_val( L       = conf,                       &
        &               thandle = restart_table,              &
        &               key     = 'read',                     &
        &               val     = me%controller%readFileName, &
        &               ErrCode = iError                      )

      if (iError == 0) then
        call aot_get_val( L       = conf,                          &
          &               thandle = restart_table,                 &
          &               key     = 'init_on_missing',             &
          &               val     = me%controller%init_on_missing, &
          &               default = .false.,                       &
          &               ErrCode = iError                         )

        ! Successfully obtained a filename to restart from, now go on and read
        ! the data from its header if it exists.
        write(logUnit(1),*) "*****************************"
        write(logUnit(1),*) "Restart read parameters: "
        write(logUnit(1),*) "  filename : "//trim(me%controller%readFileName)

        if (globProc%rank == 0) then
          inquire(file = trim(me%controller%readFileName), exist = readexists)
        end if
        call MPI_Bcast(readexists, 1, MPI_LOGICAL, 0, globProc%comm, iError)

        ! Set the restart flag if the restart file exists
        me%controller%readRestart = readexists

        if (readexists) then
          call tem_restart_readHeader( me          = me,       &
            &                          timing      = timing,   &
            &                          globProc    = globProc, &
            &                          tree        = tree      )

          if ( tree%global%nElems /= me%header%nElems ) then
            write(logUnit(0),*) 'Number of elements in restart header different ' &
              &                 // 'from mesh'
            write(logUnit(0),*) 'Stopping...'
            call tem_abort()
          end if

          call tem_time_set_clock(me = timing)
          me%lastWritten = timing

        else
          write(logUnit(1),*) ''
          write(logUnit(1),*) '!File to restart from does NOT exist!'
          write(logUnit(1),*) ''
          if (me%controller%init_on_missing) then
            write(logUnit(1),*) 'NOTE: performing initialization without'
            write(logUnit(1),*) '      reading data from restart as requested'
            write(logUnit(1),*) '      via the init_on_missing flag.'
            write(logUnit(1),*) ''
          else
            write(logUnit(1),*) '  Do not know how proceed, aborting now...'
            write(logUnit(1),*) '  If you want to perform the initialization'
            write(logUnit(1),*) '  when the restart file is missing, set'
            write(logUnit(1),*) '  the init_on_missing option in the restart'
            write(logUnit(1),*) '  table to true.'
            call tem_abort()
          end if
        end if
        write(logUnit(1),*) "*****************************"
        write(logUnit(1),*) ''
      end if ! If reading restart

      ! Now we get all the information about writing restart data.
      call aot_get_val( L       = conf,                      &
        &               thandle = restart_table,             &
        &               key     = 'write',                   &
        &               val     = me%controller%writePrefix, &
        &               ErrCode = iError                     )
      me%controller%writeRestart = (iError == 0)

      if (me%controller%writeRestart) then
        ! Read the time intervals for restart output from the Lua config.
        call tem_timeControl_load( me     = me%controller%timeControl, &
          &                        conf   = conf,                      &
          &                        parent = restart_table              )
        if (me%controller%readRestart) then
          call tem_timeControl_align_trigger(        &
            &    me     = me%controller%timeControl, &
            &    conf   = conf,                      &
            &    now    = timing,                    &
            &    parent = restart_table              )
        end if
        write(logUnit(1),*) "*****************************"
        write(logUnit(1),*) "Restart write parameters: "
        write(logUnit(1),*) "  prefix   : "//trim(me%controller%writePrefix)
        call tem_timeControl_dump(me%controller%timeControl, logUnit(2))
        write(logUnit(1),*) "*****************************"
        write(logUnit(1),*) ''
      end if

      call tem_horizontalSpacer(fUnit = logUnit(1))
    end if ! If restart table present

    call aot_table_close( L=conf, thandle=restart_table )

  end subroutine tem_load_restart
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Initialize the global restart data type and prepare for the restart output
  !!
  !! This routine is called as well for read as for write restart tasks
  !! It must be called after [[tem_load_restart]] and before the first call to
  !! any of tem_restart_*
  !! This routine is called only when restart is read from the restart table.
  !! If restart is performed from the initial condition table, then
  !! different routines (which are basically wrapped around by this one) are
  !! invoked
  !!
  subroutine tem_init_restart( me, solver, varMap, tree, subTree,   &
    &                          nDofs_write, chunkSize, solSpec_unit )
    ! -------------------------------------------------------------------- !
    !> The restart object to initialize.
    type(tem_restart_type), intent(inout)          :: me
    !> Details about the solver.
    type(tem_solveHead_type), optional, intent(in) :: solver
    !> Description of each variable system.
    !! This is ignored, if the data is provided by reading
    !! a restart.
    ! type(tem_varSys_type), intent(in)              :: varSys
    !> Contains position of variables to dump in restart file in global
    !! variable system for a scheme
    type(tem_varMap_type), intent(in)              :: varMap
    !> Mesh, provided in treelm format
    type(treelmesh_type), intent(in)               :: tree
    !> optional subTree of the given tree
    type(tem_subTree_type), optional, intent(in)   :: subTree
    !> number of degrees of freedom for each variable of the equation system
    integer, optional, intent(in)                  :: nDofs_write
    !> use predefined chunkSize
    integer, optional, intent(in)                  :: chunkSize
    !> Solver specific unit for restart header
    integer, optional, intent(in)                  :: solSpec_unit
    ! -------------------------------------------------------------------- !
    integer :: comm, rank, comm_size, locElems
    integer(kind=long_k) :: globElems, elemOff
    integer :: read_stat
    logical :: nUnitOpened
    character(len=320) :: solve_line
    ! -------------------------------------------------------------------- !

    ! set information about variables to be dumped in restart format.
    ! Do this irrespective of read or write restart
    me%varMap = varMap

    ! Get the total number of chunks necessary to write data-set through IO
    ! buffer to disk.
    if ( present( subTree ) ) then
      comm      = subTree%global%comm
      rank      = subTree%global%myPart
      comm_size = subTree%global%nParts
      globElems = subTree%global%nElems
      elemOff   = subTree%elemOffset
      locElems  = int(subTree%nElems)
    else
      comm      = tree%global%comm
      rank      = tree%global%myPart
      comm_size = tree%global%nParts
      globElems = tree%global%nElems
      elemOff   = tree%elemOffset
      locElems  = int(tree%nElems)
    end if

    ! Invoke the routine to allocate various variables
    call tem_init_restart_alloc( me          = me,         &
      &                          comm        = comm,       &
      &                          rank        = rank,       &
      &                          comm_size   = comm_size,  &
      &                          solver      = solver,     &
      &                          nDofs_write = nDofs_write )

    ! Loop over all the systems to create MPI types for reading restart
    ! Invoke the routine which creates MPI Types
    call tem_init_restart_create_types( me      = me,        &
      &                                 elemOff   = elemOff, &
      &                                 locElems  = locElems )

    call tem_restart_getTotalChunks( restart = me,         &
      &                              nElems  = locElems,   &
      &                              comm    = comm,       &
      &                              chunkSize = chunkSize )

    ! In root of this restart type,
    ! if write restart is active and solSpec_unit is present
    ! and opened then copy the content in solSpec_unit to internal
    ! restart scratch unit
    me%solSpec_unit = -1
    if (me%comm%rank == 0 .and. me%controller%writeRestart) then
      if (present(solSpec_unit)) then
        if (solSpec_unit>0) then
          write(dbgUnit(10),*) 'Writing solver specific info in restart:'
          inquire(unit=solSpec_unit, opened=nUnitOpened)
          if (nUnitOpened) then
            me%solSpec_unit = newunit()
            open(unit=me%solSpec_unit, status='scratch')
            rewind(solSpec_unit)
            do
              read(solSpec_unit,'(a)', iostat=read_stat) solve_line
              if (read_stat /= 0) EXIT
               write(dbgUnit(10),*) trim(solve_line)
              write(me%solSpec_unit,'(a)') trim(solve_line)
            end do
          end if !unitOpened
        end if  !solSpecUnit>0
      end if !present solSpecUnit
    end if !root process and active writeRestart

  end subroutine tem_init_restart
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> author: Kartik Jain
  !! This subroutine initializes the restart in case of reading initial
  !! conditions from restart file. The dependent scheme id is passed as input.
  !!
  subroutine tem_init_restart_alloc( me, comm, rank, comm_size, solver, nDofs_write )
    ! -------------------------------------------------------------------- !
    type(tem_restart_type), intent(inout) :: me
    integer, intent(in) :: comm, rank, comm_size
    type(tem_solveHead_type), optional, intent(in) :: solver
    integer, optional, intent(in) :: nDofs_write

    ! Inherit communicator information from the tree (or subTree)  we are to
    ! act on in this restart object:
    me%comm%comm      = comm
    me%comm%rank      = rank
    me%comm%comm_size = comm_size
    me%comm%root      = 0

    ! if provided reset the simulation name and solver tag
    if ( present(solver) ) then
      me%header%simName          = trim(solver%simName)
      me%header%solverTag        = tem_solverTag(solver)
      me%header%solverConfigFile = solver%configFile
    end if

    if ( present(nDofs_write) ) then
      me%write_file%nDofs = nDofs_write
    else
      me%write_file%nDofs = 1
    end if

    ! Set the prefix for the header file.
    me%header%headerPrefix = trim(me%controller%writePrefix) &
      &                      // trim(me%header%simName)
    ! Set the prefix for all binary files.
    me%header%binPrefix = trim(me%controller%writePrefix) &
      &                      // trim(me%header%simName)


    ! Reset all the values of the chunk information
    me%nChunkElems = 0
    me%write_file%chunkSize  = 0
    me%read_file%chunkSize   = 0
    me%write_file%nChunks    = 0
    me%read_file%nChunks     = 0
    me%write_file%maxnChunks = 0
    me%read_file%maxnChunks  = 0

    end subroutine tem_init_restart_alloc
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This subroutine creates MPI types for reading the restart.
  subroutine tem_init_restart_create_types( me, elemOff, locElems)
    ! -------------------------------------------------------------------- !
    type(tem_restart_type), intent(inout) :: me
    integer(kind=long_k),intent(in)       :: elemOff
    integer,intent(in)                    :: locElems
    ! -------------------------------------------------------------------- !
    ! local variables
    integer                               :: iError
    integer                               :: typesize
    ! -------------------------------------------------------------------- !

    if ( int(me%varMap%nScalars, MPI_OFFSET_KIND)             &
      &  * int(me%write_file%nDofs, MPI_OFFSET_KIND)          &
      &  * int(locelems, MPI_OFFSET_KIND) * 8_MPI_OFFSET_KIND &
      &  >= 2147483648_MPI_OFFSET_KIND                        ) then
      write(logunit(1),*) 'Error: local partition greater 2GB!'
      write(logunit(1),*) 'Most MPI implementations do not support this.'
      write(logunit(1),*) 'I will abort now, as this will likely result in'
      write(logunit(1),*) 'an error later on anyway.'
      write(logunit(1),*)
      write(logunit(1),*) 'Please make sure, to use a sufficient number of'
      write(logunit(1),*) 'processes to reduce the size of local partitions'
      write(logunit(1),*) 'to two GB.'
      write(logunit(1),*) 'Which would be less than ',                         &
        &                 int( real(2147483648_MPI_OFFSET_KIND, kind=rk)       &
        &                      / real( me%varMap%nScalars*me%write_file%nDofs  &
        &                              * 8_MPI_OFFSET_KIND, kind=rk)        ), &
        &                 ' elements for your element size.'
      call tem_abort()
    end if
    ! A contiguous type to describe the vector per element.
    ! MPI_TYPE_CONTIGUOUS(COUNT, OLDTYPE, NEWTYPE, IERROR)
    call MPI_Type_contiguous( me%varMap%nScalars*me%write_file%nDofs, &
      &                       rk_mpi,                                 &
      &                       me%write_file%vectype,                  &
      &                       iError                                  )
    call check_mpi_error( iError,                                             &
      &                   'create contiguous (write) vectype in init_restart' )

    ! Commit the type for creation
    call MPI_Type_commit( me%write_file%vectype, iError )
    call check_mpi_error( iError, 'commit (write) vectype in init_restart' )

    ! Create a MPI Contiguous as ftype for file view
    call MPI_Type_contiguous( locElems, me%write_file%vectype, &
      &                       me%write_file%ftype, iError      )
    call check_mpi_error( iError,                                           &
      &                   'create contiguous (write) ftype in init_restart' )

    ! commit the new contiguous type
    call MPI_Type_commit( me%write_file%ftype, iError )
    call check_mpi_error( iError, 'commit ftype in init_restart')

    ! get size of element
    call MPI_TYPE_SIZE(me%write_file%vectype, typesize, iError )
    call check_mpi_error(iError,'typesize in init_restart')

    ! set the start of view
    me%write_file%displacement= elemOff * typesize * 1_MPI_OFFSET_KIND

    if (me%read_file%nDofs /= me%write_file%nDofs) then
      if ( int(me%varMap%nScalars, MPI_OFFSET_KIND)             &
        &  * int(me%read_file%nDofs, MPI_OFFSET_KIND)           &
        &  * int(locelems, MPI_OFFSET_KIND) * 8_MPI_OFFSET_KIND &
        &  >= 2147483648_MPI_OFFSET_KIND                        ) then
        write(logunit(1),*) 'Error: local partition from restart greater 2GB!'
        write(logunit(1),*) 'Most MPI implementations do not support this.'
        write(logunit(1),*) 'I will abort now, as this will likely result in'
        write(logunit(1),*) 'an error later on anyway.'
        write(logunit(1),*)
        write(logunit(1),*) 'Please make sure, to use a sufficient number of'
        write(logunit(1),*) 'processes to reduce the size of local partitions'
        write(logunit(1),*) 'to two GB.'
        write(logunit(1),*) 'Which would be less than ',                       &
          &                 int( real(2147483648_MPI_OFFSET_KIND, kind=rk)     &
          &                      / real( me%varMap%nScalars*me%read_file%nDofs &
          &                              * 8_MPI_OFFSET_KIND, kind=rk)      ), &
          &                 ' elements for your elements in the restart file.'
        call tem_abort()
      end if

      ! MPI_TYPE_CONTIGUOUS(COUNT, OLDTYPE, NEWTYPE, IERROR)
      call MPI_Type_contiguous( me%varMap%nScalars*me%read_file%nDofs, &
        &                       rk_mpi,                                &
        &                       me%read_file%vectype,                  &
        &                       iError                                 )
      call check_mpi_error( iError,                                            &
        &                   'create contiguous (read) vectype in init_restart' )

      ! Commit the type for creation
      call MPI_Type_commit( me%read_file%vectype, iError )
      call check_mpi_error( iError, 'commit (read) vectype in init_restart')

      ! Create a MPI Contiguous as ftype for file view
      call MPI_Type_contiguous( locElems, me%read_file%vectype, &
        &                       me%read_file%ftype, iError      )
      call check_mpi_error( iError,                                          &
        &                   'create contiguous (read) ftype in init_restart' )

      ! commit the new contiguous type
      call MPI_Type_commit( me%read_file%ftype, iError )
      call check_mpi_error( iError, 'commit (read) ftype in init_restart' )

      ! get size of element
      call MPI_TYPE_SIZE(me%read_file%vectype, typesize, iError )
      call check_mpi_error(iError,'typesize in init_restart')

      ! set the start of view
      me%read_file%displacement = elemOff * typesize * 1_MPI_OFFSET_KIND

    else
      me%read_file%vectype = me%write_file%vectype
      me%read_file%ftype = me%write_file%ftype
      me%read_file%displacement=me%write_file%displacement
    end if

  end subroutine tem_init_restart_create_types
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This subroutine evaluated get_element and dump each chunk
  !!
  !! This routine is used in tracking to dump in data in harvester format
  !! for single variable system
  subroutine tem_restart_dump_data( restart, varSys, tree, time, subTree )
    ! &                               transientReduce                       )
    ! -------------------------------------------------------------------- !
    !> The restart object describing how and what to output.
    type(tem_restart_type), intent(inout) :: restart
    !> Description of the available variable system to get the given varnames
    !! from.
    type(tem_varSys_type), intent(in) :: varsys
    !> Mesh to write the data on.
    type(treelmesh_type), intent(in) :: tree
    !> Point in time to use for this data.
    !!
    !! Can be important for space-time function evaluations.
    type(tem_time_type), intent(in) :: time
    !> Optional restriction of the elements to output.
    type(tem_subtree_type), optional, intent(in) :: subtree
    !> transient reducution
    ! type(tem_transient_reduction_type), intent(inout) :: transientReduce(:)
    ! -------------------------------------------------------------------- !
    integer :: nVars, nElems, nScalars, elemOff, nChunkElems
    integer :: iElem, iChunk
    integer :: buf_start, buf_end
    real(kind=rk), allocatable :: res(:)
    integer, allocatable :: elemPos(:)
    integer :: ioStatus( mpi_status_size )
    integer :: iError
    ! -------------------------------------------------------------------- !
    allocate(res(io_buffer_size))

    ! Number of variables to dump
    nVars = restart%varMap%varPos%nVals

    ! Number of scalars in current output
    nScalars = restart%varMap%nScalars

    if (present(subTree)) then
      nElems = subTree%nElems
    else
      nElems = tree%nElems
    end if

    ! open transient reduction
    ! call tem_transient_reduction_open( me   = transientReduce, &
    !   &                                time = time%sim         )

    ! allocate elemPos to size of chunkSize
    allocate(elemPos(restart%write_file%chunkSize))

    ! Process all chunks to derive the quantities defined in the tracking object
    do iChunk = 1, restart%write_file%nChunks
      ! Number of elements read so far in previous chunks.
      elemOff = ((iChunk-1)*restart%write_file%chunkSize)

      ! number of elements written to THIS chunk
      nChunkElems = min(restart%write_file%chunkSize, nElems-elemOff)
      restart%nChunkElems = nChunkElems

      ! Compute the element lower and upper bound for the current chunk
      buf_start = elemOff + 1
      buf_end = elemOff + nChunkElems

      if (present(subTree)) then
        elemPos(1:nChunkElems) = subTree%map2Global(buf_start:buf_end)
      else
        elemPos(1:nChunkElems) = (/ (iElem, iElem=buf_start, buf_end) /)
      end if

      ! evaluate all variables on current chunk
      call tem_get_element_chunk(varSys  = varSys,                      &
        &                        varPos  = restart%varMap%varPos        &
        &                                                 %val(:nVars), &
        &                        elemPos = elemPos(1:nChunkElems),      &
        &                        time    = time,                        &
        &                        tree    = tree,                        &
        &                        nElems  = nChunkElems,                 &
        &                        nDofs   = restart%write_file%nDofs,    &
        &                        res     = res                          )

      ! perform transient reduction
      ! @todo KM: Check transientReduction when nDofs>1
      ! call tem_transient_reduction_apply( me          = transientReduce,    &
      !   &                                 chunk       = res,                &
      !   &                                 offset      = buf_start - 1,      &
      !   &                                 nChunkElems = nChunkElems,        &
      !   &                                 varSys      = varSys,             &
      !   &                                 varPos      = restart%varMap   &
      !   &                                               %varPos%val(:nVars) )

      ! Now write the results into the file, using the view defined in
      ! [[tem_restart_openWrite]]
      ! arguments:
      ! file handle = binary unit opened in mpi_file_open
      ! initial address of buffer = first entry to dump within the chunk
      !                             this is ad
      call mpi_file_write_all( restart%binaryUnit, res,        &
        &                      restart%nChunkElems,            &
        &                      restart%write_file%vectype,     &
        &                      iostatus, iError                )
      call check_mpi_error( iError,'File write all in tem_restart_dump_data')

    end do !iChunk

    deallocate(elemPos)
    deallocate(res)

    ! close transient reduction
    ! call tem_transient_reduction_close( me = transientReduce )

  end subroutine tem_restart_dump_data
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This subroutine dumps the given chunk to a given position in the restart
  !! file.
  !!
  subroutine tem_restart_writeData( restart, chunk )
    ! -------------------------------------------------------------------- !
    !> The restart object describing how and what to output.
    type(tem_restart_type), intent(inout) :: restart
    !> The data to output.
    !! It is organized as a serialized array of all scalar entries of all
    !! variable systems. Where first all the data for the nElems of the first
    !! variable system is provided. Within each variable system the data is
    !! organized elementwise.
    real(kind=rk), intent(in) :: chunk(:)
    ! -------------------------------------------------------------------- !
    integer :: ioStatus( mpi_status_size )
    integer :: iError
    ! -------------------------------------------------------------------- !
    ! Now write the collected data (from the state array and prob. derived
    ! quantities) into the file, using the view defined in
    ! [[tem_restart_openWrite]]
    ! arguments:
    ! file handle = binary unit opened in mpi_file_open
    ! initial address of buffer = first entry to dump within the chunk
    !                             this is ad
    call mpi_file_write_all( restart%binaryUnit, chunk,         &
      &                      restart%nChunkElems,               &
      &                      restart%write_file%vectype,        &
      &                      iostatus, iError                   )
    call check_mpi_error( iError,'File write all in tem_restart_writeData')

  end subroutine tem_restart_writeData
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Read data from a restart file.
  !!
  subroutine tem_restart_readData( restart, chunk )
    ! -------------------------------------------------------------------- !
    !> Restart description to read the data from
    type(tem_restart_type), intent(inout) :: restart
    !> Chunk of memory to put the read data into
    real(kind=rk), intent(out) :: chunk(:)
    ! -------------------------------------------------------------------- !
    ! defining local variables
    integer :: nWords
    integer :: offset        ! offset for the different var systems
    ! -------------------------------------------------------------------- !

    ! @todo:KJ: For the new restart only the requested variable system should
    ! be read.
    offset = 1

    call tem_restart_readData_single( restart, chunk, offset )
    nWords = restart%varMap%nScalars * restart%nChunkElems &
      &                              * restart%read_file%nDofs
    offset = offset + nWords

  end subroutine tem_restart_readData
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> subroutine which reads data from restart file corresponding to the input
  !! variable system
  !!
  subroutine tem_restart_readData_single( restart, chunk, offset )
    ! -------------------------------------------------------------------- !
    !> Restart object to read the data from
    type(tem_restart_type), intent(inout) :: restart
    !> Memory to chunk to put the data into
    real(kind=rk), intent(out)            :: chunk(:)
    !> Offset of the chunk of data to get in the global data
    integer, intent(in)                   :: offset
    ! -------------------------------------------------------------------- !
    ! defining local variables
    integer :: iostatus( MPI_STATUS_SIZE )
    integer :: iError
    ! -------------------------------------------------------------------- !

    ! MPI_FILE_READ_ALL(fh, buf, count, datatype, status)
    call mpi_file_read_all( restart%binaryUnit, chunk(offset), &
      &                     restart%nChunkElems,               &
      &                     restart%read_file%vectype,         &
      &                     iostatus, iError                   )

    call check_mpi_error( iError,'File write all in tem_restart_readData_single')

  end subroutine tem_restart_readData_single
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> open the restart dump file and write out the 'normal' restart header.
  !!
  subroutine tem_restart_openRead(me)
    ! -------------------------------------------------------------------- !
    !> the restart information
    type(tem_restart_type) :: me

    write(logUnit(1),*)' Open read... '

    ! read the restart header file
    ! prepare everything for the serialization of the data
    ! Invoke the routine to open single variable system from restart file
    call tem_restart_openRead_single(me)

  end subroutine tem_restart_openRead
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Open the restart dump file and set file view for the input variable system
  !!
  subroutine tem_restart_openRead_single(me)
    ! -------------------------------------------------------------------- !
    !> the restart information
    type(tem_restart_type) :: me
    ! -------------------------------------------------------------------- !
    ! variables to catch possible MPI I/O errors
    integer :: iError
    ! -------------------------------------------------------------------- !

    ! open the binary file for MPI I/O
    call MPI_FILE_OPEN( me%comm%comm,                 &
      &                 trim( me%header%binName ),    &
      &                 MPI_MODE_RDONLY,              &
      &                 MPI_INFO_NULL, me%binaryUnit, &
      &                 iError                        )
    call check_mpi_error( iError,'Open File in tem_restart_openRead_single')

    ! set the view of each process on the file opened above
    call MPI_FILE_SET_VIEW( me%binaryUnit, me%read_file%displacement,      &
      &                     me%read_file%vectype,                          &
      &                     me%read_file%ftype, "native",                  &
      &                     MPI_INFO_NULL, iError                          )
    call check_mpi_error( iError,'Set File view in tem_restart_openRead_single')

  end subroutine tem_restart_openRead_single
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This subroutine closes the restart dump file and writes the last header.
  subroutine tem_restart_closeRead( me )
    ! -------------------------------------------------------------------- !
    !> the restart information
    type(tem_restart_type) :: me
    ! -------------------------------------------------------------------- !
    ! temporary buffer array for empty reads
    real(kind=rk) :: chunk(1)
    integer :: iChunk
    ! -------------------------------------------------------------------- !

    ! If local process has performed less reads than global maximum required,
    ! we need to complete as many empty reads as the difference to maxnChunks
    ! to satisfy collective MPI-IO operations.
    me%nChunkElems = 0
    do iChunk = me%read_file%nChunks+1, me%read_file%maxnChunks
      call tem_restart_readData(restart = me, chunk = chunk)
    end do

    ! Invoke the routine to close current variable system file
    call tem_restart_closeRead_single(me)

    write(logUnit(1),*)' Closed Read.'

  end subroutine tem_restart_closeRead
  ! ************************************************************************ !

!TG: merge together
  ! ************************************************************************ !
  !> Close the restart dump file corresponding to a particular variable system
  !!
  subroutine tem_restart_closeRead_single(me)
    ! -------------------------------------------------------------------- !
    !> the restart infotmation
    type(tem_restart_type) :: me
    ! -------------------------------------------------------------------- !
    ! variables to catch possible MPI I/O errors
    integer :: iError
    ! -------------------------------------------------------------------- !

    ! now close the binary file
    call MPI_File_close(me%binaryUnit, iError)
    call check_mpi_error( iError, 'File close in tem_restart_closeRead_single' )

  end subroutine tem_restart_closeRead_single
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> open the restart dump file and write out the 'normal' restart header as
  !! well as the mesh.
  !!
  subroutine tem_restart_openWrite( me, tree, timing, varSys, subTree, label, &
    &                               suffix )
    ! -------------------------------------------------------------------- !
    !> the restart infotmation
    type(tem_restart_type)                          :: me
    !> mesh, provided in treelm format
    type(treelmesh_type)                            :: tree
    !> current simulation time information
    type(tem_time_type),intent(in)                  :: timing
    !> the used var systeme
    type(tem_varSys_type), intent(in)               :: varSys
    !> optional subTree of the given tree
    type(tem_subTree_type), optional, intent(inout) :: subTree
    !> additional label for the filename (needed for tracking in harvester
    !! format)
    character(len=*), optional, intent(in)          :: label
    !> optional suffix (if present NO timestamp will be added!!!!)
    character(len=*), optional, intent(in)          :: suffix
    ! -------------------------------------------------------------------- !
    ! variables to catch possible MPI I/O errors
    integer :: iError
    integer :: pos
    character(len=pathLen) :: prefix
    logical :: meshChange_loc
    type(tem_global_type) :: global_loc
    ! -------------------------------------------------------------------- !

    ! Update the timestamp
    me%header%timestamp = trim(tem_time_sim_stamp(timing))

    ! Set the iteration to know when the last restart file was written
    me%lastWritten = timing

    if ( present(subTree) ) then
      global_loc = subTree%global
    else
      global_loc = tree%global
    end if

    meshChange_loc = global_loc%meshChange

    ! communicate wether the mesh has changed since last time dumping it
    call MPI_ALLREDUCE( meshChange_loc, global_loc%meshChange, 1,     &
      &                 MPI_LOGICAL, MPI_LOR, global_loc%comm, iError )

    ! if the mesh has changed ...
    if (global_loc%meshChange) then
      ! ... set the meshChange to false
      global_loc%meshChange = .false.
      ! ... get the position of the last path seperator
      pos = INDEX(trim(global_loc%dirname), pathSep, .true.)
      if ( present(label) ) then
        prefix = trim(global_loc%dirname(1:pos))//trim(label)//'_'
      else
        prefix = trim(global_loc%dirname(1:pos))
      end if
      if ( present(suffix) ) then
        ! change the dirname using NO timestamp but the suffix
        write(global_loc%dirname,'(a)') trim(prefix)//trim(suffix)//'_'
      else
        ! ... change the dirname
        write(global_loc%dirname,'(a)') trim(prefix)                        &
          &                             // trim( me%header%timestamp ) // '_'
      end if
      ! ... remove a possible predefined tag
      global_loc%predefined = ''
      ! ... copy back the global information to the tree or subTree and dump it
      if ( present(subTree) ) then
        subTree%global = global_loc
        call tem_dump_subTree( subTree, tree )
      else
        tree%global = global_loc
        call dump_treelmesh( tree )
      end if
    end if

    if ( present(suffix) ) then
      ! define the name of the file to write the binary data to without
      ! timestamp but using the suffix
      write(me%header%binName,'(a)') trim( me%header%binPrefix )//'_'  &
        &                               // trim( suffix )              &
        &                               // tem_create_EndianSuffix()
    else
      ! define the name of the file to write the binary data to
      write(me%header%binName,'(a)') trim( me%header%binPrefix )//'_'  &
        &                               // trim( me%header%timestamp ) &
        &                               // tem_create_EndianSuffix()
    end if

    ! open the binary file for MPI I/O
    call MPI_FILE_OPEN( me%comm%comm,                    &
      &                 trim( me%header%binName ),       &
      &                 MPI_MODE_WRONLY+MPI_MODE_CREATE, &
      &                 MPI_INFO_NULL, me%binaryUnit,    &
      &                 iError                           )
    call check_mpi_error( iError, 'File open of '         &
      &                        // trim(me%header%binName) &
      &                        // ' for writing in '      &
      &                        // 'tem_restart_openWrite' )

    call MPI_FILE_SET_VIEW( me%binaryUnit, me%write_file%displacement,     &
      &                     me%write_file%vectype,                         &
      &                     me%write_file%ftype, "native",                 &
      &                     MPI_INFO_NULL, iError                          )
    call check_mpi_error( iError,'set File view in tem_restart_openWrite')

    ! write out a regular restart header
    ! @todo: if [[tem_restart_writeHeader]] is only called here, then it should
    !       not be public. It would be better not to call it here, but let user
    !       decide where to call it.
    call tem_restart_writeHeader( me      = me,      &
      &                           tree    = tree,    &
      &                           subTree = subTree, &
      &                           timing  = timing,  &
      &                           varSys  = varSys,  &
      &                           suffix  = suffix   )

  end subroutine tem_restart_openWrite
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This subroutine closes the restart dump file and writes the last header.
  !!
  subroutine tem_restart_writeHeader( me, tree, varSys, subTree, timing, &
    &                                 lastHeader, suffix )
    ! -------------------------------------------------------------------- !
    !> the restart header info
    type(tem_restart_type) :: me
    !> mesh, provided in treelm format
    type(treelmesh_type) :: tree
    !> global variable system defined in solver
    type(tem_varSys_type), intent(in) :: varSys
    !> optional subTree of the given tree
    type(tem_subTree_type), optional, intent(in) :: subTree
    !>
    type(tem_time_type), intent(in) :: timing
    !> is this header a last header
    logical, optional :: lastHeader
    !> optional suffix (if present NO timestamp will be added!!!!)
    character(len=*), optional, intent(in) :: suffix
    ! -------------------------------------------------------------------- !
    character(len=1024) :: filename
    character(len=320) :: solve_line
    type(aot_out_type) :: conf ! aotus lua state to write output
    logical :: isLast
    type(tem_global_type) :: global_loc
    integer(kind=long_k) :: nElems
    integer :: punit
    integer :: read_stat
    logical :: nUnitOpened
    ! -------------------------------------------------------------------- !

    if ( present(lastHeader) ) then
      isLast = lastHeader
    else
      isLast = .false.
    end if

    ! now write out the header file
    ! only one process (root) writes out the header file
    if (me%comm%rank == 0) then
      ! choose the right filename
      if (isLast) then
        filename = trim( me%header%headerPrefix )// '_lastHeader.lua'
      else
        if ( present(suffix) ) then
          filename = trim( me%header%headerPrefix )//'_header_' &
            &      // trim(suffix)//'.lua'
        else
          filename = trim( me%header%headerPrefix )//'_header_' &
            &      // trim(me%header%timestamp)//'.lua'
        end if
      end if

      ! Open up the restart header lua file, so we can write the stuff using
      ! the aotus library
      punit = newUnit()
      call tem_open( newunit = punit,          &
        &            file    = trim(filename), &
        &            action  = 'write',        &
        &            status  = 'replace',      &
        &            recl    = 360             )

      call aot_out_open(put_conf = conf, outUnit = punit)

      ! the binary names; for multiple schemes a list of binaries have to be
      ! written
      call aot_out_open_table(conf, 'binary_name')
      call aot_out_val( put_conf = conf,                         &
        &               val      = trim(me%header%binName) )
      call aot_out_close_table(conf)

      ! the solver config file names
      call aot_out_val( put_conf = conf,                             &
        &               val      = trim(me%header%solverConfigFile), &
        &               vname    = 'solver_configFile'               )

      if( present( subTree ))then
        global_loc = subTree%global
        nElems = subTree%global%nElems
      else
        global_loc = tree%global
        nElems = tree%global%nElems
      end if

      ! the mesh info
      call tem_mesh_out( me = global_loc, conf = conf )
      call aot_out_val(conf, tree%weights_file, 'weights')

      ! the time stamp
      !call aot_out_val( conf, trim(me%header%timestamp), 'time_point')
      call tem_time_out( conf = conf, me = timing, key = 'time_point' )
      ! the total number of elements
      call aot_out_val( conf, nElems, 'nElems')
      ! the number of dofs for each scalar variable of the equation system
      call aot_out_val( conf, me%write_file%nDofs, 'nDofs')
      ! the solver tag (solver name and version)
      call aot_out_val( conf, trim( me%header%solverTag), 'solver')
      ! the variable system, incl. the solver specific
      ! information from the solSpec data
      call tem_varSys_out( me         = varSys,                 &
        &                  conf       = conf,                   &
        &                  dumpVarPos = me%varMap%varPos%val(:) )
      ! close the restart header file
      call aot_out_close(conf)

      ! Append the solver specific data, stored in a scratch file to the header
      ! file, if one is provided by caller.
      if ( me%solSpec_unit>0 ) then
        inquire(unit=me%solSpec_unit, opened=nUnitOpened)
        if (nUnitOpened) then
          rewind(me%solSpec_unit)
          do
           read(me%solSpec_unit,'(a)', iostat=read_stat) solve_line
           if (read_stat /= 0) EXIT
           write(punit,'(a)') trim(solve_line)
          end do
        end if
      end if

      close(punit)

    end if

  end subroutine tem_restart_writeHeader
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> read the restart header lua file and hand the information to the
  !! required data types, re-set the time ...
  !!
  subroutine tem_restart_readHeader( me, timing, globProc, tree )
    ! -------------------------------------------------------------------- !
    !> the restart header info
    type(tem_restart_type) :: me
    !> the timing for re-setting the times
    type(tem_time_type), intent(inout) :: timing
    !> Global communicator
    type( tem_comm_env_type ), intent(in) :: globProc
    !> mesh, provided in treelm format
    type(treelmesh_type), intent(inout)    :: tree
    ! -------------------------------------------------------------------- !
    integer :: thandle
    integer :: iError
    type( flu_State ) :: conf
    character(len=labelLen) :: buffer
    ! -------------------------------------------------------------------- !
    write(logUnit(1),*) 'Opening Restart Header '         &
      &                 // trim(me%controller%readFileName)

    ! Open the restart header file
    call tem_open_distconf( L        = conf,                             &
      &                     filename = trim(me%controller%readFileName), &
      &                     proc     = globProc                          )

    ! Load the number of elements from the restart header file for sanity check
    call aot_get_val( L       = conf,             &
      &               key     = 'nElems',         &
      &               val     = me%header%nElems, &
      &               ErrCode = iError,           &
      &               default = 1                 )

    ! Load the solver name
    call aot_get_val( L       = conf,     &
      &               key     = 'solver', &
      &               val     = buffer,   &
      &               ErrCode = iError,   &
      &               default = ''        )

    me%header%solverTag = trim( buffer )
    write(logUnit(1),*) 'Solver: '// trim(me%header%solverTag )

    ! Load the solver config file
    call aot_get_val( L       = conf,                       &
      &               key     = 'solver_configFile',        &
      &               val     = me%header%solverConfigFile, &
      &               ErrCode = iError,                     &
      &               default = ''                          )

    ! Load the number of dofs for each scalar variable of the equation system
    call aot_get_val( L       = conf,               &
      &               key     = 'nDofs',            &
      &               val     = me%read_file%nDofs, &
      &               ErrCode = iError,             &
      &               default = 1                   )

    call load_tem( me      = tree,               &
      &            conf    = conf,               &
      &            myPart  = globProc%rank,      &
      &            nParts  = globProc%comm_size, &
      &            comm    = globProc%comm       )

    ! Load the timestamp from the header
    call tem_time_load( conf        = conf,              &
      &                 key         = 'time_point',      &
      &                 me          = timing,            &
      &                 clock_start = timing%clock_start )

    me%header%timestamp = trim(tem_time_sim_stamp(timing))
    write(logUnit(1),*) 'Restarting from point in time:'
    call tem_time_dump(timing, logUnit(1))

    ! Load the variable systems
    call tem_varSys_load( me = me%header%varSys, conf = conf )

    call tem_varSys_dump( me = me%header%varSys, outUnit = dbgUnit(3) )

    call aot_table_open( L=conf, thandle = thandle, key = 'binary_name' )

    ! Load the binary file names for each variable system defined
    call aot_get_val( L       = conf,                      &
      &               thandle = thandle,                   &
      &               pos     = 1,                      &
      &               val     = me%header%binName, &
      &               ErrCode = iError,                    &
      &               default = ''                         )

    call aot_table_close( L = conf, thandle = thandle )

    call close_config( conf )

  end subroutine tem_restart_readHeader
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Complete a number of empty writes (due to higher amount of mpi_file_writes
  !! from other processes to finalize the write process), close the restart dump
  !! file and write the last header.
  !!
  subroutine tem_restart_closeWrite( me, tree, timing, varSys, subTree )
    ! -------------------------------------------------------------------- !
    !> the restart information
    type(tem_restart_type) :: me
    !> mesh, provided in treelm format
    type(treelmesh_type) :: tree
    !> The timing object holding the simulation time
    type(tem_time_type), intent(in) :: timing
    !> global variable system defined in solver
    type(tem_varSys_type), intent(in) :: varSys
    !> optional subTree of given tree
    type(tem_subTree_type), optional :: subTree
    ! -------------------------------------------------------------------- !
    ! variables to catch possible MPI I/O errors
    integer :: iError
    integer :: iChunk ! chunk counter
    ! temporary buffer array for empty writes
    real(kind=rk) :: chunk(1)
    ! -------------------------------------------------------------------- !
    ! Check, if the number of calls to mpi_file_write_all corresponds to
    ! the maximum number throughout all processes
    ! if local process has performed less writes, we need to complete
    ! as many empty writes as the difference to maxnChunks
    me%nChunkElems = 0
    do iChunk = me%write_file%nChunks+1, me%write_file%maxnChunks
      call tem_restart_writeData( restart = me, chunk = chunk )
    end do

    ! now close the binary file
    call MPI_File_close(me%binaryUnit, iError)
    call check_mpi_error( iError,'File close in tem_restart_closeWrite')
    call MPI_Barrier( me%comm%comm, iError )

    ! Now write out the 'last' header which points to the last successful
    ! restart file.
    call tem_restart_writeHeader( me         = me,      &
      &                           tree       = tree,    &
      &                           subTree    = subTree, &
      &                           timing     = timing,  &
      &                           varSys     = varSys,  &
      &                           lastHeader = .true.   )

  end subroutine tem_restart_closeWrite
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> calculate the maximum number of elements which fit into the output buffer
  !! = chunk
  !!
  !! for a given set of variable systems with their nScalar values to dump
  !! Also, calculate the number of local chunks required to dump all the data
  !! = nChunks
  !! Finally find the globally largest number of nChunks
  !!
  subroutine tem_restart_getTotalChunks( restart, nElems, comm, chunkSize )
    ! -------------------------------------------------------------------- !
    !> the restart type
    type( tem_restart_type ), intent(inout) :: restart
    !> mesh, provided in treelm format
    !type(treelmesh_type), intent(in) :: tree
    !> optional subTree
    ! type(tem_subTree_type), optional, intent(in) :: subTree
    integer, intent(in) :: nElems, comm
    !> optional predefined chunksize
    integer, optional, intent(in) :: chunkSize
    ! -------------------------------------------------------------------- !
    integer :: nTotalScalars
    integer :: iError ! MPI error
    integer :: rank
    ! -------------------------------------------------------------------- !
    ! Get the number of total scalars.
    nTotalScalars = restart%varMap%nScalars
    if ( nTotalScalars ==  0 ) then
      write(logUnit(0),*) '!! Error !! No variable found in restart varSys !!'
      write(logUnit(0),*) 'May be variable label in requested system is not ' &
        &                 // 'found in global varsys'
      call tem_abort()
    endif

    if( present( chunkSize ))then
      restart%read_file%chunkSize = chunkSize
      restart%write_file%chunkSize = chunkSize
    else
      ! Get the number of elements that fit into the IO buffer.
      restart%read_file%chunkSize                            &
        &  = min( io_buffer_size                             &
        &         / (nTotalScalars*restart%read_file%nDofs), &
        &         nElems                                   )
      restart%write_file%chunkSize                            &
        &  = min( io_buffer_size                              &
        &         / (nTotalScalars*restart%write_file%nDofs), &
        &         nElems                                    )
    end if

    ! check if at least one complete element fits into the buffer
    ! if not abort since no valid output can be garanteed
    if (restart%write_file%chunkSize <= 0) then
      write(logUnit(0),*) 'The chosen io_buffer_size of ', io_buffer_size
      write(logUnit(0),*) 'is too small for outputting ', nTotalScalars
      write(logUnit(0),*) 'scalar values with ', restart%write_file%nDofs
      write(logUnit(0),*) 'degrees of freedom!'
      write(logUnit(0),*) 'Please increase the io_buffer_size to' &
        &                 // ' at least (MB) ', &
        &                 real(nTotalScalars*restart%write_file%nDofs) &
        &                 / real(131072)
      call tem_abort()
    end if

    if (restart%read_file%chunkSize <= 0) then
      write(logUnit(0),*) 'The chosen io_buffer_size of ', io_buffer_size
      write(logUnit(0),*) 'is too small for reading ', nTotalScalars
      write(logUnit(0),*) 'scalar values with ', restart%read_file%nDofs
      write(logUnit(0),*) 'degrees of freedom!'
      write(logUnit(0),*) 'Please increase the io_buffer_size to'     &
        &                 // ' at least (MB) ',                       &
        &                 real(nTotalScalars*restart%read_file%nDofs) &
        &                 / real(131072)
      call tem_abort()
    end if

    ! Get the number of chunks which are needed to dump all values to disk.
    ! This needs to be rounded up, to cover also a possible last incomplete
    ! chunk.
    restart%write_file%nChunks                           &
      &  = ceiling( real(nElems, kind=rk)              &
      &             / real(restart%write_file%chunkSize, kind=rk) )

    restart%read_file%nChunks                           &
      &  = ceiling( real(nElems, kind=rk)             &
      &             / real(restart%read_file%chunkSize, kind=rk) )

    ! identify the maximum number of chunks throughout all processes
    ! and store that into restart%maxnChunks
    call MPI_Allreduce( restart%write_file%nChunks,        &
      &                 restart%write_file%maxnChunks, 1,  &
      &                 MPI_INTEGER, MPI_MAX, comm, iError )
    call MPI_Allreduce( restart%read_file%nChunks,         &
      &                 restart%read_file%maxnChunks, 1,   &
      &                 MPI_INTEGER, MPI_MAX, comm, iError )

    call MPI_Comm_Rank( comm, rank, iError )

  end subroutine tem_restart_getTotalChunks
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Finalizing a restart object
  !!
  subroutine tem_restart_finalize( me )
    ! -------------------------------------------------------------------- !
    !> the restart type to close
    type( tem_restart_type ), intent(inout) :: me
    ! -------------------------------------------------------------------- !
    integer :: iError
    logical :: nUnitOpened
    ! -------------------------------------------------------------------- !
    ! close solverspecific scratch unit
    if (me%comm%rank == 0 .and. me%controller%writeRestart) then
      if (me%solSpec_unit>0) then
        inquire(unit=me%solSpec_unit, opened=nUnitOpened)
        if (nUnitOpened) close(me%solSpec_unit)
      end if
    end if

    ! free the contiguous type
    if (me%controller%writeRestart .or. me%controller%readRestart) then
      call MPI_Type_free(me%write_file%ftype, iError)
      call check_mpi_error( iError,'Free write-ftype in restart_finalize')
      call MPI_Type_free(me%write_file%vectype, iError)
      call check_mpi_error( iError,'Free write-vectype in restart_finalize')
      if (me%read_file%ndofs /= me%write_file%ndofs) then
        call MPI_Type_free(me%read_file%ftype, iError)
        call check_mpi_error( iError,'Free read-ftype in restart_finalize')
        call MPI_Type_free(me%read_file%vectype, iError)
        call check_mpi_error( iError,'Free read-vectype in restart_finalize')
      end if
    end if

  end subroutine tem_restart_finalize
  ! ************************************************************************ !


end module tem_restart_module
! **************************************************************************** !
