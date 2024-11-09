! Copyright (c) 2015-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2017 Raphael Haupt <Raphael.Haupt@student.uni-siegen.de>
! Copyright (c) 2017, 2019, 2021 Harald Klimach <harald.klimach@dlr.de>
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
!> This module encapsulates the data type for ascii and asciiSpatial
!! output format.
!!
!! \author Kannan Masilamani
!!
module hvs_ascii_module
  use env_module,              only: rk, PathLen, LabelLen, OutLen, newUnit,   &
    &                                io_buffer_size, long_k
  use tem_solveHead_module,    only: tem_solveHead_type, tem_solverTag
  use tem_reduction_spatial_module,                                &
    &                          only: tem_reduction_spatial_type,   &
    &                                tem_reduction_spatial_init,   &
    &                                tem_reduction_spatial_open,   &
    &                                tem_reduction_spatial_close,  &
    &                                tem_reduction_spatial_append, &
    &                                tem_reduction_spatial_toChunk
  use tem_comm_env_module,     only: tem_comm_env_type
  use tem_varsys_module,       only: tem_varsys_type, tem_varSys_out,          &
    &                                tem_get_element_chunk,                    &
    &                                tem_get_point_chunk
  use tem_subtree_type_module, only: tem_subtree_type
  use treelmesh_module,        only: treelmesh_type
  use tem_time_module,         only: tem_time_type, tem_time_sim_stamp
  use tem_timeControl_module,  only: tem_timeControl_type, tem_timeControl_out
  use tem_shape_module,        only: tem_shape_type, tem_shape_out
  use tem_logging_module,      only: logUnit
  use tem_aux_module,          only: tem_open, tem_abort
  use tem_geometry_module,     only: tem_BaryOfId

  use aot_out_module,   only: aot_out_type, aot_out_val,                       &
    &                         aot_out_open, aot_out_close

  implicit none

  !> Description of the opened files for ascii output.
  type hvs_ascii_type
    !> File handle for the ascii file with the data.
    integer :: outunit

    !> Basename of the VTK files to write
    character(len=pathLen) :: basename

    !> Timestamp to construct the filename
    character(len=labelLen) :: timestamp

    !> spatial reduction active
    !! If reduction is active, only root process of subTree dumps the data
    logical :: isReduce

    !> reduction type which saves results from reduction
    type(tem_reduction_spatial_type), allocatable  :: redSpatial(:)

    !> number of elements that fit in the buffer
    integer :: chunkSize

    !> number of chunks per output
    integer :: nChunks
  end type hvs_ascii_type

  !> Description of the opened files for ascii output.
  type hvs_asciiSpatial_type
    !> File handle for the ascii file with the data.
    integer :: outunit

    !> Basename of the VTK files to write
    character(len=pathLen) :: basename

    !> Timestamp to construct the filename
    character(len=labelLen) :: timestamp

    !> number of elements that fit in the buffer
    integer :: chunkSize

    !> number of chunks per output
    integer :: nChunks
  end type hvs_asciiSpatial_type

contains


! ****************************************************************************** !
  !> Initialize ascii output format.
  !! initialize spatial reduction if reduction is active
  subroutine hvs_ascii_init(ascii, varSys, varPos, basename, globProc,         &
    &                       outProc, solver, geometry, nElems, glob_nElems,    &
    &                       timeControl, useGetPoint, nPoints, glob_nPoints,   &
    &                       nDofs)
    ! --------------------------------------------------------------------------!
    !> Ascii output file settings
    !! It must be intent inout since ascii%redSpatial
    !! are loaded in tem_load_trackingHeader
    type(hvs_ascii_type), intent(inout) :: ascii

    !> Description of the available variable system to get the given varnames
    !! from.
    type(tem_varSys_type), intent(in) :: varsys

    !> List of variable positions that should be written in the output.
    integer, intent(in) :: varpos(:)

    !> An extension to the output basename.
    character(len=*), intent(in) :: basename

    !> Global communicator type for global rank information
    type(tem_comm_env_type ), intent(in) :: globProc

    !> Process description to use.
    type(tem_comm_env_type), intent(in) :: outProc

    !> The number of dofs for each scalar variable of the equation system
    integer, intent(in) :: nDofs

    !> Number of elements to output by local process
    integer, intent(in) :: nElems

    !> Total number of elements across process on this output
    integer(kind=long_k), intent(in) :: glob_nElems

    !> Global solver information
    type(tem_solveHead_type ),intent(in) :: solver

    !> shape defined for this ascii output
    type(tem_shape_type), optional, intent(in) :: geometry(:)

    !> output timeControl
    type(tem_timeControl_type), optional, intent(in) :: timeControl

    !> if get_point is to be used to track the point
    logical, intent(in) :: useGetPoint

    !> Number of points to output by local process
    integer, intent(in) :: nPoints

    !> Total number of points across process on this output
    integer(kind=long_k), intent(in) :: glob_nPoints
    ! ----------------------------------------------------------------------!
    integer :: nVars, nScalars, chunkSize, nChunks
    ! Define limit for max number of tracking entities in ascii file
    integer, parameter :: nElemLimit_ascii = 50
    ! ----------------------------------------------------------------------!

    nVars = size(varPos)

    ! Check if there are too many entries for an ascii line
    if ( ( (nElems*nVars*nDofs > nElemLimit_ascii) &
      & .or. (nPoints*nVars > nElemLimit_ascii) )  &
      & .and. (.not. ascii%isReduce)  ) then
      write(logUnit(1),*)'Error in ascii output: '//trim(basename)
      write(logUnit(1),*)'Reduce is: ', ascii%isReduce
      if (useGetPoint) then
        write(logUnit(1),"(A,I0)")'Entities in output:', nPoints*nVars
      else
        write(logUnit(1),"(A,I0)")'Entities in output:', nElems*nVars*nDofs
      end if
      write(logUnit(1),"(A,I0)") 'Limit: ', nElemLimit_ascii
      write(logUnit(1),*)'Error: Too many entries for the tracking '//   &
        &            'entity as it is in ascii format'
      write(logUnit(1),*)'       All this information must be packed '// &
        &                'into one line.'
      write(logUnit(1),*)'Solution1: Use the Harvester format. '//       &
        &                 'It is way more efficient.'
      write(logUnit(1),*)'Solution2: reduce number of variables or '//   &
        &                'number of elements(segments)'
      call tem_abort()
    end if

    ! Compute nChunks, Abort if nChunks>1 asciiSpatial format
    ! and for ascii format abort if nChunk>1 and no reduction defined
    nScalars = sum(varSys%method%val(varPos(:))%nComponents)

    if (useGetPoint) then
      chunkSize = min(io_buffer_size/nScalars, nPoints)
    else
      chunkSize = min(io_buffer_size/(nScalars*nDofs), nElems)
    end if
    if ( (nElems > 0 .or. nPoints > 0) .and. (chunkSize == 0) ) then
      write(logUnit(0),*)'Error in ascii output: '//trim(basename)
      write(logUnit(0),*) 'The chosen io_buffer_size of ', io_buffer_size
      write(logUnit(0),*) 'is too small for outputting ', nScalars
      write(logUnit(0),*) 'scalar values'
      write(logUnit(0),*) 'Please increase the io_buffer_size to at &
        & least ', real(nScalars) / real(131072), ' MB!'
      call tem_abort()
    end if

    if (chunkSize>0) then
      if (useGetPoint) then
        nChunks = ceiling(real(nPoints, kind=rk)/real(chunkSize, kind=rk))
      else
        nChunks = ceiling(real(nElems, kind=rk)/real(chunkSize, kind=rk))
      end if
    else
      nChunks = 0
    end if

    ! abort if nChunk> 1 and reduction is not active
    if ( nChunks > 1 .and. (.not. ascii%isReduce) ) then
      write(logUnit(0),*)'Error in ascii output: '//trim(basename)
      write(logUnit(0),*)
      write(logUnit(0),*)'Number of chunks > 1'
      write(logUnit(0),*)'Solution: Define reduction of each variable to dump ',&
        &                'data in ascii format'
      call tem_abort()
    end if

    ascii%chunkSize = chunkSize
    ascii%nChunks   = nChunks
    ascii%basename  = trim(basename)

    ! write ascii header lua
    call hvs_ascii_write_header(out_format   = 'ascii',        &
      &                         basename     = trim(basename), &
      &                         varSys       = varSys,         &
      &                         varPos       = varPos,         &
      &                         globProc     = globProc,       &
      &                         outProc      = outProc,        &
      &                         nDofs        = nDofs,          &
      &                         solver       = solver,         &
      &                         geometry     = geometry,       &
      &                         timeControl  = timeControl,    &
      &                         glob_nElems  = glob_nElems,    &
      &                         useGetPoint  = useGetPoint,    &
      &                         glob_nPoints = glob_nPoints    )

    ! If reduction is active only root of this output
    ! dumps data else all process writes their own result file
    if ( (ascii%isReduce .and. outProc%rank == 0) .or. &
        & (.not. ascii%isReduce) ) then
      ! write header for result file
      call hvs_ascii_open( ascii   = ascii,   &
        &                  outProc = outProc, &
        &                  varSys  = varSys,  &
        &                  varPos  = varPos,  &
        &                  nDofs   = nDofs    )
    end if
  end subroutine hvs_ascii_init
! ****************************************************************************** !


! ****************************************************************************** !
  !> Write the header of the ascii output files
  !!
  !! This writes a header with detailed information into the ascii file for the
  !! defined tracking object it writes into the unit, which was opened before
  !! in \ref tem_init_tracker
  !!
  subroutine hvs_ascii_write_header( out_format, basename, varSys, varPos,     &
    &                                globProc, outProc, solver, geometry,      &
    &                                timeControl, nDofs, glob_nElems,          &
    &                                useGetPoint, glob_nPoints )
    ! ---------------------------------------------------------------------------
    !> Output format
    character(len=*), intent(in) :: out_format

    !> Basename for output file. tracking%prefix//tracking%label
    character(len=*), intent(in) :: basename

    !> solver-provided variable systems
    type(tem_varSys_type), intent(in)       :: varSys

    !> List of variable positions that should be written in the output.
    integer, intent(in) :: varpos(:)

    !> Global communicator type for global rank information
    type(tem_comm_env_type ), intent(in) :: globProc

    !> Output communicator type from tracking
    type(tem_comm_env_type ), intent(in) :: outProc

    !> Global solver information
    type(tem_solveHead_type ),intent(in) :: solver

    !> shape defined for this ascii output
    type(tem_shape_type), optional, intent(in) :: geometry(:)

    !> output timeControl
    type(tem_timeControl_type), optional, intent(in) :: timeControl

    !> The number of dofs for each scalar variable of the equation system
    integer, intent(in) :: nDofs

    !> Total number of elements across process on this output
    integer(kind=long_k), intent(in) :: glob_nElems

    !> if get_point is to be used to track the point
    logical, intent(in) :: useGetPoint

    !> Total number of points across process on this output
    integer(kind=long_k), intent(in) :: glob_nPoints
    ! ---------------------------------------------------------------------------
    character(len=pathLen) :: logName
    type(aot_out_type) :: conf !< aotus lua state to write output
    character(len=pathLen)  :: buffer
    character(len=labelLen) :: solverTag
    ! ---------------------------------------------------------------------------
    if ( outProc%rank == 0 ) then
      write(logName,'(a)') trim(basename)//'.lua'

      write( buffer, '(a)' ) trim(basename)//'_p*'

      solverTag = tem_solverTag( solver )

      !! Write the header only on the root process
      !open up the mesh header lua file to dump the stuff using aotus library
      call aot_out_open( conf, logName )

      ! output format
      call aot_out_val( put_conf = conf,               &
        &               vname    = 'format',           &
        &               val      = trim(out_format)    )
      call aot_out_val( put_conf = conf,               &
        &               vname    = 'solver',           &
        &               val      = trim(solverTag)     )
      call aot_out_val( put_conf = conf,               &
        &               vname    = 'simname',          &
        &               val      = trim(solver%simname))
      call aot_out_val( put_conf = conf,               &
        &               vname    = 'basename',         &
        &               val      = trim(basename)      )
      call aot_out_val( put_conf = conf,               &
        &               vname    = 'glob_rank',        &
        &               val      = globProc%rank       )
      call aot_out_val( put_conf = conf,               &
        &               vname    = 'glob_nprocs',      &
        &               val      = globProc%comm_size  )
      call aot_out_val( put_conf = conf,               &
        &               vname    = 'sub_rank',         &
        &               val      = outProc%rank        )
      call aot_out_val( put_conf = conf,               &
        &               vname    = 'sub_nprocs',       &
        &               val      = outProc%comm_size   )
      call aot_out_val( put_conf = conf,               &
        &               vname    = 'resultfile',       &
        &               val      = trim(buffer)        )
      call aot_out_val( put_conf = conf,               &
        &               vname    = 'nDofs',            &
        &               val      = nDofs               )
      call aot_out_val( put_conf = conf,               &
        &               vname    = 'nElems',           &
        &               val      = glob_nElems         )
      if (useGetPoint) then
        call aot_out_val( put_conf = conf,               &
          &               vname    = 'use_get_point',    &
          &               val      = useGetPoint         )
        call aot_out_val( put_conf = conf,               &
          &               vname    = 'nPoints',          &
          &               val      = glob_nPoints        )
      end if
      ! write timeControl info
      if (present(timeControl)) &
        & call tem_timeControl_out( timeControl, conf )
      ! write shapes
      if (present(geometry)) call tem_shape_out( geometry, conf )
      ! write variable systems
      call tem_varSys_out( varSys, conf, varPos )
      ! close aotus file
      call aot_out_close( conf )
    end if

  end subroutine hvs_ascii_write_header
! ****************************************************************************** !


! ****************************************************************************** !
  !> open the ascii transient output unit
  !!
  subroutine hvs_ascii_open( ascii, outProc, varSys, varPos, nDofs )
    ! ---------------------------------------------------------------------------
    !> ascii output type
    type(hvs_ascii_type ), intent(inout) :: ascii
    !> Parallel environment to use for  the output.
    type(tem_comm_env_type ), intent(in)    :: outProc
    !> solver-provided variable systems
    type(tem_varSys_type), intent(in)       :: varSys
    !> Position of variables to dump in varSys
    integer, intent(in) :: varPos(:)
    !> The number of dofs for each scalar variable of the equation system
    integer, intent(in) :: nDofs
    ! ---------------------------------------------------------------------------
    character(len=labelLen)  :: logName
    character(len=1500) :: buffer
    logical             :: nUnitOpened
    logical             :: file_exists
    integer             :: UnitNumber
    ! ---------------------------------------------------------------------------
    file_exists = .false.
    nUnitOpened = .false.

    write(buffer,'(i5.5)') outProc%rank
     ! Include the rank name in the file name to avoid overwriting of files
     ! by different processes
    write(logName,'(a)') trim(ascii%basename)//'_p'//trim(buffer)//'.res'

    ! check if the file exists
    inquire( file = trim(logName), exist = file_exists )
    if( file_exists )then
      ! \todo: SZ: before simply appending check if header lua exists, open
      !            it, read the varSys and compare it with the one requested
      !            by this tracking object.
      !            if equal     -> appened file_exists
      !            if not equal -> put a new behind the filename and open
      !                            a new file

      ! in case the file exists, check wether it is already opened somewhere
      ! else (dyn load balancing)
      inquire( file=trim(logName), opened=nUnitOpened, number=UnitNumber )
      if (nUNitOpened) then
        ! if it is opened, use the corresponding unit
        ascii%outunit = UnitNumber
      else
        ! get a new unit and reopen the file
        call tem_open( newunit  = ascii%outunit, &
          &            file     = trim(logName), &
          &            action   = 'WRITE',       &
          &            position = 'APPEND',      &
          &            status   = 'OLD'          )
      end if
    else
      ! in case the file does not exist, get a new unit and open a new file
      ! and write the header
      call tem_open( newunit  = ascii%outunit, &
        &            file     = trim(logName), &
        &            action   = 'WRITE',       &
        &            status   = 'NEW'          )
      write ( ascii%outunit , '(a,i7)' ) '# Rank of the process: ',         &
        &                                 outProc%rank
    end if

    if (.not. file_exists .and. .not. nUnitOpened) then
      write( buffer, '(a)' ) '#'
      write( buffer, '(a,a23)' ) trim(buffer),'time'
      buffer = trim(buffer)                                    &
        &    //trim(getHeader( varSys, varPos, nDofs, ascii%isReduce ))
      ! ... and write the header row to the file
      write ( ascii%outunit , '(a)' ) trim(buffer)
    end if

  end subroutine hvs_ascii_open
! ****************************************************************************** !


! ****************************************************************************** !
  !> Write single log for the right scheme into its ascii file. This routine dumps
  !! the element data
  subroutine hvs_ascii_dump_elem_data( ascii, outProc, varPos, varSys, mesh,  &
    &                                   time, subTree, nDofs )
    ! ---------------------------------------------------------------------------
    !> ascii file to write data to.
    type(hvs_ascii_type), intent(inout) :: ascii

    !> Parallel environment to use for  the output.
    type(tem_comm_env_type ), intent(in)    :: outProc

    !> Position of the variable to write
    integer, intent(in) :: varpos(:)

    !> Description of the available variable system to get the given varnames
    !! from.
    type(tem_varSys_type), intent(in) :: varsys

    !> Mesh to write the data on.
    type(treelmesh_type), intent(in) :: mesh

    !> Point in time to use for this data.
    !!
    !! Can be important for space-time function evaluations.
    type(tem_time_type), intent(in) :: time

    !> Optional restriction of the elements to output.
    type(tem_subtree_type), optional, intent(in) :: subtree

    !> The number of dofs for each scalar variable of the equation system
    integer, intent(in) :: nDofs
    ! ---------------------------------------------------------------------------
    integer :: nVars, nElems, nScalars, elemOff, nChunkElems
    integer :: nDofs_out
    integer :: nScalars_out
    integer :: iElem, iChunk, iScalar, iDof
    integer :: buf_start, buf_end
    real(kind=rk), allocatable :: res(:)
    integer, allocatable :: elemPos(:)
    character(len=4000) :: log_output ! output buffer
    ! ---------------------------------------------------------------------------
    allocate(res(io_buffer_size))

    ! Number of variables to dump
    nVars = size(varPos)

    ! Number of scalars in current output
    nScalars = sum(varSys%method%val(varPos(:))%nComponents)

    if (present(subTree)) then
      nElems = subTree%nElems
    else
      nElems = mesh%nElems
    end if

    if (ascii%isReduce) then
      ! open spatial reduction
      call tem_reduction_spatial_open( me     = ascii%redSpatial, &
        &                              varSys = varSys,           &
        &                              varPos = varPos(:nVars)    )
    end if

    ! allocate elemPos to size of chunkSize
    allocate(elemPos(ascii%chunkSize))

    ! Process all chunks to derive the quantities defined in the tracking object
    do iChunk = 1, ascii%nChunks
      ! Number of elements read so far in previous chunks.
      elemOff = ((iChunk-1)*ascii%chunkSize)

      ! number of elements written to THIS chunk
      nChunkElems = min(ascii%chunkSize, nElems-elemOff)

      ! Compute the element lower and upper bound for the current chunk
      buf_start = elemOff + 1
      buf_end = elemOff + nChunkElems

      if (present(subTree)) then
        elemPos(1:nChunkElems) = subTree%map2Global(buf_start:buf_end)
      else
        elemPos(1:nChunkElems) = (/ (iElem, iElem=buf_start, buf_end) /)
      end if

      ! evaluate all variables on current chunk
      call tem_get_element_chunk(varSys  = varSys,                 &
        &                        varPos  = varPos,                 &
        &                        elemPos = elemPos(1:nChunkElems), &
        &                        time    = time,                   &
        &                        tree    = mesh,                   &
        &                        nElems  = nChunkElems,            &
        &                        nDofs   = nDofs,                  &
        &                        res     = res                     )

      ! preform spatial reduction
      if (ascii%isReduce) then
        call tem_reduction_spatial_append(                                 &
          &                        me     = ascii%redSpatial,              &
          &                        chunk  = res,                           &
          &                        nElems = nChunkElems,                   &
          &                        treeID = mesh%treeID(                   &
          &                                      elemPos(1:nChunkElems) ), &
          &                        tree   = mesh,                          &
          &                        varSys = varSys,                        &
          &                        nDofs  = nDofs,                         &
          &                        varPos = varPos                         )
      end if
    end do ! iChunk

    ndofs_out = ndofs
    nscalars_out = nscalars

    ! If a reduction is present in the tracking, then the output is
    ! changed completely to only the reduced values
    ! For each variable, a reduction has to be assigned
    if( ascii%isReduce ) then

      ! Perform the global reduction on the data which was appended
      ! inside trackVariable by tem_reduction_spatial_append
      call tem_reduction_spatial_close( me   = ascii%redSpatial, &
        &                               proc = outProc           )

      ! Re-assign the chunk here to the stuff which was produced in the
      ! reduction operation
      call tem_reduction_spatial_toChunk(me          = ascii%redSpatial, &
        &                                chunk       = res,              &
        &                                nChunkElems = nChunkElems       )

      ndofs_out = 1
      nscalars_out = sum(ascii%redSpatial(:)%nComponents)
    end if

    ! If reduction is active only root of this output
    ! dumps data else all process writes their own result file
    if ( (ascii%isReduce .and. outProc%rank == 0) .or. &
      & (.not. ascii%isReduce) ) then

      ! First assemble the complete row consisting of the time ...
      write( log_output, '(e24.16e3)' ) time%sim

      ! add all the scalars entries of the variable systems for each elements

      do iElem = 1, nChunkElems
        do iDof = 1, nDofs_out
          do iScalar = 1, nScalars_out
            write( log_output, '(a,1x,e24.16e3)' ) trim(log_output),          &
              &    res( ((iElem-1)*nDofs_out+ (iDof-1))*nScalars_out + iScalar )
          end do
        end do
      end do
      ! then write into the ascii file
      write ( ascii%outunit , '(a)' ) trim(log_output)
    end if

    deallocate(elemPos)
    deallocate(res)

  end subroutine hvs_ascii_dump_elem_data
! ****************************************************************************** !


! ****************************************************************************** !
  !> Write single log for the right scheme into its ascii file. This routine
  !! calls the get_point routine and dumps the exact point data for the point
  !! specified in the tracking table
  subroutine hvs_ascii_dump_point_data( ascii, outProc, varPos, varSys, mesh,  &
    &                                   time, subTree )
    ! ---------------------------------------------------------------------------
    !> ascii file to write data to.
    type(hvs_ascii_type), intent(inout) :: ascii

    !> Parallel environment to use for  the output.
    type(tem_comm_env_type ), intent(in)    :: outProc

    !> Position of the variable to write
    integer, intent(in) :: varpos(:)

    !> Description of the available variable system to get the given varnames
    !! from.
    type(tem_varSys_type), intent(in) :: varsys

    !> Mesh to write the data on.
    type(treelmesh_type), intent(in) :: mesh

    !> Point in time to use for this data.
    !!
    !! Can be important for space-time function evaluations.
    type(tem_time_type), intent(in) :: time

    !> Optional restriction of the points to output.
    !! Contains array of points passed in the config  to output.
    type(tem_subtree_type), optional, intent(in) :: subtree
    ! ---------------------------------------------------------------------------
    integer :: nVars, nPoints, nScalars, pointsOff, nChunkPoints
    integer :: iPoint, iChunk, iScalar, counter
    integer :: buf_start, buf_end
    real(kind=rk), allocatable :: res(:)
    character(len=4000) :: log_output ! output buffer
    real(kind=rk), allocatable :: points(:,:)
    ! ---------------------------------------------------------------------------
    if (present(subTree)) then
      nPoints = subTree%nPoints
    else
      nPoints = mesh%nElems
    end if

    allocate(res(io_buffer_size))

    ! Number of variables to dump
    nVars = size(varPos)

    ! Number of scalars in current output
    nScalars = sum(varSys%method%val(varPos(:))%nComponents)

    if (ascii%isReduce) then
      ! open spatial reduction
      call tem_reduction_spatial_open( me     = ascii%redSpatial, &
        &                              varSys = varSys,           &
        &                              varPos = varPos(:nVars)    )
    end if

    ! allocate points to size of chunkSize
    allocate(points(ascii%chunkSize,3))

    ! Process all chunks to derive the quantities defined in the tracking object
    do iChunk = 1, ascii%nChunks
      ! Number of points read so far in previous chunks.
      pointsOff = ((iChunk-1)*ascii%chunkSize)

      ! number of points written to THIS chunk
      nChunkPoints = min(ascii%chunkSize, nPoints-pointsOff)

      ! Compute the points lower and upper bound for the current chunk
      buf_start = pointsOff + 1
      buf_end = pointsOff + nChunkPoints

      if (present(subTree)) then
        points(1:nChunkPoints,:) = subTree%points(buf_start:buf_end,:)
      else
        counter = 0
        do iPoint = buf_start, buf_end
          counter = counter + 1
          points(counter, :) = tem_BaryOfId( mesh,               &
            &                                mesh%treeID(iPoint) )
        end do
      end if

      ! evaluate all variables on current chunk
      call tem_get_point_chunk(  varSys  = varSys,                   &
        &                        varPos  = varPos,                   &
        &                        point   = points(1:nChunkPoints,:), &
        &                        time    = time,                     &
        &                        tree    = mesh,                     &
        &                        nPnts   = nChunkPoints,             &
        &                        res     = res                       )

      ! preform spatial reduction
      if (ascii%isReduce) then
       ! NA : The reduction for the get points needs to be thought
       ! through as we would like at one point to have all the segment paoints
       ! stored and reduction applied. In that case, l2 norm and weighted sum
       ! requires the treeID information. Maybe we should write a new routine to
       ! do that tem_reduction_spatial_append_point
       ! KM: changed treeID argument optional because for point values
       ! no need to do volume reduction
       call tem_reduction_spatial_append( me     = ascii%redSpatial, &
         &                                chunk  = res,              &
         &                                nElems = nChunkPoints,     &
         &                                tree   = mesh,             &
         &                                varSys = varSys,           &
         &                                varPos = varPos            )
      end if
    end do ! iChunk

    ! If a reduction is present in the tracking, then the output is
    ! changed completely to only the reduced values
    ! For each variable, a reduction has to be assigned
    if( ascii%isReduce ) then

      ! Perform the global reduction on the data which was appended
      ! inside trackVariable by tem_reduction_spatial_append
      call tem_reduction_spatial_close( me   = ascii%redSpatial, &
        &                               proc = outProc           )

      ! Re-assign the chunk here to the stuff which was produced in the
      ! reduction operation
      call tem_reduction_spatial_toChunk(me          = ascii%redSpatial, &
        &                                chunk       = res,              &
        &                                nChunkElems = nChunkPoints      )
    end if

    ! If reduction is active only root of this output
    ! dumps data else all process writes their own result file
    if ( (ascii%isReduce .and. outProc%rank == 0) .or. &
      & (.not. ascii%isReduce) ) then

      ! First assemble the complete row consisting of the time ...
      write( log_output, '(e24.16e3)' ) time%sim

      ! add all the scalars entries of the variable systems for each elements

      do iPoint = 1, nChunkPoints
          do iScalar = 1, nScalars
            write( log_output, '(a,1x,e24.16e3)' ) trim(log_output),          &
              &    res( (iPoint-1)*nScalars +  iScalar )
        end do
      end do
      ! then write into the ascii file
      write ( ascii%outunit , '(a)' ) trim(log_output)
    end if

    deallocate(points)
    deallocate(res)

  end subroutine hvs_ascii_dump_point_data
! ****************************************************************************** !


! ****************************************************************************** !
  !> Initialize asciiSpatial output format.
  !! initialize reduction if reduction is active
  subroutine hvs_asciiSpatial_init(asciiSpatial, varSys, varPos, basename,     &
    &                              globProc, outProc, solver, geometry, nDofs, &
    &                              nElems, glob_nElems, useGetPoint, nPoints,  &
    &                              glob_nPoints, timeControl                   )
    ! --------------------------------------------------------------------------!
    !> AsciiSpatial output file settings
    type(hvs_asciiSpatial_type), intent(inout) :: asciiSpatial

    !> Description of the available variable system to get the given varnames
    !! from.
    type(tem_varSys_type), intent(in) :: varsys

    !> List of variable positions that should be written in the output.
    integer, intent(in) :: varpos(:)

    !> An extension to the output basename.
    character(len=*), intent(in) :: basename

    !> Global communicator type for global rank information
    type(tem_comm_env_type ), intent(in) :: globProc

    !> Process description to use.
    type(tem_comm_env_type), intent(in) :: outProc

    !> The number of dofs for each scalar variable of the equation system
    integer, intent(in) :: nDofs

    !> Number of elements to output by local process
    integer, intent(in) :: nElems

    !> Total number of elements across process on this output
    integer(kind=long_k), intent(in) :: glob_nElems

    !> Global solver information
    type(tem_solveHead_type ),intent(in) :: solver

    !> shape defined for this ascii output
    type(tem_shape_type), optional, intent(in) :: geometry(:)

    !> output timeControl
    type(tem_timeControl_type), optional, intent(in) :: timeControl

    !> if get_point is to be used to track the point
    logical, intent(in) :: useGetPoint

    !> Number of points to output by local process
    integer, intent(in) :: nPoints

    !> Total number of points across process on this output
    integer(kind=long_k), intent(in) :: glob_nPoints
    ! ----------------------------------------------------------------------!
    integer :: nScalars, chunkSize, nChunks
    ! ----------------------------------------------------------------------!

    ! Compute nChunks, Abort if nChunks>1 asciiSpatial format
    ! and for ascii format abort if nChunk>1 and no reduction defined
    nScalars = sum(varSys%method%val(varPos(:))%nComponents)

    if (useGetPoint) then
      chunkSize = min(io_buffer_size/nScalars, nPoints)
    else
      chunkSize = min(io_buffer_size/(nScalars*nDofs), nElems)
    end if
    if ( (nElems > 0) .and. (chunkSize == 0) ) then
      write(logUnit(0),*)'Error in asciiSpatial output: '//trim(basename)
      write(logUnit(0),*) 'The chosen io_buffer_size of ', io_buffer_size
      write(logUnit(0),*) 'is too small for outputting ', nScalars
      write(logUnit(0),*) 'scalar values'
      write(logUnit(0),*) 'Please increase the io_buffer_size to at &
        & least ', real(nScalars*nDofs) / real(131072), ' MB!'
      call tem_abort()
    end if

    if (chunkSize>0) then
      if (useGetPoint) then
        nChunks = ceiling(real(nPoints, kind=rk)/real(chunkSize, kind=rk))
      else
        nChunks = ceiling(real(nElems, kind=rk)/real(chunkSize, kind=rk))
      end if
    else
      nChunks = 0
    end if

    asciiSpatial%chunkSize = chunkSize
    asciiSpatial%nChunks = nChunks

    asciiSpatial%basename = trim(basename)
    ! write ascii header lua
    call hvs_ascii_write_header(out_format   = 'asciispatial', &
      &                         basename     = trim(basename), &
      &                         varSys       = varSys,         &
      &                         varPos       = varPos,         &
      &                         globProc     = globProc,       &
      &                         outProc      = outProc,        &
      &                         nDofs        = nDofs,          &
      &                         solver       = solver,         &
      &                         geometry     = geometry,       &
      &                         timeControl  = timeControl,    &
      &                         glob_nElems  = glob_nElems,    &
      &                         useGetPoint  = useGetPoint,    &
      &                         glob_nPoints = glob_nPoints    )

  end subroutine hvs_asciiSpatial_init
! ****************************************************************************** !


! ****************************************************************************** !
  !> Open the output file for AsciiSpatial format.
  !!
  !! Each time this routine is called, a new file is written
  !! Filename: {tracking_folder}{tracking_label}_spatial_{proc}_{timestamp}.res
  !!     e.g.: tracking/lineProbe_spatial_00001_01_01378.1.res
  !! Each process open its own files
  !!
  subroutine hvs_asciiSpatial_open( asciiSpatial, outProc, time, varSys,      &
    &                               varPos, nDofs )
    ! ---------------------------------------------------------------------------
    !> asciiSpatial file output
    type(hvs_asciiSpatial_type ), intent(inout) :: asciiSpatial
    !> Parallel environment to use for  the output.
    type(tem_comm_env_type ), intent(in) :: outProc
    !> current global time
    type(tem_time_type ), intent(in) :: time
    !> solver-provided variable systems
    type(tem_varSys_type), intent(in) :: varSys
    !> Position of variables to dump in varSys
    integer, intent(in) :: varPos(:)
    !> The number of dofs for each scalar variable of the equation system
    integer, intent(in) :: nDofs
    ! ---------------------------------------------------------------------------
    integer :: nScalars
    character(len=pathLen)  :: filename  ! ascii file name
    character(len=1024)  :: buffer
    ! ---------------------------------------------------------------------------
    nScalars = sum(varSys%method%val(varPos(:))%nComponents)

    ! Write the rank into the string buffer
    write(buffer,'(i5.5)') outProc%rank

    asciiSpatial%timestamp = tem_time_sim_stamp( time )

    ! Generate ascii file name
    write(filename,'(a)') trim(asciiSpatial%basename)//'_p'//trim(buffer)      &
      &                //'_t'//trim(asciiSpatial%timestamp)//'.res'

    call tem_open( newunit = asciiSpatial%outunit, &
      &            file    = trim(filename),       &
      &            recl    = 1024,                 &
      &            status  = 'replace'             )

    write ( asciiSpatial%outUnit , '(a,i7)' ) '# Rank of the process: ',       &
        &                                     outProc%rank
    write ( buffer, '(a1, 3(1x,a24))' ) '#', 'coordX', 'coordY', 'coordZ'
    buffer = trim(buffer)                                    &
      &    //trim(getHeader( varSys, varPos, nDofs, .false. ))
    ! ... and write the header row to the file
    write ( asciiSpatial%outunit , '(a)' ) trim(buffer)

  end subroutine hvs_asciiSpatial_open
! ****************************************************************************** !


! ****************************************************************************** !
  !> Write a spatial representation for elements into an ascii tracking file
  !!
  !! Each time this routine is called, a new file is written
  !! Filename: {tracking_folder}{tracking_label}_spatial_{proc}_{timestamp}.res
  !!     e.g.: tracking/lineProbe_spatial_00001_01_01378.1.res
  !! Each process writes its own files
  !!
  subroutine hvs_asciiSpatial_dump_elem_data( asciiSpatial, varPos, varSys,    &
    &                        bary, mesh, subTree, time, nDofs )
    ! ---------------------------------------------------------------------------
    !> The file description to open
    type(hvs_asciiSpatial_type), intent(inout) :: asciiSpatial

    !> solver-provided variable systems
    type(tem_varSys_type), intent(in)       :: varSys

    !> Positions of the variables to write
    integer, intent(in) :: varpos(:)

    !> Barycenter of elements
    real(kind=rk), intent(in) :: bary(:,:)

    !> Mesh to write the data on.
    type(treelmesh_type), intent(in) :: mesh

    !> Optional restriction of the elements to output.
    type(tem_subtree_type), optional, intent(in) :: subtree

    !> current global time
    type(tem_time_type ), intent(in)  :: time

    !> The number of dofs for each scalar variable of the equation system
    integer, intent(in) :: nDofs
    ! ---------------------------------------------------------------------------
    integer :: nVars, nElems, nScalars, elemOff, nChunkElems, elemSize
    integer :: iElem, iChunk, iScalar, iDof
    integer :: buf_start, buf_end
    real(kind=rk), allocatable :: res(:)
    integer, allocatable :: elemPos(:)
    character(len=1024)  :: buffer
    ! ---------------------------------------------------------------------------
    allocate(res(io_buffer_size))

    ! Number of variables to dump
    nVars = size(varPos)

    ! Number of scalars in current output
    nScalars = sum(varSys%method%val(varPos(:))%nComponents)

    ! Size of a single element
    elemsize = nScalars*nDofs

    if (present(subTree)) then
      nElems = subTree%nElems
    else
      nElems = mesh%nElems
    end if

    ! allocate elemPos to size of chunkSize
    allocate(elemPos(asciiSpatial%chunkSize))

    ! Process all chunks to derive the quantities defined in the tracking object
    do iChunk = 1, asciiSpatial%nChunks
      ! Number of elements read so far in previous chunks.
      elemOff = ((iChunk-1)*asciiSpatial%chunkSize)

      ! number of elements written to THIS chunk
      nChunkElems = min(asciiSpatial%chunkSize, nElems-elemOff)

      ! Compute the element lower and upper bound for the current chunk
      buf_start = elemOff + 1
      buf_end = elemOff + nChunkElems

      if (present(subTree)) then
        elemPos(1:nChunkElems) = subTree%map2Global(buf_start:buf_end)
      else
        elemPos(1:nChunkElems) = (/ (iElem, iElem=buf_start, buf_end) /)
      end if

      ! evaluate all variables on current chunk
      call tem_get_element_chunk(varSys  = varSys,                 &
        &                        varPos  = varPos,                 &
        &                        elemPos = elemPos(1:nChunkElems), &
        &                        time    = time,                   &
        &                        tree    = mesh,                   &
        &                        nElems  = nChunkElems,            &
        &                        nDofs   = nDofs,                  &
        &                        res     = res                     )

      ! Then gather contents into buffer, and write buffer to file
      buffer = ''
      do iElem = 1, nChunkElems
        ! write coordinates to buffer
        write( buffer, '(3(1x,e24.16e3))' ) bary(elemOff+iElem, 1:3)

        ! append values in chuck to buffer
        do iDof = 1, nDofs
          do iScalar = 1, nScalars
            write( buffer, '(a,1x,e24.16e3)' ) trim(buffer),           &
              &  res( (iElem-1)*elemSize + (iDof-1)*nScalars + iScalar )
          end do
        end do

        ! write buffer to file
        write ( asciiSpatial%outUnit , '(a)' ) trim(buffer)

      end do !nChunkElems

    end do ! iChunk

    deallocate(elemPos)
    deallocate(res)

  end subroutine hvs_asciiSpatial_dump_elem_data
! ****************************************************************************** !


! ****************************************************************************** !
  !> Write a spatial representation for list of points into an ascii
  !! tracking file
  !!
  !! Each time this routine is called, a new file is written
  !! Filename: {tracking_folder}{tracking_label}_spatial_{proc}_{timestamp}.res
  !!     e.g.: tracking/lineProbe_spatial_00001_01_01378.1.res
  !! Each process writes its own files
  !!
  subroutine hvs_asciiSpatial_dump_point_data( asciiSpatial, varPos, varSys, &
    &                             bary, mesh, subTree, time )
    ! ---------------------------------------------------------------------------
    !> The file description to open
    type(hvs_asciiSpatial_type), intent(inout) :: asciiSpatial

    !> solver-provided variable systems
    type(tem_varSys_type), intent(in)       :: varSys

    !> Positions of the variables to write
    integer, intent(in) :: varpos(:)

    !> Barycenter of elements
    real(kind=rk), intent(in) :: bary(:,:)

    !> Mesh to write the data on.
    type(treelmesh_type), intent(in) :: mesh

    !> Optional restriction of the elements to output.
    type(tem_subtree_type), optional, intent(in) :: subtree

    !> current global time
    type(tem_time_type ), intent(in)  :: time
    ! ---------------------------------------------------------------------------
    integer :: nVars, nPoints, nScalars, pointsOff, nChunkPoints
    integer :: iPoint, iChunk, iScalar, counter
    integer :: buf_start, buf_end
    real(kind=rk), allocatable :: res(:)
    real(kind=rk), allocatable :: points(:,:)
    character(len=1024)  :: buffer
    ! ---------------------------------------------------------------------------
    allocate(res(io_buffer_size))

    ! Number of variables to dump
    nVars = size(varPos)

    ! Number of scalars in current output
    nScalars = sum(varSys%method%val(varPos(:))%nComponents)

    if (present(subTree)) then
      nPoints = subTree%nPoints
    else
      nPoints = mesh%nElems
    end if

    ! allocate points to size of chunkSize
    allocate(points(asciiSpatial%chunkSize,3))

    ! Process all chunks to derive the quantities defined in the tracking object
    do iChunk = 1, asciiSpatial%nChunks
      ! Number of points read so far in previous chunks.
      pointsOff = ((iChunk-1)*asciiSpatial%chunkSize)

      ! number of points written to THIS chunk
      nChunkPoints = min(asciiSpatial%chunkSize, nPoints-pointsOff)

      ! Compute the points lower and upper bound for the current chunk
      buf_start = pointsOff + 1
      buf_end = pointsOff + nChunkPoints


      if (present(subTree)) then
        points(1:nChunkPoints,:) = subTree%points(buf_start:buf_end,:)
      else
        counter = 0
        do iPoint = buf_start, buf_end
          counter = counter + 1
          points(counter, :) = tem_BaryOfId( mesh,               &
            &                                mesh%treeID(iPoint) )
        end do
      end if

      ! evaluate all variables on current chunk
      call tem_get_point_chunk(  varSys  = varSys,                   &
        &                        varPos  = varPos,                   &
        &                        point   = points(1:nChunkPoints,:), &
        &                        time    = time,                     &
        &                        tree    = mesh,                     &
        &                        nPnts   = nChunkPoints,             &
        &                        res     = res                       )

      ! Then gather contents into buffer, and write buffer to file
      buffer = ''
      do iPoint = 1, nChunkPoints
        ! write coordinates to buffer
        write( buffer, '(3(1x,e24.16e3))' ) bary(pointsOff+iPoint, 1:3)

        ! append values in chuck to buffer
        do iScalar = 1, nScalars
          write( buffer, '(a,1x,e24.16e3)' ) trim(buffer),           &
            &  res( (iPoint-1)*nScalars + iScalar )
        end do

        ! write buffer to file
        write ( asciiSpatial%outUnit , '(a)' ) trim(buffer)

      end do !nChunkElems

    end do ! iChunk

    deallocate(points)
    deallocate(res)

  end subroutine hvs_asciiSpatial_dump_point_data
! ****************************************************************************** !


! ****************************************************************************** !
  !> close the ascii output unit
  !!
  subroutine hvs_ascii_close(ascii)
    ! ---------------------------------------------------------------------------
    !> ascii output
    type(hvs_ascii_type ), intent(inout) :: ascii
    ! ---------------------------------------------------------------------------
    !local variable
    logical :: nUnitOpened
    ! ---------------------------------------------------------------------------

    ! check output unit is open
    if ( ascii%outunit >=0 ) then
      inquire( unit = ascii%outunit, opened = nUnitOpened )
      if ( nUnitOpened ) close( ascii%outunit )
    end if

  end subroutine hvs_ascii_close
! ****************************************************************************** !


! ****************************************************************************** !
  !> close the asciiSpatial output unit
  !!
  subroutine hvs_asciiSpatial_close(asciiSpatial)
    ! ---------------------------------------------------------------------------
    !> asciiSpatial output
    type(hvs_asciiSpatial_type ), intent(inout) :: asciiSpatial
    ! ---------------------------------------------------------------------------
    !local variable
    logical :: nUnitOpened
    ! ---------------------------------------------------------------------------

    ! check output unit is open
    if ( asciiSpatial%outunit >=0 ) then
      inquire( unit = asciiSpatial%outunit, opened = nUnitOpened )
      if ( nUnitOpened ) close( asciiSpatial%outunit )
    end if

  end subroutine hvs_asciiSpatial_close
! ****************************************************************************** !

! ****************************************************************************** !
  !> This routine write variable labels into buffer and this buffer is written
  !! into the second line of ascii result file for spatial (asciiSpatial) and
  !! transient (ascii) tracking.
  !!
  function getHeader( varSys, varPos, nDofs, isReduce ) result(buffer)
    ! ---------------------------------------------------------------------------
    !> solver-provided variable systems
    type(tem_varSys_type), intent(in) :: varSys
    !> Position of variables to dump in varSys
    integer, intent(in) :: varPos(:)
    !> The number of dofs for each scalar variable of the equation system
    integer, intent(in) :: nDofs
    !> is spatial Reduction active
    logical, intent(in) :: isReduce
    !> buffer containing first column of the header
    character(len=OutLen)         :: buffer
    ! ---------------------------------------------------------------------------
    !local variable
    character(len=pathLen)  :: buffer2
    integer :: iVar, iComp, iDof
    character(len=4)  :: reduceChar
    ! ---------------------------------------------------------------------------
    buffer = ''

    ! If reduction is applied to this variable, append a ending to its label
    if( isReduce )  then
      reduceChar = '_red'
    else
      reduceChar = ''
    end if

    do iVar = 1, size(varPos)
      if (nDofs == 1 .or. isReduce) then
        if ( varSys%method%val( varPos(iVar) )%nComponents > 1 ) then
          ! when variable has more that one componenets, append icomp to the end
          ! of its label
          do iComp = 1, varSys%method%val( varPos(iVar) )%nComponents
            write( buffer2, '(a,i2.2)' )                                       &
              &  trim(varSys%varName%val(varPos(iVar)))//trim(reduceChar)//'_',&
              &  iComp
            if ( len(trim(buffer)) + 25 > OutLen ) then
              call tem_abort( errorMsg =                                  &
                &  'max char length for header exceeded, use vtk or '//   &
                &  'harvester for tracking. Restarting will result in '// &
                &  'res file without header.'                             )
            end if
            write(buffer, '(a,1x,a24)') trim(buffer), trim(buffer2)
          end do
        else
          write(buffer, '(a,1x,a24)') trim(buffer), &
            &     trim(varSys%varName%val(varPos(iVar)))//trim(reduceChar)
        end if
      else
        do iDof = 1, nDofs
          if ( varSys%method%val( varPos(iVar) )%nComponents > 1 ) then
            ! when variable has more that one componenets, append icomp to the end
            ! of its label
            do iComp = 1, varSys%method%val( varPos(iVar) )%nComponents
              write( buffer2, '(a,i3.3,a,i2.2)' )                           &
                &  trim(varSys%varName%val(varPos(iVar)))//trim(reduceChar) &
                &  //'_d', iDof, '_c', iComp
              write(buffer, '(a,1x,a24)') trim(buffer), trim(buffer2)
            end do
          else
            write(buffer, '(a,1x,a24,i3.3)') trim(buffer),                   &
              &     trim(varSys%varName%val(varPos(iVar)))//trim(reduceChar) &
              &  //'_d', iDof
          end if

        end do
      end if
    end do

  end function getHeader
! ****************************************************************************** !


end module hvs_ascii_module
