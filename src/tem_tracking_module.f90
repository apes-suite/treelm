! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011 Konstantin Kleinheinz <k.kleinheinz@grs-sim.de>
! Copyright (c) 2011-2012, 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2016,2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011 Gaurang Phadke <g.phadke@grs-sim.de>
! Copyright (c) 2011 Laura Didinger <l.didinger@grs-sim.de>
! Copyright (c) 2011 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2011-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012, 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2013, 2015 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013, 2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014, 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Rishabh Chandola <rishabh.chandola@student.uni-siegen.de>
! Copyright (c) 2014-2015 Langhammer Kay <kay.langhammer@student.uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2018 Robin Weihe <robin.weihe@student.uni-siegen.de>
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
!> Tracking enables the simulation to extract values of variables from within
!! arbitrary parts of the mesh.
!!
!! There maybe multiple tracking objects defined and each one needs to define
!! the location to track (by a `shape` subtable), the variables to track
!! (`variable` subtable). The points in time when the values should be
!! evaluated need to be defined by a `time_control` object, see
!! [[tem_timecontrol_module]] for details.
!! And finally the format of the output needs to be given in the `output`
!! subtable.
!! Each tracking object also gets a `label` to uniquely identify it, and a
!! `folder` acts as a prefix, that gets prepended to all file names.
!! Further it is possible to define some reduction operation, for example to
!! build the average for the tracked variables over the complete shape, see
!! [[tem_reduction_spatial_module]].
!!
!! The tracking table might describe a single tracking object, or it may be
!! a table of tracking objects.
!! If there are multiple tracking objects these must be simply listed, without
!! keywords in the form of:
!!
!!```lua
!!  tracking = {
!!    {label = 'atrack', ...},
!!    {label = 'btrack', ...},
!!    {label = 'ctrack', ...},
!!  }
!!```
!!
!! If the table does not have this form, it will be assumed, that a single
!! tracking object is defined.
!!
!! Overview on the settings for the tracking object:
!!
!! * `label`: A string to identify the tracking object, will be used to name
!!            files. Defaults to `'unnamed_track'`.
!! * `variable`: A bable of variable names that should be tracked. Needs to be
!!               one of the available variables and has to be provided, ther is
!!               no default.
!! * `folder`: Actually a prefix that will be prepended to filenames. If it is
!!             to be a directory, it has to end in the path seperator and
!!             the directory has to exist already! Default: empty string, which
!!             results in the tracking files being written into the current
!!             working directory.
!! * `reduction`: A spatial reduction that is to be applied across the complete
!!                shape to arrive at a single value, see
!!                [[tem_reduction_spatial_module]] for details. If reductions
!!                are defined, there needs to be one for each variable that is
!!                tracked by this tracking object. Default is to not do any
!!                reduction.
!! * `output`: The description of how the data is to be written. There are three
!!             main options to format the output: `ascii`, `vtk` and
!!             `harvester`. The `ascii` format is useful for timeseries at a
!!             single point, that provides the point in time along with the
!!             respective variable values in an ASCII text file, easily
!!             processed by plotting tools.
!!             The `vtk` format is useful for larger shapes, like a slice
!!             through the domain. It provides the mesh information along with
!!             the values of the tracked variables.
!!             The `harvester` format writes the data of the subsection of the
!!             mesh that corresponds to the tree in the same binary format as
!!             used for restarting. The resulting data can then be processed by
!!             the harvesting tool of the solver to create visualizations.
!!             See the [[hvs_output_module]] for details on this subtable and
!!             other format options.
!!             This option has to be defined for each tracking object, there is
!!             no default.
!! * `shape`: Defines the part of the domain that is to be tracked by the
!!            object.
!!            There are various shapes available to define parts of the domain,
!!            but the most basic and common ones are points, lines and boxes,
!!            which we subsumize under the name `canoND`.
!!            Please see [[tem_shape_module]] for details.
!!
!! A simple, single tracking object definition without reduction looks like
!! this:
!!
!!```lua
!!  tracking = {
!!    label = 'track_pointA',
!!    folder = './',
!!    variable = {'momentum','density','energy'},
!!    shape = {
!!      kind = 'canoND',
!!      object= {
!!        origin ={0.01, 0, 0}
!!      }
!!    },
!!    time_control = {
!!      min = 0,
!!      max = sim_control.time_control.max,
!!      interval = {iter = 10}
!!    },
!!    output = {
!!      format = 'ascii'
!!    }
!!  }
!!```
!!
!! This tracks the variables `momentum`, `density` and `energy` at a single
!! point over time and writes them every 10 iterations to disk with one line
!! per point in time.
!!
!! See the [dedicated tracking page](../page/features/tracking.html) for more
!! examples and further hints.
!!
module tem_tracking_module

  ! incude treelm modules
  use env_module,                     only: rk, labelLen, pathLen, pathSep, &
    &                                       io_buffer_size
  use treelmesh_module,               only: treelmesh_type
  use tem_bc_prop_module,             only: tem_bc_prop_type
  use tem_aux_module,                 only: tem_abort
  use tem_comm_env_module,            only: tem_comm_env_type
  use tem_reduction_spatial_module,   only: tem_load_reduction_spatial,        &
    &                                       tem_reduction_spatial_config_type, &
    &                                       tem_reduction_spatial_init
  use tem_shape_module,               only: tem_shape_type, tem_load_shape
  use tem_subTree_module,             only: tem_create_subTree_of
  use tem_subTree_type_module,        only: tem_subTree_type
  use tem_solveHead_module,           only: tem_solveHead_type
  use tem_time_module,                only: tem_time_type
  use tem_timeControl_module,         only: tem_timeControl_type, &
    &                                       tem_timeControl_load, &
    &                                       tem_timeControl_dump, &
    &                                       tem_timeControl_check
  use tem_varSys_module,              only: tem_varSys_type
  use tem_varMap_module,              only: tem_varMap_type, tem_create_varMap
  use tem_tools_module,               only: tem_horizontalSpacer
  use tem_logging_module,             only: logUnit
  use tem_stencil_module,             only: tem_stencilHeader_type
  use tem_simControl_module,          only: tem_simControl_type
  use tem_status_module,              only: tem_stat_steady_state, &
    &                                       tem_stat_stop_file,    &
    &                                       tem_status_run_terminate, &
    &                                       tem_status_run_end
  use tem_debug_module,               only: dbgUnit

  ! include aotus modules
  use aotus_module,     only: flu_State, aot_get_val, aoterr_Fatal
  use aot_table_module, only: aot_table_open, aot_table_close, &
    &                         aot_table_length, aot_get_val

  ! include libharvesting modules
  use hvs_output_module,   only: hvs_output_config_type, hvs_output_file_type, &
    &                            hvs_output_load, hvs_output_open,             &
    &                            hvs_output_write, hvs_output_close,           &
    &                            hvs_output_finalize, hvs_AsciiTransient,      &
    &                            hvs_output_init, hvs_VTK

  implicit none

  private

  public :: tem_tracking_type
  public :: tem_tracking_instance_type
  public :: tem_tracking_config_type
  public :: tem_trackingControl_type
  public :: tem_load_tracking
  public :: tem_init_tracker
  public :: tem_init_tracker_subTree
  public :: tem_tracking_finalize
  public :: tem_tracking_getData
  public :: tem_tracking_has_triggered
  public :: tem_tracker
  public :: tem_tracking_print_last_VTK_files

  !> General information about the tracking entities
  !! This data type is set in tem_load_tracking,
  !! then updated in tem_init_tracker_subTree
  !! After load balancing, it is reset in tem_reload_tracking
  type tem_trackingControl_type
    !> output status, activated?
    logical :: active = .false.
    !> number of tracking entities active on this process
    integer :: nActive = 0
    !> Total number of tracking entities defined in config file
    integer :: nDefined = 0
  end type tem_trackingControl_type

  !> Contains all config information about tracking
  !! Content in tracking config must NOT change!
  type tem_tracking_config_type

    !> log object labels
    character(len=labelLen) :: label

    !> folder to store files to
    character(len=pathLen) :: prefix

    !> array of requested variable labels
    character(len=labelLen), allocatable :: varName(:)

    !> stores time control parameters
    type(tem_timeControl_type) :: timeControl

    !> tracking shapes
    type(tem_shape_type), allocatable  :: geometry(:)

    !> originally set to true. But if false the exact polynomial is evaluated at
    ! the point tracked
    logical :: track_complete_element

    !> Data loaded from output table
    type(hvs_output_config_type) :: output_config

    !> Spatial reduction config which is loaded from disk
    type( tem_reduction_spatial_config_type ) :: redSpatial_config

  end type tem_tracking_config_type

  !> Tracking entity definition
  type tem_tracking_instance_type
    !> Contains name and position of variables to track in global varSys
    !! number of found variables can be accessed by me%varMap%varPos%nVals
    type(tem_varMap_type) :: varMap

    !> sub-tree resulting from the elements within the tracking shape
    !! The sub-tree also holds the sub-communicator
    !! This data needs to be UPDATED after balance
    type(tem_subTree_type) :: subTree

    !> Description for output file formats
    type(hvs_output_file_type) :: output_file

    !> Pointer to config array in tem_tracking_type
    integer :: pntConfig
    
  end type tem_tracking_instance_type

  type tem_tracking_type
    !> General information about the tracking entities
    type(tem_trackingControl_type) :: control
    !> tracking header for collecting the properties from the lua file
    type(tem_tracking_config_type), allocatable :: config(:)
    !> Instances of tracking type active on this process
    type(tem_tracking_instance_type), allocatable :: instance(:)
  end type tem_tracking_type

!TG  interface assignment(=)
!TG    module procedure Copy_tracking
!TG  end interface

  contains

  ! ------------------------------------------------------------------------ !
  !> Read the tracker configuration from the main lua file
  !!
  !! Setup the values for the tracking entities
  !!
  subroutine tem_load_tracking(me, conf, parent)
    ! -------------------------------------------------------------------- !
    !> list of the trackingeentities to create
    type( tem_tracking_type ), intent(out) :: me
    !> handle of the lua config file
    type( flu_state ) :: conf
    !> if the tracking table is a child-table of some other table,
    !! use the parent as a reference
    integer, optional :: parent
    ! -------------------------------------------------------------------- !
    integer :: tc_handle, sub_handle
    integer :: iTrack, nTracks
    ! -------------------------------------------------------------------- !

    ! Read the number of trackings in the lua file
    call aot_table_open( L       = conf,       &
      &                  thandle = tc_handle,  &
      &                  key     = 'tracking', &
      &                  parent  = parent      )

    if (tc_handle == 0) then
      write(logUnit(1),*) 'No Tracking entities found!'
      call aot_table_close(L=conf, thandle=tc_handle)
      call tem_horizontalSpacer(fUnit=logUnit(1))
      me%control%nActive = 0
      me%control%nDefined = 0
      allocate( me%config(0) )
      allocate( me%instance(0) )
      return
    else ! track entity exists.
      me%control%active = .true.
    end if

    write(logUnit(1),*) 'Loading tracking ...'
    ! Check whether tracking had a subtable
    ! If no, then it is a single table, load single tracking entry
    ! else load multiple tables, open tracking subtable
    call aot_table_open( L       = conf,       &
      &                  parent  = tc_handle,  &
      &                  thandle = sub_handle, &
      &                  pos     = 1           )

    ! Only single table
    if (sub_handle == 0) then
      nTracks = 1
      write(logUnit(1),*) 'Tracking is a single table'
      allocate( me%config(1) )
      call tem_load_trackingConfig( conf       = conf,         &
        &                           sub_handle = tc_handle,    &
        &                           config     = me%config(1)  )
      call aot_table_close(L=conf, thandle=sub_handle)
    else ! Multiple table
      call aot_table_close(L=conf, thandle=sub_handle)
      nTracks = aot_table_length(L=conf, thandle=tc_handle)
      ! Allocate the defined number of tracking entities
      allocate( me%config( nTracks ))
      write(logUnit(1),"(A,I0)") 'Number of Tracking entities: ', nTracks

      ! Loop over all the definitions and assign the variables from the lua
      ! file on the tem_tracking_type.
      ! Inside this routine it will open tracking subtable. Each subtable
      ! contains one or more tracking variables the stuff is done in the
      ! routine tem_load_trackingConfig
      do iTrack = 1, nTracks
        write(logUnit(1),"(A,I0)") 'Loading tracker: ', iTrack
        call aot_table_open( L       = conf,       &
          &                  parent  = tc_handle,  &
          &                  thandle = sub_handle, &
          &                  pos     = iTrack      )
        call tem_load_trackingConfig( conf       = conf,             &
          &                           sub_handle = sub_handle,       &
          &                           config     = me%config(iTrack) )
        call aot_table_close(L=conf, thandle=sub_handle)
        write(logUnit(1),"(A,I0)") 'Done tracker ', iTrack
      end do
    end if ! sub_handle

    ! me%control%nActive = nTracks
    me%control%nDefined = nTracks
    allocate(me%instance(nTracks))

    call aot_table_close(L=conf, thandle=tc_handle) ! close tracking table
    call tem_horizontalSpacer(fUnit=logUnit(1))

  end subroutine tem_load_tracking
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Read the tracker variables from tracking subtables defined in
  !! configuration from the main lua file
  !!
  !! If tracking is just a single table with single tracking entry
  !! then load only one tracking log exists with
  !! one or more variables using tem_load_trackingHeader_single.
  !! Else if tracking is table of many log then allocate log and load each
  !! log type using tem_load_trackingHeader_single
  !! Setup the values for the tracking entities
  !!
  subroutine tem_load_trackingConfig(config, conf, sub_handle)
    ! -------------------------------------------------------------------- !
    !> list of the tracking entities to create
    type( tem_tracking_config_type ), intent(out) :: config
    !> handle of the lua config file
    type( flu_state ) :: conf
    !> table sub-handle for the tracking table
    integer, intent(in) :: sub_handle
    ! -------------------------------------------------------------------- !
    integer :: iError            ! error flag handle
    integer, allocatable :: vError(:)
    !> number of requested variables
    integer :: nRequestedVars
    ! -------------------------------------------------------------------- !

    call aot_get_val( L       = conf,           &
      &               thandle = sub_handle,     &
      &               val     = config%label,   &
      &               ErrCode = iError,         &
      &               key     = 'label',        &
      &               default = 'unnamed_track' )

    write(logUnit(1),*) 'Tracking label: '//trim( config%label )

    call aot_get_val( val       = config%varName, &
      &               ErrCode   = vError,         &
      &               maxLength = 100,            &
      &               L         = conf,           &
      &               thandle   = sub_handle,     &
      &               key       = 'variable'      )

    if ( any(btest(vError, aoterr_Fatal)) ) then
      write(logUnit(1),*) 'FATAL Error occured, while retrieving'
      write(logUnit(1),*) 'list of variables to track in '//trim(config%label)
      call tem_abort()
    end if

    nRequestedVars = size(config%varName)

    ! load time control to output tracking
    call tem_timeControl_load( conf   = conf,              &
      &                        parent = sub_handle,        &
      &                        me     = config%timeControl )
    call tem_timeControl_dump(config%timeControl, logUnit(2))

    ! Where to store the tracking file?
    call aot_get_val( L       = conf,          &
      &               thandle = sub_handle,    &
      &               val     = config%prefix, &
      &               ErrCode = iError,        &
      &               key     = 'folder',      &
      &               default = ''             )

    ! Load SPATIAL reductions
    call tem_load_reduction_spatial(                                   &
      &                   conf              = conf,                    &
      &                   parent            = sub_handle,              &
      &                   redSpatial_config = config%redSpatial_config )

    if( config%redSpatial_config%active ) then
      ! Check if the number of reductions correspond to the number of variables
      ! in the system
      if( size( config%redSpatial_config%reduceType ) /= nRequestedVars ) then
        write(logUnit(1),*) 'The number of defined reductions does not ' &
          &                 //'correspond to the '
        write(logUnit(1),*)'number of variables in the system. '
        call tem_abort()
      end if
    end if

    ! Load output table for vis_kind
    call hvs_output_load( me       = config%output_config,           &
      &                   conf     = conf,                           &
      &                   parent   = sub_handle,                     &
      &                   isReduce = config%redSpatial_config%active )

    ! load tracking object shapes like point, line, plane
    call tem_load_shape( conf        = conf,                            &
      &                  parent      = sub_handle,                      &
      &                  me          = config%geometry,                 &
      &                  reqSegments = config%output_config%useGetPoint )

    if( size( config%geometry) < 1) then
      write(logUnit(1),*)'The geometrical objects for the tracker are not '//  &
        &                'defined correctly.'
      call tem_abort()
    end if

  end subroutine tem_load_trackingConfig
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Routine creates subTree for each tracking object and removes tracking
  !! objects on process which do not include any elements to track
  !!
  !! Identify, how many and which elements exist on my local process and are
  !! requested from the trackers
  !! Empty tracking entities are removed, so the track(:) might be re-allocated
  subroutine tem_init_tracker_subTree( me, tree, solver, bc_prop, stencil, &
    &                                  prefix )
    ! -------------------------------------------------------------------- !
    !> tracking entities
    type(tem_tracking_type), intent(inout)             :: me
    !> Global mesh from which the elements are identified and then stored to
    !! sub-meshes inside the trackers
    type(treelmesh_type), intent(in)                   :: tree
    !> bc property that used to identify elements of certain BCs
    type( tem_bc_prop_type ), intent(in)               :: bc_prop
    !> Global solver information
    type(tem_solveHead_type), intent(in)               :: solver
    !> stencil used to create subTree of boundary type
    type(tem_stencilHeader_type), optional, intent(in) :: stencil
    !> Prefix for output filename
    !! Usually: solver%simName
    character(len=labelLen), optional, intent(in)      :: prefix
    ! -------------------------------------------------------------------- !
    integer :: iLog, nActive
    ! temporary tracker array
    type( tem_tracking_instance_type ), allocatable :: tempTrack(:)
    ! prefix for tracking label
    character(len=pathLen) :: prefix_loc
    ! tracking%config%prefix//tracking%config%label
    character(len=pathLen) :: basename
    ! -------------------------------------------------------------------- !
    call tem_horizontalSpacer(fUnit=logUnit(1))
    write(logUnit(3),*) 'Initialize tracking subTree to remove empty objects'
    call tem_horizontalSpacer(fUnit=logUnit(1))

    nActive = 0

    if (present(prefix)) then
      prefix_loc = trim(prefix)
    else
      ! prefix for tracking label
      prefix_loc = trim(solver%simName)//'_'
    end if

    if( me%control%active ) then
      ! Allocate the temporary track
      allocate(tempTrack( me%control%nDefined ) )

      do iLog = 1, me%control%nDefined

        basename = trim(me%config(iLog)%prefix) // trim(prefix_loc) // &
          &        trim(me%config(iLog)%label)

        write(logUnit(3),*) 'Creating subTree for tracking object ' &
          &                 // trim( me%config(iLog)%label )

        !-----------------------------------------------------------------------
        ! identify tracker elements
        !-----------------------------------------------------------------------
        call tem_create_subTree_of( inTree    = tree,                       &
          &                         bc_prop   = bc_prop,                    &
          &                         stencil   = stencil,                    &
          &                         subTree   = me%instance(iLog)%subTree,  &
          &                         inShape   = me%config(iLog)%geometry,   &
          &                         storePnts = me%config(iLog)             &
          &                                     %output_config%useGetPoint, &
          &                         prefix    = trim(basename)              )

        ! get rid of the empty track in order to avoid empty writes to disk
        if ( me%instance(iLog)%subTree%useGlobalMesh .or. &
          &  ( me%instance(iLog)%subTree%nElems > 0 ) .or. &
          &  ( me%instance(iLog)%subTree%nPoints > 0) ) then
          nActive = nActive + 1
          tempTrack( nActive ) = me%instance(iLog)
          ! Pointer to array of tracking headers loaded from config file
          tempTrack( nActive )%pntConfig = iLog
        end if

      end do  ! nActive

      deallocate(me%instance)
      allocate( me%instance(nActive) )
      me%control%nActive = nActive

      do iLog = 1, nActive
        ! Copy the stuff from the temporary track
        me%instance(iLog) = temptrack(iLog)
      end do

      deallocate(temptrack)
    end if ! if tracking active

  end subroutine tem_init_tracker_subTree
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Initialize the tracker entities:
  !! * create varMap, i.e. map requested variables to global variable system
  !! * initialize spatial reduction
  !! * initialize hvs output
  !!
  subroutine tem_init_tracker( me, tree, solver, varSys, nDofs, globProc, &
    &                          solSpec_unit                               )
    ! -------------------------------------------------------------------- !
    !> tracking entities
    type(tem_tracking_type),intent(inout) :: me
    !> Global mesh from which the elements are identified and then stored to
    !! sub-meshes inside the trackers
    type(treelmesh_type), intent(in)                  :: tree
    !> Global solver information
    type(tem_solveHead_type),intent(in)               :: solver
    !> solver-provided variable systems
    type(tem_varSys_type), intent(in)                 :: varSys
    !> The number of dofs for each scalar variable of the equation system
    integer, intent(in), optional                     :: nDofs
    !> Process description to use.
    type(tem_comm_env_type), intent(in)               :: globProc
    !> Solver specific unit for restart header
    integer, optional, intent(in)                  :: solSpec_unit
    ! -------------------------------------------------------------------- !
    integer :: iLog, nVars, iVar, iConfig
    ! prefix for tracking label to differiate tracking for different scheme
    ! with same tracking label
    character(len=pathLen) :: prefix
    ! tracking%config%prefix//tracking%config%label
    character(len=pathLen) :: basename
    ! -------------------------------------------------------------------- !

    call tem_horizontalSpacer(fUnit=logUnit(1))
    write(logUnit(1),*) 'Initialize tracking objects'
    call tem_horizontalSpacer(fUnit=logUnit(1))

    ! prefix for tracking label to differiate tracking for different scheme
    ! with same tracking label
    prefix = trim(solver%simName)//'_'

    if( me%control%active ) then

      do iLog = 1, me%control%nActive
        iConfig = me%instance(iLog)%pntConfig

        write(logUnit(3),"(A,I0,A)") 'Track object: ', iLog, ', label: ' &
          &                          // trim( me%config(iConfig)%label )

        ! map variables
        ! create tracking variable position in the global varSys
        call tem_create_varMap( varname = me%config(iConfig)%varname, &
          &                     varSys  = varSys,                     &
          &                     varMap  = me%instance(iLog)%varMap    )

        nVars = me%instance(iLog)%varMap%varPos%nVals
        ! Abort if none variables of the variable defined in current
        ! tracking object are found in varSys
        if ( nVars /= size( me%config(iConfig)%varname ) ) then
          write(logUnit(1),*) ' Some of the following variables are not found:'
          do iVar = 1, size(me%config(iConfig)%varName)
            write(logUnit(1),*) trim(me%config(iConfig)%varName(iVar))
          end do
          call tem_abort()
        end if

        basename = trim(me%config(iConfig)%prefix) // trim(prefix) &
          &        // trim(me%config(iConfig)%label)

        ! Init spatial reduction
        me%instance(iLog)%output_file%ascii%isReduce = me%config(iConfig)     &
          &                                            %redSpatial_config%active
        if ( me%config(iConfig)%redSpatial_config%active ) then
          ! Initialize reduction
          call tem_reduction_spatial_init(                                     &
            &                 me = me%instance(iLog)%output_file%ascii         &
            &                        %redSpatial,                              &
            &  redSpatial_config = me%config(iConfig)%redSpatial_config,       &
            &             varSys = varSys,                                     &
            &             varPos = me%instance(iLog)%varMap%varPos%val(:nVars) )
        end if

        ! Initialize output
        if ( me%instance(iLog)%subTree%useGlobalMesh ) then
          call hvs_output_init(out_file    = me%instance(iLog)%output_file,    &
            &                  out_config  = me%config(iConfig)%output_config, &
            &                  tree        = tree,                             &
            &                  varSys      = varSys,                           &
            &                  varPos      = me%instance(iLog)%varMap%varPos   &
            &                                                 %val(:nVars),    &
            &                  basename    = trim(basename),                   &
            &                  globProc    = globProc,                         &
            &                  timeControl = me%config(iConfig)%timeControl,   &
            &                  solver      = solver,                           &
            &                  geometry    = me%config(iConfig)%geometry,      &
            &                  nDofs       = nDofs,                            &
            &                  solSpec_unit = solSpec_unit                     )
        else
          call hvs_output_init(out_file    = me%instance(iLog)%output_file,    &
            &                  out_config  = me%config(iConfig)%output_config, &
            &                  tree        = tree,                             &
            &                  varSys      = varSys,                           &
            &                  varPos      = me%instance(iLog)%varMap%varPos   &
            &                                                 %val(:nVars),    &
            &                  subTree     = me%instance(iLog)%subTree,        &
            &                  basename    = trim(basename),                   &
            &                  globProc    = globProc,                         &
            &                  timeControl = me%config(iConfig)%timeControl,   &
            &                  solver      = solver,                           &
            &                  geometry    = me%config(iConfig)%geometry,      &
            &                  nDofs       = nDofs,                            &
            &                  solspec_unit = solSpec_unit                     )
        end if
      end do

    end if ! if tracking active

  end subroutine tem_init_tracker
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> A routine to obtain tracked data.
  !!
  !! This routine will return all requested variables in the tracking object
  !! me and return it for all elements of the subtree in the res field.
  subroutine tem_tracking_getData(varMap, subTree, varSys, mesh, time, nDofs, &
    &                             res)
    !> varMap from tem_tracking_instance_type
    type(tem_varMap_type) :: varMap

    !> subTree from tem_tracking_instance_type
    type(tem_subTree_type) :: subTree

    !> Variable system describing available data.
    type(tem_varsys_type), intent(in) :: varsys

    !> Mesh definition of the input data.
    type(treelmesh_type), intent(in) :: mesh

    !> Time information for the current data.
    type(tem_time_type), intent(in) :: time

    !> Number of degrees of freedom.
    integer, intent(in) :: nDofs

    !> Tracked data, has to match the subtree definition.
    !!
    !! The memory layout is like this:
    !!  1. All variable components
    !!  2. nDofs
    !!  3. nElems (subtree%nElems)
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: maxComponents
    integer :: nComponents
    integer :: compOff
    integer :: elemOff
    integer :: nElems
    integer :: nScalars, nVars
    integer :: elemSize
    integer :: nChunks
    integer :: chunksize
    integer :: nChunkElems
    integer :: res_size
    integer :: buf_start, buf_end
    integer :: e_start, d_start, t_start
    integer :: iElem, iChunk, iDoF, iVar
    integer :: varpos
    real(kind=rk), allocatable :: tmpdat(:)
    ! -------------------------------------------------------------------- !

    nElems = subTree%nElems
    nScalars = varMap%nScalars
    nVars = varMap%varpos%nVals

    ! Need to obtain the data variable for variable, and store it in an
    ! intermediate array, because all components should be put together in the
    ! res array.
    ! The temporary array therefore needs to be sufficiently large to store the
    ! maximal number of components.
    maxComponents = maxval(varSys%method%val(varMap%varPos &
      &                          %val(:nVars))%nComponents )

    ! Number of elements to fit into a single chunk.
    chunkSize = min( io_buffer_size / (maxComponents*nDofs), nElems )

    ! Size of a single element
    elemsize = nScalars*nDofs

    if ( (nElems > 0) .and. (chunkSize == 0) ) then
      write(logUnit(0),*) 'The chosen io_buffer_size of ', io_buffer_size
      write(logUnit(0),*) 'is too small for outputting ', maxComponents
      write(logUnit(0),*) 'scalar values with ', nDofs
      write(logUnit(0),*) 'degrees of freedom!'
      write(logUnit(0),*) 'Please increase the io_buffer_size to at least ', &
        &                 real(maxComponents*nDofs) / real(131072), ' MB!'
      call tem_abort()
    end if

    ! Using a temporary array to store the variables and transfer them to res
    ! in the correct ordering afterwards.
    allocate(tmpdat(chunkSize*maxComponents*nDofs))

    nChunks = 0
    if (chunkSize > 0) then
      nChunks = ceiling( real(nElems, kind=rk)      &
        &                / real(chunkSize, kind=rk) )
    end if

    chunks: do iChunk=1,nChunks
      elemOff = ((iChunk-1)*chunkSize)
      nChunkElems = min(chunkSize, nElems - elemOff)
      buf_start = elemOff + 1
      buf_end = elemOff + nChunkElems

      compOff = 0
      vars: do iVar=1,varMap%varPos%nVals
        varpos = varMap%varPos%val(iVar)
        nComponents = varSys%method%val(varPos)%nComponents
        res_size = nChunkElems * nDofs * nComponents
        ! derive the quantities for all the elements in the current chunk
        call varSys%method%val(varpos)%get_element(                          &
          &                                varSys  = varSys,                 &
          &                                elemPos = subtree%map2global(     &
          &                                              buf_start:buf_end), &
          &                                time    = time,                   &
          &                                tree    = mesh,                   &
          &                                nElems  = nChunkElems,            &
          &                                nDofs   = nDofs,                  &
          &                                res     = tmpdat(:res_size)       )
        do iElem=1,nChunkElems
          e_start = (elemOff+iElem-1)*elemsize
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
    end do chunks

  end subroutine tem_tracking_getData
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Decision on whether the giving tracker should be written in the current
  !! iteration.
  function tem_tracking_has_triggered( timeControl, simControl, proc ) &
    &      result(triggered)
    type(tem_timeControl_type), intent(inout) :: timeControl
    type(tem_comm_env_type),    intent(inout) :: proc
    type(tem_simControl_type),  intent(in) :: simControl
    logical :: triggered

    logical :: tc_triggered

    call tem_timeControl_check( me        = timeControl,    &
      &                         now       = simControl%now, &
      &                         comm      = proc%comm,      &
      &                         triggered = tc_triggered    )

    triggered = tc_triggered                                   &
      &         .or. tem_status_run_end(simControl%status)     &
      &         .or. tem_status_run_terminate(simControl%status)
  end function tem_tracking_has_triggered
  ! ------------------------------------------------------------------------ !



  ! ------------------------------------------------------------------------ !
  !> This routine runs over each tracking object and dump requested quantities
  !! if timeControl is active on current time
  subroutine tem_tracker( track, simControl, varSys, tree )
    ! -------------------------------------------------------------------- !
    !> tracking object containing all tracking relevant information
    type(tem_tracking_type ), intent(inout)   :: track
    !> Simulation control contains current simulation time and
    !! status bits
    type(tem_simControl_type), intent(in)     :: simControl
    !> global variable system
    type(tem_varSys_type), intent(in)         :: varSys
    !> global tree
    type(treelmesh_type ), intent(in)         :: tree
    ! -------------------------------------------------------------------- !
    integer :: iLog, iConfig
    ! -------------------------------------------------------------------- !
    
    ! Run over all tracking objects
    do iLog = 1, track%control%nActive
      iConfig = track%instance(iLog)%pntConfig

      ! Skip this tracking object, if there are no entries in the
      ! variable system
      if (track%instance( iLog )%varMap%nScalars < 1) cycle

      ! dump tracking when at least one of following conditions is triggered:
      !   tracking timeControl is triggered
      !   simulation reached steady state
      !   stop file is defined
      !   simulation is terminated abruptly
      !   simulation reaches the maximum simulation time
      if ( tem_tracking_has_triggered(                           & 
        &    timeControl = track%config(iConfig)%timeControl,    &
        &    simControl  = simControl,                           &
        &    proc        = track%instance(iLog)%output_file%proc ) ) then
        
        if( track%instance( iLog )%subTree%useGlobalMesh ) then
          ! Open the output files, this also generates the vertices for the
          ! mesh, and writes the mesh data to disk. Also writes header file
          ! depends on output vis_kind
          call hvs_output_open(                                 &
            &    out_file = track%instance(iLog)%output_file,   &
            &    use_iter = track%config(iConfig)%output_config &
            &                    %vtk%iter_filename,            &
            &    mesh     = tree,                               &
            &    varSys   = varSys,                             &
            &    time     = simControl%now                      )

          ! Evaluate and write results to disk
          call hvs_output_write( out_file = track%instance(iLog)%output_file, &
            &                    varSys   = varSys,                           &
            &                    mesh     = tree                              )

          ! Close opened files
          call hvs_output_close( out_file = track%instance(iLog)%output_file, &
            &                    varSys   = varSys,                           &
            &                    mesh     = tree                              )
        else
          ! Open the output files, this also generates the vertices for the
          ! mesh, and writes the mesh data to disk. Also writes header file
          ! depends on output vis_kind
          call hvs_output_open(                                 &
            &    out_file = track%instance(iLog)%output_file,   &
            &    use_iter = track%config(iConfig)%output_config &
            &                    %vtk%iter_filename,            &
            &    mesh     = tree,                               &
            &    varSys   = varSys,                             &
            &    subTree  = track%instance(iLog)%subTree,       &
            &    time     = simControl%now                      )

          ! Evaluate and write results to disk
          call hvs_output_write( out_file = track%instance(iLog)%output_file, &
            &                    varSys   = varSys,                           &
            &                    mesh     = tree,                             &
            &                    subTree  = track%instance(iLog)%subTree      )

          ! Close opened files
          call hvs_output_close( out_file = track%instance(iLog)%output_file, &
            &                    varSys   = varSys,                           &
            &                    mesh     = tree,                             &
            &                    subTree  = track%instance(iLog)%subTree      )
        end if !Global mesh
      
      end if  ! do tracking? interval, tmin, tmax check
    end do ! iLog

  end subroutine tem_tracker
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Close all the units of the tracking objects
  !!
  subroutine tem_tracking_finalize( me )
    ! -------------------------------------------------------------------- !
    !> tracker object to close
    type(tem_tracking_type ), intent(inout) :: me
    ! -------------------------------------------------------------------- !
    integer :: iTrack !, iError
    ! -------------------------------------------------------------------- !

    do iTrack = 1, me%control%nActive
      call hvs_output_finalize( out_file = me%instance(iTrack)%output_file )
    end do

  end subroutine tem_tracking_finalize
  ! ------------------------------------------------------------------------ !

  ! ------------------------------------------------------------------------ !
  !> Print the filenames of the last written VTK files on the screen.
  !!
  !! Mainly useful for the harvesting tools to get information on the written
  !! visualization files.
  subroutine tem_tracking_print_last_VTK_files(track)
    !> List of tracking entities
    type(tem_tracking_type), intent(in) :: track

    integer :: iTrack, iConfig
    integer :: nTracks

    nTracks = track%control%nActive

    do iTrack = 1, nTracks
      iConfig = track%instance(iTrack)%pntConfig
      ! Only the root process for this output file needs to print the name
      if (track%instance(iTrack)%output_file%proc%rank == 0) then
        ! Only the VTK outputs are considered here
        if (track%config(iConfig)%output_config%vis_kind == hvs_VTK) then
          write(*,*) 'Wrote VTK file: ', &
            &        trim(track%instance(iTrack)%output_file%vtk%last_opened_file)
        end if
      end if
    end do

  end subroutine tem_tracking_print_last_VTK_files

end module tem_tracking_module
