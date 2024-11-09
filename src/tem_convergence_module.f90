! Copyright (c) 2015-2018 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016, 2019-2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2017 Sindhuja Budaraju <nagasai.budaraju@student.uni-siegen.de>
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
!
!> This module contains type definition and routines for convergence to decide
!! steady state status.
!!
!! `convergence` may be a single table with the definition of the convergence
!! parameters or it might be a table of tables with each holding the definition
!! of a convergence criterion.
!!
!! The general `convergence` definition looks like this:
!!
!!```lua
!!  convergence = {
!!    variable = {'density', 'pressure'},
!!    time_control = {
!!      min = 0,
!!      max = {sim = 10.0},
!!      interval = {iter = 1}
!!    },
!!    shape = {kind='all'},
!!    reduction = {'average', 'average'},
!!    norm = 'average',
!!    nvals = 10,
!!    absolute = true,
!!    condition = {
!!      {threshold = 1.e-15, operator = '<=' },
!!      {threshold = 1.e-12, operator = '<=' }
!!    }
!!  }
!!```
!!
!! `variable` needs to be a list of variables from the variable system, providing
!! the quantities that are to checked for convergence.
!!
!! `time_control` needs to be a time control definition, see
!! [[tem_timeControl_module]] it defines the timespan within which the convergence
!! check should be active. Its interval setting describes when to perform the
!! necessary operations to decide whether a steady state is reached.
!!
!! `shape` needs to be shape definition to describe a subsection of the overall
!! domain in which the check for convergence is to be done (might be `{kind='all'}`
!! to take the complete simulation domain into account. See [[tem_shape_module]] for
!! more details.
!!
!! `reduction` needs to be a table of spatial reduction definitions for each variable.
!! Typical reduction operations could be `average`, `sum` and `l2norm`. This reduction
!! is applied to obtain a single value for the convergence decision across all elements
!! within the domain section defined by `shape`. See [[tem_reduction_spatial_module]]
!! for more details.
!!
!! `norm` describes how to treat the spatially reduced values over time.
!! There are two options: 'simple' will just compare the current value with the one
!! from the last check. 'average' will build the average over the `nvals` last values
!! that have been obtained via the spatial reduction. Default of this setting is to
!! use the simple comparison.
!!
!! `nvals` is a setting used if `norm = 'average'` and sets the size of the sliding
!! window in the historical data across which an average will be made to obtain the
!! value to compare against in the current check.
!!
!! `absolute` denotes whether to use an absolute or relative measure to make the
!! comparisons provided in the `condition` setting. The default is `false` that
!! which results in a relative measure. Relative means here, that the difference
!! of the current value and the comparison value is checked against the threshold
!! multiplied with the current value instead of using the threshold itself directly.
!!
!! `use_get_point` indicates whether to use individual points to obtain the values
!! or whole elements. Default is `false`, so the complete state of elements in the
!! shape will be used to find the reduced values.
!!
!! `ndofs` describes how many degrees of freedom to use from elements if
!! `use_get_point` is `false`. The default is to use all degrees of freedom, which
!! can be indicated by setting `ndofs = -1`.
!!
!! `condition` needs to be table providing the condition under which convergence
!! is assumed for each variable. The condition is defined by a threshold and an
!! operator. See [[tem_condition_module]] for details.
!!
!! Convergence is assumed if all variables defined in the convergence table
!! meet their defined definition.
!! To check for convergence a reduction and condition has to be defined for each
!! variable that is to be considered for the steady state check.
!!
!!
module tem_convergence_module

  ! include treelm modules
  use env_module, only: rk, labelLen, io_buffer_size
  use tem_aux_module, only: tem_abort
  use tem_logging_module, only: logUnit
  use tem_tools_module, only: tem_horizontalSpacer, upper_to_lower
  use tem_comm_env_module, only: tem_comm_env_type
  use tem_status_module, only: tem_stat_steady_state, tem_status_type
  use tem_time_module, only: tem_time_type
  use tem_timeControl_module,  only: tem_timeControl_type,  &
    &                                tem_timeControl_load,  &
    &                                tem_timeControl_dump,  &
    &                                tem_timeControl_out,   &
    &                                tem_timeControl_check, &
    &                                tem_timeControl_reset_trigger
  use tem_varSys_module,       only: tem_varSys_type,       &
    &                                tem_get_element_chunk, &
    &                                tem_get_point_chunk
  use tem_varMap_module,       only: tem_varMap_type, &
    &                                tem_create_varMap
  use treelmesh_module,        only: treelmesh_type
  use tem_bc_prop_module,      only: tem_bc_prop_type
  use tem_reduction_spatial_module, only: tem_reduction_spatial_type,        &
    &                                     tem_reduction_spatial_close,       &
    &                                     tem_reduction_spatial_config_type, &
    &                                     tem_load_reduction_spatial,        &
    &                                     tem_reduction_spatial_init,        &
    &                                     tem_reduction_spatial_open,        &
    &                                     tem_reduction_spatial_append,      &
    &                                     tem_reduction_spatial_out
  use tem_shape_module,        only: tem_shape_type, &
    &                                tem_load_shape, &
    &                                tem_shape_out
  use tem_subTree_module,      only: tem_create_subTree_of
  use tem_subTree_type_module, only: tem_subTree_type
  use tem_solveHead_module,    only: tem_solveHead_type
  use tem_stencil_module,      only: tem_stencilHeader_type
  use tem_condition_module,    only: tem_condition_type, &
    &                                tem_load_condition, &
    &                                tem_comparator,     &
    &                                tem_condition_out
  use tem_geometry_module,     only: tem_BaryOfId

  ! include aotus modules
  use aotus_module,     only: flu_State, aot_get_val, aoterr_Fatal
  use aot_table_module, only: aot_table_open,   &
    &                         aot_table_close,  &
    &                         aot_table_length, &
    &                         aot_get_val
  use aot_out_module,   only: aot_out_type,       &
    &                         aot_out_open,       &
    &                         aot_out_close,      &
    &                         aot_out_val,        &
    &                         aot_out_open_table, &
    &                         aot_out_close_table

  implicit none
  private

  public :: tem_condition_type, tem_convergence_type
  public :: tem_convergence_load
  public :: tem_convergence_dump
  public :: tem_convergence_out
  public :: tem_convergence_check
  public :: tem_init_convergence
  public :: tem_convergence_reset

  !> Convergence description loaded from config file
  type tem_convergenceHeader_type
    !> convergence kind
    character(len=labelLen) :: norm
    !> number of defined conditions
    integer :: nConditions
    !> An instance of the condition type for each variable
    type(tem_condition_type), allocatable :: cond(:)
    !> Number of last values to check for convergence
    integer :: nLastVals
    !> absolute Error (.true.) or relative Error( .false.)
    logical :: absoluteError
    !> array of variable labels to check for convergence
    character(len=labelLen), allocatable :: varName(:)
    !> number of variables to check for convergence
    !! i.e size(varName)
    integer :: nRequestedVars
    !> stores time control parameters
    type(tem_timeControl_type) :: timeControl
    !> convergence shapes
    type(tem_shape_type), allocatable  :: geometry(:)

    !> reduction config
    type(tem_reduction_spatial_config_type) :: redSpatial_config

    !> Logic to decide to use get_point or get_element to dump data
    logical :: useGetPoint

    !> Number of dofs to check for convergence if useGetPoint = .false.
    integer :: nDofs

  end type tem_convergenceHeader_type

  !> The convergence type which contains convergence flag
  !! and an instance of the condition type
  type tem_convergence_type
    !> Convergence header info
    type(tem_convergenceHeader_type) :: header
    !> norm kind
    integer :: norm_kind
    !> state field holding the reference values for the nScalars
    !! size: nLastVals, nScalars
    real(kind=rk), allocatable :: lastState(:,:)
    !> number of performed convergence checks
    !! corresponds to the entry in the lastState array
    integer :: nChecks
    !> Process description to use for the output.
    !! Might be only a subset of the global communicator
    type(tem_comm_env_type)  :: proc
    !> Contains name and position of variables in global varSys
    type(tem_varMap_type) :: varMap
    !> sub-tree resulting from the elements within the convergence shape
    !! The sub-tree also holds the sub-communicator
    type(tem_subTree_type) :: subTree
    !> number of elements that fit in the buffer
    integer :: chunkSize
    !> number of chunks per output
    integer :: nChunks

    !> The number of dofs for each scalar variable of the equation system
    integer :: nDofs

    !> spatial reduction for each variable
    type(tem_reduction_spatial_type), allocatable :: redSpatial(:)

  end type tem_convergence_type

  !> Convergence norm kinds
  integer, parameter :: norm_simple = 1
  integer, parameter :: norm_average = 2

  interface assignment(=)
    module procedure Copy_convergence
  end interface

  interface tem_convergence_dump
    module procedure tem_convergence_dump_vector
    module procedure tem_convergence_dump_single
  end interface tem_convergence_dump

  interface tem_convergence_out
    module procedure tem_convergence_out_vector
    module procedure tem_convergence_out_single
  end interface tem_convergence_out


contains


  ! ************************************************************************ !
  !> Load the convergence definition table
  !! The convergence object must be part of a convergence object, for which the
  !! format has been set to format = 'convergence'
  !! In the convergence table, you then must define a norm:
  !!
  !! - simple: just check against the state value of the last check, and reach
  !!           convergence if below the defined threshold
  !! - average: build the average over a defined set of last checks with nvals
  !!           stop, if the difference to the current state value is below the
  !!           given threshold
  !! - nvals: define, how many last checks should be taken into account for
  !!          averaging procedure
  !!
  !! The error is by default calculated to be a relative error. If an absolute
  !! error is desired, choose absolute=true in the convergence object
  !!
  !! The stopping criterion is defined as a general condition object, where the
  !! threshold and the operator has to be given
  !!
  !!```lua
  !!  condition = { threshold = 1.E-6, operator = '<' }
  !!```
  !!
  !! A sample convergence object with a convergence definition can look as
  !! follows (within time_control table):
  !!```lua
  !!  abort_criteria = {
  !!   stop_file = 'stop',
  !!   steady_state = true,
  !!   convergence = {
  !!     variable = {'pressure','velocity'},
  !!     shape = {kind = 'all'},
  !!     time_control = {
  !!       min = {iter=0},
  !!       max = {iter=tmax},
  !!       interval = {iter=10*dt}},
  !!     reduction = {'average','average'},
  !!     norm='average', nvals = 100, absolute = true,
  !!     condition = {
  !!        { threshold = 1.e-15, operator = '<=' },
  !!        { threshold = 1.e-12, operator = '<=' }
  !!     }
  !!   }
  !!  }
  !!```
  !!
  !! Or another sample:
  !!```lua
  !!  abort_criteria = {
  !!    stop_file     = 'stop',
  !!    steady_state  = true,
  !!    convergence   = {
  !!      variable = {'pressure_phy'},
  !!      shape = {
  !!        {kind = 'canoND', object = {origin ={0.15-dx,0.2,zpos} }},
  !!        {kind = 'canoND', object = {origin ={0.25+dx,0.2,zpos} }}
  !!      },
  !!      time_control = {min = 0, max = tmax, interval = 10*dt},
  !!      reduction = {'average'},
  !!      norm      = 'average',
  !!      nvals     = 50,
  !!      absolute  = true,
  !!      condition = { threshold = 1.e-10, operator = '<=' }
  !!    }
  !!  }
  !!```
  !!
  subroutine tem_convergence_load(me, conf, parent, steady_state)
    ! ---------------------------------------------------------------------------
    !> list of the convergence entities to create
    type( tem_convergence_type ), allocatable, intent(out) :: me(:)
    !> general control parameters
    !> handle of the lua config file
    type( flu_state ) :: conf
    !> if the convergence table is a child-table of some other table,
    !! use the parent as a reference
    integer, optional :: parent
    !> Steady flag in abort_criteria to check for convergence
    logical, intent(inout) :: steady_state
    ! ---------------------------------------------------------------------------
    integer :: conv_handle, sub_handle
    integer :: iConv, nConv
    ! ---------------------------------------------------------------------------

    ! Read the number of convergences in the lua file
    call aot_table_open( L       = conf,          &
      &                  thandle = conv_handle,   &
      &                  key     = 'convergence', &
      &                  parent  = parent         )

    if (conv_handle == 0) then
      write(logUnit(1),*) 'WARNING: Abort criteria, steady state is true but'
      write(logUnit(1),*) '         No Convergence table is defined with '
      write(logUnit(1),*) '         conditions to check for steady state'
      write(logUnit(1),*) 'NOTE: Steady state is deactivated'
      steady_state = .false.
      call aot_table_close(L=conf, thandle=conv_handle)
      call tem_horizontalSpacer(fUnit=logUnit(1))
      return
    end if

    write(logUnit(1),*) 'Loading convergence for steady state...'
    ! Check whether convergence had a subtable
    ! If no, then it is a single table, load single convergence entry
    ! else load multiple tables, open convergence subtable
    call aot_table_open( L       = conf,        &
      &                  parent  = conv_handle, &
      &                  thandle = sub_handle,  &
      &                  pos     = 1            )

    ! Only single table
    if (sub_handle == 0) then
      nConv = 1
      write(logUnit(1),*) 'Convergence is a single table'
      allocate( me( nConv ) )
      call tem_load_convergenceHeader( conf       = conf,        &
        &                              sub_handle = conv_handle, &
        &                              me         = me(1)        )
      call aot_table_close(L=conf, thandle=sub_handle)
    else ! Multiple table
      call aot_table_close(L=conf, thandle=sub_handle)
      nConv = aot_table_length(L=conf, thandle=conv_handle)
      ! Allocate the defined number of convergence entities
      allocate( me( nConv ))
      write(logUnit(1),*) 'Number of Convergence entities: ', nConv

      ! Loop over all the definitions and assign the variables from the lua
      ! file on the tem_convergence_type.
      ! Inside this routine it will open convergence subtable. Each subtable
      ! contains one or more convergence variables the stuff is done in the
      ! routine tem_load_convergenceHeader
      do iConv = 1, nConv
        write(logUnit(3),*) 'Loading convergence ', iConv
        call aot_table_open( L       = conf,        &
          &                  parent  = conv_handle, &
          &                  thandle = sub_handle,  &
          &                  pos     = iConv      )
        call tem_load_convergenceHeader( conf           = conf,       &
          &                              sub_handle     = sub_handle, &
          &                              me             = me(iConv)   )
        call aot_table_close(L=conf, thandle=sub_handle)
        write(logUnit(3),*) 'Done'
      end do
    end if ! sub_handle

    call aot_table_close(L=conf, thandle=conv_handle) ! close convergence table
    call tem_horizontalSpacer(fUnit=logUnit(1))

  end subroutine tem_convergence_load
! **************************************************************************** !


! **************************************************************************** !
  !> Read the convergence variables from convergence subtables defined in
  !! configuration from the main lua file
  !!
  !! If convergence is just a single table with single convergence entry
  !! then load only one convergence log exists with
  !! one or more variables using tem_load_convergenceHeader_single.
  !! Else if convergence is table of many log then allocate log and load each
  !! log type using tem_load_convergenceHeader_single
  !! Setup the values for the convergence entities
  !!
  subroutine tem_load_convergenceHeader(me, conf, sub_handle)
    ! --------------------------------------------------------------------------
    !> list of the convergence entities to create
    type( tem_convergence_type ), intent(out) :: me
    !> handle of the lua config file
    type( flu_state ) :: conf
    !> table sub-handle for the convergence table
    integer, intent(in) :: sub_handle
    ! --------------------------------------------------------------------------
    integer :: iError            ! error flag handle
    integer, allocatable :: vError(:)
    character(len=labelLen) :: norm
    ! --------------------------------------------------------------------------

    call aot_get_val( val       = me%header%varname, &
      &               ErrCode   = vError,            &
      &               maxLength = 100,               &
      &               L         = conf,              &
      &               thandle   = sub_handle,        &
      &               key       = "variable"         )

    if ( any(btest(vError, aoterr_Fatal)) ) then
      write(logUnit(1),*) 'FATAL Error occured, while retrieving'
      write(logUnit(1),*) 'list of variables to use in convergence'
      call tem_abort()
    end if

    me%header%nRequestedVars = size(me%header%varName)

    ! load time control to output convergence
    call tem_timeControl_load( conf           = conf,                          &
      &                        parent         = sub_handle,                    &
      &                        me             = me%header%timeControl )
    call tem_timeControl_dump(me%header%timeControl, logUnit(3))

    ! load convergence object shapes like point, line, plane
    call tem_load_shape( conf = conf, parent = sub_handle,                     &
                         me = me%header%geometry )

    if( size( me%header%geometry) < 1) then
      write(logUnit(1),*)'The geometrical objects for the convergence are'//  &
        &            ' not defined correctly.'
      call tem_abort()
    end if

    ! load reductions
    call tem_load_reduction_spatial( conf   = conf,                       &
      &                              parent = sub_handle,                 &
      &                   redSpatial_config = me%header%redSpatial_config )

    if( me%header%redSpatial_config%active ) then
      ! Check if the number of reductions correspond to the number of variables
      ! in the system
      if( size( me%header%redSpatial_config%reduceType ) &
        & /= me%header%nRequestedVars ) then
        write(logUnit(1),*) 'Error: In convergence.'
        write(logUnit(1),*) 'The number of defined reductions does not '// &
          &            'correspond to the '
        write(logUnit(1),*)'number of variables in the system. '
        call tem_abort()
      end if
    else
      write(logUnit(1),*) 'Error: No spatial reduction defined.'
      write(logUnit(1),*) 'NOTE: Convergence requires spatial reduction '//&
        &                 'for each variable'
      call tem_abort()
    end if

    ! get the kind of the convergence norm
    call aot_get_val( L       = conf,                                          &
      &               thandle = sub_handle,                                    &
      &               val     = norm,                                          &
      &               ErrCode = iError,                                        &
      &               key     = 'norm',                                        &
      &               default = 'simple' )
    norm = adjustl(norm)
    norm = upper_to_lower(norm)
    me%header%norm = trim(norm)
    select case( trim( norm ))
    case( 'simple' )
      me%norm_kind = norm_simple
      ! Only need one last value to compare against in the simple case.
      me%header%nLastVals = 1
    case( 'average')
      me%norm_kind = norm_average
      ! Get number of last values to check for convergence
      call aot_get_val( L       = conf,                &
        &               thandle = sub_handle,          &
        &               val     = me%header%nLastVals, &
        &               ErrCode = iError,              &
        &               key     = 'nvals',             &
        &               default = 1                    )
    case default
      write(logUnit(1),*) 'Error: Unknown convergence norm'
      write(logUnit(1),*) 'Solution: Supported norms '
      write(logUnit(1),*) '        * simple'
      write(logUnit(1),*) '        * average'
      call tem_abort
    end select


    ! type of convergence error: relative or absolute
    call aot_get_val( L       = conf,                    &
      &               thandle = sub_handle,              &
      &               val     = me%header%absoluteError, &
      &               ErrCode = iError,                  &
      &               key     = 'absolute',              &
      &               default = .false.                  )

    ! To decide whether to use get_point or get_element
    call aot_get_val( L       = conf,                  &
      &               thandle = sub_handle,            &
      &               val     = me%header%useGetPoint, &
      &               ErrCode = iError,                &
      &               default = .false.,               &
      &               key     = 'use_get_point'        )

    ! Get the number of Dofs to be written in the output
    ! The default is set to -1. If the dofs are not specified,
    ! all the dofs should be dumped
    call aot_get_val( L       = conf,            &
      &               thandle = sub_handle,      &
      &               val     = me%header%nDofs, &
      &               ErrCode = iError,          &
      &               default = -1,              &
      &               key     = 'ndofs'          )

    ! load condition for convergence
    call tem_load_condition( me      = me%header%cond, &
      &                      conf    = conf,           &
      &                      parent  = sub_handle      )

    me%header%nConditions = size(me%header%cond)
    ! check if there is condition for each variable
    if (me%header%nConditions /= me%header%nRequestedVars) then
      write(logUnit(1),*) 'Error: Nr. of conditions \= Nr. of variables '
      write(logUnit(1),"(2(A,I0))") 'nCond: ', me%header%nConditions, &
        &                           'nVars: ', me%header%nRequestedVars
      call tem_abort()
    end if

    write(logUnit(1),"(A,I0)") ' loaded convergence with nConditions=', &
      &                 me%header%nConditions
    write(logUnit(1),*) '    Norm : '//trim(norm)
    if( me%header%absoluteError ) then
      write(logUnit(1),*)'    absolute error'
    else
      write(logUnit(1),*)'    relative error'
    end if
    write(logUnit(1),"(A,I0)")'    nVal : ', me%header%nLastVals
    write(logUnit(7),*) '  Use get_point: ', me%header%useGetPoint
    call aot_table_close(L=conf, thandle=sub_handle )

  end subroutine tem_load_convergenceHeader
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Initialize the convergence subtreee
  !!
  !! Identify, how many and which elements exist on my local process and are
  !! requested from the convergences
  !! Empty convergence  are removed, so the convergence(:) might be re-allocated
  !!
  subroutine tem_init_convergence( me, tree, varSys, bc_prop, globProc,&
    &                              stencil, nDofs )
    ! -------------------------------------------------------------------- !
    !> convergence descriptions
    type(tem_convergence_type),intent(inout), allocatable :: me(:)
    !> Global mesh from which the elements are identified and then stored to
    !! sub-meshes inside the convergences
    type(treelmesh_type), intent(in) :: tree
    !> bc property that used to identify elements of certain BCs
    type( tem_bc_prop_type ), intent(in) :: bc_prop
    !> solver-provided variable systems
    type(tem_varSys_type), intent(in) :: varSys
    !> Process description to use.
    type(tem_comm_env_type), intent(in) :: globProc
    !> stencil used to create subTree of boundary type
    type(tem_stencilHeader_type), optional, intent(in) :: stencil
    !> The number of dofs for each scalar variable of the equation system
    integer, intent(in), optional :: nDofs
    ! -------------------------------------------------------------------- !
    integer :: iConv, nConv, nVars
    integer :: nChunks, chunkSize, nElems, maxComponents, nPoints
    ! temporary convergence array
    type( tem_convergence_type ),allocatable :: tempConv(:)
    ! -------------------------------------------------------------------- !

    call tem_horizontalSpacer(fUnit=logUnit(1))
    write(logUnit(1),*) 'Setting up the convergence infrastructure'
    nConv = 0


    ! Allocate the temporary convergence
    allocate(tempConv(size(me)))

    do iConv = 1, size(me)
      write(logUnit(10),"(A,I0)") 'Initializing convergence object ', iConv

      ! map variables
      ! create convergence variable position in the global varSys
      call tem_create_varMap( varname = me(iConv)%header%varname, &
        &                     varSys  = varSys,                   &
        &                     varMap  = me(iConv)%varMap          )

      ! @todo KM: If variable not found in varSys then skip that variable
      ! reduction and condition info while copying convergence
      !
      ! Terminate if number of variables to check for convergence does not
      ! match with variables found in varMap
      if (me(iConv)%varMap%varPos%nVals /= me(iConv)%header%nRequestedVars) then
        write(logUnit(1),*) 'Error: Mapping Convergences variables'
        write(logUnit(1),*) 'Variables defined in convergence '// &
          &                 'table ', iConv, ' are not found in varSys'
        call tem_abort()
      end if

      ! -----------------------------------------------------------------------
      ! identify convergence elements
      ! -----------------------------------------------------------------------
      call tem_create_subTree_of( inTree    = tree,                           &
        &                         bc_prop   = bc_prop,                        &
        &                         stencil   = stencil,                        &
        &                         subTree   = me( iConv )%subTree,            &
        &                         storePnts = me (iConv )%header%useGetPoint, &
        &                         inShape   = me( iConv )%header%geometry     )

      ! get rid of the empty convergence in order
      ! to avoid empty writes to disk
      if ( me(iConv)%subTree%useGlobalMesh ) then

        ! set convergence communicator
        me(iConv)%proc = globProc

        nConv = nConv + 1
        tempConv( nConv ) = me( iConv )

      else if ( me(iConv)%subTree%nElems > 0 ) then

        nConv = nConv + 1 ! Increase number of log

        ! set convergence communicator
        me(iConv)%proc%comm      = me(iConv)%subTree%global%comm
        me(iConv)%proc%rank      = me(iConv)%subTree%global%myPart
        me(iConv)%proc%comm_size = me(iConv)%subTree%global%nParts
        me(iConv)%proc%root      = 0

        ! Copy all entries from the log derived type to the temporary one.
        ! this might not work on all compilers!!
        ! This assignment is realized by operator overloader Copy_convergence
        tempConv( nConv ) = me( iConv )
      end if ! useGlobalMesh
    end do  ! nConv

    deallocate(me)
    allocate( me(nConv) )

    do iConv = 1, nConv
      ! Copy the stuff from the temporary track
      me(iConv) = tempConv(iConv)

      ! number of variables in varMap
      nVars = me(iConv)%varMap%varPos%nVals

      ! nDofs is valid only for get_element
      if (me(iConv)%header%useGetPoint) then
        me(iConv)%nDofs = 1
      else
        if (present(nDofs)) then
          ! out_config%nDofs is set to -1 if unspecied
          ! in the config file. In this case all the dof's
          ! should be dumped
          if (me(iConv)%header%nDofs < 0) then
            me(iConv)%nDofs = nDofs
          else
            ! Otherwise the number of dofs dumped should
            ! be what's specified in the config
            me(iConv)%nDofs = me(iConv)%header%nDofs
          end if
        else
          me(iConv)%nDofs = 1
        end if
      end if

      if (me(iConv)%subTree%useGlobalMesh) then
        nElems = tree%nElems
        nPoints = tree%nElems
      else
        nElems = me(iConv)%subTree%nElems
        nPoints = me(iConv)%subTree%nPoints
      end if

      ! max nComponent in current convergence variables
      maxComponents = maxval(varSys%method%val(me(iConv)%varMap               &
        &                                     %varPos%val(:nVars))%nComponents)

      if (me(iConv)%header%useGetPoint) then
        chunkSize = min(io_buffer_size/maxComponents, nPoints)
      else
        chunkSize = min(io_buffer_size/(maxComponents*me(iConv)%nDofs), nElems)
      end if
      if ( (nElems > 0) .and. (chunkSize == 0) ) then
        write(logUnit(0),*)'Error in init_convergence: '
        write(logUnit(0),*) 'The chosen io_buffer_size of ', io_buffer_size
        write(logUnit(0),*) 'is too small for evaluating ', maxComponents
        write(logUnit(0),*) 'scalar values'
        write(logUnit(0),*) 'Please increase the io_buffer_size to at &
          & least ', real(maxComponents*me(iConv)%nDofs) / real(131072), ' MB!'
        call tem_abort()
      end if

      nChunks = 0
      if (chunkSize>0) then
        if (me(iConv)%header%useGetPoint) then
          nChunks = ceiling(real(nPoints, kind=rk)/real(chunkSize, kind=rk))
        else
          nChunks = ceiling(real(nElems, kind=rk)/real(chunkSize, kind=rk))
        end if
      else
        nChunks = 0
      end if

      me(iConv)%nChunks = nChunks
      me(iConv)%chunkSize = chunkSize
      me(iConv)%nChecks = 0

      ! Initialize reduction
      call tem_reduction_spatial_init(                                 &
        &      me                = me(iConv)%redSpatial,               &
        &      redSpatial_config = me(iConv)%header%redSpatial_config, &
        &      varSys            = varSys,                             &
        &      varPos            = me(iConv)%varMap%varPos%val(:nVars) )

      ! Allocate some arrays
      allocate( me(iConv)%lastState( me(iConv)%header%nLastVals, &
        &                           me(iConv)%varMap%nScalars ) )
      me(iConv)%lastState = huge( me(iConv)%lastState(1,1) )          &
        &                 / real( me(iConv)%header%nLastVals, kind=rk )

    end do

    deallocate(tempConv)

    call tem_horizontalSpacer(fUnit=logUnit(1))

  end subroutine tem_init_convergence
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine resets convergence lastState and nChecks
  subroutine tem_convergence_reset( me )
    ! -------------------------------------------------------------------- !
    !> convergence descriptions
    type(tem_convergence_type),intent(inout) :: me(:)
    ! -------------------------------------------------------------------- !
    integer :: iConv
    ! -------------------------------------------------------------------- !
    do iConv = 1, size(me)
      me(iConv)%nChecks = 0
      me(iConv)%lastState = huge( me(iConv)%lastState(1,1) )          &
        &                 / real( me(iConv)%header%nLastVals, kind=rk )
      call tem_timeControl_reset_trigger(me(iConv)%header%timeControl)
    end do

  end subroutine tem_convergence_reset
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine runs over each convergence object and check for convergence
  !! of each requested quantities timeControl is active on current time
  subroutine tem_convergence_check( me, time, status, varSys, tree )
    ! -------------------------------------------------------------------- !
    !> convergence descriptions
    type(tem_convergence_type), target, intent(inout)   :: me(:)
    !> current simulation time
    type(tem_time_type), intent(in)     :: time
    !> Status bits
    type(tem_status_type), intent(inout) :: status
    !> global variable system
    type(tem_varSys_type), intent(in)         :: varSys
    !> global tree
    type(treelmesh_type ), intent(in)         :: tree
    ! -------------------------------------------------------------------- !
    integer :: iConv
    logical :: triggered
    real(kind=rk), allocatable :: res(:)
    ! -------------------------------------------------------------------- !
    allocate(res(io_buffer_size))

    ! Run over all convergence objects
    do iConv = 1,size(me)

      ! Skip this convergence object, if there are no entries in the
      ! variable system
      if (me( iConv )%varMap%nScalars < 1) cycle

      ! Determine, if it is now the time to perform the convergence for each
      ! convergence object
      ! check time control and update it if triggered
      call tem_timeControl_check( me        = me( iConv )%header%timeControl, &
        &                         now       = time,                           &
        &                         comm      = me( iConv )%proc%comm,          &
        &                         triggered = triggered                       )

      ! check convergence if timeControl is triggered
      if( triggered )then
        if (me(iConv)%header%useGetPoint) then
          call tem_convergence_check_point( me     = me(iConv), &
            &                               time   = time,      &
            &                               status = status,    &
            &                               varSys = varSys,    &
            &                               tree   = tree,      &
            &                               res    = res        )
        else
          call tem_convergence_check_element( me     = me(iConv), &
            &                                 time   = time,      &
            &                                 status = status,    &
            &                                 varSys = varSys,    &
            &                                 tree   = tree,      &
            &                                 res    = res        )
        end if
      end if  ! do convergence? interval, tmin, tmax check
    end do ! iConv

    deallocate(res)
  end subroutine tem_convergence_check
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine runs over convergence check using get_element interface
  subroutine tem_convergence_check_element( me, time, status, varSys, tree, &
    &                                       res )
    ! -------------------------------------------------------------------- !
    !> convergence descriptions
    type(tem_convergence_type), target, intent(inout)   :: me
    !> current simulation time
    type(tem_time_type), intent(in)     :: time
    !> Status bits
    type(tem_status_type), intent(inout) :: status
    !> global variable system
    type(tem_varSys_type), intent(in)         :: varSys
    !> global tree
    type(treelmesh_type ), intent(in)         :: tree
    !> Output data
    !! size: io_buffer_size
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: nVars, nElems, nScalars, elemOff, nChunkElems
    integer :: iElem, iChunk
    integer :: buf_start, buf_end
    integer, allocatable :: elemPos(:)
    logical :: isConverged
    ! -------------------------------------------------------------------- !
    ! number of variables in varMap
    nVars = me%varMap%varPos%nVals

    ! Number of scalars in current output
    nScalars = me%varMap%nScalars

    if (me%subTree%useGlobalMesh) then
      nElems = tree%nElems
    else
      nElems = me%subTree%nElems
    end if

    ! open spatial reduction
    call tem_reduction_spatial_open( me     = me%redSpatial,               &
      &                              varSys = varSys,                      &
      &                              varPos = me%varMap%varPos%val(:nVars) )

    ! allocate elemPos to size of chunkSize
    allocate(elemPos(me%chunkSize))

    ! Process all chunks to derive the quantities defined in the tracking object
    do iChunk = 1, me%nChunks
      ! Number of elements read so far in previous chunks.
      elemOff = ((iChunk-1)*me%chunkSize)

      ! number of elements written to THIS chunk
      nChunkElems = min(me%chunkSize, nElems-elemOff)

      ! Compute the element lower and upper bound for the current chunk
      buf_start = elemOff + 1
      buf_end = elemOff + nChunkElems

      if (me%subTree%useGlobalMesh) then
        elemPos(1:nChunkElems) = (/ (iElem, iElem=buf_start, buf_end) /)
      else
        elemPos(1:nChunkElems) = me%subTree                     &
          &                        %map2Global(buf_start:buf_end)
      end if

      ! evaluate all variables on current chunk
      call tem_get_element_chunk(varSys  = varSys,                 &
        &                        varPos  = me%varMap%varPos        &
        &                                           %val(:nVars),  &
        &                        elemPos = elemPos(1:nChunkElems), &
        &                        time    = time,                   &
        &                        tree    = tree,                   &
        &                        nElems  = nChunkElems,            &
        &                        nDofs   = me%nDofs,               &
        &                        res     = res                     )

      ! preform spatial reduction
      call tem_reduction_spatial_append( me     = me%redSpatial,            &
        &                                chunk  = res,                      &
        &                                nElems = nChunkElems,              &
        &                                treeID = tree%treeID(              &
        &                                         elemPos(1:nChunkElems) ), &
        &                                tree   = tree,                     &
        &                                varSys = varSys,                   &
        &                                nDofs  = me%nDofs,                 &
        &                                varPos = me%varMap%varPos          &
        &                                                  %val(:nVars)     )

    end do ! iChunk

    deallocate(elemPos)

    ! Close reduction for current convergence
    call tem_reduction_spatial_close( me   = me%redSpatial, &
      &                               proc = me%proc        )

    ! evaluate residual and check convergence for each scalar
    ! and set steady state flag only in root process of convergence
    ! communicator
    if (me%proc%rank == 0) then
      ! Now check for convergence for each reduced scalars
      call tem_convergence_evaluate( me       = me,         &
        &                            achieved = isConverged )

      ! set steady state if any convergence object is converged
      status%bits(tem_stat_steady_state) = status%bits(tem_stat_steady_state)&
        &                                  .or. isConverged
      if (isConverged) then
        write(logUnit(5),*) 'Reached steady state ', time%iter, isConverged
      end if
    end if

  end subroutine tem_convergence_check_element
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine runs over convergence check using get_point interface
  subroutine tem_convergence_check_point( me, time, status, varSys, tree, res )
    ! -------------------------------------------------------------------- !
    !> convergence descriptions
    type(tem_convergence_type), target, intent(inout)   :: me
    !> current simulation time
    type(tem_time_type), intent(in)     :: time
    !> Status bits
    type(tem_status_type), intent(inout) :: status
    !> global variable system
    type(tem_varSys_type), intent(in)         :: varSys
    !> global tree
    type(treelmesh_type ), intent(in)         :: tree
    !> Output data
    !! size: io_buffer_size
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: nVars, nPoints, nScalars, pointsOff, nChunkPoints
    integer :: iPoint, iChunk, counter
    integer :: buf_start, buf_end
    real(kind=rk), allocatable :: points(:,:)
    logical :: isConverged
    ! -------------------------------------------------------------------- !
    ! number of variables in varMap
    nVars = me%varMap%varPos%nVals

    ! Number of scalars in current output
    nScalars = me%varMap%nScalars

    if (me%subTree%useGlobalMesh) then
      nPoints = tree%nElems
    else
      nPoints = me%subTree%nPoints
    end if

    ! open spatial reduction
    call tem_reduction_spatial_open( me     = me%redSpatial,               &
      &                              varSys = varSys,                      &
      &                              varPos = me%varMap%varPos%val(:nVars) )

    ! allocate points to size of chunkSize
    allocate(points(me%chunkSize,3))

    ! Process all chunks to derive the quantities defined in the tracking object
    do iChunk = 1, me%nChunks
      ! Number of points read so far in previous chunks.
      pointsOff = ((iChunk-1)*me%chunkSize)

      ! number of points written to THIS chunk
      nChunkPoints = min(me%chunkSize, nPoints-pointsOff)

      ! Compute the points lower and upper bound for the current chunk
      buf_start = pointsOff + 1
      buf_end = pointsOff + nChunkPoints

      if (me%subTree%useGlobalMesh) then
        counter = 0
        do iPoint = buf_start, buf_end
          counter = counter + 1
          points(counter, :) = tem_BaryOfId( tree,               &
            &                                tree%treeID(iPoint) )
        end do
      else
        points(1:nChunkPoints,:) = me%subTree%points(buf_start:buf_end,:)
      end if

      ! evaluate all variables on current chunk
      call tem_get_point_chunk(varSys  = varSys,                   &
        &                      varPos  = me%varMap%varPos          &
        &                                         %val(:nVars),    &
        &                      point   = points(1:nChunkPoints,:), &
        &                      time    = time,                     &
        &                      tree    = tree,                     &
        &                      nPnts   = nChunkPoints,             &
        &                      res     = res                       )

      ! preform spatial reduction
      call tem_reduction_spatial_append( me     = me%redSpatial,        &
        &                                chunk  = res,                  &
        &                                nElems = nChunkPoints,         &
        &                                tree   = tree,                 &
        &                                varSys = varSys,               &
        &                                varPos = me%varMap%varPos      &
        &                                                  %val(:nVars) )

    end do ! iChunk

    deallocate(points)

    ! Close reduction for current convergence
    call tem_reduction_spatial_close( me   = me%redSpatial, &
      &                               proc = me%proc        )

    ! evaluate residual and check convergence for each scalar
    ! and set steady state flag only in root process of convergence
    ! communicator
    if (me%proc%rank == 0) then
      ! Now check for convergence for each reduced scalars
      call tem_convergence_evaluate( me       = me,         &
        &                            achieved = isConverged )

      ! set steady state if any convergence object is converged
      status%bits(tem_stat_steady_state) = status%bits(tem_stat_steady_state)&
        &                                  .or. isConverged
      if (isConverged) then
        write(logUnit(10),*) 'Reached steady state ', time%iter, isConverged
      end if
    end if

  end subroutine tem_convergence_check_point
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Evaluate if the convergence was achieved
  !!
  subroutine tem_convergence_evaluate( me, achieved )
    ! -------------------------------------------------------------------- !
    !> Convergence description contians
    type(tem_convergence_type), intent(inout) :: me
    !> is all Scalars in current convergence_type are converged
    logical, intent(out)                  :: achieved
    ! -------------------------------------------------------------------- !
    integer :: iVar, iComp, iScalar
    real(kind=rk) :: residual, threshold_fac
    logical :: isConverged(me%varMap%nScalars)
    ! -------------------------------------------------------------------- !
    ! Increase the counter for the checks
    me%nChecks = me%nChecks + 1
    iScalar = 0
    do iVar = 1, me%varMap%varPos%nVals
      do iComp = 1, me%redSpatial(iVar)%nComponents
        iScalar = iScalar + 1
        ! Compare the results against threshold
        !norm = abs(state(1) - me%result_prev)
        residual = evaluate_residual(                          &
          &          me      = me,                             &
          &          state   = me%redspatial(ivar)%val(icomp), &
          &          iScalar = iScalar                         )

        if( me%header%absoluteError ) then
          ! For absolute error, just use threshold
          threshold_fac = me%header%cond(iVar)%threshold
        else
          ! For relative errors, multiply threshold with current state value
          threshold_fac = me%redSpatial(iVar)%val(iComp) &
            &           * me%header%cond(iVar)%threshold
        end if

        isConverged(iScalar) = tem_comparator(                                 &
          &                        val       = residual,                       &
          &                        operation = me%header%cond(iVar)%operation, &
          &                        threshold = threshold_fac                   )
      end do
    end do
    ! update the convergence achieved if the last compare was succesful
    achieved = all(isConverged)

  end subroutine tem_convergence_evaluate
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> evaluate the residual
  !! For relative errors (defined in convergence%absoluteError ), the result is
  !! divided by the current status value
  !!
  function evaluate_residual( me, state, iScalar ) result( res )
    ! -------------------------------------------------------------------- !
    !> Convergence description
    type(tem_convergence_type), intent(inout) :: me
    !> Spatially reduced variable value
    real(kind=rk), intent(in)              :: state
    !> Current scalar
    integer, intent(in) :: iScalar
    !> residual to check for convergence
    real(kind=rk) :: res
    ! -------------------------------------------------------------------- !
    integer :: pos_lastState
    real(kind=rk) :: average
    ! -------------------------------------------------------------------- !
    ! Reset the result
    res = huge( res )

    select case (me%norm_kind )
    case( norm_simple )
      if (me%nChecks > 1) then
        res = abs((state - me%lastState(1, iScalar)))
      end if
      ! update the result at t-1 to t as when we arrive at t+1, it will
      ! be required
      me%lastState(1, iScalar) = state
    case( norm_average )
      pos_lastState = mod( me%nChecks - 1, me%header%nLastVals ) + 1
      if ( me%nChecks <= me%header%nLastVals ) then
        average = sum( me%lastState(1:pos_lastState, iScalar) ) &
          &     / real( pos_lastState, kind=rk )
      else
        average = sum( me%lastState(:, iScalar) )    &
          &     / real( me%header%nLastVals, kind=rk )
      end if

      res = abs( (state - average ) )
      me%lastState( pos_lastState, iScalar ) = state
      !write(*,*) 'nCheck ', me%nChecks, iScalar
      !write(*,*) 'lastState', me%lastState
      !write(*,*) 'state:', state,'average', average, 'res', res, &
      !  &        'pos last', pos_lastState
    end select

  end function evaluate_residual
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Dumps array of convergence to given unit
  subroutine tem_convergence_dump_vector(me, outUnit)
    ! -------------------------------------------------------------------- !
    !> convergence to write into the lua file
    type(tem_convergence_type), intent(in) :: me(:)
    !> unit to write to
    integer, intent(in) :: outUnit
    ! -------------------------------------------------------------------- !
    ! aotus type handling the output to the file in lua format
    type(aot_out_type) :: conf
    ! -------------------------------------------------------------------- !
    call aot_out_open( put_conf = conf, outUnit = outUnit )
    call tem_convergence_out_vector( me, conf )
    call aot_out_close( put_conf = conf )

  end subroutine tem_convergence_dump_vector
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Dump single convergence to given unit
  subroutine tem_convergence_dump_single(me, outUnit)
    ! -------------------------------------------------------------------- !
    !> convergence to write into the lua file
    type(tem_convergence_type), intent(in) :: me
    !> unit to write to
    integer, intent(in) :: outUnit
    ! -------------------------------------------------------------------- !
    ! aotus type handling the output to the file in lua format
    type(aot_out_type) :: conf
    ! -------------------------------------------------------------------- !
    call aot_out_open( put_conf = conf, outUnit = outUnit )
    call tem_convergence_out_single( me, conf )
    call aot_out_close( put_conf = conf )

  end subroutine tem_convergence_dump_single
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Allows the output of array of convergence to lua out
  subroutine tem_convergence_out_vector(me, conf)
    ! -------------------------------------------------------------------- !
    !> convergence to write into the lua file
    type(tem_convergence_type), intent(in) :: me(:)
    !> aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    ! -------------------------------------------------------------------- !
    integer :: iConv
    ! -------------------------------------------------------------------- !
    call aot_out_open_table( put_conf = conf, tname='convergence' )
    do iConv = 1,size(me)
      call tem_convergence_out_single( me(iConv), conf, level=1 )
    end do
    call aot_out_close_table( put_conf = conf )

  end subroutine tem_convergence_out_vector
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Allows the output of the single convergence to lua out.
  !!
  !! The data is written into the file, the lunit is connected to.
  !! It is formatted as a Lua table.
  !!
  subroutine tem_convergence_out_single(me, conf, level)
    ! -------------------------------------------------------------------- !
    !> convergence to write into the lua file
    type(tem_convergence_type), intent(in) :: me
    !> aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    !> to dump variable with key or without key
    integer, optional, intent(in) :: level
    ! -------------------------------------------------------------------- !
    integer :: level_loc, iVar
    ! -------------------------------------------------------------------- !

    if (present(level)) then
      level_loc = level
    else
      level_loc = 0
    end if

    if( level_loc == 0) then
      call aot_out_open_table( put_conf = conf, tname = 'convergence' )
    else
      call aot_out_open_table( put_conf = conf )
    end if

    ! write variable names
    call aot_out_open_table( put_conf = conf, tname = 'variable' )
    do iVar = 1, me%header%nRequestedVars
      call aot_out_val( put_conf = conf,                         &
        &               val      = trim(me%header%varName(iVar)) )
    end do
    call aot_out_close_table( put_conf = conf )

    ! convergence norm
    call aot_out_val( put_conf = conf,           &
      &               val      = me%header%norm, &
      &               vname    = 'norm'          )
    call aot_out_val( put_conf = conf,                &
      &               val      = me%header%nLastVals, &
      &               vname    = 'nvals'              )

    call aot_out_val( put_conf = conf,                    &
      &               val      = me%header%absoluteError, &
      &               vname    = 'absolute'               )

    ! write reductions
    call tem_reduction_spatial_out( me%redSpatial, conf )

    ! write conditions
    call tem_condition_out( me%header%cond, conf )

    ! write timeControl info
    call tem_timeControl_out( me%header%timeControl, conf )
    ! write shapes
    call tem_shape_out( me%header%geometry, conf )

    call aot_out_close_table( put_conf = conf )

  end subroutine tem_convergence_out_single
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides the assignment operator of two convergence.
  !! Temporary Solution as CRAY compiler dont have = Operator
  !! Copying a convegence object (right) into other convergence (left)
  !!
  subroutine Copy_convergence(left,right)
    ! -------------------------------------------------------------------- !
    !> convegence to copy to
    type(tem_convergence_type), intent(out) :: left
    !> convegence to copy from
    type(tem_convergence_type), intent(in) :: right
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    left%header = right%header
    left%norm_kind = right%norm_kind
    left%varMap = right%varMap
    left%subTree = right%subTree
    left%proc = right%proc
    if (allocated(right%redSpatial)) then
      allocate(left%redSpatial(size(right%redSpatial)))
      left%redSpatial = right%redSpatial
    end if
    left%nDofs = right%nDofs

  end subroutine Copy_convergence
  ! ************************************************************************ !

end module tem_convergence_module
