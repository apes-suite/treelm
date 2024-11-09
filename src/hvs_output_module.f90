! Copyright (c) 2015-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
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
!> Module to encapsulate output for various formats.
module hvs_output_module
  use aotus_module,     only: flu_State, aot_get_val, aoterr_Fatal,            &
    &                         aoterr_WrongType, aoterr_NonExistent
  use aot_table_module, only: aot_table_open, aot_table_close

  use env_module,                     only: LabelLen, PathLen, long_k, rk
  use treelmesh_module,               only: treelmesh_type
  use tem_aux_module,                 only: tem_abort
  use tem_comm_env_module,            only: tem_comm_env_type
  use tem_dyn_array_module,           only: PositionOfVal
  use tem_logging_module,             only: logunit
  use tem_subtree_type_module,        only: tem_subtree_type, tem_dump_subTree
  use tem_subTree_module,             only: tem_updatePropertyBits
  use tem_time_module,                only: tem_time_type, tem_time_reset
  use tem_tools_module,               only: upper_to_lower
  use tem_varsys_module,              only: tem_varsys_type
  use tem_varMap_module,              only: tem_varMap_type, tem_create_varMap
  use tem_vrtx_module,                only: tem_vrtx_type, &
    &                                       tem_vrtx_finalize, &
    &                                       tem_calc_vrtx_coord
  use tem_geometry_module,            only: tem_BaryOfId
  ! use tem_transient_reduction_module, only: tem_transient_reduction_type,      &
  !   &                                       tem_transient_reduction_init
  use tem_restart_module,             only: tem_restart_type, tem_init_restart,&
    &                                       tem_restart_dump_data,             &
    &                                       tem_restart_finalize,              &
    &                                       tem_restart_openWrite,             &
    &                                       tem_restart_closeWrite
  use tem_timeControl_module,         only: tem_timeControl_type
  use tem_shape_module,               only: tem_shape_type
  use tem_solveHead_module,           only: tem_solveHead_type


  use hvs_vtk_module, only: hvs_vtk_config_type, hvs_vtk_config_load,          &
    &                       hvs_vtk_file_type, hvs_vtk_open,                   &
    &                       hvs_vtk_write_meshdata, hvs_vtk_close,             &
    &                       hvs_vtk_write_varsys, hvs_vtk_dump_data,           &
    &                       hvs_vtk_init, hvs_vtk_closePVD
  use hvs_ascii_module, only: hvs_ascii_type, hvs_asciiSpatial_type,           &
    &                         hvs_asciiSpatial_open, hvs_ascii_dump_elem_data, &
    &                         hvs_ascii_dump_point_data,                       &
    &                         hvs_asciiSpatial_dump_elem_data,                 &
    &                         hvs_asciiSpatial_dump_point_data,                &
    &                         hvs_ascii_close,                                 &
    &                         hvs_asciiSpatial_close, hvs_ascii_init,          &
    &                         hvs_asciiSpatial_init

  implicit none

  private

  !> parameter for ascii output, writing only one row consisting
  !! of all the elements with all variable values in that row
  !! into an initially opened file
  integer, public, parameter :: hvs_AsciiTransient = 1
  !> Parameter for gnuplot ascii, writing a complete spatial
  !! representation at each time step into a new file
  integer, public, parameter :: hvs_AsciiSpatial   = 2
  !> Parameter for harvester output
  integer, public, parameter :: hvs_Internal       = 3
  !> Parameter for precice spatial coupling
  integer, public, parameter :: hvs_PreciceSpatial = 4
  !> Parameter for VTK output
  integer, public, parameter :: hvs_VTK            = 5

  public :: hvs_output_config_type
  public :: hvs_output_file_type
  public :: hvs_output_load
  public :: hvs_output_init
  public :: hvs_output_open
  public :: hvs_output_write
  public :: hvs_output_close
  public :: hvs_output_finalize

  !> This data type contains data loaded from disk.
  type hvs_output_config_type
    !> Kind of visualization file to use for the output.
    integer :: vis_kind

    !> Description of the vtk configuration (used when vis_kind='vtk')
    type(hvs_vtk_config_type) :: vtk

    !> Logic to decide to use get_point or get_element to dump data
    logical :: useGetPoint

    !> Stores the number of output dofs to be dumped.
    !! @todo KM: 20160425 nDofs is reasonable to configure in output
    !! only for vtk so move to vtk_type.
    integer :: nDofs

  end type hvs_output_config_type

  type hvs_output_file_type
    !> Kind of visualization file to use for the output.
    integer :: vis_kind

    !> Basename to use for the output files.
    character(len=PathLen) :: basename

    !> Process description to use for the output.
    !! Might be only a subset of the global communicator
    type(tem_comm_env_type) :: proc

    !> Description for vtk output
    type(hvs_vtk_file_type) :: vtk

    !> Description for ascii output
    type(hvs_ascii_type) :: ascii

    !> Description for ascii output
    type(hvs_asciiSpatial_type) :: asciiSpatial

    !> Description for harvester output i.e restart format
    type(tem_restart_type) :: restart

    !> Point in time for which this output should be done.
    type(tem_time_type) :: time

    !> Number of variables to write to this file.
    integer :: nVars

    !> List of variable positions to write into this file.
    integer, allocatable :: varPos(:)

    !> Vertex information of elements within the tracking shape.
    !! Required for vtk output
    type(tem_vrtx_type) :: vrtx

    !> Store the barycenters for the linearized tree elements
    !! it has a size of ( nElems, 3 ).
    !! Used in dump AsciiSpatial
    !! It is set  in hvs_output_init
    !! It is used in hvs_output_write
    !! Need to UPDATE after balance
    real(kind=rk), allocatable :: bary(:,:)

    !> is there any transient reduction active?
    !! Then this tracking must be updated every timestep!
    ! logical :: isTransientReduce = .false.

    !> transient reductions which collect and reduce data over several timesteps
    ! type(tem_transient_reduction_type), allocatable :: transientReduce(:)

    !> The number of dofs for each scalar variable of the equation system
    integer :: nDofs

    !> Logic to decide to use get_point or get_element to dump data
    logical :: useGetPoint
  end type hvs_output_file_type


contains


  ! ----------------------------------------------------------------------------!
  !> Read the output configuration from a Lua script.
  subroutine hvs_output_load(me, conf, parent, isReduce)
    ! --------------------------------------------------------------------------!
    !> The output configuration settings to fill.
    type(hvs_output_config_type), intent(out) :: me

    !> Handle of the Lua script to load the configuration from.
    type(flu_state) :: conf

    !> Table handle to the table providing the output settings.
    integer, intent(in), optional :: parent

    !> true if reduction is defined
    logical, intent(in) :: isReduce
    ! ----------------------------------------------------------------------!
    character(len=labelLen) :: vkind
    integer :: thandle
    integer :: iError
    ! ----------------------------------------------------------------------!
    write(logUnit(3), *) '  Loading output table ...'

    call aot_table_open( L       = conf,    &
      &                  parent  = parent,  &
      &                  thandle = thandle, &
      &                  key     = 'output' )

    if (thandle > 0) then
      ! If reduction is active set output%vis_kind to hvs_asciiTransient
      ! else load vis_kind from output table
      if (isReduce) then
        write(logUnit(7),*) 'Spatial reduction is active, set output format to ascii'
        vkind = 'ascii'
      else
        call aot_get_val( L       = conf,     &
          &               thandle = thandle,  &
          &               key     = 'format', &
          &               val     = vkind,    &
          &               ErrCode = iError    )

        if ( btest(iError, aoterr_Fatal) ) then
          write(logUnit(0),*) 'Fatal Error: In reading format for output'
          if ( btest(iError, aoterr_NonExistent) ) then
            write(logUnit(0),*) 'NonExistent: No format for output provided.'
          end if
          if( btest(iError, aoterr_WrongType) ) then
            write(logUnit(0),*) 'WrongType: No format for output provided.'
          end if
          write(logUnit(0),*) 'STOPPING'
          call tem_abort()
        end if
      end if

      vkind = adjustl(vkind)
      vkind = upper_to_lower(vkind)
      select case(trim(vkind))
      case('vtk')
        ! Output data in paraview unstructured vtk format (.vtu)
        me%vis_kind = hvs_VTK
        call hvs_vtk_config_load( me      = me%vtk, &
          &                       conf    = conf,   &
          &                       thandle = thandle )
      case('ascii')
        me%vis_kind = hvs_AsciiTransient
      case('asciispatial')
        ! write an ascii file for each time step and write all three barycenters
        ! each line has corrdinates and variables of only ONE element
        ! #       coordX  coordY  coordZ    var1       var2        var3
        me%vis_kind = hvs_AsciiSpatial
      case('harvester')
        ! ... in the harvester format, meaning
        ! there is a mesh in treelm format and corresponding to it
        ! there are restart files which hold the state information for the
        ! elements in the mesh for the requested time step
        me%vis_kind = hvs_Internal
      case('precice')
        ! ... for preCICE
        ! there is a mesh in treelm format and corresponding to it
        ! there are restart files which hold the state information for the
        ! elements in the mesh for the requested time step
        me%vis_kind = hvs_PreciceSpatial
      case default
        write(logunit(0),*) 'ERROR in hvs_output_load: '     &
          &                 // 'unknown visualization kind ' &
          &                 // trim(vkind)
        write(logunit(0),*) 'format has to be one of:'
        write(logunit(0),*) '* vtk'
        write(logunit(0),*) '* ascii'
        write(logunit(0),*) '* asciiSpatial'
        write(logunit(0),*) '* harvester'
        write(logunit(0),*) '* precice'
        write(logunit(0),*)
        call tem_abort()
      end select
    else
      write(logUnit(0),*) 'FATAL Error: output table is not defined.'
      write(logunit(0),*) 'Please provide an output table with at least the'
      write(logunit(0),*) 'format defined in it!'
      write(logUnit(0),*) 'STOPPING'
      call tem_abort()
    end if
    write(logUnit(5),*) '  Output format: '//trim(vkind)

    ! To decide whether to use get_point or get_element
    call aot_get_val( L       = conf,           &
      &               thandle = thandle,        &
      &               val     = me%useGetPoint, &
      &               ErrCode = iError,         &
      &               default = .false.,        &
      &               key     = 'use_get_point' )
    write(logUnit(7),*) '  Use get_point: ', me%useGetPoint

    ! Get the number of Dofs to be written in the output
    ! The default is set to -1. If the dofs are not specified,
    ! all the dofs should be dumped
    call aot_get_val( L       = conf,     &
      &               thandle = thandle,  &
      &               val     = me%nDofs, &
      &               ErrCode = iError,   &
      &               default = -1,       &
      &               key     = 'ndofs'   )

    call aot_table_close( L       = conf,   &
      &                   thandle = thandle )

  end subroutine hvs_output_load
  ! ----------------------------------------------------------------------------!
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  !> Initialize the output for a given mesh.
  !!
  !! This creates vertex for a mesh and fill hvs_output_file_type.
  subroutine hvs_output_init(out_file, out_config, tree, varsys, subtree,      &
    &                        varPos, basename, timeControl, nDofs, globProc,   &
    &                        solver, geometry, solSpec_unit)
    ! --------------------------------------------------------------------------!
    !> Output file settings
    !! It must be intent inout since ascii%reduction and trasient%reduction
    !! are loaded in tem_load_trackingHeader
    type(hvs_output_file_type), intent(inout) :: out_file

    !> The output configuration settings to use.
    type(hvs_output_config_type), intent(in) :: out_config

    !> Mesh of the data to visualize.
    type(treelmesh_type), intent(in) :: tree

    !> Description of the available variable system to get the given varnames
    !! from.
    type(tem_varSys_type), intent(in) :: varsys

    !> Optional restriction of the elements to output.
    type(tem_subtree_type), optional, intent(in) :: subtree

    !> List of variable positions that should be written in the output.
    !!
    !! If this is not provided, all variables from the varsys will be written
    !! to the vtk file.
    integer, optional, intent(in) :: varPos(:)

    !> An extension to the output basename.
    !!
    !! The filename will be constructed by tracking%header%prefix//
    !! tracking%header%label
    character(len=*), intent(in) :: basename

    !> output timeControl
    type(tem_timeControl_type), optional, intent(in) :: timeControl

    !> The number of dofs for each scalar variable of the equation system
    integer, intent(in), optional :: nDofs

    !> Global communicator type for global rank information
    type(tem_comm_env_type ), intent(in) :: globProc

    !> Global solver information
    type(tem_solveHead_type ),intent(in) :: solver

    !> shape defined for this ascii output
    type(tem_shape_type), optional, intent(in) :: geometry(:)

    !> Solver specific unit for restart header
    integer, optional, intent(in) :: solSpec_unit
    ! ----------------------------------------------------------------------!
    integer :: iVar, iElem
    integer(kind=long_k) :: glob_nElems, glob_nPoints
    integer :: nElems
    integer :: nPoints
    integer(kind=long_k) :: tTreeID
    ! local varMap for restart init
    type(tem_varMap_type) :: varMap
    ! ----------------------------------------------------------------------!
    ! Copy visualization kind
    out_file%vis_kind = out_config%vis_kind

    out_file%useGetPoint = out_config%useGetPoint

    ! copy basename
    out_file%basename = trim(basename)

    ! nDofs is valid only for get_element
    if (out_file%useGetPoint) then
      out_file%nDofs = 1
    else
      if (present(nDofs)) then
        ! out_config%nDofs is set to -1 if unspecied
        ! in the config file. In this case all the dof's
        ! should be dumped
        if (out_config%nDofs < 0) then
          out_file%nDofs = nDofs
        else
          ! Otherwise the number of dofs dumped should
          ! be what's specified in the config
          out_file%nDofs = out_config%nDofs
        end if
      else
        out_file%nDofs = 1
      end if
    end if

    ! Gather global nElems, local nElems and communicator environment
    if ( present(subTree) ) then
      if (out_file%useGetPoint) then
        nPoints = subTree%nPoints
        glob_nPoints = subTree%glob_nPoints
      else
        nPoints = subTree%nElems
        glob_nPoints = subTree%global%nElems
      end if

      nElems = subTree%nElems
      glob_nElems = subTree%global%nElems
      ! set output communicator
      out_file%proc%comm      = subTree%global%comm
      out_file%proc%rank      = subTree%global%myPart
      out_file%proc%comm_size = subTree%global%nParts
      out_file%proc%root      = 0
    else
      nElems = tree%nElems
      glob_nElems = tree%global%nElems

      ! @todo KM: use nDofs to convert nElems to nPoints
      nPoints = tree%nElems
      glob_nPoints = tree%global%nElems
      ! set output communicator
      out_file%proc = globProc
    end if

    if (allocated(out_file%varpos)) deallocate(out_file%varpos)

    if (present(varPos)) then

      out_file%nVars = size(varPos)
      allocate(out_file%varpos(out_file%nVars))
      out_file%varpos = varpos

    else

      out_file%nVars = varsys%varname%nVals
      allocate(out_file%varpos(out_file%nVars))
      do iVar=1,out_file%nVars
        out_file%varpos(iVar) = iVar
      end do

    end if

    ! ! Init transient reduction
    ! if (out_file%isTransientReduce) then
    !   ! initialize transient reductions
    !   call tem_transient_reduction_init( me     = out_file%transientReduce,    &
    !     &                                nElems = nElems,                      &
    !     &                                varSys = varSys,                      &
    !     &                                varPos = out_file%varPos,             &
    !     &                                nDofs  = out_file%nDofs,              &
    !     &                                time   = timeControl%min              )
    ! end if

    select case(out_file%vis_kind)
    case(hvs_AsciiTransient)
      call hvs_ascii_init( ascii        = out_file%ascii,       &
        &                  varSys       = varSys,               &
        &                  varPos       = out_file%varPos,      &
        &                  basename     = trim(basename),       &
        &                  globProc     = globProc,             &
        &                  nDofs        = out_file%nDofs,       &
        &                  outproc      = out_file%proc,        &
        &                  nElems       = nElems,               &
        &                  glob_nElems  = glob_nElems,          &
        &                  timeControl  = timeControl,          &
        &                  solver       = solver,               &
        &                  useGetPoint  = out_file%useGetPoint, &
        &                  nPoints      = nPoints,              &
        &                  glob_nPoints = glob_nPoints,         &
        &                  geometry     = geometry              )
    case(hvs_AsciiSpatial)
      ! Store barycenter to dump at every time step
      if ( present(subTree) ) then

        if (out_file%useGetPoint) then
          allocate( out_file%bary( nPoints, 3 ) )
          out_file%bary = subTree%points
        else
          allocate( out_file%bary( nElems, 3 ) )
          do iElem = 1, nElems
            tTreeID = tree%treeID( subTree%map2global(iElem) )
            out_file%bary(iElem, :) = tem_BaryOfId(tree, tTreeID)
          end do
        end if

      else

        allocate( out_file%bary( nElems, 3 ) )
        do iElem = 1, nElems
          tTreeID = tree%treeID( iElem )
          out_file%bary(iElem, :) = tem_BaryOfId(tree, tTreeID)
        end do

      end if

      call hvs_asciiSpatial_init( asciiSpatial = out_file%asciiSpatial, &
        &                         varSys       = varSys,                &
        &                         varPos       = out_file%varPos,       &
        &                         basename     = trim(basename),        &
        &                         globProc     = globProc,              &
        &                         outproc      = out_file%proc,         &
        &                         nDofs        = out_file%nDofs,        &
        &                         nElems       = nElems,                &
        &                         glob_nElems  = glob_nElems,           &
        &                         useGetPoint  = out_file%useGetPoint,  &
        &                         nPoints      = nPoints,               &
        &                         glob_nPoints = glob_nPoints,          &
        &                         timeControl  = timeControl,           &
        &                         solver       = solver,                &
        &                         geometry     = geometry               )
    case(hvs_Internal)
      ! -----------------------------------------------------------------------
      ! Initialize restart header, communicator and chunk info
      ! -----------------------------------------------------------------------
      out_file%restart%controller%writePrefix =  trim(basename) // '_'
      out_file%restart%controller%writeRestart = .true.

      ! this restart object is not meant to read data!
      out_file%restart%controller%readRestart = .false.

      ! create varMap for restart
      call tem_create_varMap(varName = varSys%varName%val(out_file%varPos), &
        &                    varSys  = varSys,                              &
        &                    varMap  = varMap                               )

      ! init the restart typed file format
      ! should be called when the first time the mesh is dumped
      call tem_init_restart( me           = out_file%restart, &
        &                    solver       = solver,           &
        ! &                    varSys       = varsys,           &
        &                    varMap       = varMap,           &
        &                    tree         = tree,             &
        &                    subTree      = subTree,          &
        &                    solSpec_unit = solSpec_unit,     &
        &                    nDofs_write  = out_file%nDofs    )

      ! Dump tree if output is harvester format and shape is not global
      ! i.e subTree is present
      if (present(subTree)) then
        call tem_dump_subTree( subTree, tree )
      end if
    case(hvs_VTK)
      ! Calculate vertex for vtk output
      call tem_calc_vrtx_coord( tree    = tree,          &
        &                       vrtx    = out_file%vrtx, &
        &                       subtree = subtree        )

      call hvs_vtk_init( vtk_file   = out_file%vtk,   &
        &                vtk_config = out_config%vtk, &
        &                basename   = trim(basename), &
        &                proc       = out_file%proc   )
    end select

  end subroutine hvs_output_init
  ! ----------------------------------------------------------------------------!
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  !> Open the output for a given mesh.
  !!
  !! This writes the mesh data to the output files.
  !! Subsequently it can then be enriched by restart information.
  subroutine hvs_output_open(out_file, use_iter, mesh, varsys, time, subtree )
    ! --------------------------------------------------------------------------!
    type(hvs_output_file_type), intent(inout) :: out_file

    !> The output configuration settings to use.
    ! type(hvs_output_config_type), intent(in) :: out_config
    ! Use iteration in filename?
    logical, intent(in) :: use_iter

    !> Mesh of the data to visualize.
    type(treelmesh_type), intent(in) :: mesh

    !> Description of the available variable system to get the given varnames
    !! from.
    type(tem_varSys_type), intent(in) :: varsys

    !> Optional restriction of the elements to output.
    type(tem_subtree_type), optional, intent(inout) :: subtree

    !> Time information.
    !!
    !! If this is present, the filename will be built with a time stamp and
    !! the time point information is written into the vtu file.
    type(tem_time_type), intent(in), optional :: time
    ! ----------------------------------------------------------------------!
    integer :: nElems
    ! ----------------------------------------------------------------------!
    if (present(subtree)) then
      nElems = subtree%nElems
    else
      nElems = mesh%nElems
    end if

    if (present(time)) then
      out_file%time = time
    else
      call tem_time_reset(out_file%time)
    end if

    select case(out_file%vis_kind)
    case(hvs_AsciiSpatial)
      call hvs_asciiSpatial_open( asciiSpatial = out_file%asciiSpatial, &
        &                         outProc      = out_file%proc,         &
        &                         time         = out_file%time,         &
        &                         varsys       = varsys,                &
        &                         varpos       = out_file%varpos,       &
        &                         nDofs        = out_file%nDofs         )
    case(hvs_Internal)
      ! prepare the header differently for using global tree and tracking
      ! subTree. For harvester format, prepare the header.
      if( present(subTree) ) then
        ! For harvester format, prepare the header
        ! Here we pass the tracking subTree and the global tree to the header
        ! In tem_restart_openWrite the subTree will be written to disk

        ! update the property bits according to the global mesh if this has
        ! changed
        call tem_updatePropertyBits( mesh, subTree )

        call tem_restart_openWrite(me      = out_file%restart,             &
          &                        tree    = mesh,                         &
          &                        timing  = out_file%time,                &
          &                        varSys  = varSys,                     &
          &                        subTree = subTree,                      &
          &                        label   = trim(out_file%basename)//'_'  )
      else
        call tem_restart_openWrite( me     = out_file%restart, &
          &                         tree   = mesh,             &
          &                         timing = out_file%time,    &
          &                         varSys = varSys          )
      end if ! present subTree
    case(hvs_VTK)
      call hvs_vtk_open( vtk_file = out_file%vtk,   &
        &                use_iter = use_iter,       &
        &                proc     = out_file%proc,  &
        &                time     = out_file%time   )

      call hvs_vtk_write_meshdata( vtk_file = out_file%vtk,  &
        &                          vrtx     = out_file%vrtx, &
        &                          nElems   = nElems         )

      call hvs_vtk_write_varSys( vtk_file = out_file%vtk,   &
        &                        varsys   = varsys,         &
        &                        varpos   = out_file%varpos )
    end select

  end subroutine hvs_output_open
  ! ----------------------------------------------------------------------------!
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  subroutine hvs_output_write(out_file, varsys, mesh, subtree)
    ! --------------------------------------------------------------------------!
    type(hvs_output_file_type), intent(inout) :: out_file

    !> Description of the available variable system to get the given varnames
    !! from.
    type(tem_varSys_type), intent(in) :: varsys

    !> Mesh to write the data on.
    type(treelmesh_type), intent(in) :: mesh

    !> Optional restriction of the elements to output.
    type(tem_subtree_type), optional, intent(in) :: subtree
    ! --------------------------------------------------------------------------!
    integer :: iVar
    ! --------------------------------------------------------------------------!

    select case(out_file%vis_kind)
    case(hvs_AsciiTransient)
      if (out_file%useGetPoint) then
        ! Call the get point routine and dump the data at the points
        call hvs_ascii_dump_point_data(                                      &
        &                       ascii           = out_file%ascii,            &
        &                       outProc         = out_file%proc,             &
        &                       varPos          = out_file%varPos,           &
        &                       varSys          = varSys,                    &
        &                       mesh            = mesh,                      &
        &                       time            = out_file%time,             &
        &                       subTree         = subtree )
        ! &                       transientReduce = out_file%transientReduce   )
      else
        ! call the get_element routine and dump all dofs
        call hvs_ascii_dump_elem_data(                                        &
          &                     ascii           = out_file%ascii,             &
          &                     outProc         = out_file%proc,              &
          &                     varPos          = out_file%varPos,            &
          &                     varSys          = varSys,                     &
          &                     mesh            = mesh,                       &
          &                     time            = out_file%time,              &
          &                     nDofs           = out_file%nDofs,             &
          &                     subTree         = subTree )
          ! &                     transientReduce = out_file%transientReduce    )

      endif
    case(hvs_AsciiSpatial)
      if (out_file%useGetPoint) then
        call hvs_asciiSpatial_dump_point_data(                            &
          &                   asciiSpatial    = out_file%asciiSpatial,    &
          &                   varPos          = out_file%varPos,          &
          &                   varSys          = varSys,                   &
          &                   mesh            = mesh,                     &
          &                   time            = out_file%time,            &
          &                   subTree         = subTree,                  &
          ! &                   transientReduce = out_file%transientReduce, &
          &                   bary            = out_file%bary             )
      else
        call hvs_asciiSpatial_dump_elem_data(                             &
          &                   asciiSpatial    = out_file%asciiSpatial,    &
          &                   varPos          = out_file%varPos,          &
          &                   varSys          = varSys,                   &
          &                   mesh            = mesh,                     &
          &                   time            = out_file%time,            &
          &                   nDofs           = out_file%nDofs,           &
          &                   subTree         = subTree,                  &
          ! &                   transientReduce = out_file%transientReduce, &
          &                   bary            = out_file%bary             )
      end if
    case(hvs_Internal)
       call tem_restart_dump_data( restart = out_file%restart,        &
         &                         varSys  = varSys,                  &
         &                         tree    = mesh,                    &
         &                         time    = out_file%time,           &
         &                         subTree = subTree )
         ! &                         transientReduce = out_file%transientReduce )
    case(hvs_VTK)
      do iVar=1,out_file%nVars
        call hvs_vtk_dump_data( vtk_file = out_file%vtk,          &
          &                     varpos   = out_file%varpos(iVar), &
          &                     varSys   = varSys,                &
          &                     mesh     = mesh,                  &
          &                     subtree  = subtree,               &
          &                     time     = out_file%time          )
      end do
    end select

  end subroutine hvs_output_write
  ! ----------------------------------------------------------------------------!
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  !> Close all files open for current time step.
  subroutine hvs_output_close(out_file, mesh, varSys, subTree)
    ! --------------------------------------------------------------------------!
    type(hvs_output_file_type), intent(inout) :: out_file

    !> Description of the available variable system to get the given varnames
    !! from.
    type(tem_varSys_type), intent(in) :: varsys

    !> Mesh to write the data on.
    type(treelmesh_type), intent(in) :: mesh

    !> Optional restriction of the elements to output.
    type(tem_subtree_type), optional, intent(in) :: subtree
    ! --------------------------------------------------------------------------!
    select case(out_file%vis_kind)
    case(hvs_AsciiSpatial)
      call hvs_asciiSpatial_close( asciiSpatial = out_file%asciiSpatial )
    case(hvs_Internal)
      call tem_restart_closeWrite( me      = out_file%restart, &
        &                          tree    = mesh,             &
        &                          timing  = out_file%time,    &
        &                          varSys  = varSys,         &
        &                          subTree = subTree           )
    case(hvs_VTK)
      call hvs_vtk_close( vtk_file = out_file%vtk, &
        &                 proc     = out_file%proc )
    end select

  end subroutine hvs_output_close
  ! ----------------------------------------------------------------------------!
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  !> Finalize the output
  subroutine hvs_output_finalize( out_file )
    ! --------------------------------------------------------------------------!
    type(hvs_output_file_type), intent(inout) :: out_file
    ! --------------------------------------------------------------------------!
    select case(out_file%vis_kind)
    case(hvs_AsciiTransient)
      call hvs_ascii_close( ascii = out_file%ascii)
    case(hvs_VTK)
      call hvs_vtk_closePVD(vtk_file = out_file%vtk)
      call tem_vrtx_finalize(out_file%vrtx)
    case(hvs_Internal)
      call tem_restart_finalize( me = out_file%restart )
    end select
  end subroutine hvs_output_finalize
  ! ----------------------------------------------------------------------------!
  ! ----------------------------------------------------------------------------!

end module hvs_output_module
