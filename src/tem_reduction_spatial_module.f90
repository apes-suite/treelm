! Copyright (c) 2011-2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2011-2013, 2015, 2017, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Michael Gaida  <michael.gaida@student.uni-siegen.de>
! Copyright (c) 2018 Jana Gericke <jana.gericke@student.uni-siegen.de>
! Copyright (c) 2020 Raphael Haupt <Raphael.Haupt@uni-siegen.de>
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
!> author: Vyacheslav Korchagin, Manuel Hasert
!! This module applies reductions to a set of quantities
!!
!! This module takes a set of input quantities and reduces them
!! according to the chosen reduceType.
!! You can choose between the following reduceTypes:
!!
!! - average: Average of all occurring values
!! - sum: Sum of all occurring values
!! - max: Maximum of all occurring values
!! - min: Minimum of all occurring values
!! - l2norm: L2 norm of all values, i.e. L_2(x) = sqrt( sum( x_i^2 ))
!! - linfnorm: L_inf norm of all values, i.e. L_inf(x) = max( | x_i | )
!!
!!```lua
!! reduction = {
!!               'norm' -- 'average' or 'sum'
!!             }
!!```
!!
!! The table of reductions has to correspond to the number of entries in the
!! tracking variable table. For each variable, specify a reduction as
!!```lua
!! tracking = {
!!   variable = {'pressure', 'velMag', 'density'},
!!   reduction = { 'average', 'max', l2norm' }
!! }
!!```
!!
!! Usage of this module
!! --------------------
!! The tracking module handles the reductions. In the solver, the reduction
!! must then be inkoved in the compute loop
!!
!! Embedding in Treelm
!! -------------------
!! Reductions are usually defined in tracking objects. Thus,
!! the tracking module handles the loading and closing of the reductions.
!! The following steps need to be performed within the tracking module
!!
!! - Treelm: Load the configuration from the lua file
!!   with tem_load_spatial_reduction
!! - Treelm: Initialize the reductions with tem_init_reduction
!! - Solver:
!!    - Open the reduction in the current time step tem_reduction_spatial_open
!!    - invoke the append reduceType to the array on which the reduction
!!       will be applied with tem_append_reduction
!!       This can be performed for a number of chunks
!!    - Close the reduction with tem_reduction_spatial_close
!!
module tem_reduction_spatial_module

  ! include treelm modules
  use mpi
  use env_module,          only: rk, rk_mpi, labelLen, long_k
  use tem_comm_env_module, only: tem_comm_env_type
  use treelmesh_module,    only: treelmesh_type
  use tem_tools_module,    only: upper_to_lower
  use tem_varSys_module,   only: tem_varSys_type
  use tem_logging_module,  only: logUnit
  use tem_aux_module,      only: tem_abort
  use tem_topology_module, only: tem_levelOF
  use tem_geometry_module, only: tem_elemSize

  ! include aotus modules
  use aotus_module,     only: flu_State, aot_get_val, aoterr_Fatal,            &
    &                         aoterr_NonExistent, aoterr_WrongType, flu_State
  use aot_table_module, only: aot_table_open, aot_table_close,                 &
    &                         aot_table_length
  use aot_out_module,   only: aot_out_type, aot_out_open,  aot_out_close,      &
    &                         aot_out_val, aot_out_open_table,                 &
    &                         aot_out_close_table

  implicit none

  private

  public :: tem_reduction_spatial_type
  public :: tem_reduction_spatial_config_type
  public :: tem_load_reduction_spatial
  public :: tem_reduction_spatial_init
  public :: tem_reduction_spatial_open
  public :: tem_reduction_spatial_append
  public :: tem_reduction_spatial_close
  public :: tem_reduction_spatial_toChunk
  public :: tem_reduction_spatial_dump
  public :: tem_reduction_spatial_out

  !> This data type is providing the input for the reduction routines
  !! It must be filled by the solver, before the reduction is called
  !! It exists on each process
  type tem_reduction_spatial_type

    !> amount of components of the quantity to reduce
    integer :: nComponents

    !> the result from the reduction operation
    !! size: nComponents
    real(kind=rk), allocatable :: val(:)

    !> how many elements have been included into the reduction (so far)
    integer :: nElems

    !> local part of total volume of intersected elements
    real(kind=rk) :: Vloc

    !> Which operation to perform on the list of elements
    character(len=labelLen) :: reduceType = ''

  end type tem_reduction_spatial_type

  type tem_reduction_spatial_config_type
    !> Which operation to perform on the list of elements
    character(len=labelLen), allocatable :: reduceType(:)

    logical :: active = .false.
  end type tem_reduction_spatial_config_type

  interface tem_reduction_spatial_dump
    module procedure tem_reduction_spatial_dump_vector
    module procedure tem_reduction_spatial_dump_single
  end interface tem_reduction_spatial_dump

  interface tem_reduction_spatial_out
    module procedure tem_reduction_spatial_out_vector
    module procedure tem_reduction_spatial_out_single
  end interface tem_reduction_spatial_out

contains

! ****************************************************************************** !
  !> read configuration file
  !!
  subroutine tem_load_reduction_spatial(conf, redSpatial_config, parent, key)
    ! ---------------------------------------------------------------------------
    !> handle for lua file
    type(flu_State) :: conf
    !> the reduction file to fill
    type(tem_reduction_spatial_config_type),intent(out) :: redSpatial_config
    !> handle for reduce table
    integer, optional,intent(in) :: parent
    !> which key to open
    character(len=*),optional,intent(in) :: key
    ! ---------------------------------------------------------------------------
    integer :: nEntries, handle, iPos
    character(len=labelLen) :: localKey
    integer, allocatable :: vErr(:) !, vErr_NonExistent(:)
    ! ---------------------------------------------------------------------------
    if( present( key )) then
      localKey = key
    else
      localKey = 'reduction'
    endif

    allocate( redSpatial_config%reduceType(1) )
    allocate( vErr(1)   )

    redSpatial_config%reduceType(1) = ''
    ! load reduction as scalar if failed then try to load as table
    call tem_load_reduction_single( conf       = conf,             &
      &                             handle     = parent,           &
      &                             reduceType = redSpatial_config &
      &                                         %reduceType(1),    &
      &                             key        = localKey,         &
      &                             iError     = vErr(1)           )

    ! try loading it as table
    if(btest(vErr(1), aoterr_NonExistent)) then
      ! write(logUnit(1),"(A)") 'Try load reduction as a table'
      call aot_table_open( L       = conf,      &
        &                  thandle = handle,    &
        &                  parent  = parent,    &
        &                  key     = localKey   )
      ! reduction defined as table
      if ( handle /= 0 ) then
        ! load entry inside a table
        nEntries = aot_table_length( L = conf, thandle = handle )
        ! write(logUnit(1),"(A,I0)") 'Table has entries: ', nEntries
        deallocate( redSpatial_config%reduceType )
        deallocate( vErr )
        allocate( redSpatial_config%reduceType( nEntries ) )
        allocate( vErr( nEntries ) )
        do iPos = 1, nEntries
          redSpatial_config%reduceType(iPos) = ''
          call tem_load_reduction_single( conf       = conf,              &
            &                             handle     = handle,            &
            &                             reduceType = redSpatial_config  &
            &                                         %reduceType(iPos),  &
            &                             pos        = iPos,              &
            &                             iError     = vErr(iPos)         )
        enddo
        redSpatial_config%active = .true.
      end if
      call aot_table_close(conf, handle)
    else
      ! write(logUnit(1),"(A)") 'Reduction is a single entry'
      redSpatial_config%active = .true.
      nEntries = 1
    endif

    if ( redSpatial_config%active ) then
      write(logUnit(3),"(A, I0)") '  Number of reductions loaded: ', nEntries
      write(logUnit(5),"(A)")     '  their reduceTypes are:'
      do iPos = 1, nEntries
        write(logUnit(5),"(A)")   '    '//trim(redSpatial_config &
          &                                     %reduceType(iPos))
      end do
    else
      deallocate( redSpatial_config%reduceType )
      allocate( redSpatial_config%reduceType(0) )
    end if


  end subroutine tem_load_reduction_spatial
! ****************************************************************************** !


! ****************************************************************************** !
  !> Read a single entry of reductions from the lua file
  !!
  subroutine tem_load_reduction_single(conf, reduceType, handle, key, pos, &
    &                                  iError)
    ! ---------------------------------------------------------------------------
    !> handle for lua file
    type(flu_State),intent(inout) :: conf
    !> reduction type to be filled
    character(len=labelLen), intent(out) :: reduceType
    !> handle for reduce table
    integer, optional,intent(in) :: handle
    !> which key to open
    character(len=*),optional,intent(in) :: key
    !> position to load from in config file
    integer,optional,intent(in) :: pos
    !> error from aotus
    integer, intent(out) :: iError
    ! ---------------------------------------------------------------------------
    if( present( key )) then
      call aot_get_val( L       = conf,       &
        &               thandle = handle,     &
        &               val     = reduceType, &
        &               ErrCode = iError,     &
        &               key     = key         )
    elseif( present( pos ))then
      call aot_get_val( L       = conf,       &
        &               thandle = handle,     &
        &               val     = reduceType, &
        &               ErrCode = iError,     &
        &               pos     = pos         )
    else
      iError = ibset(iError,aoterr_NonExistent)
    endif

    ! if reduceType is empty. can happen when reduceType is defined as table
    if (trim(reduceType) == '' ) iError = ibset(iError,aoterr_NonExistent)

    if(.not. btest(iError, aotErr_NonExistent)) then
      ! convert into lower case
      reduceType = upper_to_lower( reduceType )
      ! Check the chosen reduceType
      select case( trim(reduceType) )
        case('sum')
        case('average')
        case('l2norm','l2_norm')
          reduceType = 'l2norm'
        case('l2normalized')
          reduceType = 'l2normalized'
        case('linfnorm', 'linf_norm', 'l_inf_norm')
          reduceType = 'linfnorm'
        case('max','maximum')
          reduceType = 'max'
        case('min','minimum')
          reduceType = 'min'
        case('weighted_sum', 'w_sum')
          ! result weighted by its volume factor
          reduceType = 'weighted_sum'
        case('none')
        case default
          write(logUnit(1),*)' Error: The chosen reduction '//  &
            &                trim(reduceType)//' is not defined.'
          call tem_abort()
      end select
    endif

  end subroutine tem_load_reduction_single
! ****************************************************************************** !


! ****************************************************************************** !
  !> Initialize reduction objects according to the variable systems
  !!
  subroutine tem_reduction_spatial_init( me, redSpatial_config, varSys, varPos )
    ! ---------------------------------------------------------------------------
    !> array of reductions to initialize
    type( tem_reduction_spatial_type ), intent(out), allocatable :: me(:)
    type( tem_reduction_spatial_config_type ), intent(in) :: redSpatial_config
    !> global variable system defined in solver
    type( tem_varSys_type ), intent(in) :: varSys
    !> position of variable to reduce in the global varSys
    integer, intent(in) :: varPos(:)
    ! ---------------------------------------------------------------------------
    ! type( tem_reduction_spatial_type ), allocatable :: tmpRdc(:)
    integer :: iRdc, nVars
    ! ---------------------------------------------------------------------------

    write(logUnit(1),*) 'Initializing reduction ...'

    nVars = size(varPos)

    ! check if the amount of reduction greater than amount of variables
    if( size( redSpatial_config%reduceType ) > nVars ) then
      write(logUnit(1),*) 'The number of reduction is more than the '// &
        &                 'number of variables!'
      write(logUnit(1),*) 'Ignore the additional reductions!'
    endif

      ! allocate( tmpRdc( nVars ))
      ! do iRdc = 1, nVars
      !   tmpRdc( iRdc ) = me( iRdc )
      ! end do
      ! deallocate( me )
      ! allocate( me( nVars ))
      ! do iRdc = 1, nVars
      !   me( iRdc ) = tmpRdc( iRdc )
      ! end do
      ! deallocate( tmpRdc )


    allocate( me( nVars ) )
    do iRdc = 1, nVars
      me( iRdc )%reduceType = trim( redSpatial_config%reduceType( iRdc ) )
      me( iRdc )%nComponents = varSys%method%val( varPos(iRdc) )%nComponents
      allocate( me( iRdc )%val( me(iRdc)%nComponents ))
      me( iRdc )%val(:) = 0._rk
    end do

  end subroutine tem_reduction_spatial_init
! ****************************************************************************** !


! ****************************************************************************** !
  !> Prepare the reduction data type
  !!
  !! @todo HK: what is going on here?
  !!           This stuff reads, as if we have a reduction per variable?
  !!           However above, we read a set of reduction objects from the
  !!           configuration, thus it seems to be independent of the number
  !!           of variables? Somehow there is a missing link between the
  !!           reduction and the variable it should act on?
  !!           How is this supposed to work? Please clarify, document or
  !!           change the code!
  !! MH: Ok, the current way it is intended to work is the following
  !!  you define a tracking object with a given set of variables
  !!  and define a number of reduction operations.
  !!  The reduction operations work on the variables starting from the first
  !!  up to the number of defined reductions
  !!  example:
  !!  tracking = { variable = {'density', 'velocity'},
  !!               reduction = 'average', ... }
  !!  Here the number of reductions is smaller than the number of variables
  !!  which results in that only the density is taken an average.
  !!
  subroutine tem_reduction_spatial_open(me, varSys, varPos)
    ! ---------------------------------------------------------------------------
    !> The reduction type to work on. All definitions should be
    !! present in here
    type( tem_reduction_spatial_type ), intent(inout)    :: me(:)
    !> global variable system defined in solver
    type( tem_varSys_type ), intent(in)          :: varSys
    integer, intent(in)                          :: varPos(:)
    ! ---------------------------------------------------------------------------
    integer :: iReduce
    ! ---------------------------------------------------------------------------

    do iReduce = 1, size( me )
      ! Reset the number of treated elements.
      ! Must be counted up during _append
      me( iReduce )%nElems = 0

      ! Take the components of the connected variable system
      me( iReduce )%nComponents = varSys%method%val(varPos(iReduce))%nComponents
      if( .not. allocated( me( iReduce )%val ))                                &
        & allocate( me( iReduce )%val( me(iReduce)%nComponents ))

      select case( trim(me( iReduce )%reduceType) )
        case('max')
          me( iReduce )%val(:) = -huge(1._rk)
        case('min')
          me( iReduce )%val(:) = huge(1._rk)
        case default ! ('sum', 'average', 'l2norm', 'none', ...
          me( iReduce )%val(:) = 0._rk
      end select
    end do

  end subroutine tem_reduction_spatial_open
! ****************************************************************************** !


! ****************************************************************************** !
  !> Local chunk-wise reduction
  !!
  !! Chunk-wise local reduction, which has to be performed for each
  !! chunk until all elements have been treated within the reduction.
  !! After all chunks are treated (=_append has been finished),
  !! the global reduction has to be performed by tem_reduction_spatial_close.
  !! NOTE: Reduction is applied only on 1st DOF
  subroutine tem_reduction_spatial_append(me, chunk, nElems, varSys, varPos, &
    &                                     tree, treeID, nDofs)
    ! ---------------------------------------------------------------------------
    !> The reduction type to work on. All definitions should be
    !! present in here
    type( tem_reduction_spatial_type ), intent(inout) :: me(:)
    !> number of elements the chunk has
    integer, intent(in) :: nElems
    !> global variable system defined in solver
    type( tem_varSys_type ), intent(in)          :: varSys
    !> position of variable to reduce in the global varSys
    integer, intent(in)                          :: varPos(:)
    !> chunk of results to reduce
    real(kind=rk), intent(in) :: chunk(:)
    !> Number of degrees of freedom.
    integer, optional, intent(in) :: nDofs
    !> the global tree
    type(treelmesh_type), intent(in) :: tree
    !> The list of treeIDs of the current chunk
    integer(kind=long_k), optional, intent(in) :: treeID( nElems )
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: temp
    real(kind=rk) :: dx, V, Vloc
    integer :: iVar, iComp, iPos, chunkVarPos, iElem, nDofs_L
    integer :: nPerElem
    integer :: minLevel, myLevel, volumeFac, nScalars
    ! ---------------------------------------------------------------------------
    minLevel = tree%global%minLevel
    nScalars = sum(varSys%method%val(varPos(:))%nComponents)

    ! spatial reduction is applied only on 1st Degrees of freedom
    if (present(nDofs)) then
      nDofs_L = nDofs
    else
      nDofs_L = 1
    end if
    nPerElem = nDofs_L*nScalars

    chunkVarPos = 0
    ! Run over all variables in the chunk
    do iVar = 1, size( me )
      me( iVar )%nElems = me( iVar )%nElems + nElems
      ! count up the me%nElems here as well
      ! update the me%val
      do iComp = 1, me( iVar )%nComponents
        chunkVarPos = chunkVarPos + 1
        temp = 0._rk
        select case( trim(me( iVar )%reduceType) )
          ! Upon extending the reductions, make sure to add them in the
          ! tem_load_spatial_reduction  routine
          case('sum', 'average')
            ! Take the values corresponding to the current component of the
            ! variable with step (nScalars)
            temp = sum( chunk( chunkVarPos : nElems*nPerElem &
              &                            : nPerElem )      )
            me( iVar )%val( iComp ) = me( iVar )%val( iComp ) + temp

          case('l2normalized')
            Vloc = 0.0_rk
            ! Sum up the squares of elemets of chunk into a temp variable
            ! Total error L2 on the domain normalize by Volum or number of points
            if (present(treeID)) then
              do iPos = chunkVarPos, nElems*nPerElem, nPerElem
                iElem  = ((iPos-1)/nPerElem) + 1
                dx = tem_Elemsize(tree, treeID(iElem))
                V = dx**3
                Vloc = Vloc + V

                temp = temp + chunk(iPos)*chunk(iPos) * V ! real(volumeFac,kind=rk)
              end do
            else
              do iPos = chunkVarPos, nElems*nPerElem, nPerElem
                temp = temp + chunk(iPos)*chunk(iPos)
                Vloc = Vloc + 1
              end do
            end if
            me( iVar )%val( iComp ) = me( iVar)%val( iComp ) + temp
            me( iVar )%Vloc = Vloc

          case('l2norm')
            ! Sum up the squares of elemets of chunk into a temp variable
            ! Total error L2 on the domain Omega is with the error function u
            ! $L2(\Omega) =  \sqrt(\int_\Omega |u|^2 ))
            !             =  \sqrt(\sum_i \int_\Omega_i |u|^2 ))
            !             =  \sqrt(\sum_i \Omega_i \cdot |u|^2 ))
            ! Value of elements on level other than minLevel is scale by a
            ! factor of 8^{myLevel-minLevel} to obtain a consistant results for
            ! multi-level.
            if (present(treeID)) then
              do iPos = chunkVarPos, nElems*nPerElem, nPerElem
                iElem  = ((iPos-1)/nPerElem) + 1
                ! get level, then get scale factor
                !myLevel = tem_levelOf( treeID(iElem) )
                !volumeFac = 8**( myLevel - minLevel )
                dx = tem_Elemsize(tree, treeID(iElem))
                V = dx**3

                temp = temp + chunk(iPos)*chunk(iPos) * V ! real(volumeFac,kind=rk)
              end do
            else
              do iPos = chunkVarPos, nElems*nPerElem, nPerElem
                temp = temp + chunk(iPos)*chunk(iPos)
              end do
            end if
            me( iVar )%val( iComp ) = me( iVar)%val( iComp ) + temp

          case('linfnorm')
            ! L_inf norm evaluates the maximum of the absolute occurring value
            me( iVar )%val( iComp ) = max( me( iVar )%val( iComp ),   &
              & maxval( abs(chunk( chunkVarPos : nElems*nPerElem      &
              &                                : nPerElem         ) )))

          case('max')
            ! maximum between previous value and current
            me( iVar )%val( iComp ) = max(me( iVar )%val( iComp ), &
              & maxval( chunk( chunkVarPos : nElems*nPerElem       &
              &                            : nScalars*nDofs_L ) )  )

          case('min')
            ! minimum between previous value and current
            me( iVar )%val( iComp ) = min(me(iVar)%val(iComp),  &
              & minval( chunk( chunkVarPos : nElems*nPerElem    &
              &                            : nPerElem         )))

          case('weighted_sum')
            ! Sum up scaled results in chunk into a temp variable
            ! Value of elements on level other than minLevel is scale by a
            ! factor of 8^{myLevel-minLevel} to obtain a consistant results for
            ! multi-level.
            if (present(treeID)) then
              do iPos = chunkVarPos, nElems*nPerElem, nPerElem
                iElem  = ((iPos-1)/nPerElem) + 1
                ! get level, then get scale factor
                myLevel = tem_levelOf( treeID(iElem) )
                volumeFac = 8**( myLevel - minLevel )
                temp = temp + chunk(iPos) / real(volumeFac,rk)
              end do
            else
              temp = sum( chunk( chunkVarPos : nElems*nPerElem &
                &                            : nPerElem )      )
            end if
            me( iVar )%val( iComp ) = me( iVar)%val( iComp ) + temp

          case('none')
            me( iVar )%val( iComp ) = 0._rk

        end select
      enddo
    enddo

  end subroutine tem_reduction_spatial_append
! ****************************************************************************** !


! ****************************************************************************** !
  !> Perform the global reduction
  !!
  !! After the local reductions have been performed (in _append),
  !! the results must be communicated between processes.
  !!
  subroutine tem_reduction_spatial_close(me, proc)
    ! ---------------------------------------------------------------------------
    !> The reduction type to work on. All definitions should be
    !! present in here
    type( tem_reduction_spatial_type ), intent(inout) :: me(:)
    !> communicator for processes participating in this reduction
    type(tem_comm_env_type), intent(in)  :: proc
    ! ---------------------------------------------------------------------------
    integer :: i, nComp, ierr, globalnElems
    real(kind=rk), allocatable :: buff(:)
    real(kind=rk)              :: Vglob
    ! ---------------------------------------------------------------------------
    globalnElems = 0

    !loop over all tracking objects
    do i = 1, size(me)

      ! get number of components
      nComp = me(i)%nComponents

      allocate(buff(nComp))
      buff = 0.0_rk

      !choose reduction operation and perform it
      select case( me( i )%reduceType )
        !sum all values
        case('sum', 'weighted_sum')
          call mpi_reduce( me(i)%val, buff,                                    &
            &              nComp, rk_mpi, mpi_sum, proc%root, proc%comm, iErr)
          me(i)%val = buff

        !sum all values and devide by number of elements
        case('average')
          globalnElems = 0
          ! get global number of elements in this reduction to rank 0
          call mpi_reduce( me(i)%nElems, globalnElems,                        &
            &              1, mpi_integer, mpi_sum, proc%root, proc%comm, iErr)
          call mpi_reduce( me(i)%val, buff,                                  &
            &              nComp, rk_mpi, mpi_sum, proc%root, proc%comm, iErr)
          if (proc%rank == proc%root) &
            & me(i)%val(:) = buff / real( globalnElems, kind = rk)

        !sum all values(sum of squares) and extract a square root
        case('l2norm')
          call mpi_reduce( me(i)%val, buff,                                    &
          &                nComp, rk_mpi, mpi_sum, proc%root, proc%comm, iErr)
          me(i)%val(:) = sqrt(buff)

        !maximium over all values
        case('max','linfnorm')
          call mpi_reduce( me(i)%val, buff,                                    &
          &                nComp, rk_mpi, mpi_max, proc%root, proc%comm, iErr)
          me(i)%val = buff

        !minimum over all values
        case('min')
          call mpi_reduce( me(i)%val, buff,                                    &
          &                nComp, rk_mpi, mpi_min, proc%root, proc%comm, iErr)
          me(i)%val = buff

        !sum all values(sum of squares), normalize and extract a square root
        case('l2normalized')
          call mpi_reduce( me(i)%Vloc, Vglob,                            &
          &                1, rk_mpi, mpi_sum, proc%root, proc%comm, iErr)
          call mpi_reduce( me(i)%val, buff,                                  &
          &                nComp, rk_mpi, mpi_sum, proc%root, proc%comm, iErr)
          me(i)%val(:) = sqrt(buff/Vglob)

        case default
      end select
      deallocate( buff )
    enddo

  end subroutine tem_reduction_spatial_close
! ****************************************************************************** !


! ****************************************************************************** !
  !> Transfer reduction results to array chunk
  !!
  subroutine tem_reduction_spatial_toChunk( me, chunk, nChunkElems )
    ! ---------------------------------------------------------------------------
    !>
    type(tem_reduction_spatial_type), intent(in) :: me(:)
    !>
    real(kind=rk), intent(inout) :: chunk(:)
    !> Number of element after spatial reduction = 1
    integer, intent(out) :: nChunkElems
    ! ---------------------------------------------------------------------------
    integer :: iReduce, iPos, iComp
    ! ---------------------------------------------------------------------------

    nChunkElems = 1
    iPos = 0
    do iReduce = 1, size( me )
      do iComp = 1, me(iReduce)%nComponents
        chunk( iPos+iComp ) = me(iReduce)%val( iComp )
      enddo
      iPos = iPos + me(iReduce)%nComponents
    enddo

  end subroutine tem_reduction_spatial_toChunk
! ****************************************************************************** !


  ! **************************************************************************** !
  !> Dumps array of reduction to given unit
  subroutine tem_reduction_spatial_dump_vector(me, outUnit)
    ! ---------------------------------------------------------------------------
    !> reduction to write into the lua file
    type(tem_reduction_spatial_type), intent(in) :: me(:)
    !> unit to write to
    integer, intent(in) :: outUnit
    ! ---------------------------------------------------------------------------
    ! aotus type handling the output to the file in lua format
    type(aot_out_type) :: conf
    ! ---------------------------------------------------------------------------
    call aot_out_open( put_conf = conf, outUnit = outUnit )
    call tem_reduction_spatial_out_vector( me, conf )
    call aot_out_close( put_conf = conf )

  end subroutine tem_reduction_spatial_dump_vector
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Dump single reduction to given unit
  subroutine tem_reduction_spatial_dump_single(me, outUnit)
    ! ---------------------------------------------------------------------------
    !> reduction to write into the lua file
    type(tem_reduction_spatial_type), intent(in) :: me
    !> unit to write to
    integer, intent(in) :: outUnit
    ! ---------------------------------------------------------------------------
    ! aotus type handling the output to the file in lua format
    type(aot_out_type) :: conf
    ! ---------------------------------------------------------------------------
    call aot_out_open( put_conf = conf, outUnit = outUnit )
    call tem_reduction_spatial_out_single( me, conf )
    call aot_out_close( put_conf = conf )

  end subroutine tem_reduction_spatial_dump_single
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Allows the output of array of reduction to lua out
  subroutine tem_reduction_spatial_out_vector(me, conf)
    ! ---------------------------------------------------------------------------
    !> reduction to write into the lua file
    type(tem_reduction_spatial_type), intent(in) :: me(:)
    !> aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    ! ---------------------------------------------------------------------------
    call aot_out_val( put_conf = conf,             &
      &               val      = me(:)%reduceType, &
      &               vname    = 'reduction'       )

  end subroutine tem_reduction_spatial_out_vector
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Allows the output of the single reduction to lua out.
  !!
  !! The data is written into the file, the lunit is connected to.
  !! It is formatted as a Lua table.
  !!
  subroutine tem_reduction_spatial_out_single(me, conf)
    ! ---------------------------------------------------------------------------
    !> reduction to write into the lua file
    type(tem_reduction_spatial_type), intent(in) :: me
    !> aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    ! ---------------------------------------------------------------------------
    call aot_out_val( put_conf = conf,          &
      &               val      = me%reduceType, &
      &               vname    = 'reduction'    )
    ! ---------------------------------------------------------------------------

  end subroutine tem_reduction_spatial_out_single
  ! **************************************************************************** !


end module tem_reduction_spatial_module
! ****************************************************************************** !
