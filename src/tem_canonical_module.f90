! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012-2013, 2015, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2013, 2016, 2019, 2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012, 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2013, 2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
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
! **************************************************************************** !
!> author: Simon Zimny
!! author: Kartik Jain
!! This module provides the user with a simple geometrical object like point,
!! line, plane and box.
!!
!! Detail description of the canonical shapes can be found
!! in the Documentation.
!!
module tem_canonicalND_module

  use mpi
  ! include treelm modules
  use env_module,            only: rk, rk_mpi, long_k, globalMaxLevels, labelLen
  use tem_float_module,      only: operator(.fne.)
  use tem_tools_module,      only: upper_to_lower
  use tem_aux_module,        only: tem_abort
  use tem_param_module,      only: qQQQ, qOffset
  use treelmesh_module,      only: treelmesh_type
  use tem_geometry_module,   only: tem_CoordOfReal, tem_PosOfId, tem_baryOfId
  use tem_topology_module,   only: tem_FirstIdAtLevel, tem_CoordOfId
  use tem_dyn_array_module,  only: dyn_intArray_type, append
  use tem_topology_module,   only: tem_levelOf, tem_IDofCoord
  use tem_logging_module,    only: logUnit, tem_toStr
  use tem_pointData_module,  only: tem_grwPoints_type, init, append
  use tem_cube_module,       only: tem_cube_type, tem_convertTreeIDtoCube
  use tem_point_module,      only: tem_point_type, tem_pointCubeOverlap
  use tem_line_module,       only: tem_line_type, tem_lineCubeOverlap
  use tem_plane_module,      only: tem_plane_type, tem_planeCubeOverlap, &
    &                              tem_createPlane
  use tem_box_module,        only: tem_box_type, tem_boxCubeOverlap, &
    &                              tem_createBox
  use tem_transformation_module, only: tem_transformation_type
  use tem_debug_module,      only: dbgUnit

  ! include aotus modules
  use aotus_module,     only: aot_get_val, aoterr_Fatal, aoterr_WrongType,     &
    &                         flu_State
  use aot_table_module, only: aot_table_open, aot_table_close, aot_table_length
  use aot_out_module,   only: aot_out_type, aot_out_val,                       &
    &                         aot_out_open_table, aot_out_close_table

  implicit none
  private

  public :: tem_canonicalND_type
  public :: tem_load_canonicalND
  public :: tem_getNextCoordOfcanonicalND
  public :: tem_canonicalND_out
  public :: tem_cano_initSubTree
  public :: tem_cano_checkNeigh
  public :: tem_cano_storePntsInSubTree
  public :: tem_transformCanoND

  !> Definition of the canonicalND
  type tem_canonicalND_type

    !> origin of the canonical shape
    real(kind=rk) :: origin(3)

    !> vector along the edge A (also defines size)
    !! 1st dimension defines x,y, z coord
    !! 2nd dimension vec number
    real(kind=rk) :: vec(3,3)

    !> how many discrete points the canonicalND is divided into
    integer       :: segments(3)

    !> spatial distribution of the points
    character(len=labellen) :: distribution

    !> kind of canonicalND (line, plane, point, box)
    character( len=labellen ) :: kind

    !> identify which vectors are active (not equal 0)
    logical :: active(3)

    !> dimension of canonical object
    !! nDim=0 - point
    !! nDim=1 - line
    !! nDim=2 - plane
    !! nDim=3 - box
    integer :: nDim

    !> To choose what to do with intersection of box
    !! if only_surface = true than the only the surface of the object
    !! is intersected
    !! if only_surface = false then the whole object is intersected
    !! default is set to false
    logical :: only_surface = .false.

    !> canonical point
    type(tem_point_type) :: point
    type(tem_line_type) :: line !< canonical line
    type(tem_plane_type) :: plane !< canonical plane
    type(tem_box_type) :: box !< canonical box
  end type tem_canonicalND_type

  !> interface to write out canonical shape(s) in lua format to a file
  interface tem_canonicalND_out
    module procedure tem_canonicalND_out_scal
    module procedure tem_canonicalND_out_vec
  end interface tem_canonicalND_out

  !> interface to load canonical objects
  interface tem_load_canonicalND
    module procedure tem_load_oneCanonicalND
    module procedure tem_load_canonicalND_vec
  end interface tem_load_canonicalND

  !> This routine apply transformations to canonical objects
  interface tem_transformCanoND
    module procedure transformCanoND
    module procedure transformCanoND_single
  end interface tem_transformCanoND

contains

! **************************************************************************** !
  !> Loading canonicalNDs from the config file valid definitions:
  !!
  subroutine tem_load_canonicalND_vec(me, transform, conf, thandle, reqSegments)
    ! --------------------------------------------------------------------------
    !> the array of canonical objects which to read in (and first allocate)
    type(tem_canonicalND_type), allocatable, intent(out) :: me(:)
    !> transformation for spatial object
    type(tem_transformation_type), intent(in) :: transform
    !> lua config handle
    type(flu_state) :: conf
    !> table handle from which to read
    integer, intent(in) :: thandle
    !> Is true if use_get_point is true in output table
    logical, optional, intent(in) :: reqSegments
    ! --------------------------------------------------------------------------
    ! lua handles
    integer :: canonicalNDs_thandle, single_thandle
    integer :: canonicalND_entries
    integer :: iCan
    integer :: ncanonicalNDs
    ! --------------------------------------------------------------------------

    write(logUnit(5),*) 'Trying to read canoND shape...'


    call aot_table_open( L       = conf,                                       &
      &                  parent  = thandle,                                    &
      &                  thandle = canonicalNDs_thandle,                       &
      &                  key     = 'object' )
    canonicalND_entries = aot_table_length( L       = conf,                    &
      &                                     thandle = canonicalNDs_thandle )

    if (canonicalNDs_thandle == 0) then
      write(logUnit(1),*) 'Error: object table is not defined for canoND'
      call tem_abort()
    end if

    ! Now check, if the first entry is a table itself:
    call aot_table_open( L       = conf,                                       &
      &                  parent  = canonicalNDs_thandle,                       &
      &                  thandle = single_thandle,                             &
      &                  pos     = 1 )

    if (single_thandle == 0) then
      ! The entries are not tables themselves, try to interpret it as
      ! single canonical object.
      call aot_table_close( L = conf, thandle = single_thandle )
      ncanonicalNDs = 1
      allocate(me( ncanonicalNDs ))
      call tem_load_oneCanonicalND( me          = me(1),                &
        &                           conf        = conf,                 &
        &                           transform   = transform,            &
        &                           thandle     = canonicalNDs_thandle, &
        &                           reqSegments = reqSegments           )
    else
      ! First entry is a table, assume all others are also tables,
      ! and properly define canonicalNDs coordinates, go on reading them.
      call aot_table_close( L = conf, thandle = single_thandle )
      ncanonicalNDs = canonicalND_entries
      allocate( me( ncanonicalNDs ))
      do iCan=1, ncanonicalNDs
        write(logUnit(5),"(A,I0,A)") 'Read Object: ', iCan, ' with '
        call aot_table_open( L       = conf,                                   &
          &                  parent  = canonicalNDs_thandle,                   &
          &                  thandle = single_thandle,                         &
          &                  pos     = iCan )
        call tem_load_oneCanonicalND( me          = me( iCan ),     &
          &                           conf        = conf,           &
          &                           transform   = transform,      &
          &                           thandle     = single_thandle, &
          &                           reqSegments = reqSegments     )
        call aot_table_close( L = conf, thandle = single_thandle )
      end do
    end if

    call aot_table_close( L = conf, thandle = canonicalNDs_thandle )

  end subroutine tem_load_canonicalND_vec
! **************************************************************************** !


! **************************************************************************** !
  !> Read one canonical object definition into a tem_canonicalND_type from a lua
  !! table.
  !!
  subroutine tem_load_oneCanonicalND(me, transform, conf, thandle, reqSegments)
    ! --------------------------------------------------------------------------
    !> contains canonicalND data
    type(tem_canonicalND_type), intent(out) :: me
    !> transformation for spatial object
    type(tem_transformation_type), intent(in) :: transform
    !> lua state
    type(flu_state) :: conf
    !> lua table identification
    integer, intent(in) :: thandle
    !> Is true if use_get_point is true in output table
    logical, optional, intent(in) :: reqSegments
    !> output messages?
    ! --------------------------------------------------------------------------
    ! error identifiers
    integer :: iError
    integer :: vError(3), errfatal(3)
    ! counter and total number of vectors
    integer :: iVec, nVecs
    ! lua handles
    integer :: vec_thandle, halfvec_thandle, single_thandle, halfsingle_thandle
    real(kind=rk), allocatable :: local_vecs(:,:)
    real(kind=rk) :: local_length
    real(kind=rk) :: vecAbsSum
    character(len=labelLen) :: buffer
    logical :: reqSegments_loc
    ! --------------------------------------------------------------------------

    ! Load segments only if reqSegments is true i.e. use_get_point in output 
    ! table is true
    if (present(reqSegments)) then
      reqSegments_loc = reqSegments
    else
      reqSegments_loc = .false.
    end if

    errfatal = aotErr_Fatal

    ! initialize the array of vectors and the array of segments
    me%vec = 0._rk

    ! open the vec and halvec table and get its size
    call aot_table_open( L       = conf,                                       &
      &                  parent  = thandle,                                    &
      &                  thandle = vec_thandle,                                &
      &                  key     = 'vec' )
    call aot_table_open( L       = conf,                                       &
      &                  parent  = thandle,                                    &
      &                  thandle = halfvec_thandle,                            &
      &                  key     = 'halfvec' )
    ! check if the vec table is given
    if( vec_thandle .ne. 0 ) then
      nVecs = aot_table_length( L = conf, thandle = vec_thandle )

      ! Now check, if the first entry is a table itself:
      call aot_table_open( L       = conf,                                     &
        &                  parent  = vec_thandle,                              &
        &                  thandle = single_thandle,                           &
        &                  pos     = 1 )

      if (single_thandle == 0) then
        ! The entries are not tables themselves, try to interpret it as
        ! single vector.
        call aot_table_close( L = conf, thandle = single_thandle )
        nVecs = 1
        allocate(local_vecs( 3, nVecs ))
        local_vecs = 0._rk
        call aot_get_val( L       = conf,                                      &
          &               thandle = thandle,                                   &
          &               val     = local_vecs(:,1),                           &
          &               key     = 'vec',                                     &
          &               ErrCode = vError )
        if( any( btest( vError, errFatal ))) then
          write(logUnit(1),*)' FATAL ERROR: in loading vec'
          write(logUnit(1),*)' You should provide 3 numbers for this vector.'
          call tem_abort()
        end if
      else
        ! First entry is a table, assume all others are also tables,
        ! and properly define canonicalNDs coordinates, go on reading them.
        call aot_table_close( L = conf, thandle = single_thandle )
        allocate( local_vecs( 3, nVecs ))
        local_vecs = 0._rk
        do iVec=1, nVecs
          call aot_get_val( L       = conf,                                    &
            &               thandle = vec_thandle,                             &
            &               val     = local_vecs(:,iVec),                      &
            &               pos     = iVec,                                    &
            &               ErrCode = vError )
        end do
      end if
      if( any( btest( vError, errFatal ))) then
        write(logUnit(1),*) ' FATAL ERROR: in loading vec'
        call tem_abort()
      endif
    ! if vec table is not set try to read halfvec
    else if( halfvec_thandle /= 0 ) then
      call aot_table_open( L       = conf,                                     &
        &                  parent  = halfvec_thandle,                          &
        &                  thandle = halfsingle_thandle,                       &
        &                  pos     = 1 )
      nVecs = aot_table_length( L = conf, thandle = halfvec_thandle )
      if( halfsingle_thandle == 0 ) then
        ! The entries are not tables themselves, try to interpret it as
        ! single vector.
        call aot_table_close( L = conf, thandle = halfsingle_thandle )
        nVecs = 1
        allocate( local_vecs( 3, nVecs ))
        local_vecs = 0._rk
        call aot_get_val( L       = conf,                                      &
          &               thandle = thandle,                                   &
          &               val     = local_vecs(:,1),                           &
          &               key     = 'halfvec',                                 &
          &               ErrCode = vError )
        local_vecs = local_vecs * 2._rk
      else
        ! First entry is a table, assume all others are also tables,
        ! and properly define canonicalNDs coordinates, go on reading them.
        call aot_table_close( L = conf, thandle = halfsingle_thandle )
        allocate( local_vecs( 3, nVecs ))
        local_vecs = 0._rk
        do iVec=1, nVecs
          call aot_get_val( L       = conf,                                    &
            &               thandle = halfvec_thandle,                         &
            &               val     = local_vecs(:,iVec),                      &
            &               pos     = iVec,                                    &
            &               ErrCode = vError )
        end do
        local_vecs = local_vecs * 2._rk
      end if
      if( any( btest( vError, errFatal ))) then
        write(logUnit(1),*) ' FATAL ERROR: in loading halfvec'
        call tem_abort()
      endif
    ! if vec and halfvec table is not set, try to read length (for a cube)
    else
      call aot_get_val( L       = conf,                                        &
        &               thandle = thandle,                                     &
        &               val     = local_length,                                &
        &               key     = 'length',                                    &
        &               default = 0.0_rk,                                      &
        &               ErrCode = iError )
      ! set the axis parallel entries of local_vecs to local_length
      nVecs = 3
      allocate( local_vecs( 3, nVecs ))
      local_vecs = 0._rk
      local_vecs(1,1) = local_length
      local_vecs(2,2) = local_length
      local_vecs(3,3) = local_length
    end if

    call aot_table_close( L = conf, thandle = vec_thandle )
    call aot_table_close( L = conf, thandle = halfvec_thandle )


    ! copy back the information from the local vector array
    me%vec(:,:nVecs) = local_vecs
    me%vec(:,nVecs+1:) = 0.0_rk

    ! read origin of the canonicalND
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               val     = me%origin,                                     &
      &               key     = 'origin',                                      &
      &               ErrCode = vError )
    if( any( btest( vError, errFatal ))) then
      call aot_get_val( L       = conf,                                        &
        &               thandle = thandle,                                     &
        &               val     = me%origin,                                   &
        &               key     = 'center',                                    &
        &               ErrCode = vError )
      if (any(btest( vError, errFatal )) ) then
        write(logUnit(1),*) 'ERROR reading the canonical objects, '//          &
          &             'origin/center is not well defined. STOPPING'
        call tem_abort()
      endif
      me%origin = me%origin - me%vec(:,1)/2.0_rk - me%vec(:,2)/2.0_rk -        &
        &         me%vec(:,3)/2.0_rk
    end if

    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               val     = me%distribution,                               &
      &               ErrCode = iError,                                        &
      &               key     = 'distribution',                                &
      &               default = 'equal' )

    me%active = .false.
    me%nDim = 0

    ! check the kind and the active vectors
    do iVec = 1, 3
      vecAbsSum = abs(me%vec(1,iVec)) + abs(me%vec(2,iVec)) &
        &       + abs(me%vec(3,iVec))
      if (vecAbsSum .fne. 0._rk) then
        me%active(iVec) = .true.
        me%nDim = me%nDim + 1
      end if
    end do

    if(me%nDim == 0) me%kind = 'point'
    if(me%nDim == 1) me%kind = 'line'
    if(me%nDim == 2) me%kind = 'plane'
    if(me%nDim == 3) me%kind = 'box'

    ! Load segments for objects other than point if reqSegments is true
    if (trim(me%kind) == 'point') then
      ! set default segments to 1
      me%segments = 1
    else
      if (reqSegments_loc) then
        ! set segments to 1 in all dimensions as default
        me%segments = 1
        ! try loading segments as integer if fails try
        ! to loading segments as vector integer
        call aot_get_val( L       = conf,           &
          &               thandle = thandle,        &
          &               val     = me%segments(1), &
          &               key     = 'segments',     &
          &               ErrCode = iError          )
        if ( btest(iError, aoterr_Fatal) ) then
          call aot_get_val( L       = conf,                &
            &               thandle = thandle,             &
            &               val     = me%segments(:nVecs), &
            &               key     = 'segments',          &
            &               ErrCode = vError(:nVecs)       )
          if ( any( btest(vError(:nVecs), errFatal(:nVecs)) ) ) then
            write(logUnit(1),*) ' FATAL ERROR: in loading segments'
            write(logUnit(1), '(A)') '  "segments" are not defined in ' &
              &                      // 'shape.object table.'
            call tem_abort('Segments are required for output.use_get_point.')
          end if
        end if
      else
        ! segments are not required setting it to -1 for error detection
        me%segments = -1
      end if
    end if

    ! if canonical object is box then check for only_surface
    if( me%kind == 'box' ) then
      call aot_get_val( L       = conf,                                        &
        &               thandle = thandle,                                     &
        &               val     = me%only_surface,                             &
        &               ErrCode = iError,                                      &
        &               key     = 'only_surface',                              &
        &               pos     = 4,                                           &
        &               default = .false. )

      if (btest(iError, aoterr_WrongType)) then
        write(logUnit(1),*) 'Error occured, while retrieving canonical box '// &
          &             'only_surface'
        write(logUnit(1),*) 'Variable has wrong type!'
        write(logUnit(1),*) 'Should be a LOGICAL!'
        call tem_abort()
      endif
    endif

    ! if tranformation is defined then tranform canoND before converting them
    ! into their respective shapes
    call tem_transformCanoND(canoND = me, transform = transform)
    ! Convert canonical shape to respective kind
    select case(upper_to_lower(trim(me%kind)))
    case('point')
      me%point%coord = me%origin
    case('line')
      me%line%origin = me%origin
      me%line%vec = me%vec(:,1)
    case('plane')
      call tem_createPlane(me     = me%plane,    &
        &                  origin = me%origin,   &
        &                  vecA   = me%vec(:,1), &
        &                  vecB   = me%vec(:,2)  )
    case('box')
      call tem_createBox(me           = me%box,        &
        &                origin       = me%origin,     &
        &                vecA         = me%vec(:,1),   &
        &                vecB         = me%vec(:,2),   &
        &                vecC         = me%vec(:,3),   &
        &                only_surface = me%only_surface)
    case default
      call tem_abort('Unknown canonical kind')
    end select

    ! @todo: DEB: Update this when converting arrays to strings is available
    write(logUnit(6),"(A)") ' kind: '//trim( me%kind )
    write(logUnit(6),*) ' origin:', me%origin(1:3)
    do iVec=1,me%nDim
      write(buffer,"(A,I0,A)") ' vec ', iVec, ': '
      write(logUnit(6),*) trim(buffer), me%vec(:,iVec)
    enddo
    if (trim(me%kind) == 'plane') then
      write(logUnit(6),*) '   unit normal: '                             &
        &             // trim(tem_toStr(me%plane%unitNormal(1)))         &
        &             // ', ' // trim(tem_toStr(me%plane%unitNormal(2))) &
        &             // ', ' // trim(tem_toStr(me%plane%unitNormal(3)))
    end if
    write(logUnit(6),*) ' only_surface: ', me%only_surface
    do iVec=1,me%nDim
      write(buffer,"(A,I0,A)") 'segments ', iVec, ': '
      write(logUnit(6),"(A,I0)") trim(buffer), me%segments(iVec)
    enddo
    write(logUnit(6),"(A)") ' distribution: '//trim( me%distribution )

  end subroutine tem_load_oneCanonicalND
! **************************************************************************** !


  ! ****************************************************************************
  !> This routine apply transformation to canonical objects.
  subroutine transformCanoND(canoND, transform)
    !--------------------------------------------------------------------------!
    !> canonical geometry object type
    type( tem_canonicalND_type ), intent(inout) :: canoND(:)
    !> transformation for spatial object
    type(tem_transformation_type), intent(in) :: transform
    !--------------------------------------------------------------------------!
    integer :: iCano
    !--------------------------------------------------------------------------!

    do iCano=1,size(canoND)
      call transformCanoND_single(canoND = canoND(iCano), &
        &                         transform = transform)
    end do

  end subroutine transformCanoND
  ! ****************************************************************************

  ! ****************************************************************************
  !> This routine apply transformation to canonical objects.
  subroutine transformCanoND_single(canoND, transform)
    !--------------------------------------------------------------------------!
    !> canonical geometry object type
    type( tem_canonicalND_type ), intent(inout) :: canoND
    !> transformation for spatial object
    type(tem_transformation_type), intent(in) :: transform
    !--------------------------------------------------------------------------!
    integer :: iVec
    !--------------------------------------------------------------------------!
    if(transform%active) then
      if(transform%deform%active) then
        canoND%origin = matmul( transform%deform%matrix, &
          &                     canoND%origin )
        do iVec=1,canoND%nDim
          canoND%vec(:,iVec) = matmul( transform%deform%matrix, &
            &                         canoND%vec(:,iVec) )
        enddo
      endif
      if(transform%translate%active) then
        canoND%origin = transform%translate%vec + canoND%origin
      endif
    endif

  end subroutine transformCanoND_single
  ! ****************************************************************************


! **************************************************************************** !
  !> Return the next coordinate of the canonical shape.
  !!
  pure function tem_getNextCoordOfcanonicalND( me, iElem ) result( coord )
    ! --------------------------------------------------------------------------
    !> current element within total amount of elements to process
    integer, intent(in) :: iElem
    !> current canonical type
    type( tem_canonicalND_type ), intent(in) :: me
    !> calulated real-world coordinate, which is returned
    real(kind=rk) :: coord(3)
    ! --------------------------------------------------------------------------

    ! Start from the origin, step through the segments on the two
    ! defining vectors successively

    coord = me%origin                                                        &
          ! offset in direction of vector A
      &   + real(mod(( iElem-1 ), me%segments( 1 )), kind=rk)                &
      &   / real( max(me%segments( 1 )-1,1), kind=rk)                        &
      &   * me%vec(:,1)                                                      &
          ! offset in direction of vector B
          + real(mod( (iElem-1)/me%segments( 1 ), me%segments( 2 )),         &
      &          kind=rk)                                                    &
      &   / real( max(me%segments( 2 )-1,1), kind=rk)                        &
      &   * me%vec(:,2)                                                      &
          ! offset in direction of vector C
          + real(mod( (iElem-1)/(me%segments( 1 )*me%segments( 2 )),         &
      &               me%segments( 3 )), kind=rk)                            &
      &   / real( max(me%segments( 3 )-1,1), kind=rk)                        &
      &   * me%vec(:,3)

  end function tem_getNextCoordOfcanonicalND
! **************************************************************************** !


! **************************************************************************** !
  !> Write out an array of canonical shapes in lua format
  !!
  subroutine tem_canonicalND_out_vec( me, conf )
    ! --------------------------------------------------------------------------
    !> canonicalND types to write out
    type( tem_canonicalND_type ), intent(in) :: me(:)
    !> Aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    ! --------------------------------------------------------------------------
    ! counter
    integer :: i
    ! --------------------------------------------------------------------------

    ! create a table with name canonicalND
    call aot_out_open_table( put_conf = conf, tname = 'object' )

    do i = 1, size(me)
      call tem_canonicalND_out_scal( me(i), conf )
    end do

    call aot_out_close_table( put_conf = conf )

  end subroutine tem_canonicalND_out_vec
! **************************************************************************** !


! **************************************************************************** !
  !> Write out a canonicalND shape in lua format
  !!
  subroutine tem_canonicalND_out_scal( me, conf )
    ! --------------------------------------------------------------------------
    !> canonicalND types to write out
    type( tem_canonicalND_type ), intent(in) :: me
    !> Aotus type handling the output to the file in lua format
    type(aot_out_type), intent(inout) :: conf
    ! --------------------------------------------------------------------------
    integer :: iVec
    ! --------------------------------------------------------------------------

    ! create a table with name canonicalND if not exist
    if( conf%level .eq. 0 ) then
      call aot_out_open_table( put_conf = conf, tname = 'object' )
    else
      call aot_out_open_table( put_conf = conf )
    end if
    !write origin
    call aot_out_val( put_conf = conf, vname = 'origin', val = me%origin )

    !write vec
    if(me%nDim >0) then
      call aot_out_open_table( put_conf = conf, tname = 'vec' )
      do iVec=1,me%nDim
        call aot_out_val( put_conf = conf, val = me%vec(:,iVec) )
      enddo
      call aot_out_close_table( put_conf = conf )

      !write segments
      call aot_out_open_table( put_conf = conf, tname = 'segments' )

      do iVec=1,me%nDim
        call aot_out_val( put_conf = conf, val = me%segments(iVec) )
      enddo
      call aot_out_close_table( put_conf = conf )
    endif

    if (trim(me%kind) == 'box') then
      call aot_out_val( put_conf = conf, vname = 'only_surface', &
        &               val = me%only_surface )
    end if

    !write distribution
    call aot_out_val( put_conf = conf,                                         &
      &               vname    = 'distribution',                               &
      &               val      = me%distribution )

    call aot_out_close_table( put_conf = conf )

  end subroutine tem_canonicalND_out_scal
! **************************************************************************** !


! **************************************************************************** !
  !> Check for neighboring elements to a given point.
  !!
  !! This routine can be used to find close by elements for points that just
  !! lie outside the domain and for which values can be extrapolated, though
  !! they do not strictly reside in the computational domain described by
  !! inTree.
  !! The closest neighboring element is identified globally and added to the
  !! map2global if any is found.
  !! This only works for points, and none of the points in the given list
  !! should have been found inside the domain.
  subroutine tem_cano_checkNeigh(me, inTree, countElems, map2global)
    ! --------------------------------------------------------------------------
    !> canonicalND objects to check, only points will be considered
    type(tem_canonicalND_type ),intent(in) :: me(:)
    !> Global tree
    type(treelmesh_type), intent(in) :: inTree
    !> How many elements there will be for each level in the track
    integer, intent( inout ) :: countElems( globalMaxLevels )
    !> growing array for the map2global
    type(dyn_intArray_type), intent(inout) :: map2global
    ! --------------------------------------------------------------------------
    integer :: iCano, tLevel, elemPos, dPos, maxLevel, iQQN
    integer :: closest_neigh(size(me))
    logical :: wasAdded
    integer :: coordOfId(4)
    integer :: nCanoNDs
    integer :: iError
    integer(kind=long_k) :: treeID, tOffset, neighID
    integer :: minproc(size(me))
    real(kind=rk) :: mindist(size(me)), dist
    real(kind=rk) :: minglobdist(size(me))
    real(kind=rk) :: elemcenter(3)
    ! --------------------------------------------------------------------------
    maxLevel = inTree%global%maxLevel
    mindist = huge(dist)
    closest_neigh = 0
    nCanoNDs = size(me)

    ! Create subTree by intersecting canoND objects with the inTree
    do iCano=1,nCanoNDs
      if (trim(me(iCano)%kind) == 'point') then
        treeID = tem_IdOfCoord( tem_CoordOfReal(inTree,                     &
          &                                     me(iCano)%point%coord(1:3), &
          &                                     maxLevel)                   )
        ! add a neighbor element of this point if any of the neighbor exist to
        ! interpolate to a point
        coordOfId = tem_CoordOfId( treeID )
        tOffset = tem_FirstIdAtLevel( coordOfId(4) )
        directionLoop: do iQQN = 1, qQQQ
          neighID = tem_IdOfCoord(                          &
            &         [ coordOfId(1) + qOffset( iQQN, 1 ),  &
            &           coordOfId(2) + qOffset( iQQN, 2 ),  &
            &           coordOfId(3) + qOffset( iQQN, 3 ),  &
            &           coordOfId(4) ], tOffset)
          elemPos = tem_PosOfId(neighID, inTree%treeID)
          if (elemPos > 0) then
            elemcenter = tem_BaryOfId(inTree, neighID)
            dist = sqrt(  (me(iCano)%point%coord(1) - elemcenter(1))**2 &
              &         + (me(iCano)%point%coord(2) - elemcenter(2))**2 &
              &         + (me(iCano)%point%coord(3) - elemcenter(3))**2 )
            if (dist < mindist(iCano)) then
              mindist(iCano) = dist
              closest_neigh(iCano) = elemPos
            else if (dist <= mindist(iCano)) then
              ! Make sure to pick the first element in the SFC ordering
              ! as the nearest neighbor, if there are multiple with the
              ! same distance (<= comparison to avoid == comparison for
              ! reals)
              closest_neigh(iCano) = min(elemPos, closest_neigh(iCano))
            end if
          end if
        end do directionLoop
      end if
    end do

    call MPI_Allreduce( mindist, minglobdist, nCanoNDs,             &
      &                 rk_mpi, MPI_MIN, inTree%global%comm, iError )

    minproc = inTree%global%nParts
    do iCano=1,nCanoNDs
      if ( minglobdist(iCano) < huge(dist) ) then
        if ( mindist(iCano) <= minglobdist(iCano) ) then
          minproc(iCano) = inTree%global%myPart
        end if
      end if
    end do

    call MPI_Allreduce( MPI_IN_PLACE, minproc, nCanoNDs,                 &
      &                 MPI_INTEGER, MPI_MIN, inTree%global%comm, iError )

    do iCano=1,nCanoNDs
      if (trim(me(iCano)%kind) == 'point') then
        found_elem: if (minproc(iCano) == inTree%global%myPart) then
          ! append the position in inTree to the map (note that already existing
          ! ones are omitted)
          call append( me       = map2global,           &
            &          pos      = dpos,                 &
            &          val      = closest_neigh(iCano), &
            &          wasAdded = wasAdded              )

          ! Count up if it was added
          if (wasAdded) then
            tLevel   = tem_levelOf( inTree%treeID(closest_neigh(iCano)) )
            countElems( tLevel ) = countElems( tLevel ) + 1
          end if ! wasAdded
        end if found_elem
      end if
    end do

  end subroutine tem_cano_checkNeigh
! **************************************************************************** !


! **************************************************************************** !
  !> Create subtree from the intersection of canonical shapes and inTree
  subroutine tem_cano_initSubTree( me, inTree, countElems, map2global, &
    &                              shapeInverted )
    ! --------------------------------------------------------------------------
    !> canonicalND objects on which to work
    type(tem_canonicalND_type ),intent(in) :: me(:)
    !> Global tree
    type(treelmesh_type), intent(in) :: inTree
    !> How many elements there will be for each level in the track
    integer, intent( inout ) :: countElems( globalMaxLevels )
    !> growing array for the map2global
    type(dyn_intArray_type), intent(inout) :: map2global
    !> If true then elements not intersected are added to subTree
    logical, intent(in) :: shapeInverted
    ! --------------------------------------------------------------------------
    integer :: iCano, tLevel, elemPos, iElem, dPos, maxLevel
    integer(kind=long_k) :: treeID
    logical :: wasAdded, intersects, addToSubTree
    type(tem_cube_type) :: cube
    ! --------------------------------------------------------------------------

    write(logUnit(4),"(A)") 'Initializing canonical shapes'

    maxLevel = inTree%global%maxLevel

    ! Create subTree by intersecting canoND objects with the inTree
    do iCano = 1, size(me)
      ! for a point there is no need to check against all elements
      ! instead just convert the point into treeID and check for its
      ! exitence in given tree
      if (trim(me(iCano)%kind) == 'point') then
        treeID = tem_IdOfCoord( tem_CoordOfReal(inTree,                     &
          &                                     me(iCano)%point%coord(1:3), &
          &                                     maxLevel) )
        ! get position of the treeID in inTree
        elemPos = tem_PosOfId( treeID, inTree%treeID)
        found_elem: if (elemPos > 0) then
          ! append the position in inTree to the map (note that already existing
          ! ones are omitted)
          call append( me       = map2global, &
            &          pos      = dpos,       &
            &          val      = elemPos,    &
            &          wasAdded = wasAdded    )

          ! Count up if it was added
          if (wasAdded) then
            tLevel   = tem_levelOf( treeID )
            countElems( tLevel ) = countElems( tLevel ) + 1
          end if ! wasAdded
        end if found_elem
      else
        do iElem = 1, inTree%nElems
          treeID = inTree%treeID(iElem)
          call tem_convertTreeIDtoCube(cube, inTree, treeID)
          intersects = .false.
          select case ( trim(me(iCano)%kind) )
          case ('line')
            intersects = tem_lineCubeOverlap(me(iCano)%line, cube)
          case ('plane')
            intersects = tem_planeCubeOverlap(me(iCano)%plane, cube)
          case ('box')
            intersects = tem_boxCubeOverlap(me(iCano)%box, cube)
          end select

          addToSubTree = .false.
          if (.not. shapeInverted .and. intersects) then
            ! Shape intersects with current element and not inverted
            addToSubTree = .true.
          else if (shapeInverted .and. .not. intersects) then
            ! shape not intersected and is inverted shape so add this to subTree
            addToSubTree = .true.
          end if

          if (addToSubTree) then
            ! append iElem in inTree to the map (note that already existing
            ! ones are omitted)
            call append( me       = map2global, &
              &          pos      = dpos,       &
              &          val      = iElem ,     &
              &          wasAdded = wasAdded    )

            ! Count up if it was added
            if( wasAdded ) then
              tLevel   = tem_levelOf( treeID )
              countElems( tLevel ) = countElems( tLevel ) + 1
            end if ! wasAdded
          end if !intersects
        end do !iElem
      end if
    end do !iCano

  end subroutine tem_cano_initSubTree
! **************************************************************************** !

! **************************************************************************** !
  !> Generate points using segments on canoND and add those points
  !! to a growing array of points if a point is found in subTree
  subroutine tem_cano_storePntsInSubTree( me, inTree, map2global, countPoints, &
    &                                     grwPnts )
    ! --------------------------------------------------------------------------
    !> canonicalND objects on which to work
    type(tem_canonicalND_type ),intent(in) :: me(:)
    !> Global tree
    type(treelmesh_type), intent(in) :: inTree
    !> How many points there will be
    integer, intent(inout) :: countPoints
    !> growing array for the map2global
    type(dyn_intArray_type), intent(in) :: map2global
    !> growing array to store tracking points
    type(tem_grwPoints_type), intent(inout) :: grwPnts
    ! --------------------------------------------------------------------------
    integer :: nElems, nPoints, maxLevel, elemPos, neighPos, coordOfId(4)
    integer :: iCano, iPnt, iQQN
    real(kind=rk) :: coord(3), offset_a, offset_b, offset_c
    real(kind=rk) :: unit_vec_a(3),unit_vec_b(3), unit_vec_c(3)
    integer(kind=long_k) :: treeID, tOffset, neighID
    integer(kind=long_k), allocatable :: subTreeID(:)
    ! --------------------------------------------------------------------------
    maxLevel = inTree%global%maxLevel
    ! Append the physical points to the growing array of points
    nElems = map2global%nVals
    allocate(subTreeID(nElems))
    subTreeID = inTree%treeID(map2global%val(map2global%sorted(1:nElems)))
    do iCano = 1, size(me)
      ! total number of elements
      nPoints = me(iCano)%segments(1) &
        &     * me(iCano)%segments(2) &
        &     * me(iCano)%segments(3)

      unit_vec_a =  me(iCano)%vec(:,1)                      &
        &        / real( max(me(iCano)%segments(1)-1,1), rk )
      unit_vec_b =  me(iCano)%vec(:,2)                      &
        &        / real( max(me(iCano)%segments(2)-1,1), rk )
      unit_vec_c =  me(iCano)%vec(:,3)                      &
        &        / real( max(me(iCano)%segments(3)-1,1), rk )

      ! Generate points and append only the points available in tree
      do iPnt = 1, nPoints
        offset_a = real(mod((iPnt-1),me(iCano)%segments(1)), kind=rk)
        offset_b = real(mod((iPnt-1)/me(iCano)%segments(1), &
          &                 me(iCano)%segments(2)), kind=rk)
        offset_c = real(mod((iPnt-1)/(me(iCano)%segments(1)  &
          &                          *me(iCano)%segments(2)), &
          &                 me(iCano)%segments(3)), kind=rk)

        coord = me(iCano)%origin + offset_a * unit_vec_a &
          &                      + offset_b * unit_vec_b &
          &                      + offset_c * unit_vec_c

        ! Get the treeID on the highest level
        treeID = tem_IdOfCoord( tem_CoordOfReal(inTree, coord(1:3), maxLevel) )
        ! get position of the treeID in subTree
        elemPos = tem_PosOfId( treeID, subTreeID )

        neighPos = 0
        if (elemPos <= 0) then
          ! Point must be outside fluid domain, its neighbor must be in subTreeID
          ! since it was added in tem_cano_initSubTree
          coordOfId = tem_CoordOfId( treeID )
          tOffset = tem_FirstIdAtLevel( coordOfId(4) )
          directionLoop: do iQQN = 1, qQQQ
            neighID = tem_IdOfCoord(                          &
              &         [ coordOfId(1) + qOffset( iQQN, 1 ),  &
              &           coordOfId(2) + qOffset( iQQN, 2 ),  &
              &           coordOfId(3) + qOffset( iQQN, 3 ),  &
              &           coordOfId(4) ], tOffset)
            neighPos = tem_PosOfId( neighID, subTreeID)
            if (neighPos > 0) exit directionLoop
          end do directionLoop
        end if

        if( elempos > 0 .or. neighPos > 0 ) then
          ! append the physical points to the growing array of points
          call append( me  = grwpnts, &
            &          val = coord    )

          countPoints = countPoints + 1
        end if
      end do !iPoint
    end do !iCano
    deallocate(subTreeID)

  end subroutine tem_cano_storePntsInSubTree

end module tem_canonicalND_module
! **************************************************************************** !
