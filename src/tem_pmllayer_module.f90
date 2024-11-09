! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> This module gathers information for the PML damping functions.
!!
module tem_pmlLayer_module

  ! include treelm modules
  use env_module,         only: rk
  use tem_param_module,   only: PI
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit

  ! include aotus modules
  use aotus_module,     only: aoterr_NonExistent, flu_state
  use aot_table_module, only: aot_table_open, aot_table_close,                 &
   &                          aot_table_length, aot_get_val

  implicit none

  public :: tem_load_pmlLayer, tem_pmlLayer_type

  !> This type contains datas to define PML layer.
  type tem_pmlLayer_type
    !> Plane origin
    real(kind=rk) :: plane_origin(3)
    !> Plane normal
    real(kind=rk) :: plane_normal(3)
    !> Damp factor for the Layer
    real(kind=rk) :: dampFactor
    !> Damping exponent for the layer
    integer :: dampExponent
  end type tem_pmlLayer_type


contains

  !> Load definition of the PML damping term.
  subroutine tem_load_pmlLayer(conf, thandle, me)
    ! ---------------------------------------------------------------------------
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: thandle
    !> Global pmlLayer data type
    type(tem_pmlLayer_type), intent(out) :: me
    ! ---------------------------------------------------------------------------
    integer :: cent_handle
    integer :: i
    integer :: iError
    ! ---------------------------------------------------------------------------
    ! Plane_origin
    call aot_table_open( L       = conf,                                       &
      &                  parent  = thandle,                                    &
      &                  thandle = cent_handle,                                &
      &                  key     = 'plane_origin' )
    if (aot_table_length(L=conf, thandle=cent_handle) == 3) then
      do i=1,3
        call aot_get_val( L       = conf,                                      &
          &               thandle = cent_handle,                               &
          &               pos     = i,                                         &
          &               val     = me%plane_origin(i),                        &
          &               ErrCode = iError )
      end do
    else
      write(*,*) 'ERROR while reading Plane origin of PML Layer,'
      write(*,*) 'should have 3 entries for each coordinate!'
      STOP
    end if
    call aot_table_close(conf, cent_handle)
    write(logUnit(1),*) ' * plane_origin =', me%plane_origin



    ! Plane_normal
    call aot_table_open( L       = conf,                                       &
      &                  parent  = thandle,                                    &
      &                  thandle = cent_handle,                                &
      &                  key     = 'plane_normal' )
    if (aot_table_length(L=conf, thandle=cent_handle) == 3) then
      do i=1,3
        call aot_get_val( L       = conf,                                      &
          &               thandle = cent_handle,                               &
          &               pos     = i,                                         &
          &               val     = me%plane_normal(i),                        &
          &               ErrCode = iError )
      end do
    else
      write(*,*) 'ERROR while reading Plane origin of PML Layer,'
      write(*,*) 'should have 3 entries for each coordinate!'
      STOP
    end if
    call aot_table_close(conf, cent_handle)
    write(logUnit(1),*) ' * plane_normal =', me%plane_normal

    !damp_factor
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'damp_factor',                                 &
      &               val     = me%dampFactor,                                 &
      &               ErrCode = iError )
    write(logUnit(1),*) ' * Damping_factor =', me%dampFactor

    !damp_exponent
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'damp_exponent',                               &
      &               val     = me%dampExponent,                               &
      &               ErrCode = iError )
    if (iError .ne. 0) then
      me%dampExponent = 1
    endif

    write(logUnit(1),*) ' * Damping_exponent =', me%dampExponent

  end subroutine tem_load_pmlLayer


  !> Calculate damping functions times normal and derivatives
  !! times normal for the PML evaluation.
  function tem_evaluate_pml(me, nComp, coord, n)  &
    &                           result(res)
    !> Spacetime function to evaluate
    type(tem_pmlLayer_type) :: me
    !> Number of arrays to return
    integer, intent(in) :: n
    !> Number of entrys in each array
    integer, intent(in) :: ncomp
    !> barycentric Ids of an elements.
    !! 1st index goes over number of elements and
    !! 2nd index goes over x,y,z coordinates
    real(kind=rk), intent( in ) :: coord(n,3)
    !> return value
    real(kind=rk) :: res(n,ncomp)
    ! ---------------------------------------------------------------------------
    integer :: i
    real(kind=rk) :: origin(3), normal(3), vec(3), projection
    ! ---------------------------------------------------------------------------
    origin(:) = me%plane_origin
    normal(:) = me%plane_normal

    select case(ncomp)
    ! 2D case
    case(4)
      do i = 1,n
        vec (:) = coord(i,:) - origin(:)
        projection = vec(1)*normal(1)+vec(2)*normal(2)
        if(projection > 0.0_rk) then
          res(i,1:2) = me%dampFactor*abs(normal(1:2))*((coord(i,1:2) - origin(1:2))**me%dampExponent)/sqrt(sum(normal**2))
          res(i,3:4) = me%dampExponent*me%dampFactor*abs(normal(1:2))*((coord(i,1:2)&
                     & - origin(1:2))**(me%dampExponent-1))/sqrt(sum(normal**2))
        else
          res(i,1:4) = 0.0_rk
        end if
      enddo
    ! 3D case
    case(6)
      do i = 1,n
        vec (:) = coord(i,:) - origin(:)
        projection = vec(1)*normal(1)+vec(2)*normal(2)+vec(3)*normal(3)
        if(projection > 0.0_rk) then
          res(i,1:3) = me%dampFactor*abs(normal(1:3))*(((coord(i,1:3) - origin(1:3))/sqrt(sum(normal**2)))**me%dampExponent)
          res(i,4:6) = me%dampExponent*me%dampFactor*abs(normal(1:3))*(((coord(i,1:3)&
                     & - origin(1:3))/sqrt(sum(normal**2)))**(me%dampExponent-1)) / sqrt(sum(normal**2))
        else
          res(i,1:6) = 0.0_rk
        end if
      enddo
    case default
      write(*,*) 'ERROR in tem_evaluate_pml: Unknown number of components, stopping ...'
      stop
    end select

  end function tem_evaluate_pml


end module tem_pmlLayer_module
