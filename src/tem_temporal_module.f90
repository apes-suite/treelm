! Copyright (c) 2011-2013,2016,2019,2021 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2011-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
!> author: Kannan Masilamani
!! author: Simon Zimny
!! This module gathers different temporal effect treatments.
!! like: constant, lua function and predefined Fortran function. \n
!! Currently supports:
!!
!! * const: data is a constant value
!! * lua_fun: data is a lua function
!! * linear: data is ramped linearly
!! * smooth: data is ramped with a smooth function
!! * cos: oscillation as a cosine function
!! * from_file: data is read from file (supports periodicity)
!!
!! @note Boundaries can be defined with superposed conditions in
!! space and time. Therefore the function defined here acts as a
!! factor applied to the overall boundary value.
!!
!!
module tem_temporal_module

  ! include treelm modules
  use env_module,            only: rk, LabelLen, PathLen, eps
  use tem_aux_module,        only: tem_open, tem_abort
  use tem_param_module,      only: PI
  use tem_time_module,       only: tem_time_type
  use tem_grow_array_module, only: grw_real2darray_type, &
    &                              init, append, expand, destroy
  use tem_logging_module,    only: logUnit

  ! include aotus modules
  use aotus_module,     only: aoterr_Fatal, aoterr_NonExistent, flu_State, &
    &                         aoterr_WrongType, aot_top_get_val, aot_get_val
  use aot_table_module, only: aot_table_open, aot_table_close
  use aot_fun_module,   only: aot_fun_type, aot_fun_open, aot_fun_close, &
    &                         aot_fun_put, aot_fun_do
  use aot_references_module, only: aot_reference_for, aot_reference_to_top

  implicit none

  private

  public :: tem_temporal_type
  public :: tem_load_temporal, tem_temporal_for

  !> contains information for predefined temporal functions
  type tem_linear_type
    !> start time
    real(kind=rk) :: from_time
    !> end time
    real(kind=rk) :: to_time
    !> minimum factor at start time
    real(kind=rk) :: min_factor
    !> maximum factor at end time
    real(kind=rk) :: max_factor
  end type tem_linear_type

  !> contains information for loading inlet velocities from a datafile
  !! The data has to be stored as tuples (time,velocity) columnwise.
  !!         t1 v1
  !!         t2 v2
  !!         ...
  !!         tn vn
  !! and v1 .eq. vn has to be fullfilled.
  !! The data has to be provided in the format '(e15.8)'
  type tem_from_file_temporal_type
    !> filename of the data
    character(len=PathLen) :: datafile
    !> interpolate linearly between the data
    character(len=LabelLen) :: intp
    !> growing array of tuples (time, velocity)
    type( grw_real2darray_type ) :: signal
    !> ramping active?
    logical :: ramp = .false.
    !> ramping value at the end of rampT
    real(kind=rk) :: rampVal
    !> ramping time
    real(kind=rk) :: rampT
    !> factor to multiply data with
    real(kind=rk) :: fac
    !> is the data periodic?
    logical :: periodic
  end type tem_from_file_temporal_type

  !> defines different temporal types like const, lua func or predefined func
  type tem_temporal_type
    !> temporal kind
    character(len=LabelLen) :: kind
    !> constant value
    real(kind=rk) :: const

    !> Reference to the Lua function if the temporal function is defined
    !! as a Lua function.
    integer :: lua_fun_ref = 0

    !> Handle to the Lua script for the Lua function
    type(flu_state) :: conf

    !> contains information for predefined functions
    type( tem_linear_type ) :: linear
    !> contains information for reading the data from file
    type( tem_from_file_temporal_type) :: from_file

    !> frequency of oscillation
    !! Load in routine: load_temporal_cos
    real(kind=rk) :: freq
    !> initial phase
    real(kind=rk) :: phi
    !> offset
    real(kind=rk) :: offset
  end type tem_temporal_type

  contains

! ****************************************************************************** !
  !> This subroutine load temporal table defined for a boundary.\n
  !!
  !! If temporal is defined as lua function then set kind = temporal_lua
  !! else if temporal block is defined then load temporal table for predefined
  !! Fortran function variables and set kind = temporal_\a function_name
  !! else temporal is a constant value and set kind = temporal_const.\n
  !! \n
  !! Valid definitions:
  !! \li Constant
  !! ~~~~~~~~~~~~~~~~~~~~~{.lua}
  !! temporal = 1.0
  !! ~~~~~~~~~~~~~~~~~~~~~
  !! \li lua_function
  !! ~~~~~~~~~~~~~~~~~~~~~
  !! temporal = 'linear'
  !! ~~~~~~~~~~~~~~~~~~~~~
  !! Example: \a linear lua function
  !! \verbatim
  !! function linear(iTime)
  !!  local to_time = 1000
  !!  if iTime < to_time then
  !!    return iTime/to_time
  !!  else
  !!    return 1.0
  !!  end
  !!end
  !! \endverbatim
  !! \li Predefined Fortran function
  !! ~~~~~~~~~~~~~~~~~~~~~{.lua}
  !! temporal = {predefined="Fortranfun_name",
  !!              min_factor = 0.0, max_factor = 1.0,
  !!              from_time = 0.0, to_time = 1000.0}
  !! ~~~~~~~~~~~~~~~~~~~~~
  !! \li Data from a file (periodic data supported)
  !! ~~~~~~~~~~~~~~~~~~~~~{.lua}
  !! temporal = {predefined="datafile",
  !!              filename='data.dat',      -- path/name of the datafile
  !!              intp='linear',            -- interpolation between the time tics ('linear','none')
  !!              periodic= true}           -- is the data periodic?
  !! ~~~~~~~~~~~~~~~~~~~~~
  !! \n
  subroutine tem_load_temporal( me, conf, parent, key )
    ! ---------------------------------------------------------------------------
    !> boundary temporal type
    type(tem_temporal_type), intent(out) :: me
    !> lua state
    type(flu_State) :: conf
    !> parent handle contains temporal table
    integer, intent(in) :: parent
    !> state variable key string defined in lua
    character(len=*), intent(in), optional :: key
    ! ---------------------------------------------------------------------------
    type(aot_fun_type) :: fun
    integer :: thandle
    integer :: iError
    character(len=labelLen) :: local_key
    ! ---------------------------------------------------------------------------

    if(present(key)) then
      local_key = trim(key)
    else
      local_key = 'temporal'
    endif

    ! First test for a lua function
    call aot_fun_open( L      = conf,                                          &
      &                parent = parent,                                        &
      &                fun    = fun,                                           &
      &                key    = trim( local_key ))

    me%conf = conf

    if (fun%handle /= 0) then
      ! This temporal modifier is defined as a Lua function
      me%kind = 'lua_fun'
      me%lua_fun_ref = aot_reference_for(conf)
      write(logUnit(1),*)'    Defined lua temporal function'
      call aot_fun_close( L = conf, fun = fun )
    else
      ! It is not defined as a function, try to interpret it as a table
      call aot_table_open( L       = conf,                                     &
        &                  thandle = thandle,                                  &
        &                  parent  = parent,                                   &
        &                  key     = trim( local_key ))
      if (thandle /= 0) then
        call aot_get_val( L       = conf,                                      &
          &               thandle = thandle,                                   &
          &               key     = 'predefined',                              &
          &               val     = me%kind,                                   &
          &               default = 'unknown',                                 &
          &               ErrCode = iError)

        write(logUnit(1),*)'   A predefined temporal function is chosen: '&
          &                //trim(me%kind)
        select case(trim(me%kind))
        case('linear', 'smooth')
          ! Load the standard parameters necessary to describe the
          ! temporal behavior of the BC.
          call load_temporal_linear( me      = me%linear, &
            &                         conf    = conf,          &
            &                         thandle = thandle )
        case('datafile')
          ! Load the filename for the datafile
          call load_temporal_from_file( me      = me%from_file,               &
            &                            conf    = conf,                       &
            &                            thandle = thandle )
        case('cos')
          call load_temporal_cos( freq    = me%freq,&
            &                      phi     = me%phi, &
            &                      offset  = me%offset, &
            &                      conf    = conf,   &
            &                      thandle = thandle )
        case default
          write(logUnit(1),*)'ERROR in definition of the temporal '//         &
            &            ' boundary Conditions:'
          write(logUnit(1),*)'Selected an unknown temporal boundary'
          write(logUnit(1),*)trim(me%kind)
          call tem_abort()
        end select
      else
        ! As the entry for the variable is neither a function nor a table, try
        ! to interpret it as a constant variable
        call aot_get_val( L       = conf,                                      &
          &               thandle = parent,                                    &
          &               key     = trim(local_key),                           &
          &               val     = me%const,                                  &
          &               ErrCode = iError )
        if (btest(iError, aoterr_WrongType)) then
          write(logUnit(1),*)'FATAL Error occured in definition of the '//     &
            &                'temporal boundary conditions'
          write(logUnit(1),*)'while retrieving temporal constant:'
          write(logUnit(1),*)'Variable has wrong type (should be a '//         &
            &                'real number)!'
          write(logUnit(1),*)'STOPPING'
          call tem_abort()
        end if
        if (btest(iError, aoterr_NonExistent)) then
          ! "temporal" variable not specified
          me%kind = 'none'
          me%const = 1.0_rk
        else
          me%kind = 'const'
        end if
      end if
    end if

  end subroutine tem_load_temporal
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine load standard temporal function variables from LUA file.
  !!
  !! Default values:
  !! \li min_factor = 0.0
  !! \li max_factor = 1.0
  !! \li from_time = 0
  !! \li to_time = tmax/2
  !!
  !! Valid definition:
  !! \li linear function
  !! ~~~~~~~~~~~~~~~~~~~~~{.lua}
  !! temporal = {predefined="linear", min_factor = 0.0, max_factor = 1.0,
  !! from_time = 0, to_time = 1000}
  !! ~~~~~~~~~~~~~~~~~~~~~
  !! Example: Transient inlet velocity which starts from 0 to 1000 with maximum
  !! value 0.08 is shown below for linear and smooth with definition and image.
  !! ~~~~~~~~~~~~~~~~~~~~~{.lua}
  !! boundary_condition =
  !!         {
  !!           { label = 'inlet',
  !!             kind = 'inlet_ubb',
  !!             velocityX = { spatial = 1.0
  !!                ,temporal = {predefined="linear", min_factor = 0.0,
  !!                      max_factor = 1.0, from_time = 0, to_time = 1000}}
  !! --                ,temporal = {predefined="smooth", min_factor = 0.0,
  !! --                      max_factor = 1.0, from_time = 0, to_time = 1000}}
  !!             velocityY = 0.0
  !!             velocityZ = 0.0
  !!           }
  !!         }
  !! ~~~~~~~~~~~~~~~~~~~~~
  !! \image html transient.png
  !!
  subroutine load_temporal_linear( me, conf, thandle )
    ! ---------------------------------------------------------------------------
    !> temporal predefined fun type
    type(tem_linear_type),intent(inout) :: me
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: thandle
    ! ---------------------------------------------------------------------------
    ! local variables
    integer :: iError
    ! ---------------------------------------------------------------------------

    !min_factor
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'min_factor',                                  &
      &               val     = me%min_factor,                                 &
      &               ErrCode = iError,                                        &
      &               default = 0.0_rk )
    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*)'FATAL Error occured, while retrieving min_factor :'
      if (btest(iError, aoterr_NonExistent))                                   &
        & write(logUnit(1),*)'Variable not existent!'
      if (btest(iError, aoterr_WrongType))                                     &
        & write(logUnit(1),*)'Variable has wrong type!'
      write(logUnit(1),*)'STOPPING'
      call tem_abort()
    end if

    !max_factor
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'max_factor',                                  &
      &               val     = me%max_factor,                                 &
      &               ErrCode = iError,                                        &
      &               default = 1.0_rk )
    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*)'FATAL Error occured, while retrieving max_factor :'
      if (btest(iError, aoterr_NonExistent))                                   &
        & write(logUnit(1),*)'Variable not existent!'
      if (btest(iError, aoterr_WrongType))                                     &
        & write(logUnit(1),*)'Variable has wrong type!'
      write(logUnit(1),*)'STOPPING'
      call tem_abort()
    end if

    !from_time
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'from_time',                                   &
      &               val     = me%from_time,                                  &
      &               ErrCode = iError,                                        &
      &               default = 0.0_rk )
    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*)'FATAL Error occured, while retrieving from_time :'
      if (btest(iError, aoterr_NonExistent))                                   &
        & write(logUnit(1),*)'Variable not existent!'
      if (btest(iError, aoterr_WrongType))                                     &
        & write(logUnit(1),*)'Variable has wrong type!'
      write(logUnit(1),*)'STOPPING'
      call tem_abort()
    end if

    ! to_time
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'to_time',                                     &
      &               val     = me%to_time,                                    &
      &               ErrCode = iError,                                        &
      &               default = 1.0_rk )
    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*)'FATAL Error occured, while retrieving to_time :'
      if (btest(iError, aoterr_NonExistent))                                   &
        & write(logUnit(1),*)'Variable not existent!'
      if (btest(iError, aoterr_WrongType))                                     &
        & write(logUnit(1),*)'Variable has wrong type!'
      write(logUnit(1),*)'STOPPING'
      call tem_abort()
    end if

    write(logUnit(3), "(A, F10.5)") "     min factor: ", me%min_factor
    write(logUnit(3), "(A, F10.5)") "     max factor: ", me%max_factor
    write(logUnit(3), "(A, E10.5)") "      from time: ", me%from_time
    write(logUnit(3), "(A, E10.5)") "        to time: ", me%to_time

  end subroutine load_temporal_linear
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine loads the information needed to read data from a file.
  !!
  !! Example: Read time periodic velocity signals from file 'data.dat' as inlet
  !!          condition without interpolation between the velocity signals.
  !!
  !! ~~~~~~~~~~~~~~~~~~~~~{.lua}
  !! boundary_condition =
  !!         {
  !!           { label = 'inlet',
  !!             kind = 'inlet_ubb',
  !!             velocityX = { kind = 'combined',
  !!                           temporal = { predefined = 'datafile',
  !!                                         filename   = 'data.dat',
  !!                                         intp       = 'none',
  !!                                         periodic   = true
  !!                                       }
  !!                         },
  !!             velocityY = 0.0
  !!             velocityZ = 0.0
  !!           }
  !!         }
  !! ~~~~~~~~~~~~~~~~~~~~~
  !!
  !! The options for interpolating the values are
  !! \li \a 'none': No interpolation is used. The signals will be chosen by:
  !!                [t1,t2[ = v1, [t2,t3[ = v2, ...
  !! \li \a 'linear': A linear interpolation is done between neighboring
  !!                  values according to the following equation:
  !!                  \f[
  !!                         \frac{v_{i+1}-v_{i}}{t_{i+1}-t_{i}} \cdot
  !!                              (t_{sim} \% t_{period} - t_{i}) + v_{i}
  !!                  \f]
  !!
  subroutine load_temporal_from_file( me, conf, thandle )
    ! ---------------------------------------------------------------------------
    !> temporal predefined from_file type
    type(tem_from_file_temporal_type),intent(inout) :: me
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: thandle
    ! ---------------------------------------------------------------------------
    ! local variables
    integer :: iError
    ! aotus handle for ramping table
    integer :: rampHandle
    ! ---------------------------------------------------------------------------

    ! read the filename of the datafile
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'filename',                                    &
      &               val     = me%datafile,                                   &
      &               ErrCode = iError )
    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*)'FATAL Error occured, while retrieving datafile:'
      if (btest(iError, aoterr_NonExistent))                                   &
        & write(logUnit(1),*)'Variable not existent!'
      if (btest(iError, aoterr_WrongType))                                     &
        & write(logUnit(1),*)'Variable has wrong type!'
      write(logUnit(1),*)'STOPPING'
      call tem_abort()
    end if

    ! read the interpolation kind for the data
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'intp',                                        &
      &               val     = me%intp,                                       &
      &               ErrCode = iError,                                        &
      &               default = 'none' )
    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*)'FATAL Error occured, while retrieving interpolation:'
      if (btest(iError, aoterr_NonExistent))                                   &
        & write(logUnit(1),*)'Variable not existent!'
      if (btest(iError, aoterr_WrongType))                                     &
        & write(logUnit(1),*)'Variable has wrong type!'
      write(logUnit(1),*)'STOPPING'
      call tem_abort()
    end if

    ! open the ramping table
    call aot_table_open( L       = conf,                                       &
      &                  thandle = rampHandle,                                 &
      &                  parent  = thandle,                                    &
      &                  key     = 'ramping')
    if (rampHandle /= 0) then
      me%ramp = .true.
      ! read the ramping time
      call aot_get_val( L       = conf,                                        &
        &               thandle = rampHandle,                                  &
        &               key     = 'rampVal',                                   &
        &               val     = me%rampVal,                                  &
        &               ErrCode = iError)
      if (btest(iError, aoterr_Fatal)) then
        write(logUnit(1),*)'FATAL Error occured, while retrieving ramping '//  &
          &                'value:'
        if (btest(iError, aoterr_NonExistent))                                 &
          & write(logUnit(1),*)'Variable not existent!'
        if (btest(iError, aoterr_WrongType))                                   &
          & write(logUnit(1),*)'Variable has wrong type!'
        write(logUnit(1),*)'STOPPING'
        call tem_abort()
      end if

      ! read the ramping time
      call aot_get_val( L       = conf,                                        &
        &               thandle = rampHandle,                                  &
        &               key     = 'rampT',                                     &
        &               val     = me%rampT,                                    &
        &               ErrCode = iError)
      if (btest(iError, aoterr_Fatal)) then
        write(logUnit(1),*)'FATAL Error occured, while retrieving ramping time:'
        if (btest(iError, aoterr_NonExistent))                                 &
          & write(logUnit(1),*)'Variable not existent!'
        if (btest(iError, aoterr_WrongType))                                   &
          & write(logUnit(1),*)'Variable has wrong type!'
        write(logUnit(1),*)'STOPPING'
        call tem_abort()
      end if
    end if

    ! read the factor
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'fac',                                         &
      &               val     = me%fac,                                        &
      &               ErrCode = iError,                                        &
      &               default = 1.0_rk )
    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*)'FATAL Error occured, while retrieving factor: '
      if (btest(iError, aoterr_NonExistent))                                   &
        & write(logUnit(1),*)'Variable not existent!'
      if (btest(iError, aoterr_WrongType))                                     &
        & write(logUnit(1),*)'Variable has wrong type!'
      write(logUnit(1),*)'STOPPING'
      call tem_abort()
    end if

    ! read in wether the data should be treated periodic
    call aot_get_val( L       = conf,                                          &
      &               thandle = thandle,                                       &
      &               key     = 'periodic',                                    &
      &               val     = me%periodic,                                   &
      &               ErrCode = iError,                                        &
      &               default = .false. )
    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*)'FATAL Error occured, while retrieving periodic:'
      if (btest(iError, aoterr_NonExistent))                                   &
        & write(logUnit(1),*)'Variable not existent!'
      if (btest(iError, aoterr_WrongType))                                     &
        & write(logUnit(1),*)'Variable has wrong type!'
      write(logUnit(1),*)'STOPPING'
      call tem_abort()
    end if

    write(logUnit(3),*)"Using data from file: "//trim(me%datafile)
    write(logUnit(3),*)"  interpolation:      "//trim(me%intp)
    if (rampHandle /= 0) then
      write(logUnit(3),*)"  ramping time:       ",me%rampT
      write(logUnit(3),*)"  ramping value:      ",me%rampVal
    end if
    write(logUnit(3),*)"  periodic:           ",me%periodic

    ! load the data from file
    call load_datafile( me )

  end subroutine load_temporal_from_file
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine load_temporal_cos( freq, phi, offset, conf, thandle )
    ! ---------------------------------------------------------------------------
    !> temporal predefined fun type
    real(kind=rk), intent(inout) :: freq, phi, offset
    !> lua state type
    type(flu_State) :: conf
    !> aotus parent handle
    integer, intent(in) :: thandle
    ! ---------------------------------------------------------------------------
    ! local variables
    integer :: iError
    ! ---------------------------------------------------------------------------

    ! load angular velocity
    call aot_get_val( L       = conf,      &
      &               thandle = thandle,   &
      &               key     = 'frequency',       &
      &               val     = freq,         &
      &               ErrCode = iError,                                        &
      &               default = 1.0_rk)
    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*)'FATAL Error occured, while retrieving min_factor :'
      if (btest(iError, aoterr_NonExistent))                                   &
        & write(logUnit(1),*)'Variable not existent!'
      if (btest(iError, aoterr_WrongType))                                     &
        & write(logUnit(1),*)'Variable has wrong type!'
      write(logUnit(1),*)'STOPPING'
      call tem_abort()
    end if

    ! load initial phase
    call aot_get_val( L       = conf,    &
      &               thandle = thandle, &
      &               key     = 'phase', &
      &               val     = phi,     &
      &               ErrCode = iError,                                        &
      &               default = 0.0_rk   )
    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*)'FATAL Error occured, while retrieving max_factor :'
      if (btest(iError, aoterr_NonExistent))                                   &
        & write(logUnit(1),*)'Variable not existent!'
      if (btest(iError, aoterr_WrongType))                                     &
        & write(logUnit(1),*)'Variable has wrong type!'
      write(logUnit(1),*)'STOPPING'
      call tem_abort()
    end if

    ! load initial phase
    call aot_get_val( L       = conf,    &
      &               thandle = thandle, &
      &               key     = 'offset', &
      &               val     = offset,   &
      &               ErrCode = iError,                                        &
      &               default = 0.0_rk   )
    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*)'FATAL Error occured, while retrieving max_factor :'
      if (btest(iError, aoterr_NonExistent))                                   &
        & write(logUnit(1),*)'Variable not existent!'
      if (btest(iError, aoterr_WrongType))                                     &
        & write(logUnit(1),*)'Variable has wrong type!'
      write(logUnit(1),*)'STOPPING'
      call tem_abort()
    end if

    write(logUnit(3), "(A, F10.5)") " frequency: ", freq
    write(logUnit(3), "(A, F10.5)") "     phase: ", phi
    write(logUnit(3), "(A, F10.5)") "    offset: ", offset

  end subroutine load_temporal_cos
! **************************************************************************** !

! **************************************************************************** !
  !> This subroutine reads the data from disc and stores it in the
  !! tem_from_file_temporal_type.
  !!
  subroutine load_datafile( me )
    ! --------------------------------------------------------------------------
    !> temporal predefined from_file type
    type(tem_from_file_temporal_type),intent(inout) :: me
    ! --------------------------------------------------------------------------
    ! unit of the datafile
    integer :: datafile_unit
    ! status of IO
    integer :: stat = 0
    ! tmp time value
    real(kind=rk) :: tmp_time
    ! tmp vel value
    real(kind=rk) :: tmp_vel
    ! --------------------------------------------------------------------------

    ! initialize the data array
    call init( me = me%signal, width = 2 )

    ! open the datafile
    call tem_open( newunit = datafile_unit,     &
      &            file    = trim(me%datafile), &
      &            action  = 'read',            &
      &            status  = 'old',             &
      &            form    = 'formatted',       &
      &            access  = 'sequential'       )

    ! as long as there is no io error (end of file is not reached)
    do
      ! ... read the data
      read( unit = datafile_unit, fmt='(2e15.8)', iostat=stat ) tmp_time, &
        &                                                       tmp_vel
      ! ... check wether the file has reached its end to prevent storing
      !     the last item twice
      if ( stat .ne. 0) then
        exit
      end if
      ! ... and append them to the growing array of tuples
      call append( me = me%signal, val = (/ tmp_time, tmp_vel/) )
    end do

    ! multiply the data with the provided factor
    me%signal%val(2,:)=me%fac*me%signal%val(2,:)

    if( stat .gt. 0 ) then
      write(logUnit(1),*)'Error when reading the inlet data from file!'
      write(logUnit(1),*)'STOPPING'
      call tem_abort()
    end if

    ! close the file
    close( unit = datafile_unit )

  end subroutine load_datafile
! **************************************************************************** !


! **************************************************************************** !
  !> This function returns value of linear function which is defined by
  !! from_time, to_time, min_factor and max_factor
  !!
  pure function temporal_linear_for( me, t ) result(res)
    ! --------------------------------------------------------------------------
    !> temporal linear type
    type(tem_linear_type), intent(in) :: me
    !> current time
    real(kind=rk), intent(in) :: t
    !> return value of a function
    real(kind=rk) :: res
    ! --------------------------------------------------------------------------

    if (t <= me%from_time) then
      res = me%min_factor
    else if (t >= me%to_time) then
      res = me%max_factor
    else
      res = me%min_factor + (me%max_factor - me%min_factor)                    &
        & * (  (t - me%from_time) / (me%to_time - me%from_time) )
    end if
  end function temporal_linear_for
! **************************************************************************** !

! **************************************************************************** !
  !> This function returns value of smooth sin function which is defined by
  !! from_time, to_time, min_factor and max_factor
  !!
  pure function temporal_smooth_for( me, t ) result(res)
    ! --------------------------------------------------------------------------
    !> temporal smooth type
    type(tem_linear_type), intent(in) :: me
    !> current time
    real(kind=rk), intent(in) :: t
    !> return value of a function
    real(kind=rk) :: res
    ! --------------------------------------------------------------------------

    if (t <= me%from_time) then
      res = me%min_factor
    else if (t >= me%to_time) then
      res = me%max_factor
    else
      res = me%min_factor                                                      &
        & + (me%max_factor - me%min_factor)                                    &
        & * sin( PI * 0.5_rk * (t          - me%from_time)&
        &                    / (me%to_time - me%from_time))
    end if
  end function temporal_smooth_for
! **************************************************************************** !

! **************************************************************************** !
  !> This function evaluate lua function and return its result
  !!
  function temporal_lua_for(fun_ref, time, conf) result(res)
    ! --------------------------------------------------------------------------
    !> Lua reference to the function to evaluate.
    integer, intent(in) :: fun_ref
    !> timer object incl. the current time information
    type(tem_time_type), intent( in )  :: time
    !> optional lua state
    type(flu_State) :: conf
    !> return value
    real(kind=rk) :: res
    ! --------------------------------------------------------------------------
    ! local variables
    type(aot_fun_type) :: fun
    integer :: iError
    ! --------------------------------------------------------------------------

    call aot_fun_open( L=conf, fun=fun, ref=fun_ref )

    call aot_fun_put( L=conf, fun=fun, arg=time%sim )
    call aot_fun_do( L=conf, fun=fun, nresults=1 )
    call aot_top_get_val( L=conf, val=res, ErrCode=iError )

    call aot_fun_close( L=conf, fun=fun )

  end function temporal_lua_for
! **************************************************************************** !


! **************************************************************************** !
  !> This function searches for the right values in the periodic data read
  !! from file and interpolates them if needed.
  !!
  !! Two different kinds of interpolation between the data tuples are available.
  !! \li \a 'none': No interpolation between the data points is performed.
  !!                The data is evaluated as ]t1,t2] -> v1, ]t2,t3] -> v2, ...
  !! \li \a 'linear': The value is interpolated in a linear fashion from the
  !!                  provided data.
  !!
  function temporal_from_file_periodic_for( me, time ) result(res)
    ! --------------------------------------------------------------------------
    !> datatype incl. the data read from file
    type(tem_from_file_temporal_type) :: me
    !> timer object incl. the current time information
    type(tem_time_type), intent(in)  :: time
    !> return value
    real(kind=rk) :: res
    ! --------------------------------------------------------------------------
    ! local variables
    integer :: iData
    ! period in physical units
    real(kind=rk) :: rPeriod
    real(kind=rk) :: ratio
    ! remainder in physical units
    real(kind=rk) :: rRem
    ! --------------------------------------------------------------------------

    ! calculate the period in physical units (last time signal
    ! - first time signal)
    rPeriod = me%signal%val(1,me%signal%nVals) - me%signal%val(1,1)
    ratio = time%sim / rPeriod
    rRem =  (ratio - floor(ratio))*rPeriod
    res = 1.0_rk

    ! search for the correct entry in the data read from file
    do iData = 1, me%signal%nVals-1
      if ( ( (rRem+eps) >= me%signal%val(1, iData) ) .and.                   &
        &  ( (rRem) <= me%signal%val(1, iData+1) ) )then
        exit
      end if
    end do

    if( .not. me%ramp .or. time%sim >= me%rampT)then
      select case( trim(adjustl(me%intp)) )
        case( 'none' )
          ! result is the data stored in the position directly
          res = me%signal%val(2, iData)
        case( 'linear' )
          ! result is the linear interpolation between the data from
          ! iData and iData+1
          res = ( me%signal%val(2, iData+1) - me%signal%val(2, iData) ) &
            & / ( me%signal%val(1, iData+1) - me%signal%val(1, iData) ) &
            & * (rRem - me%signal%val(1, iData) )                       &
            & + me%signal%val(2, iData)
      end select
    else
      res = me%rampVal * sin( PI * 0.5_rk * (time%sim/me%rampT))
    end if

  end function temporal_from_file_periodic_for
! **************************************************************************** !


! **************************************************************************** !
  !> This function invokes the type of the boundary such as constant, lua or
  !! predefined Fortran function.
  !!
  !! If temporal block is not defined than it returns value = 1.0.
  function tem_temporal_for( temporal, time ) result( res )
    ! --------------------------------------------------------------------------
    !> boundary state
    type( tem_temporal_type ) :: temporal
    !> timer object incl. the current time information
    type(tem_time_type), intent( in )  :: time
    !> return value of a function
    real(kind=rk) :: res
    ! --------------------------------------------------------------------------

    res = 0.0_rk
    select case( trim(adjustl(temporal%kind)) )
    case( 'none' )
      res = temporal%const
    case( 'const' )
      res = temporal%const
    case( 'lua_fun' )
      res = temporal_lua_for( fun_ref = temporal%lua_fun_ref, &
        &                      time    = time,                &
        &                      conf    = temporal%conf        )
    case( 'linear' )
      res = temporal_linear_for( temporal%linear, time%sim )
    case( 'smooth' )
      res = temporal_smooth_for( temporal%linear, time%sim )
    case( 'datafile' )
      if ( temporal%from_file%periodic )then
        res = temporal_from_file_periodic_for( me     = temporal%from_file,  &
          &                                     time   = time )
      else
        write(logUnit(1),*)'Only Periodic datasets are currently implemented!'
        write(logUnit(1),*)'STOPPING'
        call tem_abort()
      end if
    case( 'cos' )
      res = cos( 2.0_rk * PI * temporal%freq * time%sim + temporal%phi ) &
        & + temporal%offset
    end select

  end function tem_temporal_for
! **************************************************************************** !


end module tem_temporal_module
! **************************************************************************** !
