! Copyright (c) 2011-2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2011, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2012 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> Definitions and routines to help deployment
!! in the environment of the hosting system.
!! @todo: Should generate MPI infos during the build phase in waf.
module env_module

  ! include treelm modules
  use mpi
  use iso_fortran_env, only: output_unit, input_unit

  ! include aotus modules
  use flu_binding,  only: flu_State
  use aotus_module, only: aot_get_val

  implicit none
  private

  public :: output_unit

  !> Minimal number of entries in growing/dynamic array
  integer, parameter, public :: minLength = 1
  integer, parameter, public :: zeroLength = 0

  !> Number of digits for single precision numbers
  integer, parameter, public :: single_prec = 6
  !> Number of digits for double precision numbers
  integer, parameter, public :: double_prec = 15
  !> Number of digits for extended double precision numbers
  integer, parameter, public :: xdble_prec = 18
  !> Number of digits for quadruple precision numbers
  integer, parameter, public :: quad_prec = 33

  !> Number of digits for long integer numbers
  integer, parameter, public :: long_prec = 15
  !> Number of digits for normal long integer numbers
  integer, parameter, public :: int_prec = 6

  !> If the default precision of reals is to be changed, use this
  !! parameter. The real kind (rk will then be selected accordingly).
  !! For your convenience there are two predefined parameters available,
  !! to select single or double precision: double_prec and single_prec.
  integer, parameter, public :: rk_prec = double_prec

  !> The kind to select for default reals
  integer, parameter, public :: rk = selected_real_kind(rk_prec)

  !> The kind to select for reals using quadruple-precision
  integer, parameter, public :: quad_k   = selected_real_kind(quad_prec)
  !> The kind to select for reals using double-precision
  integer, parameter, public :: double_k = selected_real_kind(double_prec)
  !> The kind to select for reals using single-precision
  integer, parameter, public :: single_k = selected_real_kind(single_prec)

  !> The kind to select for default long integers
  integer, parameter, public :: long_k = selected_int_kind(long_prec)
  !> The kind to select for default normal length integers
  integer, parameter, public :: int_k = selected_int_kind(int_prec)

  !> machine error for selected real kind, rk
  real(kind=rk), parameter, public :: eps = epsilon(1._rk)

  !> close value to 0, rk
  real(kind=rk), parameter, public :: zero_rk = 1024*tiny(1._rk)

  !> machine error for selected real kind, single_k
  real(kind=rk), parameter, public :: eps_single_k = epsilon(1.0_single_k)

  !> Length for labels describing selectable options, like names of
  !! boundary or initial conditions.
  integer, parameter, public :: LabelLen = 80

  !> Length for strings describing file system paths.
  integer, parameter, public :: PathLen = 800

  !> Length for data strings dumped to file in output routines
  !! e.g. ASCII output (1kb).
  integer, parameter, public :: OutLen = 1024

  !> Length for strings describing the solver specific character (64kb).
  integer, parameter, public :: SolSpecLen = 7198

  !> Length for MPI I/O errors.
  integer, parameter, public :: ErrorLen = 100

  !> The default output unit, inherited from the iso_fortran_env module
  !! @todo if the iso_fortran_env module is not available, this should
  !! be set to the usual 6 during compilation.
  integer, parameter, public :: stdOutUnit = output_unit

  !> Maximal level of refinements in the Octree
  integer, parameter, public :: globalMaxLevels = 20

  !> Endianness of the system, is going to be determined at runtime in init_env.
  logical, public :: isLittleEndian

  !> path seperator of the system
  character,parameter, public :: pathSep = '/'

  !> Name of the null device to use for output to be discarded.
  character(len=pathLen), public :: null_device = '/dev/null'

  !> MPI Datatype for the used real kind, set in init_env.
  integer, public :: rk_mpi
  !> MPI Datatype for the used long kind, set in init_env.
  integer, public :: long_k_mpi

  !> Size of Buffer for IO in 8 Byte words. (Use a sane default of 8 MB)
  integer, save, public :: io_buffer_size = 1048576

  !> Should the information on the runtime obtained from /proc/self/status
  !! be printed during finalize?
  !!
  !! If this Linux pseudo file is not accessible on the current system, a
  !! reporting message will be printed instead.
  !! The data will only be printed by the first process and thus only shows
  !! its state in parallel runs.
  logical, save, public :: printRuntimeInfo = .true.

  !> Some parameters for the error handling.
  !!
  !! They indicate the bits to set in case of
  !! the corresponding error, to allow appropiate
  !! reactions of the calling application.
  integer, parameter, public :: tem_err_fatal = 0
  integer, parameter, public :: tem_err_warning = 1

  public :: init_env
  public :: fin_env
  public :: newunit
  public :: Sys_is_Little_Endian
  public :: tem_load_env_params
  public :: print_self_status
  public :: tem_create_EndianSuffix
  public :: my_status_string, my_status_int
  public :: my_status_string_vec, my_status_int_vec
  public :: tem_connect_toNull
  public :: init_random_seed

contains

! ****************************************************************************** !
  !> Initialize the environment, should be the very first call in the
  !! program and includes the call of MPI_Init.
  !!
  subroutine init_env()
    ! ---------------------------------------------------------------------------
    integer :: iError
    integer :: rank
    ! ---------------------------------------------------------------------------
    call MPI_Init(iError)
    !call mpi_type_create_f90_real(rk_prec, MPI_UNDEFINED, rk_mpi, iError)
    !call mpi_type_create_f90_integer(long_prec, long_k_mpi, iError)
    ! Workaround shitty MPI-Implementations.
    if( rk_prec == double_prec ) then
      rk_mpi = mpi_double_precision
    elseif( rk_prec == quad_prec ) then
      rk_mpi = MPI_REAL16
    else
      write(*,*) 'unknown real type specified for mpi_type creation'
    end if
    long_k_mpi = MPI_INTEGER8
    isLittleEndian = Sys_is_Little_Endian()

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, iError)

    call init_random_seed( rank )
  end subroutine init_env
! ****************************************************************************** !

! ****************************************************************************** !
  !> Initialized random seed with idx
  subroutine init_random_seed( idx )
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: idx
    integer :: nSeeds
    integer, allocatable :: seed(:)
    integer :: i
    ! ---------------------------------------------------------------------------

    call random_seed(size = nSeeds)
    allocate(seed(nSeeds))
    seed = 37 * [ (idx*nSeeds + i, i=1,nSeeds) ]
    call random_seed(put = seed)
    deallocate(seed)

  end subroutine init_random_seed
! ****************************************************************************** !

! ****************************************************************************** !
  !> Finalize the environment, should be the very last call in the program, and
  !! includes the call of MPI_Finalize.
  !!
  subroutine fin_env()
    ! ---------------------------------------------------------------------------
    integer :: iError
    ! ---------------------------------------------------------------------------
    call MPI_Finalize(iError)
  end subroutine fin_env
! ****************************************************************************** !


! ****************************************************************************** !
  !> Helper function to provide new unit, as long as F2008 newunit argument
  !! in open statement is not commonly available.
  !! To be used right in front of the open statement like this:
  !!  myUnit = newunit()
  !!  open(myUnit, ...)
  !!
  function newunit() result(nu)
    ! ---------------------------------------------------------------------------
    ! Melven Zoellner:
    ! for the solspecunit and in other places the unit is not kept open,
    ! so this does only work correctly if we count up over different calls
    ! to this function, currently not thread save...
    !HK: Sorry this is not the task of newunit, it is supposed to serve
    !HK: as an replacement for the new_unit keyword in the open statement
    !HK: of Fortran, which will definitely not keep track of the units
    !HK: opened so far. Also for iterative calls to newunit, you could
    !HK: reach the limit of integers easily!
    !HK: If you need to keep track of a file over closes and subsequent
    !HK: reopenings, you need to do this in some other way!
    integer :: nu
    ! ---------------------------------------------------------------------------
    integer, parameter :: nu_start = 22
    logical :: connected
    ! ---------------------------------------------------------------------------

    nu = nu_start
    inquire(unit=nu, opened=connected)
    do while(connected)
      nu = nu + 1
      inquire(unit=nu, opened=connected)
    end do
  end function newunit
! ****************************************************************************** !


! ****************************************************************************** !
  !> Determine if the system is little or big endian
  !!
  function Sys_is_Little_Endian() result(endian)
    ! ---------------------------------------------------------------------------
    logical :: endian
    ! ---------------------------------------------------------------------------
    integer :: defInt
    character :: dummy(4)
    ! ---------------------------------------------------------------------------

    defInt = ichar('a')

    dummy = transfer(defInt, dummy(1), 4)

    endian = (dummy(1) == 'a')
  end function sys_is_little_endian
! ****************************************************************************** !


! ****************************************************************************** !
  !> Load globally configurable environment parameters from the config script
  !!
  subroutine tem_load_env_params( conf )
    ! ---------------------------------------------------------------------------
    !> lua state to load information from
    type(flu_state) :: conf
    ! ---------------------------------------------------------------------------
    real    :: sizeMB
    integer :: iError
    ! ---------------------------------------------------------------------------
    ! Get the io_buffer_size from a given lua file
    call aot_get_val( L       = conf,                                          &
      &               key     = 'io_buffer_size',                              &
      &               val     = sizeMB,                                        &
      &               ErrCode = iError,                                        &
      &               default = 8. )

    ! Get the null device which might be used to do outputs without effect.
    call aot_get_val( L       = conf,                                          &
      &               key     = 'null_device',                                 &
      &               val     = null_device,                                   &
      &               ErrCode = iError,                                        &
      &               default = '/dev/null')

    ! Configure, if we are to write the Runtime info in tem_finalize
    call aot_get_val( L       = conf,                                          &
      &               key     = 'printRuntimeInfo',                            &
      &               val     = printRuntimeInfo,                              &
      &               ErrCode = iError,                                        &
      &               default = .true. )

    ! Calculate the size of the buffer array and store it in io_buffer_size.
    ! The buffer is supposed to be built up by 8 byte reals, thus for a total
    ! size given in MB, the number entries is given by sizeMB * 1024*1024/8.
    io_buffer_size = int( sizeMB * 131072. )

  end subroutine tem_load_env_params
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine reads out the status of the process
  !! if it is available in /proc/self/status, which is
  !! provided by the Linux operating system.
  !!
  subroutine print_self_status(unit, info, text)
    ! ---------------------------------------------------------------------------
    !> A preconnected unit, to write the content
    !! of /proc/self/status to.
    integer, optional, intent(in) :: unit
    !> optional filename prefix
    character(len=*), optional, intent(in) :: info
    !> optional header text to be dumped to file
    character(len=*), optional, intent(in) :: text
    ! ---------------------------------------------------------------------------
    character(len=128) :: cInfo
    character(len=128) :: line
    integer :: stat
    integer :: inUnit
    integer :: outUnit
    ! ---------------------------------------------------------------------------
    if( present( info )) then
      cInfo = info
    else
      cInfo = '/proc/self/status'
    end if
    if( present( unit )) then
      outUnit = unit
    else
      outUnit = output_unit
    end if
    inUnit = newUnit()
    open( file   = trim(cInfo), &
      &   action = 'read',      &
      &   iostat = stat,        &
      &   unit   = inUnit       )

    if( present( text )) then
      write(outUnit, '(a)')                                                    &
        &  '-------------------------------------------------------'
      write(outUnit, '(2a)') '        ', trim(text)
    end if
    if (stat==0) then
      read(inUnit, '(a)', iostat=stat) line
      do while (stat==0)
        if (present(unit)) then
          write(unit,'(a)') trim(line)
        else
          write(*,'(a)') trim(line)
        end if
        read(inUnit, '(a)', iostat=stat) line
      end do

      close(inUnit)
    else
      write(outUnit,*) '/proc/self/'//trim(cInfo)//' not available'
    end if
    if( present( text )) then
      write(outUnit, '(a)')                                                    &
        &  '-------------------------------------------------------'
    end if

  end subroutine print_self_status
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function returns the string in the line which starts with the
  !! specified key string in the text returned from print_self_status.
  !! The key string itself is excluded from the returned string.
  !! If the specified key is not found in the text an empty string will be
  !! returned.
  !!
  function tem_create_EndianSuffix()  result(suffix)
    ! ---------------------------------------------------------------------------
    character(len=4) :: suffix
    ! ---------------------------------------------------------------------------

    if (isLittleEndian) then
      suffix = '.lsb'
    else
      suffix = '.msb'
    end if

  end function tem_create_endianSuffix
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function returns the string in the line which starts with the
  !! specified key string in the text returned from print_self_status.
  !! The key string itself is excluded from the returned string.
  !! If the specified key is not found in the text an empty string will be
  !! returned.
  !!
  !! info can be used to change the input that should be read, it defaults to
  !! '/proc/self/status'.
  function my_status_string(key, info) result(val)
    ! ---------------------------------------------------------------------------
    character(len=*), intent(in) :: key !< Case sensitive!
    character(len=*), intent(in), optional :: info
    character(len=80) :: val
    ! ---------------------------------------------------------------------------
    integer :: scratchUnit
    character(len=80) :: bufline
    integer :: ios
    integer :: keylen
    logical :: match
    ! ---------------------------------------------------------------------------

    scratchUnit = newUnit()

    val = ''
    keylen = len_trim(key)
    open(unit = scratchUnit, status = 'scratch')

    call print_self_status(scratchUnit, info=info)
    rewind(scratchUnit)
    read(scratchUnit, '(a)', iostat=ios) bufline
    match = (bufline(1:keylen) == trim(key))
    do while ((ios == 0) .and. (.not.match))
      read(scratchUnit, '(a)', iostat=ios) bufline
      match = (bufline(1:keylen) == trim(key))
    end do
    close(scratchUnit)
    if (match) val = trim(bufline(len_trim(key)+1:))
    val = adjustl(val)
  end function my_status_string
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function returns the strings in the first line which starts with the
  !! specified key strings in the text returned from print_self_status.
  !!
  !! The key strings themselves are excluded from the returned string.
  !! If the specified key is not found in the text an empty string will be
  !! returned.
  !!
  !! info can be used to change the input that should be read, it defaults to
  !! '/proc/self/status'.
  function my_status_string_vec(key, info) result(val)
    ! ---------------------------------------------------------------------------
    character(len=*), intent(in) :: key(:) !< Case sensitive!
    character(len=*), intent(in), optional :: info
    character(len=80) :: val(size(key))
    ! ---------------------------------------------------------------------------
    integer :: scratchUnit
    character(len=80) :: bufline
    integer :: ios
    integer :: n, i
    integer :: keylen(size(key))
    logical :: match
    ! ---------------------------------------------------------------------------

    scratchUnit = newUnit()
    n = size(key)

    val = ''
    keylen = len_trim(key)
    open(unit = scratchUnit, status = 'scratch')

    call print_self_status(scratchUnit, info=info)
    do i=1,n
      rewind(scratchUnit)
      do
        read(scratchUnit, '(a)', iostat=ios) bufline
        if (ios /= 0) EXIT  ! end of file reached or reading error occurred
        match = ( bufline(1:keylen(i)) == trim(key(i)) )
        if (match) then
          val(i) = trim(bufline(keylen(i)+1:))
          EXIT
        end if
      end do
    end do
    close(scratchUnit)
    val = adjustl(val)
  end function my_status_string_vec
! ****************************************************************************** !


! ****************************************************************************** !
  !> Read out the first integer in the line matched by the specified key
  !! in the text provided by print_self_status.
  !! The key should include all text in front of the integer.
  !! For example the peak memory consumption up to now could be extracted with:
  !!   peakval = my_status_int('VmPeak:')
  !! If the key is not found 0 or the optionally defined value in def
  !! will be returned.
  !!
  !! info can be used to change the input that should be read, it defaults to
  !! '/proc/self/status'.
  function my_status_int(key, info, def) result(val)
    ! ---------------------------------------------------------------------------
    character(len=*), intent(in) :: key !< Case sensitive!
    character(len=*), intent(in), optional :: info
    integer, intent(in), optional :: def !< Result for non-existing keys
    integer :: val
    ! ---------------------------------------------------------------------------
    character(len=labelLen) :: vstring
    ! ---------------------------------------------------------------------------

    val = 0
    if (present(def)) val = def

    vstring = my_status_string(key, info=info)

    if (scan(vstring, set='0123456789') > 0) then
      read(vstring,*) val
    end if

  end function my_status_int
! ****************************************************************************** !


! ****************************************************************************** !
  !> Reads out the first integer in the line matched by the specified key
  !! in the text provided by print_self_status.
  !! The key should include all text in front of the integer.
  !! For example the peak memory consumption up to now could be extracted with:
  !! `peakval = my_status_int('VmPeak:')`
  !! If the key is not found 0 or the optionally defined value in def
  !! will be returned.
  !!
  !! info can be used to change the input that should be read, it defaults to
  !! '/proc/self/status'.
  function my_status_int_vec(key, info, def) result(val)
    ! ---------------------------------------------------------------------------
    character(len=*), intent(in) :: key(:) !< Case sensitive!
    character(len=*), intent(in), optional :: info
    integer, intent(in), optional :: def(:) !< Result for non-existing keys
    integer :: val(size(key))
    ! ---------------------------------------------------------------------------
    character(len=labelLen) :: vstring(size(key))
    integer :: i,n
    ! ---------------------------------------------------------------------------

    val = 0
    if (present(def)) val = def
    n = size(key)

    vstring = my_status_string_vec(key, info=info)

    do i=1,n
      if (scan(vstring(i), set='0123456789') > 0) then
        read(vstring(i),*) val(i)
      end if
    end do

  end function my_status_int_vec
! ****************************************************************************** !


! ****************************************************************************** !
  subroutine tem_connect_toNull(null_unit)
    integer, intent(out) :: null_unit

    logical :: null_connected

    inquire(file=trim(null_device), opened=null_connected)
    if (null_connected) then
      inquire(file=trim(null_device), number=null_unit)
      ! Some MPI implementations might connect stdin to the null device.
      ! Need to check if this is the case and attempt to open it again
      ! with a different unit to make it writable.
      if (null_unit == input_unit) then
        null_unit = newUnit()
        open(file=trim(null_device), unit=null_unit)
      end if
    else
      null_unit = newUnit()
      open(file=trim(null_device), unit=null_unit)
    end if

  end subroutine tem_connect_toNull
! ****************************************************************************** !


end module env_module
! ****************************************************************************** !
