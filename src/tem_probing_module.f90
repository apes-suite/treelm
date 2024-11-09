! Copyright (c) 2017, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> This module provides some little functionality to easily write
!! data into files for debugging.
!!
!! It is related to [[tem_debug_module]], but intended for even smaller
!! tasks and specifically for tracking the evolution of values over the
!! runtime.
module tem_probing_module
  use env_module, only: newunit, rk

  implicit none


contains


  ! ------------------------------------------------------------------------ !
  !> Delete the file with the given name.
  !!
  !! Use this in the beginning of the application, to clear old logs.
  subroutine tem_probing_delete(filename)
    ! -------------------------------------------------------------------- !
    !> Name of the file to delete.
    character(len=*) :: filename
    ! -------------------------------------------------------------------- !
    logical :: is_oldfile
    integer :: probeunit
    ! -------------------------------------------------------------------- !

    inquire(file=trim(filename), exist = is_oldfile)
    probeunit = newunit()
    if (is_oldfile) then
      open( probeunit, file=filename, status='old', position='append', &
        &   action='write' )
    else
      open( probeunit, file=filename, status='new', action='write' )
    end if
    close(probeunit, status='delete')
  end subroutine tem_probing_delete
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Write an array of reals into a file of the given name.
  subroutine tem_probing_write(filename, label, dat)
    ! -------------------------------------------------------------------- !
    !> Name of the file to write the data to.
    !! If the file already exists, the data will be appended.
    character(len=*), intent(in) :: filename

    !> Label to identify the point in the program execution, where the
    !! data was written.
    character(len=*), intent(in) :: label

    !> Data to write to the file
    real(kind=rk), intent(in) :: dat(:)
    ! -------------------------------------------------------------------- !
    logical :: is_oldfile
    integer :: probeunit
    ! -------------------------------------------------------------------- !

    inquire(file=trim(filename), exist = is_oldfile)
    probeunit = newunit()
    if (is_oldfile) then
      open( probeunit, file=filename, status='old', position='append', &
        &   action='write' )
    else
      open( probeunit, file=filename, status='new', action='write' )
    end if
    write(probeunit,*) dat, ' !', trim(label)

    close(probeunit)
  end subroutine tem_probing_write
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !

end module tem_probing_module
