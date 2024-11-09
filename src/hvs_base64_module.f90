! Copyright (c) 2015 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
! ****************************************************************************** !
!> This module provides the methods for base64Encoding.
!!
module hvs_base64_module

  ! include treelm modules
  use env_module, only: rk

  ! include harvester modules
  use hvs_sizeof_module, only: c_sizeof

  use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_float,      &
    &                                    c_int_least8_t, c_int_least64_t,      &
    &                                    c_char, c_loc

  implicit none

  integer, parameter :: EncoderLen = 4
  integer, parameter :: EncoderBlockLen = 3

  !> Convert one block of input (3 Bytes) into its base64
  !! representation (4 Bytes).
  interface
    function EncodeBlock( input, output, iplen, oplen )                        &
      & bind(c, name='encodeblock')
      use, intrinsic :: iso_c_binding
      !> input character stream in ascii
      type(c_ptr), value :: input

      !> output character stream in base64
      type(c_ptr), value :: output

      !> length of the input stream
      integer(kind=c_int), value, intent(in) :: iplen

      !> length of the output stream
      integer(kind=c_int), value :: oplen

      !> Result, indicating the status of encode
      integer(kind=c_int) :: EncodeBlock
    end function EncodeBlock
  end interface

  !> Interface to convert ascii to binary base64 encoder
  interface
    function EncodeIndex( input, output, iplen, ind, ipindex )                 &
      & bind(c, name='EncodeIndex')
      use, intrinsic :: iso_c_binding
      !> input character stream in ascii
      type(c_ptr), value :: input

      !> output character stream in base64
      type(c_ptr), value :: output

      !> length of the input stream
      integer(kind=c_int), value, intent(in) :: iplen

      !> current length
      integer(kind=c_int), value :: ind

      !> current length
      integer(kind=c_int), value :: ipindex

      !> Result, indicating the status of encode
      integer(kind=c_int) :: EncodeIndex
    end function EncodeIndex
  end interface

  !> Interface to convert ascii to binary base64 encoder
  interface
    function Base64Encode( input, output, iplen, oplen, ipindex )              &
      & bind(c, name='Base64Encode')
      use, intrinsic :: iso_c_binding
      !> input character stream in ascii
      type(c_ptr), value :: input

      !> output character stream in base64
      type(c_ptr), value :: output

      !> length of the input stream
      integer(kind=c_int), value, intent(in) :: iplen

      !> length of the output stream
      integer(kind=c_int), value :: oplen

      !> current length
      integer(kind=c_int), value :: ipindex

      !> Result, indicating the status of encode
      integer(kind=c_int) :: Base64Encode
    end function Base64Encode
  end interface

  interface convert_to_Base64
    module procedure real32_to_base64
    module procedure real64_to_base64
    module procedure char_to_base64
    module procedure int8_to_base64
    module procedure int32_to_base64
    module procedure int64_to_base64
  end interface convert_to_Base64

  interface convert_to_Base64_single
    module procedure real32_to_base64_single
    module procedure real64_to_base64_single
    module procedure char_to_base64_single
    module procedure int8_to_base64_single
    module procedure int32_to_base64_single
    module procedure int64_to_base64_single
  end interface convert_to_Base64_single


contains



! ****************************************************************************** !
  !> this routine encodes data of type real64 to base64 format
  !!
  subroutine real64_to_base64( indata, iplen, outfile )
    ! ---------------------------------------------------------------------------
    !> size of data to be encoded
    integer, intent(in) :: iplen
    !> data to be encoded
    real(kind=c_double), target, intent(in) :: indata(iplen)
    !> output file unit
    integer, intent(in) :: outfile
    ! ---------------------------------------------------------------------------
    integer(kind=c_int) :: baserc
    integer(kind=c_int) :: insize, outsize, ipindex, min_iplen
    integer :: ind
    type(c_ptr) :: base64_out
    type(c_ptr) :: base64_in
    type(c_ptr) :: encoder_in
    character, target :: base64_string(encoderlen), encoder_str(encoderblocklen)
    ! ---------------------------------------------------------------------------

    base64_in = c_loc(indata)

    insize = int(c_sizeof(indata(1))*iplen, kind=4)

    ! write insize i.e bit size in the beginning of the string
    call convert_to_base64_single( insize, outfile )

    outsize = ceiling(insize/3._rk)*4

    base64_out = c_loc(base64_string)

    encoder_in = c_loc(encoder_str)

    ipindex = 0
    do
      do ind = 0, 2
        baserc = encodeindex( base64_in, encoder_in, insize, ind, ipindex )
        ipindex = ipindex + 1
      end do
      min_iplen = min(insize - ipindex + 3, 3)
      baserc = encodeblock( encoder_in, base64_out, min_iplen, outsize )
      write(outfile) base64_string
      if (ipindex >= insize) exit
    end do

  end subroutine real64_to_base64
! ****************************************************************************** !


! ****************************************************************************** !
  !> this routine encodes a single variable of type real64 into base64 format
  !!
  subroutine real64_to_base64_single( indata, outfile )
    ! ---------------------------------------------------------------------------
    !> data to be encoded
    real(kind=c_double), target, intent(in) :: indata
    !> output file unit
    integer, intent(in) :: outfile
    ! ---------------------------------------------------------------------------
    integer(kind=c_int) :: baserc
    integer(kind=c_int) :: insize, outsize, ipindex, min_iplen
    integer :: ind
    type(c_ptr) :: base64_out
    type(c_ptr) :: base64_in
    type(c_ptr) :: encoder_in
    character, target :: base64_string(encoderlen), encoder_str(encoderblocklen)
    ! ---------------------------------------------------------------------------

    base64_in = c_loc(indata)

    insize = int(c_sizeof(indata), kind=4)

    outsize = ceiling(insize/3._rk)*4

    base64_out = c_loc(base64_string)

    encoder_in = c_loc(encoder_str)

    ipindex = 0
    do
      do ind = 0, 2
        baserc = encodeindex( base64_in, encoder_in, insize, ind, ipindex )
        ipindex = ipindex + 1
      end do
      min_iplen = min(insize - ipindex + 3, 3)
      baserc = encodeblock( encoder_in, base64_out, min_iplen, outsize )
      write(outfile) base64_string
      if (ipindex >= insize) exit
    end do

  end subroutine real64_to_base64_single
! ****************************************************************************** !

! ****************************************************************************** !
  !> this routine encodes data of type real32 to base64 format
  !!
  subroutine real32_to_base64( indata, iplen, outfile )
    ! ---------------------------------------------------------------------------
    !> size of data to be encoded
    integer, intent(in) :: iplen
    !> data to be encoded
    real(kind=c_float), target, intent(in) :: indata(iplen)
    !> output file unit
    integer, intent(in) :: outfile
    ! ---------------------------------------------------------------------------
    integer(kind=c_int) :: baserc
    integer(kind=c_int) :: insize, outsize, ipindex, min_iplen
    integer :: ind
    type(c_ptr) :: base64_out
    type(c_ptr) :: base64_in
    type(c_ptr) :: encoder_in
    character, target :: base64_string(encoderlen), encoder_str(encoderblocklen)
    ! ---------------------------------------------------------------------------

    base64_in = c_loc(indata)

    insize = int(c_sizeof(indata(1))*iplen, kind=4)

    ! write insize i.e bit size in the beginning of the string
    call convert_to_base64_single( insize, outfile )

    outsize = ceiling(insize/3._rk)*4

    base64_out = c_loc(base64_string)

    encoder_in = c_loc(encoder_str)

    ipindex = 0
    do
      do ind = 0, 2
        baserc = encodeindex( base64_in, encoder_in, insize, ind, ipindex )
        ipindex = ipindex + 1
      end do
      min_iplen = min(insize - ipindex + 3, 3)
      baserc = encodeblock( encoder_in, base64_out, min_iplen, outsize )
      write(outfile) base64_string
      if (ipindex >= insize) exit
    end do

  end subroutine real32_to_base64
! ****************************************************************************** !


! ****************************************************************************** !
  !> this routine encodes a single variable of type real32 into base64 format
  !!
  subroutine real32_to_base64_single( indata, outfile )
    ! ---------------------------------------------------------------------------
    !> data to be encoded
    real(kind=c_float), target, intent(in) :: indata
    !> output file unit
    integer, intent(in) :: outfile
    ! ---------------------------------------------------------------------------
    integer(kind=c_int) :: baserc
    integer(kind=c_int) :: insize, outsize, ipindex, min_iplen
    integer :: ind
    type(c_ptr) :: base64_out
    type(c_ptr) :: base64_in
    type(c_ptr) :: encoder_in
    character, target :: base64_string(encoderlen), encoder_str(encoderblocklen)
    ! ---------------------------------------------------------------------------

    base64_in = c_loc(indata)

    insize = int(c_sizeof(indata), kind=4)

    outsize = ceiling(insize/3._rk)*4

    base64_out = c_loc(base64_string)

    encoder_in = c_loc(encoder_str)

    ipindex = 0
    do
      do ind = 0, 2
        baserc = encodeindex( base64_in, encoder_in, insize, ind, ipindex )
        ipindex = ipindex + 1
      end do
      min_iplen = min(insize - ipindex + 3, 3)
      baserc = encodeblock( encoder_in, base64_out, min_iplen, outsize )
      write(outfile) base64_string
      if (ipindex >= insize) exit
    end do

  end subroutine real32_to_base64_single
! ****************************************************************************** !

! ****************************************************************************** !
  !> this routine encodes data of type int8 to base64 format
  !!
  subroutine int8_to_base64( indata, iplen, outfile )
    ! ---------------------------------------------------------------------------
    !> size of data to be encoded
    integer, intent(in) :: iplen
    !> data to be encoded
    integer(kind=c_int_least8_t), target, intent(in) :: indata(iplen)
    !> output file unit
    integer, intent(in) :: outfile
    ! ---------------------------------------------------------------------------
    integer(kind=c_int) :: baserc
    integer(kind=c_int) :: insize, outsize, ipindex, min_iplen
    integer :: ind
    type(c_ptr) :: base64_out
    type(c_ptr) :: base64_in
    type(c_ptr) :: encoder_in
    character, target :: base64_string(encoderlen), encoder_str(encoderblocklen)
    ! ---------------------------------------------------------------------------

    base64_in = c_loc(indata)

    insize = int(c_sizeof(indata(1))*iplen, kind=4)

    ! write insize i.e bit size in the beginning of the string
    call convert_to_base64_single( insize, outfile )

    outsize = ceiling(insize/3._rk)*4

    base64_out = c_loc(base64_string)

    encoder_in = c_loc(encoder_str)

    ipindex = 0
    do
      do ind = 0, 2
        baserc = encodeindex( base64_in, encoder_in, insize, ind, ipindex )
        ipindex = ipindex + 1
      end do
      min_iplen = min(insize - ipindex + 3, 3)
      baserc = encodeblock( encoder_in, base64_out, min_iplen, outsize )
      write(outfile) base64_string
      if (ipindex >= insize) exit
    end do

  end subroutine int8_to_base64
! ****************************************************************************** !


! ****************************************************************************** !
  !> this routine encodes a single variable of type int8 into base64 format
  !!
  subroutine int8_to_base64_single( indata, outfile )
    ! ---------------------------------------------------------------------------
    !> data to be encoded
    integer(kind=c_int_least8_t), target, intent(in) :: indata
    !> output file unit
    integer, intent(in) :: outfile
    ! ---------------------------------------------------------------------------
    integer(kind=c_int) :: baserc
    integer(kind=c_int) :: insize, outsize, ipindex, min_iplen
    integer :: ind
    type(c_ptr) :: base64_out
    type(c_ptr) :: base64_in
    type(c_ptr) :: encoder_in
    character, target :: base64_string(encoderlen), encoder_str(encoderblocklen)
    ! ---------------------------------------------------------------------------

    base64_in = c_loc(indata)

    insize = int(c_sizeof(indata), kind=4)

    outsize = ceiling(insize/3._rk)*4

    base64_out = c_loc(base64_string)

    encoder_in = c_loc(encoder_str)

    ipindex = 0
    do
      do ind = 0, 2
        baserc = encodeindex( base64_in, encoder_in, insize, ind, ipindex )
        ipindex = ipindex + 1
      end do
      min_iplen = min(insize - ipindex + 3, 3)
      baserc = encodeblock( encoder_in, base64_out, min_iplen, outsize )
      write(outfile) base64_string
      if (ipindex >= insize) exit
    end do

  end subroutine int8_to_base64_single
! ****************************************************************************** !

! ****************************************************************************** !
  !> this routine encodes data of type int32 to base64 format
  !!
  subroutine int32_to_base64( indata, iplen, outfile )
    ! ---------------------------------------------------------------------------
    !> size of data to be encoded
    integer, intent(in) :: iplen
    !> data to be encoded
    integer(kind=c_int), target, intent(in) :: indata(iplen)
    !> output file unit
    integer, intent(in) :: outfile
    ! ---------------------------------------------------------------------------
    integer(kind=c_int) :: baserc
    integer(kind=c_int) :: insize, outsize, ipindex, min_iplen
    integer :: ind
    type(c_ptr) :: base64_out
    type(c_ptr) :: base64_in
    type(c_ptr) :: encoder_in
    character, target :: base64_string(encoderlen), encoder_str(encoderblocklen)
    ! ---------------------------------------------------------------------------

    base64_in = c_loc(indata)

    insize = int(c_sizeof(indata(1))*iplen, kind=4)

    ! write insize i.e bit size in the beginning of the string
    call convert_to_base64_single( insize, outfile )

    outsize = ceiling(insize/3._rk)*4

    base64_out = c_loc(base64_string)

    encoder_in = c_loc(encoder_str)

    ipindex = 0
    do
      do ind = 0, 2
        baserc = encodeindex( base64_in, encoder_in, insize, ind, ipindex )
        ipindex = ipindex + 1
      end do
      min_iplen = min(insize - ipindex + 3, 3)
      baserc = encodeblock( encoder_in, base64_out, min_iplen, outsize )
      write(outfile) base64_string
      if (ipindex >= insize) exit
    end do

  end subroutine int32_to_base64
! ****************************************************************************** !


! ****************************************************************************** !
  !> this routine encodes a single variable of type int32 into base64 format
  !!
  subroutine int32_to_base64_single( indata, outfile )
    ! ---------------------------------------------------------------------------
    !> data to be encoded
    integer(kind=c_int), target, intent(in) :: indata
    !> output file unit
    integer, intent(in) :: outfile
    ! ---------------------------------------------------------------------------
    integer(kind=c_int) :: baserc
    integer(kind=c_int) :: insize, outsize, ipindex, min_iplen
    integer :: ind
    type(c_ptr) :: base64_out
    type(c_ptr) :: base64_in
    type(c_ptr) :: encoder_in
    character, target :: base64_string(encoderlen), encoder_str(encoderblocklen)
    ! ---------------------------------------------------------------------------

    base64_in = c_loc(indata)

    insize = int(c_sizeof(indata), kind=4)

    outsize = ceiling(insize/3._rk)*4

    base64_out = c_loc(base64_string)

    encoder_in = c_loc(encoder_str)

    ipindex = 0
    do
      do ind = 0, 2
        baserc = encodeindex( base64_in, encoder_in, insize, ind, ipindex )
        ipindex = ipindex + 1
      end do
      min_iplen = min(insize - ipindex + 3, 3)
      baserc = encodeblock( encoder_in, base64_out, min_iplen, outsize )
      write(outfile) base64_string
      if (ipindex >= insize) exit
    end do

  end subroutine int32_to_base64_single
! ****************************************************************************** !

! ****************************************************************************** !
  !> this routine encodes data of type int64 to base64 format
  !!
  subroutine int64_to_base64( indata, iplen, outfile )
    ! ---------------------------------------------------------------------------
    !> size of data to be encoded
    integer, intent(in) :: iplen
    !> data to be encoded
    integer(kind=c_int_least64_t), target, intent(in) :: indata(iplen)
    !> output file unit
    integer, intent(in) :: outfile
    ! ---------------------------------------------------------------------------
    integer(kind=c_int) :: baserc
    integer(kind=c_int) :: insize, outsize, ipindex, min_iplen
    integer :: ind
    type(c_ptr) :: base64_out
    type(c_ptr) :: base64_in
    type(c_ptr) :: encoder_in
    character, target :: base64_string(encoderlen), encoder_str(encoderblocklen)
    ! ---------------------------------------------------------------------------

    base64_in = c_loc(indata)

    insize = int(c_sizeof(indata(1))*iplen, kind=4)

    ! write insize i.e bit size in the beginning of the string
    call convert_to_base64_single( insize, outfile )

    outsize = ceiling(insize/3._rk)*4

    base64_out = c_loc(base64_string)

    encoder_in = c_loc(encoder_str)

    ipindex = 0
    do
      do ind = 0, 2
        baserc = encodeindex( base64_in, encoder_in, insize, ind, ipindex )
        ipindex = ipindex + 1
      end do
      min_iplen = min(insize - ipindex + 3, 3)
      baserc = encodeblock( encoder_in, base64_out, min_iplen, outsize )
      write(outfile) base64_string
      if (ipindex >= insize) exit
    end do

  end subroutine int64_to_base64
! ****************************************************************************** !


! ****************************************************************************** !
  !> this routine encodes a single variable of type int64 into base64 format
  !!
  subroutine int64_to_base64_single( indata, outfile )
    ! ---------------------------------------------------------------------------
    !> data to be encoded
    integer(kind=c_int_least64_t), target, intent(in) :: indata
    !> output file unit
    integer, intent(in) :: outfile
    ! ---------------------------------------------------------------------------
    integer(kind=c_int) :: baserc
    integer(kind=c_int) :: insize, outsize, ipindex, min_iplen
    integer :: ind
    type(c_ptr) :: base64_out
    type(c_ptr) :: base64_in
    type(c_ptr) :: encoder_in
    character, target :: base64_string(encoderlen), encoder_str(encoderblocklen)
    ! ---------------------------------------------------------------------------

    base64_in = c_loc(indata)

    insize = int(c_sizeof(indata), kind=4)

    outsize = ceiling(insize/3._rk)*4

    base64_out = c_loc(base64_string)

    encoder_in = c_loc(encoder_str)

    ipindex = 0
    do
      do ind = 0, 2
        baserc = encodeindex( base64_in, encoder_in, insize, ind, ipindex )
        ipindex = ipindex + 1
      end do
      min_iplen = min(insize - ipindex + 3, 3)
      baserc = encodeblock( encoder_in, base64_out, min_iplen, outsize )
      write(outfile) base64_string
      if (ipindex >= insize) exit
    end do

  end subroutine int64_to_base64_single
! ****************************************************************************** !

! ****************************************************************************** !
  !> this routine encodes data of type char to base64 format
  !!
  subroutine char_to_base64( indata, iplen, outfile )
    ! ---------------------------------------------------------------------------
    !> size of data to be encoded
    integer, intent(in) :: iplen
    !> data to be encoded
    character(kind=c_char), target, intent(in) :: indata(iplen)
    !> output file unit
    integer, intent(in) :: outfile
    ! ---------------------------------------------------------------------------
    integer(kind=c_int) :: baserc
    integer(kind=c_int) :: insize, outsize, ipindex, min_iplen
    integer :: ind
    type(c_ptr) :: base64_out
    type(c_ptr) :: base64_in
    type(c_ptr) :: encoder_in
    character, target :: base64_string(encoderlen), encoder_str(encoderblocklen)
    ! ---------------------------------------------------------------------------

    base64_in = c_loc(indata)

    insize = int(c_sizeof(indata(1))*iplen, kind=4)

    ! write insize i.e bit size in the beginning of the string
    call convert_to_base64_single( insize, outfile )

    outsize = ceiling(insize/3._rk)*4

    base64_out = c_loc(base64_string)

    encoder_in = c_loc(encoder_str)

    ipindex = 0
    do
      do ind = 0, 2
        baserc = encodeindex( base64_in, encoder_in, insize, ind, ipindex )
        ipindex = ipindex + 1
      end do
      min_iplen = min(insize - ipindex + 3, 3)
      baserc = encodeblock( encoder_in, base64_out, min_iplen, outsize )
      write(outfile) base64_string
      if (ipindex >= insize) exit
    end do

  end subroutine char_to_base64
! ****************************************************************************** !


! ****************************************************************************** !
  !> this routine encodes a single variable of type char into base64 format
  !!
  subroutine char_to_base64_single( indata, outfile )
    ! ---------------------------------------------------------------------------
    !> data to be encoded
    character(kind=c_char), target, intent(in) :: indata
    !> output file unit
    integer, intent(in) :: outfile
    ! ---------------------------------------------------------------------------
    integer(kind=c_int) :: baserc
    integer(kind=c_int) :: insize, outsize, ipindex, min_iplen
    integer :: ind
    type(c_ptr) :: base64_out
    type(c_ptr) :: base64_in
    type(c_ptr) :: encoder_in
    character, target :: base64_string(encoderlen), encoder_str(encoderblocklen)
    ! ---------------------------------------------------------------------------

    base64_in = c_loc(indata)

    insize = int(c_sizeof(indata), kind=4)

    outsize = ceiling(insize/3._rk)*4

    base64_out = c_loc(base64_string)

    encoder_in = c_loc(encoder_str)

    ipindex = 0
    do
      do ind = 0, 2
        baserc = encodeindex( base64_in, encoder_in, insize, ind, ipindex )
        ipindex = ipindex + 1
      end do
      min_iplen = min(insize - ipindex + 3, 3)
      baserc = encodeblock( encoder_in, base64_out, min_iplen, outsize )
      write(outfile) base64_string
      if (ipindex >= insize) exit
    end do

  end subroutine char_to_base64_single
! ****************************************************************************** !


end module hvs_base64_module
! ****************************************************************************** !
