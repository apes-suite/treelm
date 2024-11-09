! Copyright (c) 2012-2013, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
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
!> This module is used to track the requires in loaded Lua scripts.
!!
!! The purpose of this functionality is it to enable the loading of all Lua
!! scripts on a single process and broadcasting them to all others, without the
!! need to access the files by all processes.
module tem_lua_requires_module

  ! include treelm modules
  use, intrinsic :: iso_c_binding, only: C_NEW_LINE
  use env_module,                  only: labelLen, pathLen

  ! include aotus modules
  use flu_binding,      only: flu_State, cbuf_type, fluL_loadbuffer, flu_dump, &
    &                         flu_pushnil, flu_pop
  use aotus_module,     only: aot_get_val, open_config_chunk
  use aot_table_module, only: aot_table_open, aot_table_close,                 &
    &                         aot_table_set_top, aot_table_set_val,            &
    &                         aot_table_push

  implicit none

  private

  public :: tem_get_required_Lua, tem_require_track_rq, tem_pop_from_track_rq
  public :: tem_require_rq_store, tem_push_to_rq_store

  character, parameter :: nl = C_NEW_LINE

  !> This script amends the require and the Lua searcher, such that they report
  !! the loaded Lua files in the track_rq table.
  !!
  !! The filled table can then be used to fetch all required files for the
  !! given script.
  !! It therefore should be executed before the actual configuration, and
  !! evaluated after the configuration script is executed.
  character(len=1300), parameter :: track_rq_Lua =                             &
    &   "track_rq = {}" // nl                                                  &
    & //"track_rq.nesting = 0" // nl                                           &
    & //"track_rq.code = {}" // nl                                             &
    & //"track_rq.open_files = {}" // nl                                       &
    & //"track_rq.open_files.count = 0" // nl                                  &
    & //"track_rq.open_files.name = {}" // nl                                  &
    & //"track_rq.open_files.file = {}" // nl                                  &
    & //"track_rq.open_files.head = {}" // nl                                  &
    & // nl                                                                    &
    & //"track_rq.rq_orig = require" // nl                                     &
    & //"function require(name)" // nl                                         &
    & //"  track_rq.nesting = track_rq.nesting+1" // nl                        &
    & //"  local orig_ret = track_rq.rq_orig(name)" // nl                      &
    & //"  track_rq.nesting = track_rq.nesting-1" // nl                        &
    & //"  return orig_ret" // nl                                              &
    & //"end" // nl                                                            &
    & // nl                                                                    &
    & //"track_rq.search_orig = package.searchers[2]" // nl                    &
    & //"function track_rq.searcher(name)" // nl                               &
    & //"  local modfun, fname = track_rq.search_orig(name)" // nl             &
    & //"  local nest = track_rq.nesting" // nl                                &
    & //"  if fname then" // nl                                                &
    & //"    track_rq.open_files.count = track_rq.open_files.count + 1" // nl  &
    ! The heads for each nesting level are used to keep the required files in
    ! order. Deeper nested modules are put in front of those which require them.
    ! By this we can later on simply work down the list loading the various
    ! modules with all requires satisfied by the preloaded searcher.
    & //"    if nest > 1 then" // nl                                           &
    & //"      track_rq.open_files.head[nest] "                                &
    &       //"= track_rq.open_files.head[nest-1]" // nl                       &
    & //"      for i = 1,nest-1 do" // nl                                      &
    & //"        track_rq.open_files.head[i] "                                 &
    &         //"= track_rq.open_files.head[i] + 1" // nl                      &
    & //"      end" // nl                                                      &
    & //"    else" // nl                                                       &
    & //"      track_rq.open_files.head[nest] = track_rq.open_files.count"     &
    &          // nl                                                           &
    & //"    end" // nl                                                        &
    ! Now insert the found file at the sorted position.
    & //"    local head = track_rq.open_files.head[nest]" // nl                &
    & //"    table.insert(track_rq.open_files.name, head, name)" // nl         &
    & //"    table.insert(track_rq.open_files.file, head, fname)" // nl        &
    & //"    track_rq.code[name] = modfun" // nl                               &
    & //"    return modfun, fname" // nl                                       &
    & //"  else" // nl                                                         &
    ! An error was encountered.
    & //"    return modfun" // nl                                              &
    & //"  end" // nl                                                          &
    & //"end" // nl                                                            &
    & //"package.searchers[2] = track_rq.searcher"

  character(len=500), parameter :: rq_store_Lua =                              &
    &    "rq_store = {}" // nl                                                 &
    & // "rq_store.file = {}" // nl                                            &
    & // "rq_store.code = {}" // nl                                            &
    & // nl                                                                    &
    & // "function rq_store.searcher(name)" // nl                              &
    & // "  if rq_store.code[name] then" // nl                                 &
    & // "    local mycode = rq_store.code[name]" // nl                        &
    & // "    local myfile = rq_store.file[name]" // nl                        &
    & // "    rq_store.code[name] = nil" // nl                                 &
    & // "    rq_store.file[name] = nil" // nl                                 &
    & // "    return mycode, myfile" // nl                                     &
    & // "  else" // nl                                                        &
    & // "    return 'Module '..name..' not found in the rq_store!'" // nl     &
    & // "  end" // nl                                                         &
    & // "end" // nl                                                           &
    & // "package.searchers[2] = rq_store.searcher"

  contains

! ****************************************************************************** !
  !> Load the track_rq script into the Lua state.
  !!
  !! After calling this routine, all following requires in Lua will record
  !! the opened files in the track_rq table.
  !!
  subroutine tem_require_track_rq(L)
    ! ---------------------------------------------------------------------------
    !>
    type(flu_State) :: L
    ! ---------------------------------------------------------------------------

    call open_config_chunk(L, trim(track_rq_Lua))

  end subroutine tem_require_track_rq
! ****************************************************************************** !


! ****************************************************************************** !
  !> Subroutine to pop the code from the tracked required files.
  !!
  !! This routine reads the given module from the previously stored list of
  !! required codes into the buffer and removes it from the table in Lua.
  !!
  subroutine tem_pop_from_track_rq(L, modname, buffer)
    ! ---------------------------------------------------------------------------
    !>
    type(flu_State) :: L
    !> Name of the module to pop.
    character(len=*), intent(in) :: modname
    !> Buffer to hold the code.
    type(cbuf_type) :: buffer
    ! ---------------------------------------------------------------------------
    integer :: track_handle
    integer :: code_handle
    integer :: buflen
    integer :: iError
    ! ---------------------------------------------------------------------------

    call aot_table_open( L       = L,                                          &
      &                  thandle = track_handle,                               &
      &                  key     = 'track_rq' )
    call aot_table_open( L       = L,                                          &
      &                  parent  = track_handle,                               &
      &                  thandle = code_handle,                                &
      &                  key     = 'code' )

    ! Get field track_rq.code[modname].
    call aot_table_push(L, thandle = code_handle, key=trim(modname))
    ! Dump to buffer
    call flu_dump(L, buffer, buflen, iError)
    ! Pop it from the stack again.
    call flu_pop(L)
    ! Set field track_rq.code[modname] to nil.
    call flu_pushNil(L)
    call aot_table_set_top(L, code_handle, key=trim(modname))

    call aot_table_close(L, code_handle)
    call aot_table_close(L, track_handle)

  end subroutine tem_pop_from_track_rq
! ****************************************************************************** !


! ****************************************************************************** !
  !> Load the rq_store script into the Lua state.
  !!
  !! After calling this routine, following requires will not look up Lua scripts
  !! from the file system but instead try to fetch them from the rq_store.
  !!
  subroutine tem_require_rq_store(L)
    ! ---------------------------------------------------------------------------
    type(flu_State) :: L
    ! ---------------------------------------------------------------------------

    call open_config_chunk(L, trim(rq_store_Lua))

  end subroutine tem_require_rq_store
! ****************************************************************************** !


! ****************************************************************************** !
  !> Push a script in a buffer into the rq_store.
  !!
  !! This has to be done after tem_require_rq_store, any values put into the
  !! store before the call to tem_require_rq_store will be discarded.
  !!
  subroutine tem_push_to_rq_store(L, modname, filename, buffer)
    ! ---------------------------------------------------------------------------
    !>
    type(flu_State) :: L
    !> Name of the module to store.
    character(len=*), intent(in) :: modname
    !> Original file name.
    character(len=*), intent(in) :: filename
    !> Actual buffer holding the code.
    character, intent(in) :: buffer(:)
    ! ---------------------------------------------------------------------------
    integer :: store_handle
    integer :: list_handle
    integer :: iError
    ! ---------------------------------------------------------------------------

    call aot_table_open( L       = L,                                          &
      &                  thandle = store_handle,                               &
      &                  key     = "rq_store" )

    ! Set the filename.
    call aot_table_open( L       = L,                                          &
      &                  parent  = store_handle,                               &
      &                  thandle = list_handle,                                &
      &                  key     = "file" )
    call aot_table_set_val( val     = trim(filename),                          &
      &                     L       = L,                                       &
      &                     thandle = list_handle,                             &
      &                     key     = trim(modname))
    call aot_table_close(L, list_handle)

    ! Set the execution code.
    call aot_table_open( L       = L,                                          &
      &                  parent  = store_handle,                               &
      &                  thandle = list_handle,                                &
      &                  key     = "code")
    iError =  fluL_loadbuffer(L, buffer, trim(modname))
    call aot_table_set_top( L       = L,                                       &
      &                     thandle = list_handle,                             &
      &                     key     = trim(modname))
    call aot_table_close(L, list_handle)

    call aot_table_close(L, store_handle)

  end subroutine tem_push_to_rq_store
! ****************************************************************************** !


! ****************************************************************************** !
  !> Obtain the recorded open files by track_rq from the Lua state.
  !!
  !! The track_rq function tracks all requires that are made after its call
  !! and stores the corresponding filename along with the given module name.
  !! These two lists can be retrieved with this routine.
  !!
  subroutine tem_get_required_Lua(L, fileList, labelList)
    ! ---------------------------------------------------------------------------
    !> Lua state to get the requires from
    type(flu_State) :: L
    !> List of required filenames in the Lua state.
    character(len=PathLen), allocatable, intent(out) :: fileList(:)
    !> List of the module names associated to the corresponding files.
    character(len=LabelLen), allocatable, intent(out) :: labelList(:)
    ! ---------------------------------------------------------------------------
    integer :: track_handle
    integer :: of_handle
    integer :: list_handle
    integer :: nFiles
    integer :: iError
    integer :: i
    ! ---------------------------------------------------------------------------

    call aot_table_open( L       = L,                                          &
      &                  thandle = track_handle,                               &
      &                  key     = "track_rq" )
    call aot_table_open( L       = L,                                          &
      &                  parent  = track_handle,                               &
      &                  thandle = of_handle,                                  &
      &                  key     = "open_files" )
    call aot_get_val( L       = L,                                             &
      &               val     = nFiles,                                        &
      &               ErrCode = iError,                                        &
      &               thandle = of_handle,                                     &
      &               key     = "count" )
    allocate(fileList(nFiles))
    allocate(labelList(nFiles))

    call aot_table_open( L       = L,                                          &
      &                  parent  = of_handle,                                  &
      &                  thandle = list_handle,                                &
      &                  key     = "name" )
    do i=1,nFiles
      call aot_get_val( L       = L,                                           &
        &               val     = labelList(i),                                &
        &               ErrCode = iError,                                      &
        &               thandle = list_handle,                                 &
        &               pos     = i )
    end do
    call aot_table_close(L, list_handle)

    call aot_table_open( L       = L,                                          &
      &                  parent  = of_handle,                                  &
      &                  thandle = list_handle,                                &
      &                  key     = "file" )
    do i=1,nFiles
      call aot_get_val( L       = L,                                           &
        &               val     = fileList(i),                                 &
        &               ErrCode = iError,                                      &
        &               thandle = list_handle,                                 &
        &               pos     = i )
    end do
    call aot_table_close(L, list_handle)
    call aot_table_close(L, of_handle)
    call aot_table_close(L, track_handle)

  end subroutine tem_get_required_Lua
! ****************************************************************************** !


end module tem_lua_requires_module
! ****************************************************************************** !
