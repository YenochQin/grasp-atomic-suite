program system_paths
    use suite_system
    use sys_chdir_I
    use sys_getwd_I
    use sys_mkdir_I
    implicit none
    character(4096) :: before, after
    character(1) :: short
    character(8) :: mode
    character(*), parameter :: nested = "parent with space/quote'and;$literal"
    integer :: status, length, unit

    if (command_argument_count() > 0) then
        call get_command_argument(1,mode)
        open(newunit=unit,file='exit-flush',status='replace')
        write(unit,'(A)') 'buffered output'
        if (mode == 'zero') call terminate_process()
        call terminate_process(7)
        error stop 'exit returned'
    end if

    call sys_getwd(before,length)
    if (length /= len_trim(before)) error stop 'cwd length'
    if (run_shell('exit 7') /= 1792) error stop 'shell wait status changed'
    call sys_mkdir(nested,len(nested),status)
    if (status /= 0) error stop 'quoted mkdir failed'
    call sys_mkdir(nested,len(nested),status)
    if (status /= 0) error stop 'mkdir not idempotent'
    call sys_chdir(nested,len(nested),status)
    if (status /= 0) error stop 'chdir failed'
    call sys_getwd(after,length)
    if (trim(after) /= trim(before)//'/'//nested) error stop 'cwd changed incorrectly'
    status = working_directory(short)
    if (status == 0) error stop 'short cwd buffer accepted'
    open(newunit=unit,file='regular-file',status='new')
    close(unit)
    if (path_access('regular-file','rw') /= 0) error stop 'access existing file'
    if (path_access('missing-file',' ') == 0) error stop 'access missing file'
    call sys_mkdir('regular-file/child',len('regular-file/child'),status)
    if (status == 0) error stop 'mkdir through regular file accepted'
    call sys_chdir('missing-directory',len('missing-directory'),status)
    if (status == 0) error stop 'missing directory accepted'
    call sys_getwd(after,length)
    if (trim(after) /= trim(before)//'/'//nested) error stop 'failed chdir changed cwd'
    call sys_chdir(before,len_trim(before),status)
    if (status /= 0) error stop 'restore cwd failed'
end program system_paths
