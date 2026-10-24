module npol_initialization
    use vast_kind_param, only: double
    use parameter_def, only: NNNP
    use npol_C
    use npol_input, only: npol_read_config, npol_validate_config
    use npol_kernel, only: npol_check_grid, npol_build_grid, npol_compose
    implicit none
    private
    public :: npol_init, npol_write_summary, npol_summary_path
contains
    ! Caller supplies the freshly initialized Coulomb ZZ, never a prior total.
    subroutine npol_init(path, z, n, r, zz, status, message)
        character(*), intent(in) :: path
        real(double), intent(in) :: z, r(:)
        integer, intent(in) :: n
        real(double), intent(inout) :: zz(:)
        integer, intent(out) :: status
        character(*), intent(out) :: message
        type(npol_config) :: cfg
        real(double) :: a(NNNP), b(NNNP), v(NNNP), candidate(NNNP)
        call npol_reset()
        call npol_check_grid(n, r, status, message)
        if (status /= 0) return
        status = 1
        message = 'npol: Coulomb array shorter than N'
        if (size(zz) < n) return
        call npol_capture(n, zz(:n))
        call npol_read_config(path, cfg, status, message)
        if (status /= 0) return
        call npol_validate_config(cfg, z, status, message)
        if (status /= 0) return
        call npol_build_grid(cfg, n, r, a, b, v, status, message)
        if (status /= 0) return
        call npol_compose(n, r, zz_coul, v, candidate, status, message)
        if (status /= 0) return
        call npol_write_summary(npol_summary_path(path), cfg, n, status, message)
        if (status /= 0) return
        zz(:n) = candidate(:n)
        call npol_publish(cfg, n, a(:n), b(:n), v(:n))
    end subroutine

    function npol_summary_path(path) result(summary)
        character(*), intent(in) :: path
        character(:), allocatable :: summary
        integer :: slash
        slash = scan(path, '/'//achar(92), back=.true.)
        summary = path(:slash)//'npol.summary'
    end function

    subroutine npol_write_summary(path, cfg, n, status, message)
        character(*), intent(in) :: path
        type(npol_config), intent(in) :: cfg
        integer, intent(in) :: n
        integer, intent(out) :: status
        character(*), intent(out) :: message
        integer :: unit, ios, close_status, i, bytes
        character(512) :: detail
        character(1200) :: records(15)
        character(18000) :: expected, observed
        integer :: used, width
        status = 0
        message = ''
        if (.not.cfg%enabled) return
        ! Format in memory first, then write a bounded stream and verify it.
        ! Some Fortran runtimes lose device/flush errors; IOSTAT alone must
        ! not turn failed output into a successful initialization.
        write(records, '(A,I0,/,3(A,L1,/),2(A,I0,/),3(A,A,/),4(A,ES25.17,/),A,I0,/,A)', &
            iostat=ios, iomsg=detail) &
            'schema_version = ', cfg%schema_version, 'enabled = ', cfg%enabled, &
            'use_e1 = ', cfg%use_e1, 'use_e2 = ', cfg%use_e2, &
            'isotope_z = ', cfg%isotope_z, 'isotope_a = ', cfg%isotope_a, &
            'units = ', trim(cfg%units), 'model_id = ', trim(cfg%model_id), &
            'parameter_source = ', trim(cfg%parameter_source), &
            'k1 [Eh a0^4] = ', cfg%k1, 'b1 [a0] = ', cfg%b1, &
            'k2 [Eh a0^6] = ', cfg%k2, 'b2 [a0] = ', cfg%b2, 'N = ', n, 'composed = T'
        if (ios /= 0) then
            status = 1
            message = 'npol: cannot format summary: '//trim(detail)
            return
        end if
        used = 0
        do i = 1, size(records)
            width = len_trim(records(i))
            expected(used+1:used+width+1) = trim(records(i))//new_line('a')
            used = used+width+1
        end do
        open(newunit=unit, file=path, status='unknown', position='rewind', access='stream', &
            form='unformatted', action='write', iostat=ios, iomsg=detail)
        if (ios /= 0) then
            status = 1
            message = 'npol: cannot open summary: '//trim(detail)
            return
        end if
        write(unit, iostat=ios, iomsg=detail) expected(:used)
        if (ios == 0) endfile(unit, iostat=ios, iomsg=detail)
        if (ios == 0) flush(unit, iostat=ios, iomsg=detail)
        if (ios /= 0) then
            close(unit, status='delete', iostat=close_status)
        else
            close(unit, iostat=ios, iomsg=detail)
        end if
        if (ios /= 0) then
            status = 1
            message = 'npol: cannot write/close summary: '//trim(detail)
            return
        end if
        open(newunit=unit, file=path, status='old', access='stream', form='unformatted', &
            action='read', iostat=ios, iomsg=detail)
        if (ios == 0) then
            bytes = 0
            inquire(unit=unit, size=bytes, iostat=ios)
            if (bytes /= used) ios = 1
            if (ios == 0) read(unit, iostat=ios, iomsg=detail) observed(:used)
            if (ios == 0) then
                if (observed(:used) /= expected(:used)) ios = 1
            end if
            if (ios /= 0) then
                close(unit, status='delete', iostat=close_status)
            else
                close(unit, iostat=ios, iomsg=detail)
            end if
        end if
        if (ios /= 0) then
            status = 1
            message = 'npol: cannot write/close summary: read-back verification failed'
        end if
    end subroutine
end module npol_initialization
