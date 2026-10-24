module npol_input
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use vast_kind_param, only: double
    use npol_C, only: npol_config, NPOL_MODEL
    implicit none
    private
    public :: npol_read_config, npol_validate_config
contains
    subroutine npol_read_config(path, cfg, status, message)
        character(*), intent(in) :: path
        type(npol_config), intent(out) :: cfg
        integer, intent(out) :: status
        character(*), intent(out) :: message
        integer :: schema_version, isotope_z, isotope_a, unit, ios, bytes
        logical :: enabled, use_e1, use_e2, exists
        real(double) :: k1, b1, k2, b2
        ! Match the maximum accepted file size so namelist strings cannot be
        ! truncated before length validation. Fixed scratch buffers also avoid
        ! saved allocatable descriptors under the suite's -fno-automatic flag.
        character(65536) :: units, model_id, parameter_source, text
        character(512) :: detail
        namelist /npol/ schema_version, enabled, use_e1, use_e2, units, k1, b1, k2, b2, &
            isotope_z, isotope_a, model_id, parameter_source

        cfg = npol_config()
        status = 1
        message = 'npol: cannot inspect input file'
        if (len_trim(path) == 0 .or. index(path, achar(0)) /= 0) then
            message = 'npol: invalid empty or NUL-containing input path'
            return
        end if
        inquire(file=path, exist=exists, iostat=ios)
        if (ios /= 0) return
        if (.not.exists) then
            status = 0
            message = ''
            return
        end if
        open(newunit=unit, file=path, status='old', access='stream', form='unformatted', &
            action='read', iostat=ios, iomsg=detail)
        if (ios /= 0) then
            message = 'npol: cannot open input: '//trim(detail)
            return
        end if
        bytes = 0
        inquire(unit=unit, size=bytes, iostat=ios)
        if (ios /= 0 .or. bytes < 1 .or. bytes > 65536) then
            message = 'npol: input size must be between 1 and 65536 bytes'
            close(unit, iostat=ios)
            return
        end if
        detail = ''
        read(unit, iostat=ios, iomsg=detail) text(:bytes)
        close(unit, iostat=status)
        if (ios /= 0 .or. status /= 0) then
            status = 1
            message = 'npol: cannot read/close input: '//trim(detail)
            return
        end if
        call check_envelope(text(:bytes), status, message)
        if (status /= 0) return
        schema_version = -1
        enabled = .false.
        use_e1 = .false.
        use_e2 = .false.
        units(:) = ''
        model_id(:) = ''
        parameter_source(:) = ''
        k1 = -1
        b1 = -1
        k2 = -1
        b2 = -1
        isotope_z = 0
        isotope_a = 0
        read(text(:bytes), nml=npol, iostat=ios, iomsg=detail)
        status = 1
        if (ios /= 0) then
            message = 'npol: invalid namelist (syntax, field or value type): '//trim(detail)
            return
        end if
        if (schema_version /= 1) then
            message = 'npol: schema_version must be explicitly set to 1'
            return
        end if
        if (len_trim(units) > len(cfg%units)) then
            message = 'npol: units exceeds 16 characters'
            return
        end if
        if (len_trim(model_id) > len(cfg%model_id)) then
            message = 'npol: model_id exceeds 64 characters'
            return
        end if
        if (len_trim(parameter_source) > len(cfg%parameter_source)) then
            message = 'npol: parameter_source exceeds 1024 characters'
            return
        end if
        do ios = 1, len_trim(parameter_source)
            if (iachar(parameter_source(ios:ios)) < 32 .or. iachar(parameter_source(ios:ios)) == 127) then
                message = 'npol: parameter_source must not contain control characters'
                return
            end if
        end do
        cfg%schema_version = schema_version
        cfg%enabled = enabled
        cfg%use_e1 = use_e1
        cfg%use_e2 = use_e2
        cfg%units = units
        cfg%model_id = model_id
        cfg%parameter_source = parameter_source
        cfg%k1 = k1
        cfg%b1 = b1
        cfg%k2 = k2
        cfg%b2 = b2
        cfg%isotope_z = isotope_z
        cfg%isotope_a = isotope_a
        status = 0
        message = ''
    end subroutine

    subroutine check_envelope(text, status, message)
        character(*), intent(in) :: text
        integer, intent(out) :: status
        character(*), intent(out) :: message
        character(len(text)) :: clean
        character :: quote, c
        integer :: i, last, finish, code
        ! The intrinsic namelist reader searches for its group and ignores
        ! trailing data. Reject that otherwise silent acceptance of bad files.
        clean = text
        quote = ' '
        i = 1
        last = len(text)
        finish = 0
        do while (i <= last)
            c = text(i:i)
            if (quote /= ' ') then
                if (c == quote) then
                    if (i < last) then
                        if (text(i+1:i+1) == quote) then
                            i = i+2
                            cycle
                        end if
                    end if
                    quote = ' '
                end if
            else if (c == "'" .or. c == '"') then
                quote = c
            else if (c == '!') then
                do while (i <= last)
                    if (text(i:i) == new_line('a')) exit
                    clean(i:i) = ' '
                    i = i+1
                end do
                cycle
            else
                code = iachar(c)
                if (code == 9 .or. code == 10 .or. code == 13) clean(i:i) = ' '
                if (code >= iachar('A') .and. code <= iachar('Z')) clean(i:i) = achar(code+32)
                if (c == '/' .and. finish == 0) finish = i
            end if
            i = i+1
        end do
        status = 1
        message = 'npol: invalid namelist; expected exactly one &npol ... / group'
        if (quote /= ' ' .or. finish == 0) return
        if (len_trim(clean(finish+1:)) /= 0) return
        clean = adjustl(clean)
        if (len_trim(clean) < 7) return
        if (clean(:6) /= '&npol ') return
        status = 0
        message = ''
    end subroutine

    subroutine npol_validate_config(cfg, z, status, message)
        type(npol_config), intent(in) :: cfg
        real(double), intent(in) :: z
        integer, intent(out) :: status
        character(*), intent(out) :: message
        status = 1
        message = 'npol: schema_version must be 1'
        if (cfg%schema_version /= 1) return
        if (.not.cfg%enabled) then
            status = 0
            message = ''
            return
        end if
        message = 'npol: enabled requires use_e1 or use_e2'
        if (.not.(cfg%use_e1 .or. cfg%use_e2)) return
        message = 'npol: units must be atomic'
        if (cfg%units /= 'atomic') return
        message = 'npol: model_id must be '//NPOL_MODEL
        if (cfg%model_id /= NPOL_MODEL) return
        message = 'npol: parameter_source is required'
        if (len_trim(cfg%parameter_source) == 0) return
        message = 'npol: isotope_z must be positive and match the running nuclear charge'
        if (.not.ieee_is_finite(z)) return
        if (cfg%isotope_z <= 0 .or. real(cfg%isotope_z,double) /= z) return
        message = 'npol: isotope_a must be supplied and at least isotope_z'
        if (cfg%isotope_a < cfg%isotope_z) return
        if (cfg%use_e1) then
            message = 'npol: k1 is required, finite and nonnegative'
            if (.not.ieee_is_finite(cfg%k1)) return
            if (cfg%k1 < 0) return
            message = 'npol: b1 is required, finite and strictly positive'
            if (.not.ieee_is_finite(cfg%b1)) return
            if (cfg%b1 <= 0) return
        end if
        if (cfg%use_e2) then
            message = 'npol: k2 is required, finite and nonnegative'
            if (.not.ieee_is_finite(cfg%k2)) return
            if (cfg%k2 < 0) return
            message = 'npol: b2 is required, finite and strictly positive'
            if (.not.ieee_is_finite(cfg%b2)) return
            if (cfg%b2 <= 0) return
        end if
        status = 0
        message = ''
    end subroutine
end module npol_input
