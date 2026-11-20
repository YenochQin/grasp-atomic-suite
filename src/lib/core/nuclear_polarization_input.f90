! Read and validate the optional nuclear-polarization namelist for the selected isotope.
module nuclear_polarization_input
    ! Module dependencies
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use nuclear_polarization_state, only: npol_config, npol_model
    implicit none (type, external)

    ! Module interface
    private
    public :: npol_read_config, npol_validate_config
contains
    subroutine npol_read_config(config_path, config, status, message)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64

        ! Arguments
        character(*), intent(in) :: config_path
        type(npol_config), intent(out) :: config
        integer, intent(out) :: status
        character(*), intent(out) :: message

        ! Local declarations
        integer :: schema_version, isotope_z, isotope_a, input_unit, io_status, input_size_bytes
        integer :: source_index
        logical :: enabled, use_e1, use_e2, file_exists
        real(real64) :: k1, b1, k2, b2
        ! The string buffers cover the input size limit so validation sees the full supplied values.
        character(65536) :: units, model_id, parameter_source, input_contents
        character(512) :: io_message
        namelist /npol/ schema_version, enabled, use_e1, use_e2, units, k1, b1, k2, b2, &
            isotope_z, isotope_a, model_id, parameter_source

        config = npol_config()
        status = 1
        message = "npol: cannot inspect input file"
        if (len_trim(config_path) == 0 .or. index(config_path, achar(0)) /= 0) then
            message = "npol: invalid empty or NUL-containing input path"
            return
        end if
        inquire(file=config_path, exist=file_exists, iostat=io_status)
        if (io_status /= 0) return
        if (.not.file_exists) then
            status = 0
            message = ""
            return
        end if
        open(newunit=input_unit, file=config_path, status="old", access="stream", form="unformatted", &
            action="read", iostat=io_status, iomsg=io_message)
        if (io_status /= 0) then
            message = "npol: cannot open input: "//trim(io_message)
            return
        end if
        input_size_bytes = 0
        inquire(unit=input_unit, size=input_size_bytes, iostat=io_status)
        if (io_status /= 0 .or. input_size_bytes < 1 .or. input_size_bytes > 65536) then
            message = "npol: input size must be between 1 and 65536 bytes"
            close(input_unit, iostat=io_status)
            return
        end if
        io_message = ""
        read(input_unit, iostat=io_status, iomsg=io_message) input_contents(:input_size_bytes)
        close(input_unit, iostat=status)
        if (io_status /= 0 .or. status /= 0) then
            status = 1
            message = "npol: cannot read/close input: "//trim(io_message)
            return
        end if
        call check_envelope(input_contents(:input_size_bytes), status, message)
        if (status /= 0) return
        schema_version = -1
        enabled = .false.
        use_e1 = .false.
        use_e2 = .false.
        units(:) = ""
        model_id(:) = ""
        parameter_source(:) = ""
        k1 = -1
        b1 = -1
        k2 = -1
        b2 = -1
        isotope_z = 0
        isotope_a = 0
        read(input_contents(:input_size_bytes), nml=npol, iostat=io_status, iomsg=io_message)
        status = 1
        if (io_status /= 0) then
            message = "npol: invalid namelist (syntax, field or value type): "//trim(io_message)
            return
        end if
        if (schema_version /= 1) then
            message = "npol: schema_version must be explicitly set to 1"
            return
        end if
        if (len_trim(units) > len(config%units)) then
            message = "npol: units exceeds 16 characters"
            return
        end if
        if (len_trim(model_id) > len(config%model_id)) then
            message = "npol: model_id exceeds 64 characters"
            return
        end if
        if (len_trim(parameter_source) > len(config%parameter_source)) then
            message = "npol: parameter_source exceeds 1024 characters"
            return
        end if
        do source_index = 1, len_trim(parameter_source)
            if (iachar(parameter_source(source_index:source_index)) < 32 .or. &
                iachar(parameter_source(source_index:source_index)) == 127) then
                message = "npol: parameter_source must not contain control characters"
                return
            end if
        end do
        config%schema_version = schema_version
        config%enabled = enabled
        config%use_e1 = use_e1
        config%use_e2 = use_e2
        config%units = units
        config%model_id = model_id
        config%parameter_source = parameter_source
        config%k1 = k1
        config%b1 = b1
        config%k2 = k2
        config%b2 = b2
        config%isotope_z = isotope_z
        config%isotope_a = isotope_a
        status = 0
        message = ""
    end subroutine npol_read_config

    ! Require exactly one namelist group, accounting for quotes and comments.
    subroutine check_envelope(input_contents, status, message)
        ! Arguments
        character(*), intent(in) :: input_contents
        integer, intent(out) :: status
        character(*), intent(out) :: message

        ! Local declarations
        character(len(input_contents)) :: normalized_contents
        character :: open_quote, current_character
        integer :: character_index, contents_length, group_end_index, character_code
        normalized_contents = input_contents
        open_quote = " "
        character_index = 1
        contents_length = len(input_contents)
        group_end_index = 0
        do while (character_index <= contents_length)
            current_character = input_contents(character_index:character_index)
            if (open_quote /= " ") then
                if (current_character == open_quote) then
                    if (character_index < contents_length) then
                        if (input_contents(character_index+1:character_index+1) == open_quote) then
                            character_index = character_index+2
                            cycle
                        end if
                    end if
                    open_quote = " "
                end if
            else if (current_character == "'" .or. current_character == '"') then
                open_quote = current_character
            else if (current_character == "!") then
                do while (character_index <= contents_length)
                    if (input_contents(character_index:character_index) == new_line("a")) exit
                    normalized_contents(character_index:character_index) = " "
                    character_index = character_index+1
                end do
                cycle
            else
                character_code = iachar(current_character)
                if (character_code == 9 .or. character_code == 10 .or. character_code == 13) &
                    & normalized_contents(character_index:character_index) = " "
                if (character_code >= iachar("A") .and. character_code <= iachar("Z")) &
                    & normalized_contents(character_index:character_index) = achar(character_code+32)
                if (current_character == "/" .and. group_end_index == 0) group_end_index = character_index
            end if
            character_index = character_index+1
        end do
        status = 1
        message = "npol: invalid namelist; expected exactly one &npol ... / group"
        if (open_quote /= " " .or. group_end_index == 0) return
        if (len_trim(normalized_contents(group_end_index+1:)) /= 0) return
        normalized_contents = adjustl(normalized_contents)
        if (len_trim(normalized_contents) < 7) return
        if (normalized_contents(:6) /= "&npol ") return
        status = 0
        message = ""
    end subroutine check_envelope

    ! Check active channel parameters, units, and isotope identity.
    subroutine npol_validate_config(config, nuclear_charge, status, message)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64

        ! Arguments
        type(npol_config), intent(in) :: config
        real(real64), intent(in) :: nuclear_charge
        integer, intent(out) :: status
        character(*), intent(out) :: message
        status = 1
        message = "npol: schema_version must be 1"
        if (config%schema_version /= 1) return
        if (.not.config%enabled) then
            status = 0
            message = ""
            return
        end if
        message = "npol: enabled requires use_e1 or use_e2"
        if (.not.(config%use_e1 .or. config%use_e2)) return
        message = "npol: units must be atomic"
        if (config%units /= "atomic") return
        message = "npol: model_id must be "//npol_model
        if (config%model_id /= npol_model) return
        message = "npol: parameter_source is required"
        if (len_trim(config%parameter_source) == 0) return
        message = "npol: isotope_z must be positive and match the running nuclear charge"
        if (.not.ieee_is_finite(nuclear_charge)) return
        if (config%isotope_z <= 0 .or. real(config%isotope_z,real64) /= nuclear_charge) return
        message = "npol: isotope_a must be supplied and at least isotope_z"
        if (config%isotope_a < config%isotope_z) return
        if (config%use_e1) then
            message = "npol: k1 is required, finite and nonnegative"
            if (.not.ieee_is_finite(config%k1)) return
            if (config%k1 < 0) return
            message = "npol: b1 is required, finite and strictly positive"
            if (.not.ieee_is_finite(config%b1)) return
            if (config%b1 <= 0) return
        end if
        if (config%use_e2) then
            message = "npol: k2 is required, finite and nonnegative"
            if (.not.ieee_is_finite(config%k2)) return
            if (config%k2 < 0) return
            message = "npol: b2 is required, finite and strictly positive"
            if (.not.ieee_is_finite(config%b2)) return
            if (config%b2 <= 0) return
        end if
        status = 0
        message = ""
    end subroutine npol_validate_config
end module nuclear_polarization_input
