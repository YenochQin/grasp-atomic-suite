! Load optional nuclear-polarization data and publish the combined potential after writing its summary.
module nuclear_polarization_initialization
    ! Module dependencies
    use grasp_parameters, only: max_radial_points
    use nuclear_polarization_state, only: npol_capture, npol_config, npol_publish, npol_reset, coulomb_potential_baseline
    use nuclear_polarization_input, only: npol_read_config, npol_validate_config
    use nuclear_polarization_potential, only: npol_check_grid, npol_build_grid, npol_compose
    implicit none (type, external)

    ! Module interface
    private
    public :: npol_init, npol_write_summary, npol_summary_path
contains
    ! Supply the freshly computed Coulomb grid so polarization is composed exactly once.
    subroutine npol_init(config_path, nuclear_charge, grid_point_count, radial_grid, negative_radius_times_potential, &
        & status, message)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64

        ! Arguments
        character(*), intent(in) :: config_path
        real(real64), intent(in) :: nuclear_charge, radial_grid(:)
        integer, intent(in) :: grid_point_count
        real(real64), intent(inout) :: negative_radius_times_potential(:)
        integer, intent(out) :: status
        character(*), intent(out) :: message

        ! Local declarations
        type(npol_config) :: config
        real(real64) :: dipole_potential(max_radial_points), quadrupole_potential(max_radial_points), &
            & polarization_potential(max_radial_points), composed_potential(max_radial_points)
        call npol_reset()
        call npol_check_grid(grid_point_count, radial_grid, status, message)
        if (status /= 0) return
        status = 1
        message = "npol: Coulomb array shorter than N"
        if (size(negative_radius_times_potential) < grid_point_count) return
        call npol_capture(grid_point_count, negative_radius_times_potential(:grid_point_count))
        call npol_read_config(config_path, config, status, message)
        if (status /= 0) return
        call npol_validate_config(config, nuclear_charge, status, message)
        if (status /= 0) return
        call npol_build_grid(config, grid_point_count, radial_grid, dipole_potential, quadrupole_potential, &
            & polarization_potential, status, message)
        if (status /= 0) return
        call npol_compose(grid_point_count, radial_grid, coulomb_potential_baseline, polarization_potential, &
            composed_potential, status, message)
        if (status /= 0) return
        call npol_write_summary(npol_summary_path(config_path), config, grid_point_count, status, message)
        if (status /= 0) return
        negative_radius_times_potential(:grid_point_count) = composed_potential(:grid_point_count)
        call npol_publish(config, grid_point_count, dipole_potential(:grid_point_count), &
            & quadrupole_potential(:grid_point_count), polarization_potential(:grid_point_count))
    end subroutine npol_init

    ! Place the summary beside the input file.
    function npol_summary_path(path) result(summary)
        ! Arguments
        character(*), intent(in) :: path
        character(len=scan(path, "/"//achar(92), back=.true.)+len("npol.summary")) :: summary

        ! Local declarations
        integer :: last_directory_separator
        last_directory_separator = scan(path, "/"//achar(92), back=.true.)
        summary = path(:last_directory_separator)//"npol.summary"
    end function npol_summary_path

    ! Verify the written stream by reading back its length and contents before reporting success.
    subroutine npol_write_summary(path, config, grid_point_count, status, message)
        ! Arguments
        character(*), intent(in) :: path
        type(npol_config), intent(in) :: config
        integer, intent(in) :: grid_point_count
        integer, intent(out) :: status
        character(*), intent(out) :: message

        ! Local declarations
        integer :: summary_unit, io_status, close_status, record_index, file_size_bytes
        character(512) :: io_message
        character(1200) :: summary_records(15)
        character(18000) :: expected_contents, read_back_contents
        integer :: contents_length, record_length
        status = 0
        message = ""
        if (.not.config%enabled) return
        write(summary_records, "(A,I0,/,3(A,L1,/),2(A,I0,/),3(A,A,/),4(A,ES25.17,/),A,I0,/,A)", &
            iostat=io_status, iomsg=io_message) &
            "schema_version = ", config%schema_version, "enabled = ", config%enabled, &
            "use_e1 = ", config%use_e1, "use_e2 = ", config%use_e2, &
            "isotope_z = ", config%isotope_z, "isotope_a = ", config%isotope_a, &
            "units = ", trim(config%units), "model_id = ", trim(config%model_id), &
            "parameter_source = ", trim(config%parameter_source), &
            "k1 [Eh a0^4] = ", config%k1, "b1 [a0] = ", config%b1, &
            "k2 [Eh a0^6] = ", config%k2, "b2 [a0] = ", config%b2, "N = ", grid_point_count, "composed = T"
        if (io_status /= 0) then
            status = 1
            message = "npol: cannot format summary: "//trim(io_message)
            return
        end if
        contents_length = 0
        do record_index = 1, size(summary_records)
            record_length = len_trim(summary_records(record_index))
            expected_contents(contents_length+1:contents_length+record_length+1) = &
                & trim(summary_records(record_index))//new_line("a")
            contents_length = contents_length+record_length+1
        end do
        open(newunit=summary_unit, file=path, status="unknown", position="rewind", access="stream", &
            form="unformatted", action="write", iostat=io_status, iomsg=io_message)
        if (io_status /= 0) then
            status = 1
            message = "npol: cannot open summary: "//trim(io_message)
            return
        end if
        write(summary_unit, iostat=io_status, iomsg=io_message) expected_contents(:contents_length)
        if (io_status == 0) endfile(summary_unit, iostat=io_status, iomsg=io_message)
        if (io_status == 0) flush(summary_unit, iostat=io_status, iomsg=io_message)
        if (io_status /= 0) then
            close(summary_unit, status="delete", iostat=close_status)
        else
            close(summary_unit, iostat=io_status, iomsg=io_message)
        end if
        if (io_status /= 0) then
            status = 1
            message = "npol: cannot write/close summary: "//trim(io_message)
            return
        end if
        open(newunit=summary_unit, file=path, status="old", access="stream", form="unformatted", &
            action="read", iostat=io_status, iomsg=io_message)
        if (io_status == 0) then
            file_size_bytes = 0
            inquire(unit=summary_unit, size=file_size_bytes, iostat=io_status)
            if (file_size_bytes /= contents_length) io_status = 1
            if (io_status == 0) read(summary_unit, iostat=io_status, iomsg=io_message) read_back_contents(:contents_length)
            if (io_status == 0) then
                if (read_back_contents(:contents_length) /= expected_contents(:contents_length)) io_status = 1
            end if
            if (io_status /= 0) then
                close(summary_unit, status="delete", iostat=close_status)
            else
                close(summary_unit, iostat=io_status, iomsg=io_message)
            end if
        end if
        if (io_status /= 0) then
            status = 1
            message = "npol: cannot write/close summary: read-back verification failed"
        end if
    end subroutine npol_write_summary
end module nuclear_polarization_initialization
