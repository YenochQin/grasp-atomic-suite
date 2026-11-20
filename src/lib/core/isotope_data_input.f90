! Check the isotope-file header, load its records, and close the file.
module isotope_data_input
    implicit none (type, external)

    ! Module interface
    private
    public :: read_isotope_data
contains
    subroutine read_isotope_data(isotope_file_path)
        ! Module dependencies
        use standard_io_units, only: stderr_unit
        use file_open, only: open_calculation_file
        use isotope_data_records, only: read_isotope_records
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: isotope_input_unit
        character (len = *), intent(in) :: isotope_file_path
        character(len=3), parameter  :: requested_status = "OLD"
        character(len=6), parameter  :: routine_name = "SETISO"
        character(len=9), parameter  :: record_format = "FORMATTED"
        character(len=14), parameter :: expected_file_header = "Atomic number:"
        integer           :: path_length, open_status, read_status
        logical           :: file_exists
        character(len=14) :: file_header
        inquire(file=isotope_file_path, exist=file_exists)
        if (.not.file_exists) then
            path_length = len_trim(isotope_file_path)
            write (stderr_unit,*)                                              &
                routine_name,"- file: ",isotope_file_path(1:path_length)," does not exist"
            stop
        end if
        call open_calculation_file (isotope_input_unit, isotope_file_path, record_format, requested_status, open_status)
        if (open_status /= 0) then
            write (stderr_unit, *) "Error opening isodata file: ", isotope_file_path(1:path_length)
            stop
        end if
        read (isotope_input_unit, "(A)", iostat=read_status) file_header
        if (read_status/=0 .or. file_header/=expected_file_header) then
            write (stderr_unit, *) "Not an ISOtope Data File;"
            close(isotope_input_unit)
            stop
        end if
        call read_isotope_records(isotope_input_unit)
        close(isotope_input_unit)

        return
    end subroutine read_isotope_data
end module isotope_data_input
