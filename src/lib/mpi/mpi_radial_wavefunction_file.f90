module mpi_radial_wavefunction_file
    implicit none (type, external)

    ! Module interface
    private
    public :: read_mpi_radial_wavefunction_file
contains
    ! Rank 0 opens the unformatted radial-wavefunction file and checks the G92RWF header.
    ! All ranks call the reader to receive the orbital records.
    subroutine read_mpi_radial_wavefunction_file (wavefunction_path)

        ! Module dependencies
        use mpi_state, only: mpi_rank
        use standard_io_units, only: stderr_unit
        use file_open, only: open_calculation_file
        use mpi_radial_wavefunction_reader, only: read_mpi_radial_wavefunction_records
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: wavefunction_input_unit
        character (len = *) :: wavefunction_path
        integer :: header_read_status, setup_status
        character (len = 6) :: file_header
        wavefunction_input_unit = -1
        if (mpi_rank == 0) then
            call open_calculation_file (wavefunction_input_unit, wavefunction_path, "UNFORMATTED", "OLD", setup_status)
            if (setup_status == 1) then
                write (stderr_unit,*) "Error opening", wavefunction_path(1:len_trim (wavefunction_path))
                stop
            end if

            read (wavefunction_input_unit,iostat = header_read_status) file_header
            if ((header_read_status /= 0) .or. (file_header /= "G92RWF")) then
                write (stderr_unit,*) "This is not a Radial WaveFunction File;"
                close (wavefunction_input_unit)
                stop
            end if
        end if

        call read_mpi_radial_wavefunction_records (wavefunction_input_unit, setup_status)

        if (setup_status /= 0) then
            if (mpi_rank == 0) then
                write (stderr_unit,*) "Radial wavefunctions defined in CSL file" &
                    , " not found in Radial WaveFunction File"
                close (wavefunction_input_unit)
            end if
            stop
        end if

        if (mpi_rank == 0) close (wavefunction_input_unit)

        return
    end subroutine read_mpi_radial_wavefunction_file
end module mpi_radial_wavefunction_file
