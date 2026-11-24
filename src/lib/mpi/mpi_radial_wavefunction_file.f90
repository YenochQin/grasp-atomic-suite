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
        use mpi_io_errors, only: fail_input_io_mpi
        use mpi_radial_wavefunction_reader, only: read_mpi_radial_wavefunction_records
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: wavefunction_input_unit
        character (len = *) :: wavefunction_path
        integer :: header_read_status, setup_status
        character (len = 6) :: file_header
        wavefunction_input_unit = -1
        if (mpi_rank == 0) then
            open(newunit=wavefunction_input_unit, file=wavefunction_path, form="UNFORMATTED", &
                status="OLD", action="READ", iostat=setup_status)
            if (setup_status /= 0) call fail_input_io_mpi(setup_status, wavefunction_input_unit, &
                "Cannot open radial wavefunction file: "//trim(wavefunction_path))

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
