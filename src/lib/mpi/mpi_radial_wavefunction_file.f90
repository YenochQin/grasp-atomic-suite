module mpi_radial_wavefunction_file
    implicit none (type, external)

    ! Module interface
    private
    public :: read_mpi_radial_wavefunction_file
contains
    ! Rank 0 opens the unformatted radial-wavefunction file.
    ! All ranks call the reader to receive the orbital records.
    subroutine read_mpi_radial_wavefunction_file (wavefunction_path)

        ! Module dependencies
        use radial_wavefunction_format, only: open_wavefunction_input
        use mpi_state, only: mpi_rank
        use standard_io_units, only: stderr_unit
        use mpi_io_errors, only: fail_input_io_mpi
        use mpi_radial_wavefunction_reader, only: read_mpi_radial_wavefunction_records
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: wavefunction_input_unit
        character (len = *) :: wavefunction_path
        integer :: setup_status
        wavefunction_input_unit = -1
        if (mpi_rank == 0) then
            call open_wavefunction_input(wavefunction_input_unit, wavefunction_path, setup_status)
            if (setup_status /= 0) call fail_input_io_mpi(setup_status, wavefunction_input_unit, &
                "Cannot open radial wavefunction file: "//trim(wavefunction_path))
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
