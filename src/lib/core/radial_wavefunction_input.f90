! Check the radial wavefunction file signature, load its orbital records, and close the file.
module radial_wavefunction_input
    implicit none (type, external)

    ! Module interface
    private
    public :: read_radial_wavefunction_file
contains
    subroutine read_radial_wavefunction_file(radial_file_path)
        ! Module dependencies
        use standard_io_units, only: stderr_unit
        use radial_wavefunction_records, only: read_radial_wavefunction_records
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: wavefunction_input_unit
        character (len = *), intent(in) :: radial_file_path
        integer :: load_status, read_status
        character :: file_signature*6

        open(newunit=wavefunction_input_unit, file=radial_file_path, form="UNFORMATTED", &
            status="OLD", action="READ", iostat=load_status)
        if (load_status /= 0) error stop "Cannot open radial wavefunction file: "//trim(radial_file_path)
        read (wavefunction_input_unit, iostat=read_status) file_signature
        if (read_status/=0 .or. file_signature/="G92RWF") then
            write (stderr_unit, *) "This is not a Radial WaveFunction File;"
            close(wavefunction_input_unit)
            stop
        end if
        call read_radial_wavefunction_records (wavefunction_input_unit, load_status)

        if (load_status /= 0) then
            write (stderr_unit, *) "Radial wavefunctions defined in CSL file", &
                " not found in Radial WaveFunction File"
            close(wavefunction_input_unit)
            stop
        end if

        close(wavefunction_input_unit)

        return
    end subroutine read_radial_wavefunction_file
end module radial_wavefunction_input
