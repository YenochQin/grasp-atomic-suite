! Save the optimized radial orbitals for subsequent CI and property calculations.
module write_scf_radial_wavefunctions_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_scf_radial_wavefunctions
contains

    subroutine write_scf_radial_wavefunctions(output_wavefunction_path)
        ! Module dependencies
        use radial_wavefunction_format, only: write_radial_orbitals
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: wavefunction_output_unit
        character , intent(in) :: output_wavefunction_path*(*)

        open(newunit=wavefunction_output_unit, file=output_wavefunction_path, status="replace", &
            form="unformatted", access="stream", action="write")
        call write_radial_orbitals(wavefunction_output_unit)
        close(wavefunction_output_unit)

        return
    end subroutine write_scf_radial_wavefunctions
end module write_scf_radial_wavefunctions_interface
