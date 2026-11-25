! Write the initial radial orbital estimates in the Suite wavefunction format.
module write_initial_radial_wavefunctions_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_initial_radial_wavefunctions
contains

    !   Open, write a header and all subshell radial wavefunctions, and
    !   close the  .rwf  file.
    !
    ! Written by Farid A. Parpia. Last revision: 18 Dec 1992
    subroutine write_initial_radial_wavefunctions(file_path)
        ! Module dependencies
        use array_memory_management, only: dalloc
        use radial_wavefunction_format, only: write_radial_orbitals
        use radial_wavefunctions, only: radial_large_components, radial_small_components
        implicit none (type, external)

        ! Arguments
        character(len=*), intent(in) :: file_path

        ! Local declarations
        integer :: wavefunction_output_unit, open_status

        open(newunit=wavefunction_output_unit, file=file_path, status="replace", form="unformatted", &
            access="stream", action="write", iostat=open_status)
        if (open_status /= 0) error stop "Cannot create initial radial wavefunction: "//trim(file_path)
        call write_radial_orbitals(wavefunction_output_unit)

        close(wavefunction_output_unit)
        ! Deallocate the storage for the radial wavefunctions
        call dalloc (radial_large_components, "radial_large_components", "write_initial_radial_wavefunctions")
        call dalloc (radial_small_components, "radial_small_components", "write_initial_radial_wavefunctions")
        return
    end subroutine write_initial_radial_wavefunctions
end module write_initial_radial_wavefunctions_interface
