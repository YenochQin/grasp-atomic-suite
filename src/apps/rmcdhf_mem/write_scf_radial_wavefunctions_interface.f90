module write_scf_radial_wavefunctions_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_scf_radial_wavefunctions (output_wavefunction_path)
            implicit none (type, external)

            ! Arguments
            character (len = *), intent(in) :: output_wavefunction_path
        end subroutine write_scf_radial_wavefunctions
    end interface
end module write_scf_radial_wavefunctions_interface
