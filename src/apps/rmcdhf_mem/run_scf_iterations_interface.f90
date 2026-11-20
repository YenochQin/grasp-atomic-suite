module run_scf_iterations_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine run_scf_iterations (extended_optimal_level, output_wavefunction_path)
            implicit none (type, external)

            ! Arguments
            logical, intent(in) :: extended_optimal_level
            character (len = *) :: output_wavefunction_path
        end subroutine run_scf_iterations
    end interface
end module run_scf_iterations_interface
