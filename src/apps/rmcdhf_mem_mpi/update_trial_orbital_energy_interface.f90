module update_trial_orbital_energy_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine update_trial_orbital_energy (orbital_index, first_oscillation_sign, effective_principal_n, &
            & node_count_difference, energy_correction, &
            & solve_failed, orbital_sign_reversal_count)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: orbital_index
            real(real64), intent(in) :: first_oscillation_sign
            integer, intent(in) :: effective_principal_n
            integer, intent(in) :: node_count_difference
            real(real64), intent(inout) :: energy_correction
            logical, intent(out) :: solve_failed
            integer, intent(out) :: orbital_sign_reversal_count
        end subroutine update_trial_orbital_energy
    end interface
end module update_trial_orbital_energy_interface
