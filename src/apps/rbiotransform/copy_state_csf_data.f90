module copy_state_csf_data_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: copy_state_csf_data
contains

    !***********************************************************************
    !                                                                      *
    subroutine copy_state_csf_data(n)
        ! Module dependencies
        !                                                                      *
        !   This subroutine transfers data to the initial and final state      *
        !   common blocks                                                      *
        !                                                                      *
        !   Written by Per Jonsson                Last revision: June   1996   *
        !                                                                      *
        !***********************************************************************
        use biorthogonal_orbital_sets, only: final_orbital_symmetry_labels, initial_orbital_symmetry_labels, final_csf_count, &
            final_orbital_count, initial_csf_count, initial_orbital_count, final_orbital_principal_n, final_orbital_kappa, &
            initial_orbital_principal_n, initial_orbital_kappa, final_orbital_l, final_orbital_two_j, initial_orbital_l, &
            initial_orbital_two_j
        use orbital_csf_basis, only: csf_count, orbital_count, orbital_kappa, orbital_l, orbital_principal_n, &
            orbital_symmetry_labels, orbital_two_j
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: n

        ! Local declarations
        !-----------------------------------------------
        !
        !  Initial state commons
        !
        !
        !  Final state commons
        !
        !
        if (n == 1) then
            initial_csf_count = csf_count
            initial_orbital_count = orbital_count
            initial_orbital_principal_n(:orbital_count) = orbital_principal_n(:orbital_count)
            initial_orbital_kappa(:orbital_count) = orbital_kappa(:orbital_count)
            initial_orbital_l(:orbital_count) = orbital_l(:orbital_count)
            initial_orbital_two_j(:orbital_count) = orbital_two_j(:orbital_count)
            initial_orbital_symmetry_labels(:orbital_count) = orbital_symmetry_labels(:orbital_count)
        else
            final_csf_count = csf_count
            final_orbital_count = orbital_count
            final_orbital_principal_n(:orbital_count) = orbital_principal_n(:orbital_count)
            final_orbital_kappa(:orbital_count) = orbital_kappa(:orbital_count)
            final_orbital_l(:orbital_count) = orbital_l(:orbital_count)
            final_orbital_two_j(:orbital_count) = orbital_two_j(:orbital_count)
            final_orbital_symmetry_labels(:orbital_count) = orbital_symmetry_labels(:orbital_count)
        end if

        return
    end subroutine copy_state_csf_data
end module copy_state_csf_data_interface
