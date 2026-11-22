! Load the final-state descriptor basis and retain its original orbital and block metadata.
module load_final_transition_basis_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: load_final_transition_basis
contains

    ! Written by Farid A. Parpia. Last revision: 23 Dec 1992
    !
    ! Translated by Charlotte Froese Fischer
    !                        Gediminas Gaigalas  10/05/17
    subroutine load_final_transition_basis(final_core_subshell_count, state_name)
        !  Module dependencies
        use atomic_calculation_state, only: electron_count, final_electron_count
        use orbital_csf_basis, only: csf_count, orbital_count, orbital_kappa, orbital_l, orbital_principal_n, &
            orbital_symmetry_labels, orbital_two_j
        use biorthogonal_orbital_sets, only: final_orbital_symmetry_labels, final_csf_count, final_orbital_count, &
            final_orbital_principal_n, final_orbital_kappa, final_orbital_l, final_orbital_two_j

        use mpi_csf_basis_input, only: load_csf_basis_mpi
        use transition_basis_blocks, only: transition_block_count, cumulative_csf_counts, final_transition_block_count, &
            final_cumulative_csf_counts
        implicit none (type, external)

        ! Arguments
        integer , intent(out) :: final_core_subshell_count
        character , intent(in) :: state_name*128

        ! Local declarations
        integer :: io_status, ncorer
        !
        ! Load the descriptor pair for this state.
        call load_csf_basis_mpi(trim(state_name),ncorer)
        io_status=0
        if (io_status /= 0) error stop "Transition descriptor load failed"
        if (transition_block_count > size(final_cumulative_csf_counts)) error stop "Transition block capacity exceeded"
        final_transition_block_count=transition_block_count
        final_cumulative_csf_counts=0
        final_cumulative_csf_counts(:transition_block_count)=cumulative_csf_counts(:transition_block_count)
        ! Keep the closed-core count from the descriptor basis.
        final_core_subshell_count = ncorer
        final_electron_count = electron_count
        final_orbital_count = orbital_count
        final_csf_count = csf_count
        final_orbital_symmetry_labels(:orbital_count) = orbital_symmetry_labels(:orbital_count)
        final_orbital_principal_n(:orbital_count) = orbital_principal_n(:orbital_count)
        final_orbital_kappa(:orbital_count) = orbital_kappa(:orbital_count)
        final_orbital_l(:orbital_count) = orbital_l(:orbital_count)
        final_orbital_two_j(:orbital_count) = orbital_two_j(:orbital_count)
        return
    end subroutine load_final_transition_basis
end module load_final_transition_basis_interface
