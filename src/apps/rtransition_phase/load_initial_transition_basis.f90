module load_initial_transition_basis_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: load_initial_transition_basis
contains

    !***********************************************************************
    !                                                                      *
    subroutine load_initial_transition_basis (initial_core_subshell_count,state_name)
        ! Module dependencies
        !                                                                      *
        !   Open, check, load data from and close the  .csl  file. This file   *
        !   is always attached to stream 21.                                   *
        !                                                                      *
        !   Call(s) to: [LIB92]: ALLOC, IQ, LENGTH, LODCSL, OPENFL.            *
        !                                                                      *
        !                                                                      *
        !***********************************************************************
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        use grasp_parameters,   only: max_orbitals
        use array_memory_management, only: alloc, dalloc
        use atomic_calculation_state, only: electron_count, initial_electron_count, left_electron_count
        use biorthogonal_orbital_sets, only: max_orbitals, initial_orbital_symmetry_labels, initial_csf_count, &
            initial_orbital_count, initial_orbital_principal_n, initial_orbital_kappa, initial_orbital_l, initial_orbital_two_j
        use orbital_csf_basis, only: csf_count, csf_subshell_occupations, left_occupations, max_orbitals, left_kappa, &
            left_csf_count, left_orbital_labels, left_principal_n, left_orbital_count, orbital_count, orbital_kappa, &
            orbital_l, orbital_principal_n, orbital_symmetry_labels, orbital_two_j
        use csf_coupling_data, only: csf_coupling_quantum_codes, csf_subshell_quantum_codes, left_coupling_codes, &
            left_subshell_quantum_codes
        use csf_quantum_queries, only: csf_occupation
        use csf_full_basis_reader, only: load_csf_basis
        use transition_basis_blocks, only: transition_block_count, cumulative_csf_counts, initial_transition_block_count, &
            initial_cumulative_csf_counts
        use file_open, only: open_calculation_file
        implicit none (type, external)

        ! Arguments
        integer :: initial_core_subshell_count
        character(len=24) :: state_name

        ! Local declarations
        integer :: i, j, io_status, core_subshell_count
        !
        !   Check the first record of the file; if not as expected, try again
        !
        call load_csf_basis(trim(state_name),core_subshell_count,io_status)
        if (io_status /= 0) error stop "Transition descriptor load failed"
        if (transition_block_count > size(initial_cumulative_csf_counts)) error stop "Transition block capacity exceeded"
        initial_transition_block_count=transition_block_count
        initial_cumulative_csf_counts=0
        initial_cumulative_csf_counts(:transition_block_count)=cumulative_csf_counts(:transition_block_count)
        !
        !   Check if the core should be redefined
        !
        initial_core_subshell_count = core_subshell_count
        !
        !   Copy data to reference arrays
        !
        left_electron_count = electron_count
        initial_electron_count = electron_count
        left_csf_count = csf_count
        left_orbital_count = orbital_count
        !
        call alloc (left_occupations, max_orbitals, csf_count, "IQAR", "LDCSL1" )
        call alloc (left_subshell_quantum_codes,max_orbitals,3,csf_count, "JQSAR", "LDCSL1")
        call alloc (left_coupling_codes,max_orbitals, csf_count, "JCUPAR", "LDCSL1")
        !
        do i = 1,csf_count
            do j = 1,max_orbitals
                left_occupations(j,i) = csf_subshell_occupations(j,i)
                left_subshell_quantum_codes(j,1,i) = csf_subshell_quantum_codes(j,1,i)
                left_subshell_quantum_codes(j,2,i) = csf_subshell_quantum_codes(j,2,i)
                left_subshell_quantum_codes(j,3,i) = csf_subshell_quantum_codes(j,3,i)
                left_coupling_codes(j,i) = csf_coupling_quantum_codes(j,i)
            end do
        end do
        !
        initial_orbital_count=orbital_count
        initial_csf_count=csf_count
        do i = 1,orbital_count
            initial_orbital_symmetry_labels(i) = orbital_symmetry_labels(i)
            initial_orbital_principal_n(i) = orbital_principal_n(i)
            initial_orbital_kappa(i) = orbital_kappa(i)
            initial_orbital_l(i) = orbital_l(i)
            initial_orbital_two_j(i) = orbital_two_j(i)
            left_principal_n(i) = orbital_principal_n(i)
            left_kappa(i) = orbital_kappa(i)
            left_orbital_labels(i) = orbital_symmetry_labels(i)
        end do
        !
        !   Deallocate storage so that LODCSL can restart correctly
        !
        call dalloc (csf_subshell_occupations, "IQA", "LDCSL1")
        call dalloc (csf_subshell_quantum_codes, "JQSA", "LDCSL1")
        call dalloc (csf_coupling_quantum_codes, "JCUPA", "LDCSL1")
        !


        return
    end subroutine load_initial_transition_basis
end module load_initial_transition_basis_interface
