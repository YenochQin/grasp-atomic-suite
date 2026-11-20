module read_eol_average_options_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_eol_average_options
contains

    !***********************************************************************
    !                                                                      *
    subroutine read_eol_average_options(idblk)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,    only: max_orbitals
        use array_memory_management, only: alloc
        use transition_block_indices, only: level_block_indices
        use scf_damping_factors, only: mixing_damping_factors, orbital_damping_factors
        use atomic_calculation_state, only: level_average_weights, selected_level_count, maximum_scf_iterations, &
            orbital_updates_per_scf_iteration, minimum_radial_solver_iterations, requested_level_weights
        use calculation_input_defaults, only: nondefault_options_flag
        use fixed_orbital_flags, only: is_orbital_fixed, fixed_orbital_count
        use hamiltonian_symmetry_blocks, only: eigenvector_counts_by_block, highest_selected_level_by_block, &
            symmetry_block_count
        use standard_io_units, only: stderr_unit
        use orbital_inversion_flags, only: suppress_radial_inversion
        use scf_orthogonality_control, only: enforce_orbital_orthogonality
        use orbital_csf_basis, only: max_orbitals, orbital_count
        use scf_orbital_order, only: orbital_update_order
        use correlation_orbital_flags, only: is_correlation_orbital
        use scf_iteration_state, only: orbital_self_consistency_errors, radial_solution_methods
        use block_level_selection, only: read_block_level_selection

        ! Procedure interfaces
        use read_eol_level_weights_interface, only: read_eol_level_weights
        use print_scf_results_interface, only: print_scf_results

        ! Module dependencies
        use subshell_label_selection, only: select_subshell_labels
        implicit none (type, external)

        ! Arguments
        character  :: idblk(*)*8

        ! Local declarations
        integer  :: indx(max_orbitals)
        integer :: i, nsubs, norder, loc, noffset, symmetry_block_index, level_index
        !-----------------------------------------------

        !WRITE (istde,*) 'EOL type calculation;'

        ! lodstate fills
        !    nevblk(), ncmaxblk()
        !    ncmin, iccmin(1:ncmin) -- via items (memories allocated there)

        call alloc (highest_selected_level_by_block, symmetry_block_count, "NCMAXBLK", "GETOLD")
        call alloc (eigenvector_counts_by_block, symmetry_block_count, "NEVBLK", "GETOLD" )

        !cjb  ncmaxblk & nevblk initialised in lodstate
        !     NEVBLK = 0
        !     NCMAXBLK = 0
        !cjb

        !cjb LODSTATE(NBLOCK, NCFBLK, IDBLK, NEVBLK, NCMAXBLK) -> (IDBLK)
        !     CALL LODSTATE (NBLOCK, NCFBLK(1), IDBLK, NEVBLK, NCMAXBLK)
        call read_block_level_selection (idblk)
        !cjb
        !
        !   Allocate the storage for and set the weights
        !
        call alloc (requested_level_weights, selected_level_count, "WEIGHT", "GETOLD")

        call read_eol_level_weights (nondefault_options_flag, selected_level_count, requested_level_weights)
        !
        !   Eigenvector damping
        !
        call alloc (mixing_damping_factors, selected_level_count, "CDAMP", "GETOLD")
        !
        mixing_damping_factors(:selected_level_count) = 0.0_real64
        !
        !   Print the list of all subshells
        !
        write (stderr_unit, *) "Radial functions"
        call print_scf_results
        !
        !   Determine which orbitals are to be varied, which are fixed.
        !   Quantities determined here:
        !     nfix, lfix(1:nw), iorder(1:nw), scnsty(1:nw)
        !   Instead of broadcasting these quantities, we broadcast
        !   the intermediate result from GETRSL (see below)
        !
        is_orbital_fixed(:orbital_count) = .true.

        write (stderr_unit, *) "Enter orbitals to be varied (Updating order)"
        call select_subshell_labels (indx, nsubs)

        is_orbital_fixed(indx(:nsubs)) = .false.
        !XHH      give a big value, rather than zero to scnsty()
        orbital_self_consistency_errors(indx(:nsubs)) = 1.0e20_real64
        fixed_orbital_count = orbital_count - nsubs
        if (fixed_orbital_count == orbital_count) write (stderr_unit, *) &
            "All subshell radial wavefunctions are fixed;", &
            " performing CI calculation."

        !   Determine orbital updating order

        norder = 0
        do i = 1, orbital_count
            orbital_update_order(i) = i
            if (is_orbital_fixed(i)) cycle
            norder = norder + 1
            orbital_update_order(i) = indx(norder)
        end do
        !
        !XHH added a array to store the index of the correlation functions
        !
        is_correlation_orbital(:orbital_count) = .true.

        write (stderr_unit, *) "Which of these are spectroscopic orbitals?"
        call select_subshell_labels (indx, nsubs)
        if (nsubs > 0) then
            do i = 1, nsubs
                loc = indx(i)
                if (is_orbital_fixed(loc)) cycle
                radial_solution_methods(loc) = 1
                suppress_radial_inversion(loc) = .false.
                orbital_damping_factors(loc) = 0.0_real64
                is_correlation_orbital(loc) = .false.
            end do
        end if

        ! Set NSIC. It will be non-zero if all orbitals to be varied are
        ! spectroscopic orbitals

        orbital_updates_per_scf_iteration = (orbital_count - fixed_orbital_count)/4
        do i = 1, orbital_count
            if (.not.(.not.is_orbital_fixed(i) .and. is_correlation_orbital(i))) cycle
            orbital_updates_per_scf_iteration = 0
            exit
        end do
        !
        maximum_scf_iterations = 24
        minimum_radial_solver_iterations = 3
        enforce_orbital_orthogonality = .true.
        !
        !   Make the allocation for the auxiliary vector required
        !   by SUBROUTINE NEWCO
        !
        call alloc (level_average_weights, selected_level_count, "RWTDUM", "GETOLD")
        !
        !   Place the block numbers of the all ncmin eigenstate(wanted)
        !   in array idxblk
        !
        call alloc (level_block_indices, selected_level_count, "IDXBLK", "GETOLD")
        noffset = 0
        do symmetry_block_index = 1, symmetry_block_count
            do level_index = 1, eigenvector_counts_by_block(symmetry_block_index)
                level_block_indices(level_index+noffset) = symmetry_block_index
            end do
            noffset = noffset + eigenvector_counts_by_block(symmetry_block_index)
        end do
        if (noffset /= selected_level_count) stop "getold: ncmin trouble"

        return
    end subroutine read_eol_average_options
end module read_eol_average_options_interface
