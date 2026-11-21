module read_biorthogonal_mcp_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_biorthogonal_mcp
contains

    !***********************************************************************
    !                                                                      *
    subroutine read_biorthogonal_mcp (name,ik,global_trace_level,mixing_file_kind)
        ! Module dependencies
        use biorthogonal_subshell_layout, only: biorthogonal_angular_unit
        use calculation_file_units, only: mixing_input_unit
        !                                                                      *
        !   This routine controls the computation  and storage of the values   *
        !   and all indices of the angular coefficients                        *
        !                                                                      *
        !                                                                      *
        !                   T  (ab)                                            *
        !                    rs                                                *
        !                                                                      *
        !   k is the multipolarity of a two-particle Coulomb integral. a, b,   *
        !   are orbital sequence numbers.  r and s are configuration           *
        !   state function indices.                                            *
        !                                                                      *
        !   Call(s) to: [LIB92]: ALLOC, DALLOC, TNSRJJ                         *
        !               [GENMCP]: QSORT.                                       *
        !                                                                      *
        !   Written by Per Jonsson                Last revision: JUne 1996     *
        !   Bug corrected 2005-10-18                                           *
        !                                                                      *
        !***********************************************************************
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: orbital_index_radix, max_orbitals
        use array_memory_management, only: alloc, dalloc
        use atomic_calculation_state, only: electron_count
        use eigenstate_expansions, only: hamiltonian_energy_shift, shifted_eigenvalues, csf_mixing_coefficients
        use orbital_csf_basis,           only: csf_count, orbital_count
        use selected_level_indices, only: active_level_indices, active_level_count
        use level_symmetry_metadata, only: level_parity_signs, level_two_j_plus_one
        use mcp_coefficient_storage, only: mcp_coefficient_values, mcp_packed_orbital_pairs, mcp_coefficient_end_offsets, &
            mcp_bra_csf_indices, mcp_ket_csf_indices, mcp_coefficient_count, mcp_integral_count
        use biorthogonal_subshell_layout,         only: kappa_group_count, initial_orbital_counts_by_group, &
            final_orbital_counts_by_group, initial_orbital_positions_by_group, final_orbital_positions_by_group, &
            initial_core_counts_by_group, final_core_counts_by_group, initial_orbitals_by_group, final_orbitals_by_group, &
            initial_kappa_group_by_orbital, final_kappa_group_by_orbital, max_biorthogonal_kappa_groups
        use biorthogonal_subshell_mapping,        only: active_orbitals_by_group, active_orbital_positions_by_group
        use biorthogonal_orbital_rotations,         only: initial_orbital_rotation_values, final_orbital_rotation_values
        use transition_basis_blocks, only: transition_block_count, mixing_file_csf_count, mixing_file_level_count

        ! Procedure interfaces
        !      USE tnsrjj_I
        !      USE cor_I
        !      USE cord_I
        !      USE itjpo_I
        use transform_ci_coefficients_interface, only: transform_ci_coefficients
        use read_state_mixing_coefficients_interface, only: read_state_mixing_coefficients
        implicit none (type, external)

        ! Local declarations
        integer :: biorthogonal_mixing_output_unit

        ! Arguments
        character(len=24), intent(in) :: name
        integer, intent(in) :: ik, global_trace_level, mixing_file_kind

        ! Local constants
        !   Set the encoding key and its square
        integer, parameter :: orbital_label_radix_squared = orbital_index_radix*orbital_index_radix, nvmax=100
        real(real64), parameter :: cutoff = 1.0e-10_real64

        ! Local declarations
        !
        integer :: kappa_group_by_orbital(max_orbitals)
        integer :: shell_counts(max_biorthogonal_kappa_groups), inactive_shell_counts(max_biorthogonal_kappa_groups)
        real(real64) :: ci_rotation_matrices(20*max_biorthogonal_kappa_groups*max_biorthogonal_kappa_groups)
        !
        real(real64), pointer :: transfer_power_workspace(:), transformed_ci_workspace(:)
        integer :: trace_level, local_trace_level, i, j, symmetry_block_index, cached_csf_count, cached_orbital_count, &
                & cached_symmetry_count,  &
            symmetry_index, ci_rotation_offset, nvecsize
        !
        local_trace_level = 00
        trace_level = max(global_trace_level,local_trace_level)
        trace_level = 00
        !
        !   Set up data for the initial and final state case respectively
        !
        if (ik==1) then
            do i = 1,max_biorthogonal_kappa_groups
                shell_counts(i) = initial_orbital_counts_by_group(i)
                inactive_shell_counts(i) = initial_core_counts_by_group(i)
            end do

            do i = 1,max_orbitals
                kappa_group_by_orbital(i) = initial_kappa_group_by_orbital(i)
            end do

            do j = 1,max_biorthogonal_kappa_groups
                do i = 1,max_biorthogonal_kappa_groups
                    active_orbitals_by_group(i,j) = initial_orbitals_by_group(i,j)
                end do
            end do

            do j = 1,max_orbitals
                do i = 1,max_biorthogonal_kappa_groups
                    active_orbital_positions_by_group(i,j) = initial_orbital_positions_by_group(i,j)
                end do
            end do

            do i = 1,20*max_biorthogonal_kappa_groups*max_biorthogonal_kappa_groups
                ci_rotation_matrices(i) = initial_orbital_rotation_values(i)
            end do

        else if (ik==2) then

            do i = 1,max_biorthogonal_kappa_groups
                shell_counts(i) = final_orbital_counts_by_group(i)
                inactive_shell_counts(i) = final_core_counts_by_group(i)
            end do

            do i = 1,max_orbitals
                kappa_group_by_orbital(i) = final_kappa_group_by_orbital(i)
            end do

            do j = 1,max_biorthogonal_kappa_groups
                do i = 1,max_biorthogonal_kappa_groups
                    active_orbitals_by_group(i,j) = final_orbitals_by_group(i,j)
                end do
            end do

            do j = 1,max_orbitals
                do i = 1,max_biorthogonal_kappa_groups
                    active_orbital_positions_by_group(i,j) = final_orbital_positions_by_group(i,j)
                end do
            end do

            do i = 1,20*max_biorthogonal_kappa_groups*max_biorthogonal_kappa_groups
                ci_rotation_matrices(i) = final_orbital_rotation_values(i)
            end do
        end if

        rewind (biorthogonal_angular_unit)
        do symmetry_block_index = 1, transition_block_count
            !
            !   Read the CI vectors
            !
            call read_state_mixing_coefficients(name,mixing_file_kind,symmetry_block_index)
            !
            !   Allocate memory for use in CITRAG
            !
            call alloc (transfer_power_workspace, csf_count*active_level_count, "transfer_power_workspace", "read_biorthogonal_mcp")
            call alloc (transformed_ci_workspace, csf_count*active_level_count, "transformed_ci_workspace", "read_biorthogonal_mcp")
            !
            !   if data available on file read the datafile. Else sort the
            !   MCP data into inegral based lists.
            !
            read(biorthogonal_angular_unit) cached_csf_count,cached_orbital_count,cached_symmetry_count
            do symmetry_index = 1,kappa_group_count

                !************
                !
                !. Offset for given L in shell matrices
                !
                if (symmetry_index==1) then
                    ci_rotation_offset = 1
                else
                    ci_rotation_offset = ci_rotation_offset + shell_counts(symmetry_index-1)** 2
                end if
                !  Corrected PER J

                !**************

                read(biorthogonal_angular_unit) mcp_integral_count,mcp_coefficient_count
                if(mcp_coefficient_count*mcp_integral_count==0) cycle
                !
                !   Allocate memory. Note that in this case we can,
                !   since the data is sorted, allocate less
                !   memory than is done in qqsort.
                !
                call alloc (mcp_bra_csf_indices,mcp_coefficient_count, "mcp_bra_csf_indices", "read_biorthogonal_mcp")
                call alloc (mcp_ket_csf_indices,mcp_coefficient_count, "mcp_ket_csf_indices", "read_biorthogonal_mcp")
                call alloc (mcp_packed_orbital_pairs,mcp_integral_count, "mcp_packed_orbital_pairs", "read_biorthogonal_mcp")
                call alloc (mcp_coefficient_values,mcp_coefficient_count, "mcp_coefficient_values", "read_biorthogonal_mcp")
                call alloc (mcp_coefficient_end_offsets,mcp_integral_count, "mcp_coefficient_end_offsets", "read_biorthogonal_mcp")

                do i = 1,mcp_integral_count
                    read(biorthogonal_angular_unit) mcp_packed_orbital_pairs(i),mcp_coefficient_end_offsets(i)
                end do
                do i = 1,mcp_coefficient_count
                    read(biorthogonal_angular_unit) mcp_coefficient_values(i),mcp_bra_csf_indices(i),mcp_ket_csf_indices(i)
                end do
                !*
                !*. Offset for given L in shell matrices
                !*
                !        IF (L.EQ.1) THEN
                !          IIOFF = 1
                !        ELSE
                !          IIOFF = IIOFF + NSHL(L-1)** 2
                !        END IF
                !   BUG THE OFFSET SHOULD BE BEFORE THE CYCLE STATEMENT!!
                !
                ! Transform for given L
                !
                if (shell_counts(symmetry_index) > 0) &
                    call transform_ci_coefficients(csf_mixing_coefficients(:), csf_count, active_level_count, symmetry_index, &
                        shell_counts(symmetry_index), ci_rotation_matrices(ci_rotation_offset), &
                & inactive_shell_counts(symmetry_index), local_trace_level, transformed_ci_workspace, transfer_power_workspace)
                !
                !   Deallocate storage for this kappa.
                !
                call dalloc (mcp_bra_csf_indices, "mcp_bra_csf_indices", "read_biorthogonal_mcp")
                call dalloc (mcp_ket_csf_indices, "mcp_ket_csf_indices", "read_biorthogonal_mcp")
                call dalloc (mcp_packed_orbital_pairs, "mcp_packed_orbital_pairs", "read_biorthogonal_mcp")
                call dalloc (mcp_coefficient_values, "mcp_coefficient_values", "read_biorthogonal_mcp")
                call dalloc (mcp_coefficient_end_offsets, "mcp_coefficient_end_offsets", "read_biorthogonal_mcp")

            end do

            call dalloc (transfer_power_workspace, "transfer_power_workspace", "read_biorthogonal_mcp")
            call dalloc (transformed_ci_workspace, "transformed_ci_workspace", "read_biorthogonal_mcp")
            !
            !   Write the rotated CI vectors on file
            !
            if(symmetry_block_index==1) then
                j = index(name," ")
                if (mixing_file_kind==0) then
                    open (newunit=biorthogonal_mixing_output_unit,file=name(1:j-1)//".cbm",form="UNFORMATTED",  &
                        &    status="UNKNOWN")
                else
                    open (newunit=biorthogonal_mixing_output_unit,file=name(1:j-1)//".bm",form="UNFORMATTED",   &
                        &    status="UNKNOWN")
                end if

                write(biorthogonal_mixing_output_unit) "G92MIX"
                write(biorthogonal_mixing_output_unit) electron_count,mixing_file_csf_count,orbital_count, &
                    mixing_file_level_count,nvecsize,transition_block_count
            end if

            write(biorthogonal_mixing_output_unit) symmetry_block_index,csf_count,active_level_count,level_two_j_plus_one(1), &
                & level_parity_signs(1)
            write(biorthogonal_mixing_output_unit) (active_level_indices(i),i = 1,active_level_count)
            !     WRITE(31) (IATJPO(I),IASPAR(I),I = 1,NVEC)
            write(biorthogonal_mixing_output_unit) hamiltonian_energy_shift,(shifted_eigenvalues(i),i = 1,active_level_count)
            write(biorthogonal_mixing_output_unit) ((csf_mixing_coefficients(i+(j-1)*csf_count),i = 1,csf_count),j = 1, &
                active_level_count)

            call dalloc (shifted_eigenvalues, "shifted_eigenvalues", "read_biorthogonal_mcp")
            call dalloc (csf_mixing_coefficients, "csf_mixing_coefficients", "read_biorthogonal_mcp")
            call dalloc (active_level_indices, "active_level_indices", "read_biorthogonal_mcp")
            call dalloc (level_two_j_plus_one, "level_two_j_plus_one", "read_biorthogonal_mcp")
            call dalloc (level_parity_signs, "level_parity_signs", "read_biorthogonal_mcp")
        end do
        close (biorthogonal_angular_unit)
        !
        !   Close the  mixing  file
        !
        close (mixing_input_unit)
        close (biorthogonal_mixing_output_unit)
        !
        return
    end subroutine read_biorthogonal_mcp
end module read_biorthogonal_mcp_interface
