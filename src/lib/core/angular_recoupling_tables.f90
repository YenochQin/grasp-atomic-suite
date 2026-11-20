! Build or restore the coupling trees and angular labels used to recouple a CSF pair.
module angular_recoupling_tables
    implicit none (type, external)

    ! Module interface
    private
    public :: initialize_angular_recoupling_tables
contains
    subroutine initialize_angular_recoupling_tables(orbital_indices, coupling_node_indices, active_two_j_plus_one, &
        & active_shell_count, coupling_tables_ready)
        ! Module dependencies
        use csf_pair_subshell_couplings, only: bra_background_subshell_codes, bra_active_pair_codes, ket_active_pair_codes
        use angular_coupling_triads, only: saved_bra_coupling_triads, saved_ket_coupling_triads
        use csf_pair_angular_state,      only: bra_intermediate_couplings, ket_intermediate_couplings, &
            bra_subshell_quantum_codes, ket_subshell_quantum_codes, peel_orbital_indices
        use angular_coupling_workspace, only: coupling_angular_value_count, coupling_triad_count, coupling_two_j_plus_one, &
            bra_coupling_triads, ket_coupling_triads
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: active_shell_count
        integer, intent(inout) :: coupling_tables_ready
        integer, intent(in) :: orbital_indices(2, 2)
        integer, intent(in) :: coupling_node_indices(2, 2)
        integer, intent(in) :: active_two_j_plus_one(2, 2)

        ! Local declarations
        integer :: angular_list_end, orbital_index, table_index, intermediate_coupling_count, coupling_table_size, &
            three_times_shell_count, first_coupling_node, second_coupling_node
        angular_list_end = 0
        if (active_shell_count > 0) then
            coupling_two_j_plus_one(:active_shell_count) = bra_background_subshell_codes(3, &
                peel_orbital_indices(:active_shell_count))
            angular_list_end = active_shell_count
        end if
        if (active_shell_count /= 1) then
            intermediate_coupling_count = active_shell_count - 1
            if (intermediate_coupling_count > 0) then
                coupling_two_j_plus_one(angular_list_end+1:intermediate_coupling_count+angular_list_end) = &
                    & bra_intermediate_couplings(:intermediate_coupling_count)
                angular_list_end = intermediate_coupling_count + angular_list_end
            end if
            if (intermediate_coupling_count > 0) then
                coupling_two_j_plus_one(angular_list_end+1:intermediate_coupling_count+angular_list_end) = &
                    & ket_intermediate_couplings(:intermediate_coupling_count)
                angular_list_end = intermediate_coupling_count + angular_list_end
            end if
        end if
        do table_index = 1, 2
            angular_list_end = angular_list_end + 1
            orbital_index = orbital_indices(table_index,1)
            coupling_two_j_plus_one(angular_list_end) = bra_subshell_quantum_codes(3,orbital_index)
            if (table_index==1 .and. orbital_indices(1,1)==orbital_indices(2,1)) coupling_two_j_plus_one(angular_list_end) = &
                bra_active_pair_codes(3)
            coupling_two_j_plus_one(angular_list_end+4) = active_two_j_plus_one(table_index,1)
        end do
        do table_index = 1, 2
            angular_list_end = angular_list_end + 1
            orbital_index = orbital_indices(table_index,2)
            coupling_two_j_plus_one(angular_list_end) = ket_subshell_quantum_codes(3,orbital_index)
            if (table_index==1 .and. orbital_indices(1,2)==orbital_indices(2,2)) coupling_two_j_plus_one(angular_list_end) = &
                ket_active_pair_codes(3)
            coupling_two_j_plus_one(angular_list_end+4) = active_two_j_plus_one(table_index,2)
        end do
        coupling_table_size = max(4,active_shell_count + 2)
        ! Build both coupling trees when the cache is not marked ready.
        if (coupling_tables_ready <= 0) then
            do table_index = 4, coupling_table_size
                bra_coupling_triads(table_index,1) = active_shell_count + table_index - 4
                bra_coupling_triads(table_index,2) = table_index - 2
                bra_coupling_triads(table_index,3) = active_shell_count + table_index - 3
                ket_coupling_triads(table_index,1) = bra_coupling_triads(table_index,1) + active_shell_count - 1
                ket_coupling_triads(table_index,2) = table_index - 2
                ket_coupling_triads(table_index,3) = bra_coupling_triads(table_index,3) + active_shell_count - 1
            end do
            bra_coupling_triads(4,1) = 1
            ket_coupling_triads(4,1) = 1
            three_times_shell_count = 3*active_shell_count
            bra_coupling_triads(1,1) = three_times_shell_count + 5
            bra_coupling_triads(1,2) = three_times_shell_count + 7
            bra_coupling_triads(1,3) = three_times_shell_count + 3
            bra_coupling_triads(2,1) = coupling_node_indices(1,1)
            bra_coupling_triads(2,2) = three_times_shell_count + 3
            bra_coupling_triads(2,3) = three_times_shell_count - 1
            bra_coupling_triads(3,1) = coupling_node_indices(2,1)
            bra_coupling_triads(3,2) = three_times_shell_count + 4
            bra_coupling_triads(3,3) = three_times_shell_count
            ket_coupling_triads(1,1) = three_times_shell_count + 7
            ket_coupling_triads(1,2) = three_times_shell_count + 4
            ket_coupling_triads(1,3) = three_times_shell_count + 6
            ket_coupling_triads(2,1) = coupling_node_indices(1,2)
            ket_coupling_triads(2,2) = three_times_shell_count + 5
            ket_coupling_triads(2,3) = three_times_shell_count + 1
            ket_coupling_triads(3,1) = coupling_node_indices(2,2)
            ket_coupling_triads(3,2) = three_times_shell_count + 6
            ket_coupling_triads(3,3) = three_times_shell_count + 2
            first_coupling_node = coupling_node_indices(1,1)
            second_coupling_node = coupling_node_indices(2,1)
            if (second_coupling_node > 1) bra_coupling_triads(second_coupling_node+2,2) = bra_coupling_triads(3,3)
            if (second_coupling_node == 1) bra_coupling_triads(4,1) = bra_coupling_triads(3,3)
            if (first_coupling_node == second_coupling_node) then
                bra_coupling_triads(3,1) = bra_coupling_triads(2,3)
            else
                if (first_coupling_node > 1) bra_coupling_triads(first_coupling_node+2,2) = bra_coupling_triads(2,3)
                if (first_coupling_node == 1) bra_coupling_triads(4,1) = bra_coupling_triads(2,3)
            end if
            first_coupling_node = coupling_node_indices(1,2)
            second_coupling_node = coupling_node_indices(2,2)
            if (second_coupling_node > 1) ket_coupling_triads(second_coupling_node+2,2) = ket_coupling_triads(3,3)
            if (second_coupling_node == 1) ket_coupling_triads(4,1) = ket_coupling_triads(3,3)
            if (first_coupling_node == second_coupling_node) then
                ket_coupling_triads(3,1) = ket_coupling_triads(2,3)
            else
                if (first_coupling_node > 1) ket_coupling_triads(first_coupling_node+2,2) = ket_coupling_triads(2,3)
                if (first_coupling_node == 1) ket_coupling_triads(4,1) = ket_coupling_triads(2,3)
            end if
            coupling_tables_ready = 1
            coupling_angular_value_count = three_times_shell_count + 7
            coupling_triad_count = active_shell_count + 3
            saved_bra_coupling_triads(:coupling_table_size,:) = bra_coupling_triads(:coupling_table_size,:)
            saved_ket_coupling_triads(:coupling_table_size,:) = ket_coupling_triads(:coupling_table_size,:)
            return
        end if
        bra_coupling_triads(:coupling_table_size,:) = saved_bra_coupling_triads(:coupling_table_size,:)
        ket_coupling_triads(:coupling_table_size,:) = saved_ket_coupling_triads(:coupling_table_size,:)
        return
    end subroutine initialize_angular_recoupling_tables
end module angular_recoupling_tables
