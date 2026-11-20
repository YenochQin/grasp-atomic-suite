! Build or restore the coupling trees and angular labels used to recouple a CSF pair.
module angular_recoupling_tables
    implicit none (type, external)

    ! Module interface
    private
    public :: setj
contains
    subroutine setj(orbital_indices, coupling_node_indices, active_two_j_plus_one, active_shell_count, coupling_tables_ready)
        ! Module dependencies
        use csf_pair_subshell_couplings, only: jbq1, jbq2, jtq1, jtq2, max_orbitals
        use angular_coupling_triads, only: j2s, j3s, mtriad
        use csf_pair_angular_state,      only: bra_intermediate_couplings => jjc1, ket_intermediate_couplings => jjc2, &
            & bra_subshell_quantum_codes => jjq1, ket_subshell_quantum_codes => jjq2, peel_orbital_indices => jlist
        use angular_coupling_workspace, only: mja, nja, j1, j2, j3, mangm
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: active_shell_count
        integer, intent(inout) :: coupling_tables_ready
        integer, intent(in) :: orbital_indices(2, 2)
        integer, intent(in) :: coupling_node_indices(2, 2)
        integer, intent(in) :: active_two_j_plus_one(2, 2)

        ! Local declarations
        integer :: angular_list_end, orbital_index, table_index, intermediate_coupling_count, coupling_table_size, &
            & three_times_shell_count, first_coupling_node, second_coupling_node, j
        angular_list_end = 0
        if (active_shell_count > 0) then
            j1(:active_shell_count) = jbq1(3,peel_orbital_indices(:active_shell_count))
            angular_list_end = active_shell_count
        end if
        if (active_shell_count /= 1) then
            intermediate_coupling_count = active_shell_count - 1
            if (intermediate_coupling_count > 0) then
                j1(angular_list_end+1:intermediate_coupling_count+angular_list_end) = &
                    & bra_intermediate_couplings(:intermediate_coupling_count)
                angular_list_end = intermediate_coupling_count + angular_list_end
            end if
            if (intermediate_coupling_count > 0) then
                j1(angular_list_end+1:intermediate_coupling_count+angular_list_end) = &
                    & ket_intermediate_couplings(:intermediate_coupling_count)
                angular_list_end = intermediate_coupling_count + angular_list_end
            end if
        end if
        do table_index = 1, 2
            angular_list_end = angular_list_end + 1
            orbital_index = orbital_indices(table_index,1)
            j1(angular_list_end) = bra_subshell_quantum_codes(3,orbital_index)
            if (table_index==1 .and. orbital_indices(1,1)==orbital_indices(2,1)) j1(angular_list_end) = jtq1(3)
            j1(angular_list_end+4) = active_two_j_plus_one(table_index,1)
        end do
        do table_index = 1, 2
            angular_list_end = angular_list_end + 1
            orbital_index = orbital_indices(table_index,2)
            j1(angular_list_end) = ket_subshell_quantum_codes(3,orbital_index)
            if (table_index==1 .and. orbital_indices(1,2)==orbital_indices(2,2)) j1(angular_list_end) = jtq2(3)
            j1(angular_list_end+4) = active_two_j_plus_one(table_index,2)
        end do
        coupling_table_size = max(4,active_shell_count + 2)
        ! Build both coupling trees when the cache is not marked ready.
        if (coupling_tables_ready <= 0) then
            do table_index = 4, coupling_table_size
                j2(table_index,1) = active_shell_count + table_index - 4
                j2(table_index,2) = table_index - 2
                j2(table_index,3) = active_shell_count + table_index - 3
                j3(table_index,1) = j2(table_index,1) + active_shell_count - 1
                j3(table_index,2) = table_index - 2
                j3(table_index,3) = j2(table_index,3) + active_shell_count - 1
            end do
            j2(4,1) = 1
            j3(4,1) = 1
            three_times_shell_count = 3*active_shell_count
            j2(1,1) = three_times_shell_count + 5
            j2(1,2) = three_times_shell_count + 7
            j2(1,3) = three_times_shell_count + 3
            j2(2,1) = coupling_node_indices(1,1)
            j2(2,2) = three_times_shell_count + 3
            j2(2,3) = three_times_shell_count - 1
            j2(3,1) = coupling_node_indices(2,1)
            j2(3,2) = three_times_shell_count + 4
            j2(3,3) = three_times_shell_count
            j3(1,1) = three_times_shell_count + 7
            j3(1,2) = three_times_shell_count + 4
            j3(1,3) = three_times_shell_count + 6
            j3(2,1) = coupling_node_indices(1,2)
            j3(2,2) = three_times_shell_count + 5
            j3(2,3) = three_times_shell_count + 1
            j3(3,1) = coupling_node_indices(2,2)
            j3(3,2) = three_times_shell_count + 6
            j3(3,3) = three_times_shell_count + 2
            first_coupling_node = coupling_node_indices(1,1)
            second_coupling_node = coupling_node_indices(2,1)
            if (second_coupling_node > 1) j2(second_coupling_node+2,2) = j2(3,3)
            if (second_coupling_node == 1) j2(4,1) = j2(3,3)
            if (first_coupling_node == second_coupling_node) then
                j2(3,1) = j2(2,3)
            else
                if (first_coupling_node > 1) j2(first_coupling_node+2,2) = j2(2,3)
                if (first_coupling_node == 1) j2(4,1) = j2(2,3)
            end if
            first_coupling_node = coupling_node_indices(1,2)
            second_coupling_node = coupling_node_indices(2,2)
            if (second_coupling_node > 1) j3(second_coupling_node+2,2) = j3(3,3)
            if (second_coupling_node == 1) j3(4,1) = j3(3,3)
            if (first_coupling_node == second_coupling_node) then
                j3(3,1) = j3(2,3)
            else
                if (first_coupling_node > 1) j3(first_coupling_node+2,2) = j3(2,3)
                if (first_coupling_node == 1) j3(4,1) = j3(2,3)
            end if
            coupling_tables_ready = 1
            mja = three_times_shell_count + 7
            nja = active_shell_count + 3
            j2s(:coupling_table_size,:) = j2(:coupling_table_size,:)
            j3s(:coupling_table_size,:) = j3(:coupling_table_size,:)
            return
        end if
        j2(:coupling_table_size,:) = j2s(:coupling_table_size,:)
        j3(:coupling_table_size,:) = j3s(:coupling_table_size,:)
        return
    end subroutine setj
end module angular_recoupling_tables
