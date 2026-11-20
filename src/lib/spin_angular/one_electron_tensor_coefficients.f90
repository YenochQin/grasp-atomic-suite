! Generate one-electron angular coefficients for a specified rank and parity.
module one_electron_tensor_coefficients
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_one_electron_tensor_coefficients
contains
    subroutine calculate_one_electron_tensor_coefficients(tensor_rank,operator_parity,bra_csf_index,ket_csf_index, &
        & bra_orbital_index, &
        ket_orbital_index,shell_coefficients)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only:  max_orbitals
        use numeric_constants, only: zero
        use csf_pair_angular_state, only: bra_intermediate_couplings, ket_intermediate_couplings, bra_subshell_quantum_codes, &
            ket_subshell_quantum_codes, peel_orbital_indices, max_orbitals, peel_orbital_count, bra_occupations, ket_occupations
        use orbital_csf_basis,   only: orbital_count, orbital_kappa
        use angular_recoupling_backup,  only: saved_peel_orbitals, saved_bra_j_weights, saved_ket_j_weights
        use angular_momentum_triangle, only: satisfies_angular_triangle
        use one_electron_csf_filter, only: check_one_electron_csf_pair
        use subshell_occupation_class, only: classify_subshell_occupation
        use csf_quantum_queries, only: csf_parity_code, csf_total_j_code
        use csf_pair_quantum_numbers, only: prepare_csf_pair_quantum_numbers
        use one_electron_tensor_same_subshell, only: calculate_tensor_same_subshell_coefficient
        use one_electron_tensor_between_subshells, only: calculate_tensor_transfer_coefficient
        implicit none (type, external)

        ! Arguments
        integer, intent(in)  :: tensor_rank,operator_parity,bra_csf_index,ket_csf_index
        integer, intent(out) :: bra_orbital_index,ket_orbital_index
        real(real64), intent(out) :: shell_coefficients(max_orbitals)

        ! Local declarations
        integer :: two_rank_plus_one,orbital_index,total_occupation_difference,bra_shell_position,ket_shell_position, &
            peel_position,occupation_difference,active_shell_count,diagonal_orbital,shell_position,shift_step, &
            first_orbital,shift_count,  &
            two_bra_j_plus_one,two_ket_j_plus_one,orbital_difference,peel_coupling_count,subshell_occupation
        real(real64) :: angular_coefficient

        ! A zero bra orbital index signals that no coefficient channel has been identified.
        bra_orbital_index = 0
        two_rank_plus_one = tensor_rank+tensor_rank+1
        if(satisfies_angular_triangle(csf_total_j_code(bra_csf_index),csf_total_j_code(ket_csf_index),two_rank_plus_one)==0)return
        if((operator_parity/=0).and.(csf_parity_code(bra_csf_index)*csf_parity_code(ket_csf_index)*operator_parity/=1))return
        if(check_one_electron_csf_pair(bra_csf_index,ket_csf_index)==0)return
        call prepare_csf_pair_quantum_numbers (bra_csf_index,ket_csf_index)
        do orbital_index = 1,orbital_count
            shell_coefficients(orbital_index) = zero
        end do
        total_occupation_difference = 0
        bra_shell_position = 0
        ket_shell_position = 0
        if (peel_orbital_count /= 0) then
            do peel_position = 1,peel_orbital_count
                orbital_index = peel_orbital_indices(peel_position)
                occupation_difference = bra_occupations(orbital_index)-ket_occupations(orbital_index)
                if (abs (occupation_difference) > 1) return
                if (occupation_difference > 0) then
                    bra_shell_position = peel_position
                    total_occupation_difference = total_occupation_difference+1
                else if (occupation_difference < 0) then
                    ket_shell_position = peel_position
                    total_occupation_difference = total_occupation_difference+1
                end if
            end do
            if (total_occupation_difference > 2) return
            active_shell_count = peel_orbital_count
        end if

        ! A one-electron transfer fixes both active subshells. Equal occupations require a scan.
        diagonal_orbital = 0
        if (total_occupation_difference /= 2) then
            do shell_position = 1,peel_orbital_count
                saved_peel_orbitals(shell_position) = peel_orbital_indices(shell_position)
            end do
            peel_coupling_count = peel_orbital_count-1
            do shell_position = 1,peel_coupling_count
                saved_bra_j_weights(shell_position) = bra_intermediate_couplings(shell_position)
                saved_ket_j_weights(shell_position) = ket_intermediate_couplings(shell_position)
            end do
        end if

        do
            if (total_occupation_difference /= 2) then
                diagonal_orbital = diagonal_orbital+1
                if (diagonal_orbital > orbital_count) return
                subshell_occupation = classify_subshell_occupation(diagonal_orbital,bra_csf_index)
                if (subshell_occupation == -1) cycle
                if (subshell_occupation == 0) then
                    active_shell_count = peel_orbital_count
                    do peel_position = 1,peel_orbital_count
                        orbital_difference = diagonal_orbital-peel_orbital_indices(peel_position)
                        if (orbital_difference == 0) exit
                    end do
                    bra_shell_position = peel_position
                    ket_shell_position = peel_position
                else
                    shell_position = 1
                    if (peel_orbital_count /= 0) then
                        do shell_position = 1,peel_orbital_count
                            orbital_index = peel_orbital_indices(shell_position)
                            if (diagonal_orbital < orbital_index) exit
                        end do

                        ! Open a slot for a closed subshell in the ordered active coupling chain.
                        shift_count = peel_orbital_count-shell_position+1
                        do shift_step = 1,shift_count
                            peel_orbital_indices(peel_orbital_count+2-shift_step) = peel_orbital_indices( &
                                peel_orbital_count+1-shift_step)
                            if (peel_orbital_count == shift_step) exit
                            bra_intermediate_couplings(peel_orbital_count+1-shift_step) = bra_intermediate_couplings( &
                                peel_orbital_count-shift_step)
                            ket_intermediate_couplings(peel_orbital_count+1-shift_step) = ket_intermediate_couplings( &
                                peel_orbital_count-shift_step)
                        end do
                        if (shell_position < 3) then
                            first_orbital = peel_orbital_indices(1)
                            bra_intermediate_couplings(1) = bra_subshell_quantum_codes(3,first_orbital)
                            ket_intermediate_couplings(1) = ket_subshell_quantum_codes(3,first_orbital)
                        else
                            bra_intermediate_couplings(shell_position-1) = bra_intermediate_couplings(shell_position-2)
                            ket_intermediate_couplings(shell_position-1) = ket_intermediate_couplings(shell_position-2)
                        end if
                    end if
                    peel_orbital_indices(shell_position) = diagonal_orbital
                    bra_shell_position = shell_position
                    ket_shell_position = shell_position
                    active_shell_count = peel_orbital_count+1
                end if
            end if

            bra_orbital_index = peel_orbital_indices(bra_shell_position)
            ket_orbital_index = peel_orbital_indices(ket_shell_position)
            two_bra_j_plus_one = 2*abs (orbital_kappa(bra_orbital_index))
            two_ket_j_plus_one = 2*abs (orbital_kappa(ket_orbital_index))
            if (satisfies_angular_triangle(two_bra_j_plus_one,two_ket_j_plus_one,two_rank_plus_one) == 1) then
                if (total_occupation_difference == 0) then
                    call calculate_tensor_same_subshell_coefficient(active_shell_count,tensor_rank,bra_csf_index,ket_csf_index, &
                        & bra_shell_position, &
                        ket_shell_position,angular_coefficient)
                else if (total_occupation_difference == 2) then
                    call calculate_tensor_transfer_coefficient(active_shell_count,tensor_rank,bra_shell_position, &
                        & ket_shell_position,angular_coefficient)

                    ! Store a transfer coefficient in slot 1 and return the orbital indices separately.
                    shell_coefficients(1) = angular_coefficient
                    return
                end if

                ! An equal-orbital contribution is indexed by its full orbital label.
                shell_coefficients(diagonal_orbital) = angular_coefficient
            end if
            if (total_occupation_difference == 2) return

            do shell_position = 1,peel_orbital_count
                peel_orbital_indices(shell_position) = saved_peel_orbitals(shell_position)
            end do
            peel_coupling_count = peel_orbital_count-1
            do shell_position = 1,peel_coupling_count
                bra_intermediate_couplings(shell_position) = saved_bra_j_weights(shell_position)
                ket_intermediate_couplings(shell_position) = saved_ket_j_weights(shell_position)
            end do
        end do
    end subroutine calculate_one_electron_tensor_coefficients
end module one_electron_tensor_coefficients
