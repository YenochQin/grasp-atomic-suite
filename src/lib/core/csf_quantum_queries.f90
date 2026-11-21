! Query occupations, subshell and intermediate angular momenta, and total CSF symmetry.
module csf_quantum_queries
    ! Module dependencies
    use standard_io_units, only: stderr_unit
    use orbital_csf_basis, only: csf_count, csf_subshell_occupations, orbital_count, orbital_two_j, &
        left_orbital_count, left_csf_count, left_occupations, left_kappa
    use grasp_parameters, only: max_orbitals
    use csf_coupling_data, only: csf_angular_codes, csf_subshell_quantum_codes, left_angular_codes, &
        csf_total_j_slot, csf_parity_slot
    implicit none (type, external)

    ! Module interface
    private
    public :: csf_occupation, csf_subshell_quantum_code, csf_cumulative_coupling
    public :: csf_intermediate_coupling, left_intermediate_coupling, csf_parity_code, csf_total_j_code
contains

    ! Return the stored occupation as a default-kind integer.
    integer function csf_occupation(subshell_index, csf_index)
        ! Arguments
        integer, intent(in) :: subshell_index, csf_index

        csf_occupation = csf_subshell_occupations(subshell_index,csf_index)
    end function csf_occupation

    ! Selectors 1, 2, and 3 return seniority, w, and subshell 2J+1, respectively.
    integer function csf_subshell_quantum_code(quantum_number_selector, subshell_index, csf_index)
        ! Arguments
        integer, intent(in) :: quantum_number_selector, subshell_index, csf_index

        csf_subshell_quantum_code = csf_subshell_quantum_codes(subshell_index,quantum_number_selector,csf_index)
    end function csf_subshell_quantum_code

    integer function csf_cumulative_coupling(orbital_index, csf_index)
        integer, intent(in) :: orbital_index, csf_index
        csf_cumulative_coupling = csf_angular_codes(orbital_index,csf_index)
    end function csf_cumulative_coupling

    ! Step t ends at the (t+1)th open shell. Older callers also query unused steps.
    integer function csf_intermediate_coupling(coupling_index, csf_index) result(code)
        integer, intent(in) :: coupling_index, csf_index
        integer :: orbital_index, open_shells, occupation

        if (csf_index < 1 .or. csf_index > csf_count) error stop "CSF index is out of range"
        if (coupling_index < 1 .or. coupling_index >= max_orbitals) error stop "Coupling step is out of range"
        code = 0
        open_shells = 0
        do orbital_index = 1, orbital_count
            occupation = csf_subshell_occupations(orbital_index,csf_index)
            if (occupation == 0 .or. occupation == orbital_two_j(orbital_index)+1) cycle
            open_shells = open_shells + 1
            if (open_shells /= coupling_index+1) cycle
            code = csf_cumulative_coupling(orbital_index,csf_index)
            return
        end do
    end function csf_intermediate_coupling

    integer function left_intermediate_coupling(coupling_index, csf_index) result(code)
        integer, intent(in) :: coupling_index, csf_index
        integer :: orbital_index, open_shells, occupation

        if (csf_index < 1 .or. csf_index > left_csf_count) error stop "Left CSF index is out of range"
        if (coupling_index < 1 .or. coupling_index >= left_orbital_count) error stop "Left coupling step is out of range"
        code = 0
        open_shells = 0
        do orbital_index = 1, left_orbital_count
            occupation = left_occupations(orbital_index,csf_index)
            if (occupation == 0 .or. occupation == 2*abs(left_kappa(orbital_index))) cycle
            open_shells = open_shells + 1
            if (open_shells /= coupling_index+1) cycle
            code = left_angular_codes(orbital_index,csf_index)
            return
        end do
    end function left_intermediate_coupling

    integer function csf_parity_code(csf_index)
        ! Arguments
        integer, intent(in) :: csf_index

        if (csf_index>=1 .and. csf_index<=csf_count) then
            csf_parity_code = csf_angular_codes(csf_parity_slot,csf_index)
        else
            write (stderr_unit, *) "csf_parity_code: Argument csf_index is out of range."
            write (stderr_unit, *) "ICSF =",csf_index," NCF =",csf_count
            stop
        end if
    end function csf_parity_code

    ! Return the positive total 2J+1.
    integer function csf_total_j_code(csf_index)
        ! Arguments
        integer, intent(in) :: csf_index

        if (csf_index>=1 .and. csf_index<=csf_count) then
            csf_total_j_code = csf_angular_codes(csf_total_j_slot,csf_index)
        else
            write (stderr_unit, *) "csf_total_j_code: CSF index is out of range."
            stop
        end if
    end function csf_total_j_code
end module csf_quantum_queries
