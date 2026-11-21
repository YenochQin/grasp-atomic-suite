! Query occupations, subshell and intermediate angular momenta, and total CSF symmetry.
module csf_quantum_queries
    ! Module dependencies
    use standard_io_units, only: stderr_unit
    use orbital_csf_basis, only: csf_count, csf_subshell_occupations
    use grasp_parameters, only: max_orbitals
    use csf_coupling_data, only: csf_coupling_quantum_codes, csf_subshell_quantum_codes
    implicit none (type, external)

    ! Module interface
    private
    public :: csf_occupation, csf_subshell_quantum_code, csf_intermediate_coupling, csf_parity_code, csf_total_j_code
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

    ! Read 2J+1 for the nontrivial intermediate coupling at coupling_index.
    integer function csf_intermediate_coupling(coupling_index, csf_index)
        ! Arguments
        integer, intent(in) :: coupling_index, csf_index

        if (coupling_index>=1 .and. coupling_index<=max_orbitals-1) then
            if (csf_index>=1 .and. csf_index<=csf_count) then
                csf_intermediate_coupling = csf_coupling_quantum_codes(coupling_index,csf_index)
            else
                write (stderr_unit, *) "csf_intermediate_coupling: Argument csf_index is out of range."
                stop
            end if
        else
            write (stderr_unit, *) "csf_intermediate_coupling: Argument coupling index is out of range."
            stop
        end if
    end function csf_intermediate_coupling

    ! Decode parity from the sign of the packed total symmetry.
    integer function csf_parity_code(csf_index)
        ! Arguments
        integer, intent(in) :: csf_index

        if (csf_index>=1 .and. csf_index<=csf_count) then
            csf_parity_code = csf_coupling_quantum_codes(max_orbitals,csf_index)
            if (csf_parity_code > 127) csf_parity_code = csf_parity_code - 256
            csf_parity_code = sign(1,csf_parity_code)
        else
            write (stderr_unit, *) "csf_parity_code: Argument csf_index is out of range."
            write (stderr_unit, *) "ICSF =",csf_index," NCF =",csf_count
            stop
        end if
    end function csf_parity_code

    ! Decode 2J+1 from the magnitude of the packed total symmetry.
    integer function csf_total_j_code(csf_index)
        ! Arguments
        integer, intent(in) :: csf_index

        if (csf_index>=1 .and. csf_index<=csf_count) then
            csf_total_j_code = csf_coupling_quantum_codes(max_orbitals,csf_index)
            if (csf_total_j_code > 127) csf_total_j_code = 256 - csf_total_j_code
            csf_total_j_code = iabs(csf_total_j_code)
        else
            write (stderr_unit, *) "ITJPO: Argument ICSF is out of range."
            stop
        end if
    end function csf_total_j_code
end module csf_quantum_queries
