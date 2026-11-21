! Step-indexed coupling arrays used only by the historical text readers.
module csf_text_oracle_coupling_data
    use iso_fortran_env, only: int8
    use grasp_parameters, only: max_orbitals
    implicit none (type, external)
    private :: int8, max_orbitals
    integer, pointer :: left_subshell_quantum_codes(:,:,:) => null(), left_coupling_codes(:,:) => null()
    integer(int8), pointer :: csf_subshell_quantum_codes(:,:,:) => null(), csf_coupling_quantum_codes(:,:) => null()
contains
    integer function csf_subshell_quantum_code(selector, orbital_index, csf_index)
        integer, intent(in) :: selector, orbital_index, csf_index
        csf_subshell_quantum_code = csf_subshell_quantum_codes(orbital_index,selector,csf_index)
    end function csf_subshell_quantum_code

    integer function csf_intermediate_coupling(step, csf_index)
        integer, intent(in) :: step, csf_index
        csf_intermediate_coupling = csf_coupling_quantum_codes(step,csf_index)
    end function csf_intermediate_coupling

    integer function csf_total_j_code(csf_index)
        integer, intent(in) :: csf_index
        csf_total_j_code = abs(int(csf_coupling_quantum_codes(max_orbitals,csf_index)))
    end function csf_total_j_code

    integer function csf_parity_code(csf_index)
        integer, intent(in) :: csf_index
        csf_parity_code = sign(1,int(csf_coupling_quantum_codes(max_orbitals,csf_index)))
    end function csf_parity_code
end module csf_text_oracle_coupling_data
