! Encode one subshell's bra and ket term, occupation and quasispin data.
module jj_subshell_state_packing
    implicit none (type, external)

    ! Module interface
    private
    public :: pack_jj_subshell_state
contains
    subroutine pack_jj_subshell_state(shell_position,bra_shell_values,bra_shell_codes,ket_shell_values,ket_shell_codes)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: half
        use csf_pair_angular_state, only: peel_orbital_indices, bra_occupations, ket_occupations, bra_subshell_quantum_codes, &
            ket_subshell_quantum_codes
        use orbital_csf_basis,           only: orbital_principal_n, orbital_kappa
        use jj_subshell_term_lookup, only: lookup_jj_subshell_term
        implicit none (type, external)

        ! Arguments
        integer,      intent(in)                :: shell_position
        integer,      intent(out) :: bra_shell_codes(7), ket_shell_codes(7)
        real(real64), intent(out) :: bra_shell_values(3), ket_shell_values(3)

        ! Local declarations
        integer :: orbital_index
        orbital_index=peel_orbital_indices(shell_position)

        ! Integer slots 1:7 store the term index, n, 2j, occupation, l, 2J and 2Q.
        bra_shell_codes(2)=orbital_principal_n(orbital_index)
        ket_shell_codes(2)=bra_shell_codes(2)
        bra_shell_codes(3)=(iabs(orbital_kappa(orbital_index))*2)-1
        ket_shell_codes(3)=bra_shell_codes(3)
        bra_shell_codes(4)=bra_occupations(orbital_index)
        ket_shell_codes(4)=ket_occupations(orbital_index)
        bra_shell_codes(5)=(bra_shell_codes(3)+orbital_kappa(orbital_index)/iabs(orbital_kappa(orbital_index)))/2
        ket_shell_codes(5)=bra_shell_codes(5)
        bra_shell_codes(6)=bra_subshell_quantum_codes(3,orbital_index)-1
        ket_shell_codes(6)=ket_subshell_quantum_codes(3,orbital_index)-1
        ! With v the seniority and N the occupation, 2Q = |kappa|-v and 2M_Q = N-|kappa|.
        bra_shell_codes(7)=iabs(orbital_kappa(orbital_index))-bra_subshell_quantum_codes(1,orbital_index)
        ket_shell_codes(7)=iabs(orbital_kappa(orbital_index))-ket_subshell_quantum_codes(1,orbital_index)

        ! Real slots 1:3 store Q, J and the quasispin projection.
        bra_shell_values(1)=half*dble(bra_shell_codes(7))
        ket_shell_values(1)=half*dble(ket_shell_codes(7))
        bra_shell_values(2)=half*dble(bra_shell_codes(6))
        ket_shell_values(2)=half*dble(ket_shell_codes(6))
        bra_shell_values(3)=-half*dble(iabs(orbital_kappa(orbital_index))-bra_shell_codes(4))
        ket_shell_values(3)=-half*dble(iabs(orbital_kappa(orbital_index))-ket_shell_codes(4))
        bra_shell_codes(1)=lookup_jj_subshell_term(bra_shell_codes(7),bra_shell_codes(6),bra_shell_codes(3),ket_shell_codes(4), &
            bra_shell_codes(4))
        ket_shell_codes(1)=lookup_jj_subshell_term(ket_shell_codes(7),ket_shell_codes(6),ket_shell_codes(3),ket_shell_codes(4), &
            bra_shell_codes(4))
        return
    end subroutine pack_jj_subshell_state
end module jj_subshell_state_packing
