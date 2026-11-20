! Query occupations, subshell and intermediate angular momenta, and total CSF symmetry.
module csf_quantum_queries
    use standard_io_units, only: stderr_unit
    use orbital_csf_basis, only: csf_count, csf_subshell_occupations
    use grasp_parameters, only: max_orbitals
    use csf_coupling_data, only: csf_coupling_quantum_codes, csf_subshell_quantum_codes
    implicit none (type, external)
    private
    public :: iq, jqs, jcup, ispar, itjpo
contains

    ! Return the stored occupation as a default-kind integer.
    integer function iq(subshell_index, csf_index)
        integer, intent(in) :: subshell_index, csf_index

        iq = csf_subshell_occupations(subshell_index,csf_index)
    end function iq

    ! Selectors 1, 2, and 3 return seniority, w, and subshell 2J+1, respectively.
    integer function jqs(quantum_number_selector, subshell_index, csf_index)
        integer, intent(in) :: quantum_number_selector, subshell_index, csf_index

        jqs = csf_subshell_quantum_codes(subshell_index,quantum_number_selector,csf_index)
    end function jqs

    ! Read 2J+1 for the nontrivial intermediate coupling at coupling_index.
    integer function jcup(coupling_index, csf_index)
        integer, intent(in) :: coupling_index, csf_index

        if (coupling_index>=1 .and. coupling_index<=max_orbitals-1) then
            if (csf_index>=1 .and. csf_index<=csf_count) then
                jcup = csf_coupling_quantum_codes(coupling_index,csf_index)
            else
                write (stderr_unit, *) "JCUP: Argument ICSF is out of range."
                stop
            end if
        else
            write (stderr_unit, *) "JCUP: Argument LOC is out of range."
            stop
        end if
    end function jcup

    ! Decode parity from the sign of the packed total symmetry.
    integer function ispar(csf_index)
        integer, intent(in) :: csf_index

        if (csf_index>=1 .and. csf_index<=csf_count) then
            ispar = csf_coupling_quantum_codes(max_orbitals,csf_index)
            if (ispar > 127) ispar = ispar - 256
            ispar = sign(1,ispar)
        else
            write (stderr_unit, *) "ISPAR: Argument ICSF is out of range."
            write (stderr_unit, *) "ICSF =",csf_index," NCF =",csf_count
            stop
        end if
    end function ispar

    ! Decode 2J+1 from the magnitude of the packed total symmetry.
    integer function itjpo(csf_index)
        integer, intent(in) :: csf_index

        if (csf_index>=1 .and. csf_index<=csf_count) then
            itjpo = csf_coupling_quantum_codes(max_orbitals,csf_index)
            if (itjpo > 127) itjpo = 256 - itjpo
            itjpo = iabs(itjpo)
        else
            write (stderr_unit, *) "ITJPO: Argument ICSF is out of range."
            stop
        end if
    end function itjpo
end module csf_quantum_queries
