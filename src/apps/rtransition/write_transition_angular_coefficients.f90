module write_transition_angular_coefficients_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_transition_angular_coefficients
contains

    ! This routine loads  coefficients with parity and  rank specified
    ! by KP(JKP) into the arrays ISLDR and XSLDR.  IOPAR is the parity
    ! (+/- 1) and is determined from the sign of  KP(JKP).
    !
    ! Last revision: 28 Dec 1992
    ! Updated by Jacek Bieron. Last revision: 10 Mar 1994
    subroutine write_transition_angular_coefficients(operator_parity, transition_type_index, name, transition_angular_unit)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: orbital_index_radix, max_orbitals
        use array_memory_management, only: alloc, dalloc, ralloc
        use transition_basis_blocks, only: initial_transition_block_count, final_transition_block_count, &
            initial_cumulative_csf_counts, final_cumulative_csf_counts
        use debug_output_settings,         only: ldbpa
        use orbital_csf_basis, only: max_orbitals, orbital_count
        use transition_off_diagonal_flags, only: different_j_only
        use transition_matrix_workspace, only: transition_multipole_rank, requested_multipole_codes

        ! Procedure interfaces
        use transition_angular_cache, only: transition_cache_key
        use open_transition_angular_cache_interface, only: open_transition_angular_cache

        ! Module dependencies
        use csf_quantum_queries, only: csf_total_j_code
        use one_electron_tensor_coefficients, only: calculate_one_electron_tensor_coefficients

        ! Procedure interfaces
        use sort_transition_angular_coefficients_interface, only: sort_transition_angular_coefficients
        implicit none (type, external)

        ! Arguments
        integer, intent(out) :: transition_angular_unit
        integer, intent(out) :: operator_parity
        integer, intent(in) :: transition_type_index
        character  :: name(2)*24

        ! Local constants
        real(real64), parameter :: cutoff = 1.0e-10_real64

        ! Local declarations
        integer :: coefficient_scratch_unit
        integer :: initial_block_index, final_block_index, nmct, nlabel, first_initial_csf_index, first_final_csf_index, &
            & bra_csf_index, ket_csf_index, ncr&
            , orbital_index, ket_orbital_index, newsiz, i
        type(transition_cache_key) :: key
        logical :: cache_exists, unequal_j
        integer, pointer :: label(:)
        real(real64), pointer:: coeff(:)
        real(real64) :: shell_coefficients(max_orbitals)
        ! NCFI(I): the end position of the Ith block for the initial states in the globle CSF list
        ! NCFF(I): the end position of the Ith block for the final states in the globle CSF list
        !
        !   Check if angular data is available on file
        call open_transition_angular_cache(name, cache_exists, transition_type_index, transition_angular_unit, key)

        ! If angular data is not available open the scratch file to store the
        ! coefficients; position file
        ! to beginning
        transition_multipole_rank = abs(requested_multipole_codes(transition_type_index))
        operator_parity = isign(1,requested_multipole_codes(transition_type_index))
        if (cache_exists) return
        unequal_j = different_j_only(requested_multipole_codes(transition_type_index))
        write (output_unit, *) "LK,IOPAR,from write_transition_angular_coefficients"
        write (output_unit, *) transition_multipole_rank, operator_parity
        ! Start of the block loops
        do initial_block_index = 1, initial_transition_block_count
            do final_block_index = 1, final_transition_block_count
                ! Sometimes, when there has been an error, status need to be "unknown"
                open(newunit=coefficient_scratch_unit, file="fort.93", status="unknown", form="UNFORMATTED")
                nmct = 0
                ! If angular data is not available
                ! Generate MCT coefficients for given rank/parity combination and
                ! store them by CSF on coefficient_scratch_unit
                !
                ! Allocate storage to buffer arrays
                nlabel = 32
                call alloc (label, nlabel, "LABEL", "write_transition_angular_coefficients")
                call alloc (coeff, nlabel, "COEFF", "write_transition_angular_coefficients")

                if (initial_block_index == 1) then
                    first_initial_csf_index = 1
                else
                    first_initial_csf_index = initial_cumulative_csf_counts(initial_block_index-1) + 1
                end if
                if (final_block_index == 1) then
                    first_final_csf_index = initial_cumulative_csf_counts(initial_transition_block_count) + 1
                else
                    first_final_csf_index = initial_cumulative_csf_counts(initial_transition_block_count) + &
                        final_cumulative_csf_counts(final_block_index-1) + 1
                end if
                do bra_csf_index = first_initial_csf_index, initial_cumulative_csf_counts(initial_block_index)
                    do ket_csf_index = first_final_csf_index, initial_cumulative_csf_counts(initial_transition_block_count) + &
                        final_cumulative_csf_counts(final_block_index)
                        ! IR = IC
                        ncr = 0

                        ! In many case one is interested only in M1 and E2 transitions between
                        ! levels with different J values. If this is the case then the do check
                        ! on the J quantum numbers of the CSFs before calling TNSRJJ.
                        if (unequal_j) then
                            if (csf_total_j_code(bra_csf_index) == csf_total_j_code(ket_csf_index)) cycle
                        end if
                        call calculate_one_electron_tensor_coefficients(transition_multipole_rank,operator_parity, &
                            bra_csf_index,ket_csf_index, &
                            & orbital_index,ket_orbital_index, &
                            & shell_coefficients)
                        if (orbital_index /= 0) then
                            if (orbital_index == ket_orbital_index) then
                                do orbital_index = 1, orbital_count
                                    if (abs(shell_coefficients(orbital_index)) <= cutoff) cycle
                                    ncr = ncr + 1
                                    if (ncr > nlabel) then
                                        newsiz = 2*nlabel
                                        call ralloc (label,  newsiz, "LABEL", "write_transition_angular_coefficients")
                                        call ralloc (coeff,  newsiz, "COEFF", "write_transition_angular_coefficients")
                                        nlabel = newsiz
                                    end if
                                    label(ncr) = orbital_index*orbital_index_radix + orbital_index
                                    coeff(ncr) = shell_coefficients(orbital_index)
                                end do
                            else
                                if (abs(shell_coefficients(1)) > cutoff) then
                                    ncr = ncr + 1
                                    if (ncr > nlabel) then
                                        newsiz = 2*nlabel
                                        call ralloc (label,  newsiz, "LABEL", "write_transition_angular_coefficients")
                                        call ralloc (coeff,  newsiz, "COEFF", "write_transition_angular_coefficients")
                                        nlabel = newsiz
                                    end if
                                    label(ncr) = orbital_index*orbital_index_radix + ket_orbital_index
                                    coeff(ncr) = shell_coefficients(1)
                                end if
                            end if
                        end if
                        if (ncr <= 0) cycle
                        write (coefficient_scratch_unit) bra_csf_index - first_initial_csf_index + 1, ket_csf_index - &
                            & first_final_csf_index + 1, ncr
                        write (coefficient_scratch_unit) (label(i),coeff(i),i=1,ncr)
                        nmct = nmct + ncr

                    end do
                end do
                ! Deallocate storage for buffer arrays
                call dalloc (label,  "LABEL", "write_transition_angular_coefficients")
                call dalloc (coeff,  "COEFF", "write_transition_angular_coefficients")
                write (*, '(/,/,/,1x,i8," MCT coefficients generated for rank ",i2," and parity ",i2,/,/)') &
                    nmct, transition_multipole_rank, operator_parity
                ! Sort the MCT coefficients by integral labels
                call sort_transition_angular_coefficients(coefficient_scratch_unit, transition_angular_unit, ldbpa(2), &
                    initial_block_index, final_block_index)
                close(coefficient_scratch_unit, status="delete")
                ! end of the loops for blocks
            end do
        end do
        ! Read the data back as required by OSCL conventions
        write (transition_angular_unit) key
        rewind (transition_angular_unit)
        return
        return
    end subroutine write_transition_angular_coefficients
end module write_transition_angular_coefficients_interface
