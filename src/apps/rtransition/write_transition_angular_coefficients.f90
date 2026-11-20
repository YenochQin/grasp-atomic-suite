!***********************************************************************
!                                                                      *
subroutine write_transition_angular_coefficients(operator_parity, transition_type_index, name, transition_angular_unit)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    !                                                                      *
    !   This routine loads  coefficients with parity and  rank specified   *
    !   by KP(JKP) into the arrays ISLDR and XSLDR.  IOPAR is the parity   *
    !   (+/- 1) and is determined from the sign of  KP(JKP).               *
    !                                                                      *
    !                                         Last revision: 28 Dec 1992   *
    !   Updated by Jacek Bieron               Last revision: 10 Mar 1994   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: orbital_index_radix, max_orbitals
    use array_memory_management, only: alloc, dalloc, ralloc
    use transition_basis_blocks, only: ncfblk, nevinblk, ncfinblk, two_j, nblock, idxblk, nelectot, ncftot, nwtot, nvectot, &
        nvecsiztot, nblock1, nblocki, nblockf, ncfi, ncff
    use debug_output_settings,         only: ldbpa
    use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
    use first_order_csf_cutoffs, only: iccut
    use orbital_csf_basis, only: max_orbitals, orbital_count
    use transition_off_diagonal_flags, only: noffd1, noffd2
    use transition_matrix_workspace, only: nsdim, hb1, hb2, hc1, hc2, hm1, hm2, jja, jjb, lk, kk, ntdim, &
        xsldr, totc, totb, isldr, isldr1, nint, ninteg, nptr, lab, nkp, kp, ltc

    ! Procedure interfaces
    use open_transition_angular_cache_interface, only: open_transition_angular_cache

    ! Module dependencies
    use csf_quantum_queries, only: itjpo
    use one_electron_tensor_coefficients, only: oneparticlejj

    ! Procedure interfaces
    use sort_transition_angular_coefficients_interface, only: sort_transition_angular_coefficients
    implicit none (type, external)

    ! Arguments
    integer, intent(out) :: transition_angular_unit
    integer  :: operator_parity
    integer, intent(in) :: transition_type_index
    character  :: name(2)*24

    ! Local constants
    real(real64), parameter :: cutoff = 1.0e-10_real64

    ! Local declarations
    integer :: coefficient_scratch_unit
    integer :: iblki, iblkf, nmct, nlabel, ncfi0, ncff0, bra_csf_index, ket_csf_index, ncr&
        , orbital_index, ket_orbital_index, newsiz, i
    logical :: cache_exists
    integer, pointer :: label(:)
    real(real64), pointer:: coeff(:)
    real(real64) :: shell_coefficients(max_orbitals)
    !-----------------------------------------------
    !
    !
    !
    ! NCFI(I): the end position of the Ith block for the initial states in the globle CSF list
    ! NCFF(I): the end position of the Ith block for the final states in the globle CSF list
    !
    !   Check if angular data is available on file
    !
    call open_transition_angular_cache (name, cache_exists, transition_type_index, transition_angular_unit)

    !
    !   If angular data is not available open the scratch file to store the
    !   coefficients; position file
    !   to beginning
    !
    !
    lk = abs(kp(transition_type_index))
    operator_parity = isign(1,kp(transition_type_index))
    if (cache_exists) return
    write (output_unit, *) "LK,IOPAR,from MCTOUT"
    write (output_unit, *) lk, operator_parity
    !
    ! Start of the block loops
    do iblki = 1, nblocki
        do iblkf = 1, nblockf
            !        OPEN (coefficient_scratch_unit,STATUS = 'new', FORM = 'UNFORMATTED')
            !        Sometimes, when there has been an error, status need to be "unknown"
            open(newunit=coefficient_scratch_unit, file="fort.93", status="unknown", form="UNFORMATTED")
            nmct = 0
            !
            !
            !   If angular data is not available
            !   Generate MCT coefficients for given rank/parity combination and
            !   store them by CSF on coefficient_scratch_unit
            !
            !
            !   Allocate storage to buffer arrays
            !
            nlabel = 32
            call alloc (label, nlabel, "LABEL", "MCTOUT")
            call alloc (coeff, nlabel, "COEFF", "MCTOUT")

            if (iblki == 1) then
                ncfi0 = 1
            else
                ncfi0 = ncfi(iblki-1) + 1
            end if
            !
            if (iblkf == 1) then
                ncff0 = ncfi(nblocki) + 1
            else
                ncff0 = ncfi(nblocki) + ncff(iblkf-1) + 1
            end if
            do bra_csf_index = ncfi0, ncfi(iblki)
                do ket_csf_index = ncff0, ncfi(nblocki) + ncff(iblkf)
                    !
                    !           IR = IC
                    !
                    ncr = 0

                    !
                    !   In many case one is interested only in M1 and E2 transitions between
                    !   levels with different J values. If this is the case then the do check
                    !   on the J quantum numbers of the CSFs before calling TNSRJJ.
                    !
                    if (kp(transition_type_index)==1 .and. noffd1==1) then
                        if (itjpo(bra_csf_index) == itjpo(ket_csf_index)) cycle
                    end if
                    if (kp(transition_type_index)==2 .and. noffd2==1) then
                        if (itjpo(bra_csf_index) == itjpo(ket_csf_index)) cycle
                    end if
                    call oneparticlejj(lk,operator_parity,bra_csf_index,ket_csf_index,orbital_index,ket_orbital_index, &
                        & shell_coefficients)
                    if (orbital_index /= 0) then
                        if (orbital_index == ket_orbital_index) then
                            do orbital_index = 1, orbital_count
                                if (abs(shell_coefficients(orbital_index)) <= cutoff) cycle
                                ncr = ncr + 1
                                if (ncr > nlabel) then
                                    newsiz = 2*nlabel
                                    call ralloc (label,  newsiz, "LABEL", "MCTOUT")
                                    call ralloc (coeff,  newsiz, "COEFF", "MCTOUT")
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
                                    call ralloc (label,  newsiz, "LABEL", "MCTOUT")
                                    call ralloc (coeff,  newsiz, "COEFF", "MCTOUT")
                                    nlabel = newsiz
                                end if
                                label(ncr) = orbital_index*orbital_index_radix + ket_orbital_index
                                coeff(ncr) = shell_coefficients(1)
                            end if
                        end if
                    end if
                    if (ncr <= 0) cycle
                    write (coefficient_scratch_unit) bra_csf_index - ncfi0 + 1, ket_csf_index - ncff0 + 1, ncr
                    !             WRITE (coefficient_scratch_unit) IC,IR,NCR
                    write (coefficient_scratch_unit) (label(i),coeff(i),i=1,ncr)
                    nmct = nmct + ncr

                    !
                end do
            end do
            !
            !   Deallocate storage for buffer arrays
            !
            call dalloc (label,  "LABEL", "MCTOUT")
            call dalloc (coeff,  "COEFF", "MCTOUT")
            !
            write (*, '(/,/,/,1x,i8," MCT coefficients generated for rank ",i2," and parity ",i2,/,/)') &
                nmct, lk, operator_parity
            !
            !   Sort the MCT coefficients by integral labels
            !
            call sort_transition_angular_coefficients (name, coefficient_scratch_unit, transition_angular_unit, ldbpa(2), &
                & transition_type_index, iblki, iblkf)
            close(coefficient_scratch_unit, status="delete")
            ! end of the loops for blocks
        end do
    end do
    !
    !   Read the data back as required by OSCL conventions
    !
    rewind (transition_angular_unit)
    return
    !
    return
    !
end subroutine write_transition_angular_coefficients
