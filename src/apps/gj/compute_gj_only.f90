!***********************************************************************
!                                                                      *
subroutine compute_gj_only(gjc_diag, dgjc_diag)
    ! Module dependencies
    !                                                                      *
    !   This routine evaluates only the diagonal ASF contributions needed  *
    !   for g_J, delta g_J, and total g_J.                                 *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: max_orbitals
    use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
    use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
        csf_mixing_coefficients, evecff, evecii
    use first_order_csf_cutoffs,        only: iccut
    use standard_io_units,        only: stderr_unit
    use orbital_csf_basis,           only: csf_count, orbital_count
    use selected_level_indices,          only: nvec
    use integer_format, only: convrt
    use csf_quantum_queries, only: itjpo

    ! Procedure interfaces
    use hyperfine_one_electron_element_interface, only: hyperfine_one_electron_element

    ! Module dependencies
    use one_electron_tensor_coefficients, only: oneparticlejj
    use weighted_radial_overlap, only: rint

    ! Procedure interfaces
    use hyperfine_radial_integral_interface, only: hyperfine_radial_integral
    implicit none (type, external)

    ! Arguments
    real(real64), intent(out) :: gjc_diag(nvec), dgjc_diag(nvec)

    ! Local constants
    real(real64), parameter :: cutoff = 1.0e-10_real64

    ! Local declarations
    integer :: orbital_index, j, operator_parity, csf_index, text_length, ir, itjpoc, itjpor, idiff, bra_orbital_index, &
        & ket_orbital_index, &
        k, loc
    real(real64) :: shell_coefficients(max_orbitals)
    real(real64) :: rintgj(max_orbitals, max_orbitals), rintdgj(max_orbitals, max_orbitals), &
        gjmelt(max_orbitals, max_orbitals), dgjmelt(max_orbitals, max_orbitals)
    real(real64) :: apart, gjpart, dgjpart, elemntgj, elemntdgj, contrgj, &
        contrdgj
    character :: number_text*11
    !-----------------------------------------------
    !
    gjc_diag(:nvec) = 0.0_real64
    dgjc_diag(:nvec) = 0.0_real64
    !
    !   Precompute the radial integrals and angular factors required by
    !   the two one-body operators contributing to the Landé factor.
    !
    do orbital_index = 1, orbital_count
        do j = 1, orbital_count
            rintgj(orbital_index,j) = hyperfine_radial_integral(orbital_index,j,1)
            rintdgj(orbital_index,j) = rint(orbital_index,j,0)
            call hyperfine_one_electron_element(orbital_index, 1, j, apart, gjpart, dgjpart)
            gjmelt(orbital_index,j) = gjpart
            dgjmelt(orbital_index,j) = dgjpart
        end do
    end do
    !
    operator_parity = 1
    !
    !   Sweep over the CSF pairs. Only the diagonal-in-J branch from HFSGG
    !   is required, and only ASF diagonal matrix elements are accumulated.
    !
    do csf_index = 1, csf_count
        if (mod(csf_index,100) == 0) then
            call convrt(csf_index, number_text, text_length)
            write (stderr_unit, *) "Column "//number_text(1:text_length)//" complete;"
        end if
        !
        do ir = 1, csf_count
            if (lfordr .and. csf_index>iccut .and. csf_index/=ir) cycle
            !
            itjpoc = itjpo(csf_index)
            itjpor = itjpo(ir)
            idiff = itjpoc - itjpor
            if (.not.(idiff == 0 .and. ir >= csf_index)) cycle
            !
            elemntgj = 0.0_real64
            elemntdgj = 0.0_real64
            !
            call oneparticlejj(1, operator_parity, csf_index, ir, bra_orbital_index, ket_orbital_index, shell_coefficients)
            if (bra_orbital_index == 0) cycle
            !
            if (bra_orbital_index == ket_orbital_index) then
                do orbital_index = 1, orbital_count
                    if (abs(shell_coefficients(orbital_index)) <= cutoff) cycle
                    elemntgj = elemntgj + gjmelt(orbital_index,orbital_index)*rintgj(orbital_index, &
                        & orbital_index)*shell_coefficients(orbital_index)
                    elemntdgj = elemntdgj + dgjmelt(orbital_index,orbital_index)*rintdgj(orbital_index, &
                        & orbital_index)*shell_coefficients(orbital_index)
                end do
            else
                if (abs(shell_coefficients(1)) <= cutoff) cycle
                elemntgj = elemntgj + gjmelt(bra_orbital_index,ket_orbital_index)*rintgj(bra_orbital_index, &
                    & ket_orbital_index)*shell_coefficients(1)
                elemntdgj = elemntdgj + dgjmelt(bra_orbital_index,ket_orbital_index)*rintdgj(bra_orbital_index, &
                    & ket_orbital_index)*shell_coefficients(1)
            end if
            !
            if (abs(elemntgj) <= cutoff .and. abs(elemntdgj) <= cutoff) cycle
            !
            do k = 1, nvec
                loc = (k - 1)*csf_count
                if (ir /= csf_index) then
                    contrgj = elemntgj*(csf_mixing_coefficients(csf_index + loc)*csf_mixing_coefficients(ir + loc) + &
                        & csf_mixing_coefficients(ir &
                        + loc)*csf_mixing_coefficients(csf_index + loc))
                    contrdgj = elemntdgj*(csf_mixing_coefficients(csf_index + loc)*csf_mixing_coefficients(ir + loc) + &
                        & csf_mixing_coefficients(ir&
                        + loc)*csf_mixing_coefficients(csf_index + loc))
                else
                    contrgj = elemntgj*csf_mixing_coefficients(csf_index + loc)*csf_mixing_coefficients(ir + loc)
                    contrdgj = elemntdgj*csf_mixing_coefficients(csf_index + loc)*csf_mixing_coefficients(ir + loc)
                end if
                gjc_diag(k) = gjc_diag(k) + contrgj
                dgjc_diag(k) = dgjc_diag(k) + contrdgj
            end do
        end do
    end do
    !
    return
end subroutine compute_gj_only
