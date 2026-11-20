!***********************************************************************
!                                                                      *
subroutine accumulate_closed_shell_breit(bra_csf_index, ket_csf_index, ja1, ipca, jb1)
    ! Module dependencies
    use debug_output_settings, only: angular_debug_unit, ensure_debug_output
    !                                                                      *
    !   Computes closed shell contributions - aaaa and exchange only.      *
    !                                                                      *
    !   Call(s) to: [LIB92]: CLRX, CXK, ITRIG, TALK, SNRC.                 *
    !                                                                      *
    !                                           LAST UPDATE: 09 Oct 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use breit_core_occupations, only: icore, max_orbitals
    use numeric_constants, only: eleven, eps, four, half, one, seven, ten, tenth, three, two, zero
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use csf_pair_angular_state, only: bra_intermediate_couplings => jjc1, ket_intermediate_couplings => jjc2, &
        & bra_subshell_quantum_codes => jjq1, ket_subshell_quantum_codes => jjq2, peel_orbital_indices => jlist, &
        & core_orbital_indices => klist, max_orbitals, core_orbital_count => ncore, peel_orbital_count => npeel, &
        & bra_occupations => nq1, ket_occupations => nq2
    use orbital_csf_basis, only: orbital_principal_n, orbital_kappa
    use breit_radial_rank_ranges, only: snrc
    use relativistic_three_j, only: clrx
    use breit_angular_coefficient_buffer, only: talk
    use angular_momentum_triangle, only: itrig
    use breit_radial_integral_coefficients, only: cxk
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: bra_csf_index
    integer, intent(in) :: ket_csf_index
    integer, intent(in) :: ja1
    integer, intent(in) :: ipca
    integer, intent(in) :: jb1

    ! Local constants
    integer, parameter :: max_radial_order_count = 20

    ! Local declarations
    integer :: orbital_indices(4), twice_orbital_kappa(4), two_j_plus_one(4)
    integer :: bra_orbital_index, ket_orbital_index, isg, bra_orbital_occupation, ket_orbital_occupation, i, first_direct_nu, &
        & direct_nu_count, first_exchange_nu, exchange_nu_count, direct_integral_type, exchange_integral_type, &
        n, radial_integral_order, k, kap1, itype, mu, ip, ipp, kk, nup1
    real(real64) :: cone(7, 20)
    real(real64) :: radial_integral_coefficients(12)
    real(real64) :: const, gam, dksks, dnunu1, coef, proc, prod
    !-----------------------------------------------
    !
    !   1.0  Initialization
    !
    if (ipca == 2) then
        bra_orbital_index = core_orbital_indices(ja1)
    else
        bra_orbital_index = peel_orbital_indices(ja1)
    end if
    ket_orbital_index = core_orbital_indices(jb1)
    !
    isg = 1
    if (bra_csf_index == ket_csf_index) then
        if (icore(bra_orbital_index)/=0 .and. icore(ket_orbital_index)/=0) then
            if (bra_csf_index > 1) return
            isg = -1
        end if
    end if
    !
    orbital_indices(1) = bra_orbital_index
    orbital_indices(2) = ket_orbital_index
    orbital_indices(3) = bra_orbital_index
    orbital_indices(4) = ket_orbital_index
    bra_orbital_occupation = bra_occupations(bra_orbital_index)
    ket_orbital_occupation = ket_occupations(ket_orbital_index)
    do i = 1, 4
        twice_orbital_kappa(i) = 2*orbital_kappa(orbital_indices(i))
        two_j_plus_one(i) = iabs(twice_orbital_kappa(i))
    end do
    const = bra_orbital_occupation*ket_orbital_occupation
    if (ibug2 /= 0) then
        call ensure_debug_output()
        write (angular_debug_unit, '("BREID: orbitals ",2i3)') bra_orbital_index, ket_orbital_index
    end if
    !
    !   2.0  Set range of tensor indices
    !
    call snrc (orbital_indices, twice_orbital_kappa, two_j_plus_one, first_direct_nu, direct_nu_count, first_exchange_nu, &
        & exchange_nu_count, direct_integral_type, exchange_integral_type)
    if (ibug2 /= 0) then
        call ensure_debug_output()
        write (angular_debug_unit, '(2x,"ND1 ND2 NE1 NE2 IBRD IBRE ",6i5)') first_direct_nu, direct_nu_count, &
            & first_exchange_nu, exchange_nu_count, direct_integral_type, exchange_integral_type
    end if
    if (bra_orbital_index == ket_orbital_index) then
        !
        !   3.0 Calculate aaaa interaction
        !
        do n = 1, direct_nu_count
            radial_integral_order = first_direct_nu + 2*(n - 1)
            k = radial_integral_order
            if (mod(k,2) /= 1) return
            kap1 = twice_orbital_kappa(1)/2
            gam = clrx(kap1,radial_integral_order,kap1)
            dksks = two_j_plus_one(1)*two_j_plus_one(1)
            dnunu1 = radial_integral_order*(radial_integral_order + 1)
            coef = const*two*dksks*gam*gam/dnunu1
            if (ibug2 /= 0) then
                call ensure_debug_output()
                write (angular_debug_unit, '(2x,"aaaa contribution: NU,GAM,COEF",i5,2(3x,1p,d15.8))') radial_integral_order, &
                    & gam, coef
            end if
            itype = isg*4
            call talk (bra_csf_index, ket_csf_index, radial_integral_order, bra_orbital_index, bra_orbital_index, &
                & bra_orbital_index, bra_orbital_index, itype, coef)
        end do
        return
    end if
    !
    !   Calculate exchange interactions
    !
    if (exchange_integral_type < 0) return
    if (exchange_nu_count > max_radial_order_count) then
        write (*, '("BREID: Dimension error for NUMAX.")')
        stop
    end if
    !
    cone(:,:exchange_nu_count) = zero
    !
    proc = -const/dble(two_j_plus_one(1)*two_j_plus_one(2))
    !
    !   Negative sign arises from Pauli phase factor
    !
    do n = 1, exchange_nu_count
        radial_integral_order = first_exchange_nu + 2*(n - 1)
        k = radial_integral_order
        ip = (two_j_plus_one(1)-two_j_plus_one(2))/2 + k
        ipp = ip + 1
        if (radial_integral_order /= 0) then
            kk = k + k + 1
            if (itrig(two_j_plus_one(1),two_j_plus_one(2),kk) /= 0) then
                prod = proc
                if (mod(ip,2) /= 0) prod = -prod
                call cxk (radial_integral_coefficients, orbital_indices, twice_orbital_kappa, radial_integral_order, k, &
                    & exchange_integral_type, 2)
                if (ibug2 /= 0) then
                    call ensure_debug_output()
                    write (angular_debug_unit, '(2x,"PROD = ",1p,d15.8,/," S",7d15.8)') prod, (radial_integral_coefficients(mu), &
                        & mu=1,3)
                end if
                cone(:3,n) = cone(:3,n) + prod*radial_integral_coefficients(:3)
            end if
            !
            k = radial_integral_order - 1
            kk = k + k + 1
            if (itrig(two_j_plus_one(1),two_j_plus_one(2),kk) /= 0) then
                prod = proc
                if (mod(ipp,2) /= 0) prod = -prod
                call cxk (radial_integral_coefficients, orbital_indices, twice_orbital_kappa, radial_integral_order, k, &
                    & exchange_integral_type, 2)
                if (ibug2 /= 0) then
                    call ensure_debug_output()
                    write (angular_debug_unit, '(2x,"PROD = ",1p,d15.8,/," S",7d15.8)') prod, (radial_integral_coefficients(mu), &
                        & mu=1,3)
                end if
                cone(:3,n) = cone(:3,n) + prod*radial_integral_coefficients(:3)
                !
            end if
        end if
        if (n == exchange_nu_count) exit
        k = radial_integral_order + 1
        kk = k + k + 1
        prod = proc
        if (mod(ipp,2) /= 0) prod = -prod
        call cxk (radial_integral_coefficients, orbital_indices, twice_orbital_kappa, radial_integral_order, k, &
            & exchange_integral_type, 2)
        if (ibug2 /= 0) then
            call ensure_debug_output()
            write (angular_debug_unit, '(2x,"PROD = ",1p,d15.8,/," S",7d15.8)') prod, (radial_integral_coefficients(mu),mu=1,7)
        end if
        cone(:,n) = cone(:,n) + prod*radial_integral_coefficients(:7)
    end do
    !
    !   4.0  Output results
    !
    do n = 1, exchange_nu_count
        radial_integral_order = first_exchange_nu + 2*(n - 1)
        itype = isg*5
        call talk (bra_csf_index, ket_csf_index, radial_integral_order, ket_orbital_index, bra_orbital_index, &
            & ket_orbital_index, bra_orbital_index, itype, cone(1,n))
        call talk (bra_csf_index, ket_csf_index, radial_integral_order, bra_orbital_index, ket_orbital_index, &
            & ket_orbital_index, bra_orbital_index, itype, cone(2,n))
        call talk (bra_csf_index, ket_csf_index, radial_integral_order, bra_orbital_index, ket_orbital_index, &
            & bra_orbital_index, ket_orbital_index, itype, cone(3,n))
        if (n == exchange_nu_count) cycle
        nup1 = radial_integral_order + 1
        itype = isg*6
        call talk (bra_csf_index, ket_csf_index, nup1, bra_orbital_index, ket_orbital_index, bra_orbital_index, &
            & ket_orbital_index, itype, cone(4,n))
        call talk (bra_csf_index, ket_csf_index, nup1, ket_orbital_index, bra_orbital_index, ket_orbital_index, &
            & bra_orbital_index, itype, cone(5,n))
        call talk (bra_csf_index, ket_csf_index, nup1, bra_orbital_index, ket_orbital_index, ket_orbital_index, &
            & bra_orbital_index, itype, cone(6,n))
        call talk (bra_csf_index, ket_csf_index, nup1, ket_orbital_index, bra_orbital_index, bra_orbital_index, &
            & ket_orbital_index, itype, cone(7,n))
    end do
    return
    !
    return
    !
end subroutine accumulate_closed_shell_breit
