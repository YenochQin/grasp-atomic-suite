! Classify Breit orbital pairings and find their admissible radial ranks.
module breit_radial_rank_ranges
    implicit none (type, external)

    ! Module interface
    private
    public :: snrc
contains
    subroutine snrc(integral_orbitals, two_kappas, two_j_plus_one, first_direct_nu, direct_nu_count, &
        first_exchange_nu, exchange_nu_count, direct_integral_type, exchange_integral_type)
        implicit none (type, external)

        ! Arguments
        integer, intent(out) :: first_direct_nu
        integer, intent(out) :: direct_nu_count
        integer, intent(out) :: first_exchange_nu
        integer, intent(out) :: exchange_nu_count
        integer, intent(out) :: direct_integral_type
        integer, intent(out) :: exchange_integral_type
        integer, intent(in) :: integral_orbitals(4)
        integer, intent(in) :: two_kappas(4)
        integer, intent(in) :: two_j_plus_one(4)

        ! Local declarations
        integer :: relative_kappa_sign_1, relative_kappa_sign_2, second_direct_lower_nu, second_direct_upper_nu, &
            second_exchange_lower_nu, second_exchange_upper_nu

        ! Direct terms pair orbitals (1,3) and (2,4). Exchange terms pair (1,4) and (2,3).
        ! The allowed orders are first_nu + 2*(i-1), where i runs through the returned count.
        direct_nu_count = 0
        exchange_nu_count = 0
        relative_kappa_sign_1 = 1
        if (two_kappas(1)*two_kappas(3) < 0) relative_kappa_sign_1 = -1
        relative_kappa_sign_2 = 1
        if (two_kappas(2)*two_kappas(4) < 0) relative_kappa_sign_2 = -1
        first_direct_nu = abs(two_j_plus_one(1)-two_j_plus_one(3))/2 - 1
        if (relative_kappa_sign_1 == (-1)) first_direct_nu = first_direct_nu + 1
        if (first_direct_nu == (-1)) first_direct_nu = 1
        second_direct_lower_nu = abs(two_j_plus_one(2)-two_j_plus_one(4))/2 - 1
        if (relative_kappa_sign_2 == (-1)) second_direct_lower_nu = second_direct_lower_nu + 1
        if (second_direct_lower_nu == (-1)) second_direct_lower_nu = 1

        ! Two rank ranges intersect only if their lower bounds have equal parity.
        if (mod(first_direct_nu - second_direct_lower_nu,2) /= 0) then
            direct_integral_type = -1
        else
            direct_nu_count = abs(two_j_plus_one(1)+two_j_plus_one(3))/2
            if (relative_kappa_sign_1 == (-1)) direct_nu_count = direct_nu_count + 1
            second_direct_upper_nu = abs(two_j_plus_one(2)+two_j_plus_one(4))/2
            if (relative_kappa_sign_2 == (-1)) second_direct_upper_nu = second_direct_upper_nu + 1
            first_direct_nu = max(first_direct_nu,second_direct_lower_nu)
            direct_nu_count = min(direct_nu_count,second_direct_upper_nu)
            direct_nu_count = (direct_nu_count - first_direct_nu)/2 + 1
            direct_integral_type = 1
            if (integral_orbitals(1)==integral_orbitals(3) .and. integral_orbitals(2)/=integral_orbitals(4) .or. &
                integral_orbitals(1)/=integral_orbitals(3) .and. integral_orbitals(2)==integral_orbitals(&
                4)) direct_integral_type = 2
            if (integral_orbitals(1)==integral_orbitals(3) .and. integral_orbitals(2)==integral_orbitals(4)) &
                direct_integral_type = 3
        end if

        ! The exchange channel is unnecessary when either same-side orbital pair coincides.
        if (integral_orbitals(1)==integral_orbitals(2) .or. integral_orbitals(3)==integral_orbitals(4)) then
            exchange_integral_type = -1
            return
        end if
        relative_kappa_sign_1 = 1
        if (two_kappas(1)*two_kappas(4) < 0) relative_kappa_sign_1 = -1
        relative_kappa_sign_2 = 1
        if (two_kappas(2)*two_kappas(3) < 0) relative_kappa_sign_2 = -1
        first_exchange_nu = iabs(two_j_plus_one(1)-two_j_plus_one(4))/2 - 1
        if (relative_kappa_sign_1 == (-1)) first_exchange_nu = first_exchange_nu + 1
        if (first_exchange_nu == (-1)) first_exchange_nu = 1
        second_exchange_lower_nu = abs(two_j_plus_one(2)-two_j_plus_one(3))/2 - 1
        if (relative_kappa_sign_2 == (-1)) second_exchange_lower_nu = second_exchange_lower_nu + 1
        if (second_exchange_lower_nu == (-1)) second_exchange_lower_nu = 1
        if (mod(first_exchange_nu - second_exchange_lower_nu,2) /= 0) then
            exchange_integral_type = -1
            return
        end if
        exchange_nu_count = abs(two_j_plus_one(1)+two_j_plus_one(4))/2
        if (relative_kappa_sign_1 == (-1)) exchange_nu_count = exchange_nu_count + 1
        second_exchange_upper_nu = abs(two_j_plus_one(2)+two_j_plus_one(3))/2
        if (relative_kappa_sign_2 == (-1)) second_exchange_upper_nu = second_exchange_upper_nu + 1
        first_exchange_nu = max(first_exchange_nu,second_exchange_lower_nu)
        exchange_nu_count = min(exchange_nu_count,second_exchange_upper_nu)
        exchange_nu_count = (exchange_nu_count - first_exchange_nu)/2 + 1
        exchange_integral_type = 1
        if (integral_orbitals(1)==integral_orbitals(4) .and. integral_orbitals(2)/=integral_orbitals(3) .or. &
            integral_orbitals(1)/=integral_orbitals(4) .and. integral_orbitals(2)==integral_orbitals(3)&
            ) exchange_integral_type = 2
        if (integral_orbitals(1)==integral_orbitals(3) .and. integral_orbitals(2)==integral_orbitals(4)) &
            exchange_integral_type = 4
        return
    end subroutine snrc
end module breit_radial_rank_ranges
