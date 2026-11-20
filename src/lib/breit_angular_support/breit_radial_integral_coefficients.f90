module breit_radial_integral_coefficients
    implicit none (type, external)

    private
    public :: calculate_breit_radial_integral_coefficients

contains

    !***********************************************************************
    !                                                                      *
    subroutine calculate_breit_radial_integral_coefficients(radial_integral_coefficients, orbital_indices, &
        twice_orbital_kappa, radial_integral_order, tensor_rank, radial_integral_case, direct_exchange_code)
        ! Module dependencies
        !                                                                      *
        !   Computes  the  coefficients of radial integrals in the expansion   *
        !   of the effective interaction strength: X(K,IA1,IB1,IA2,IB2).       *
        !                                                                      *
        !   Input variables:                                                   *
        !                                                                      *
        !      orbital_indices       Four orbital labels                      *
        !      twice_orbital_kappa   Values of 2*kappa                         *
        !      radial_integral_order Order of the radial integral (NU)         *
        !      tensor_rank           Rank of the tensor operator (K)           *
        !      direct_exchange_code 1 for direct, 2 for exchange terms        *
        !      radial_integral_case  Type of radial integral (IBR below)       *
        !            cases:                                                    *
        !            IBR = 1 A. All states distinct                            *
        !                    B. ((IA .EQ. IB) .AND. (IC .NE. ID)), or          *
        !                       ((IA .NE. IB) .AND. (IC .EQ. ID))              *
        !                    These give 12 distinct  radial  integrals, with   *
        !                    values of K and NU limited only by angular mom-   *
        !                    entum and parity                                  *
        !            IBR = 2 ((IA .EQ. IC) .AND. (IB .NE. ID)) or              *
        !                    ((IA .NE. IC) .AND. (IB .EQ. ID))                 *
        !                    This case gives one non-zero integral when K =    *
        !                    NU is ODD                                         *
        !            IBR = 3 ((IA .EQ. IC) .AND. (IB .EQ. ID)) AND             *
        !                    (IA .NE. IB)                                      *
        !                    Integrals of magnetic F-type when K = NU is odd   *
        !            IBR = 4 ((IA .EQ. ID) .AND. (IB .EQ. IC)) gives 3  mag-   *
        !                    netic G-type integrals and  four  H-TYPE  inte-   *
        !                    grals                                             *
        !                                                                      *
        !   Output:                                                            *
        !                                                                      *
        !      radial_integral_coefficients   12 coefficient slots            *
        !                                                                      *
        !                                                                      *
        !   Call(s) to: [LIB92] CRE.                                           *
        !                                                                      *
        !                                           Last update: 09 Oct 1992   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use reduced_angular_matrix_element, only: reduced_angular_element
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: radial_integral_order
        integer, intent(in) :: tensor_rank
        integer, intent(in) :: radial_integral_case
        integer, intent(in) :: direct_exchange_code
        integer, intent(in) :: orbital_indices(4)
        integer, intent(in) :: twice_orbital_kappa(4)
        real(real64), intent(inout) :: radial_integral_coefficients(12)

        ! Local declarations
        integer :: orbital_a, orbital_b, orbital_c, orbital_d, kappa_a, kappa_b, kappa_c, kappa_d, kappa_swap, orbital_swap, &
            & phase_exponent, two_rank_plus_one
        real(real64) :: rank_norm, reduced_angular_product, first_kappa_difference, second_kappa_difference, &
            rank, next_rank, g1, g2, g3, g4, leading_prefactor, f1, f2, f3, f4, mixed_prefactor, exchange_kappa_difference
        !-----------------------------------------------
        !
        !
        !   1.0  Initialization
        !
        radial_integral_coefficients = 0.0_real64
        !
        orbital_a = orbital_indices(1)
        orbital_b = orbital_indices(2)
        orbital_c = orbital_indices(3)
        orbital_d = orbital_indices(4)
        kappa_a = twice_orbital_kappa(1)/2
        kappa_b = twice_orbital_kappa(2)/2
        kappa_c = twice_orbital_kappa(3)/2
        kappa_d = twice_orbital_kappa(4)/2
        if (direct_exchange_code == 2) then
            kappa_swap = kappa_d
            orbital_swap = orbital_d
            kappa_d = kappa_c
            orbital_d = orbital_c
            kappa_c = kappa_swap
            orbital_c = orbital_swap
        end if
        select case (radial_integral_case)
        case default
            call stop_invalid_orbital_combination
        !
        !   2.0  IBR = 1 --- The general case
        !
        case (1)
            if (radial_integral_order - tensor_rank >= 0) then
                if (radial_integral_order - tensor_rank <= 0) then
                    radial_integral_coefficients(1) = -(kappa_a + kappa_c)*(kappa_d + kappa_b)
                    if (tensor_rank == 0) then
                        write (*, '("CXK: Illegal value K = 0 -"/,1x,4i3,2x,i3,2x,2i2)') &
                            orbital_indices(1), orbital_indices(2), orbital_indices(3), orbital_indices(4), &
                            radial_integral_order, radial_integral_case, direct_exchange_code
                        stop
                    end if
                    rank_norm = tensor_rank*(tensor_rank + 1)
                    reduced_angular_product = reduced_angular_element(kappa_a,tensor_rank, &
                        & kappa_c)*reduced_angular_element(kappa_b,tensor_rank,kappa_d)
                    if (mod(tensor_rank,2) /= 0) reduced_angular_product = -reduced_angular_product
                    radial_integral_coefficients(1) = radial_integral_coefficients(1)*reduced_angular_product/rank_norm
                    radial_integral_coefficients(2:4) = radial_integral_coefficients(1)
                    return
                end if
                !
                !   2.2  NU = K+1
                !
                first_kappa_difference = kappa_c - kappa_a
                second_kappa_difference = kappa_d - kappa_b
                rank = tensor_rank
                next_rank = tensor_rank + 1
                g1 = first_kappa_difference - next_rank
                g2 = first_kappa_difference + next_rank
                g3 = second_kappa_difference - next_rank
                g4 = second_kappa_difference + next_rank
                two_rank_plus_one = tensor_rank + tensor_rank + 1
                reduced_angular_product = reduced_angular_element(kappa_a,tensor_rank,kappa_c)*reduced_angular_element(kappa_b, &
                    & tensor_rank,kappa_d)
                if (mod(tensor_rank,2) /= 0) reduced_angular_product = -reduced_angular_product
                leading_prefactor = reduced_angular_product*rank/next_rank/dble(two_rank_plus_one*(two_rank_plus_one + 2))
                radial_integral_coefficients(1) = leading_prefactor*g1*g3
                radial_integral_coefficients(2) = leading_prefactor*g2*g4
                radial_integral_coefficients(3) = leading_prefactor*g1*g4
                radial_integral_coefficients(4) = leading_prefactor*g2*g3
                return
            end if
            !
            !   2.2  NU = K-1
            !
            first_kappa_difference = kappa_c - kappa_a
            second_kappa_difference = kappa_d - kappa_b
            rank = tensor_rank
            next_rank = tensor_rank + 1
            f1 = first_kappa_difference - rank
            f2 = first_kappa_difference + rank
            f3 = second_kappa_difference - rank
            f4 = second_kappa_difference + rank
            g1 = first_kappa_difference - next_rank
            g2 = first_kappa_difference + next_rank
            g3 = second_kappa_difference - next_rank
            g4 = second_kappa_difference + next_rank
            two_rank_plus_one = tensor_rank + tensor_rank + 1
            reduced_angular_product = reduced_angular_element(kappa_a,tensor_rank,kappa_c)*reduced_angular_element(kappa_b, &
                & tensor_rank,kappa_d)
            if (mod(tensor_rank,2) /= 0) reduced_angular_product = -reduced_angular_product
            leading_prefactor = reduced_angular_product*next_rank/rank/dble(two_rank_plus_one*(two_rank_plus_one - 2))
            radial_integral_coefficients(1) = leading_prefactor*f2*f4
            radial_integral_coefficients(2) = leading_prefactor*f1*f3
            radial_integral_coefficients(3) = leading_prefactor*f2*f3
            radial_integral_coefficients(4) = leading_prefactor*f1*f4
            mixed_prefactor = reduced_angular_product/dble(two_rank_plus_one*two_rank_plus_one)
            radial_integral_coefficients(5) = mixed_prefactor*f2*g3
            radial_integral_coefficients(6) = mixed_prefactor*f4*g1
            radial_integral_coefficients(7) = mixed_prefactor*f1*g4
            radial_integral_coefficients(8) = mixed_prefactor*f3*g2
            radial_integral_coefficients(9) = mixed_prefactor*f2*g4
            radial_integral_coefficients(10) = mixed_prefactor*f3*g1
            radial_integral_coefficients(11) = mixed_prefactor*f1*g3
            radial_integral_coefficients(12) = mixed_prefactor*f4*g2
            return
        !
        !   3.0  IBR = 2  Degenerate case: only one non-zero R-integral
        !
        case (2)
            if (orbital_a/=orbital_c .or. orbital_b==orbital_d) then
                if (orbital_a==orbital_c .or. orbital_b/=orbital_d) call stop_invalid_orbital_combination
                !
                orbital_swap = orbital_b
                orbital_b = orbital_a
                orbital_a = orbital_swap
                orbital_swap = orbital_d
                orbital_d = orbital_c
                orbital_c = orbital_swap
                !
                kappa_swap = kappa_b
                kappa_b = kappa_a
                kappa_a = kappa_swap
                kappa_swap = kappa_d
                kappa_d = kappa_c
                kappa_c = kappa_swap
            end if
            !
            if (mod(tensor_rank,2) /= 1) return
            rank_norm = tensor_rank*(tensor_rank + 1)
            reduced_angular_product = reduced_angular_element(kappa_a,tensor_rank,kappa_c)*reduced_angular_element(kappa_b, &
                & tensor_rank,kappa_d)/rank_norm
            radial_integral_coefficients(1) = reduced_angular_product*dble(4*kappa_a*(kappa_b + kappa_d))
            return
        !
        !   4.0  IBR = 3. Direct magnetic F-integrals
        !
        case (3)
            if (orbital_a/=orbital_c .or. orbital_b/=orbital_d) call stop_invalid_orbital_combination
            if (mod(tensor_rank,2) /= 1) return
            rank_norm = tensor_rank*(tensor_rank + 1)
            reduced_angular_product = reduced_angular_element(kappa_a,tensor_rank,kappa_a)*reduced_angular_element(kappa_b, &
                & tensor_rank,kappa_b)/rank_norm
            radial_integral_coefficients(1) = reduced_angular_product*dble(16*kappa_a*kappa_b)
            return
        !
        !   5.0   IBR = 4. Exchange magnetic G- and H-integrals
        !
        case (4)
            if (orbital_a/=orbital_d .or. orbital_b/=orbital_c) call stop_invalid_orbital_combination
            if (radial_integral_order - tensor_rank >= 0) then
                if (radial_integral_order - tensor_rank <= 0) then
                    radial_integral_coefficients(1) = dble(kappa_a + kappa_b)*reduced_angular_element(kappa_a,tensor_rank,kappa_b)
                    phase_exponent = abs(kappa_a) - abs(kappa_b) + tensor_rank + 1
                    radial_integral_coefficients(1) = radial_integral_coefficients(1)*radial_integral_coefficients(1) &
                        /dble(tensor_rank*(tensor_rank + 1))
                    if (mod(phase_exponent,2) /= 0) radial_integral_coefficients(1) = -radial_integral_coefficients(1)
                    radial_integral_coefficients(3) = radial_integral_coefficients(1)
                    radial_integral_coefficients(2) = radial_integral_coefficients(1) + radial_integral_coefficients(1)
                    return
                end if
                !
                !   5.2  NU = K+1
                !
                exchange_kappa_difference = kappa_b - kappa_a
                next_rank = tensor_rank + 1
                rank = tensor_rank
                g1 = exchange_kappa_difference + next_rank
                g2 = exchange_kappa_difference - next_rank
                two_rank_plus_one = tensor_rank + tensor_rank + 1
                reduced_angular_product = reduced_angular_element(kappa_a,tensor_rank,kappa_b)**2
                if (kappa_a*kappa_b < 0) reduced_angular_product = -reduced_angular_product
                leading_prefactor = reduced_angular_product*rank/next_rank/dble(two_rank_plus_one*(two_rank_plus_one + 2))
                radial_integral_coefficients(1) = -leading_prefactor*g1*g1
                radial_integral_coefficients(2) = -2.0_real64*leading_prefactor*g1*g2
                radial_integral_coefficients(3) = -leading_prefactor*g2*g2
                return
            end if
            !
            !   5.3  NU = K-1
            !
            exchange_kappa_difference = kappa_b - kappa_a
            rank = tensor_rank
            next_rank = tensor_rank + 1
            f1 = exchange_kappa_difference + rank
            f2 = exchange_kappa_difference - rank
            g1 = exchange_kappa_difference + next_rank
            g2 = exchange_kappa_difference - next_rank
            two_rank_plus_one = tensor_rank + tensor_rank + 1
            reduced_angular_product = reduced_angular_element(kappa_a,tensor_rank,kappa_b)**2
            if (kappa_a*kappa_b < 0) reduced_angular_product = -reduced_angular_product
            leading_prefactor = reduced_angular_product*next_rank/rank/dble(two_rank_plus_one*(two_rank_plus_one - 2))
            radial_integral_coefficients(1) = -leading_prefactor*f2*f2
            radial_integral_coefficients(2) = -2.0_real64*leading_prefactor*f1*f2
            radial_integral_coefficients(3) = -leading_prefactor*f1*f1
            mixed_prefactor = reduced_angular_product/dble(two_rank_plus_one*two_rank_plus_one)
            mixed_prefactor = mixed_prefactor + mixed_prefactor
            radial_integral_coefficients(4) = -mixed_prefactor*f1*g2
            !     S(5) = S(4)
            radial_integral_coefficients(5) = -mixed_prefactor*f2*g1
            radial_integral_coefficients(6) = -mixed_prefactor*f1*g1
            radial_integral_coefficients(7) = -mixed_prefactor*f2*g2
            return
        end select
    !
    contains

        ! Stop on an illegal combination of orbital indices.
        subroutine stop_invalid_orbital_combination
            write (*, '("CXK: Type ",i2,"-"/,1x,i2,3x,4i3,2x,2i3,2x,i2)') &
                radial_integral_case, orbital_indices(1), orbital_indices(2), orbital_indices(3), orbital_indices(4), &
                radial_integral_order, tensor_rank, direct_exchange_code
            stop
        end subroutine stop_invalid_orbital_combination
    end subroutine calculate_breit_radial_integral_coefficients

end module breit_radial_integral_coefficients
