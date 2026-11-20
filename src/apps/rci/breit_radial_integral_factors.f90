!***********************************************************************
!                                                                      *
subroutine cxk(radial_integral_coefficients, orbital_indices, twice_orbital_kappa, radial_integral_order, &
    tensor_rank, radial_integral_case, direct_exchange_code)
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
    use reduced_angular_matrix_element, only: cre
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
    integer :: mu, ia, ib, ic, id, ka, kb, kc, kd, kk, ik, ip
    real(real64) :: d, h, dk1, dk2, fk, gk, g1, g2, g3, g4, a, f1, f2, f3, f4&
        , b, dk
    !-----------------------------------------------
    !
    !
    !   1.0  Initialization
    !
    radial_integral_coefficients = 0.0_real64
    !
    ia = orbital_indices(1)
    ib = orbital_indices(2)
    ic = orbital_indices(3)
    id = orbital_indices(4)
    ka = twice_orbital_kappa(1)/2
    kb = twice_orbital_kappa(2)/2
    kc = twice_orbital_kappa(3)/2
    kd = twice_orbital_kappa(4)/2
    if (direct_exchange_code == 2) then
        kk = kd
        ik = id
        kd = kc
        id = ic
        kc = kk
        ic = ik
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
                radial_integral_coefficients(1) = -(ka + kc)*(kd + kb)
                if (tensor_rank == 0) then
                    write (*, '("CXK: Illegal value K = 0 -"/,1x,4i3,2x,i3,2x,2i2)') &
                        orbital_indices(1), orbital_indices(2), orbital_indices(3), orbital_indices(4), &
                        radial_integral_order, radial_integral_case, direct_exchange_code
                    stop
                end if
                d = tensor_rank*(tensor_rank + 1)
                h = cre(ka,tensor_rank,kc)*cre(kb,tensor_rank,kd)
                if (mod(tensor_rank,2) /= 0) h = -h
                radial_integral_coefficients(1) = radial_integral_coefficients(1)*h/d
                radial_integral_coefficients(2:4) = radial_integral_coefficients(1)
                return
            end if
            !
            !   2.2  NU = K+1
            !
            dk1 = kc - ka
            dk2 = kd - kb
            fk = tensor_rank
            gk = tensor_rank + 1
            g1 = dk1 - gk
            g2 = dk1 + gk
            g3 = dk2 - gk
            g4 = dk2 + gk
            kk = tensor_rank + tensor_rank + 1
            h = cre(ka,tensor_rank,kc)*cre(kb,tensor_rank,kd)
            if (mod(tensor_rank,2) /= 0) h = -h
            a = h*fk/gk/dble(kk*(kk + 2))
            radial_integral_coefficients(1) = a*g1*g3
            radial_integral_coefficients(2) = a*g2*g4
            radial_integral_coefficients(3) = a*g1*g4
            radial_integral_coefficients(4) = a*g2*g3
            return
        end if
        !
        !   2.2  NU = K-1
        !
        dk1 = kc - ka
        dk2 = kd - kb
        fk = tensor_rank
        gk = tensor_rank + 1
        f1 = dk1 - fk
        f2 = dk1 + fk
        f3 = dk2 - fk
        f4 = dk2 + fk
        g1 = dk1 - gk
        g2 = dk1 + gk
        g3 = dk2 - gk
        g4 = dk2 + gk
        kk = tensor_rank + tensor_rank + 1
        h = cre(ka,tensor_rank,kc)*cre(kb,tensor_rank,kd)
        if (mod(tensor_rank,2) /= 0) h = -h
        a = h*gk/fk/dble(kk*(kk - 2))
        radial_integral_coefficients(1) = a*f2*f4
        radial_integral_coefficients(2) = a*f1*f3
        radial_integral_coefficients(3) = a*f2*f3
        radial_integral_coefficients(4) = a*f1*f4
        b = h/dble(kk*kk)
        radial_integral_coefficients(5) = b*f2*g3
        radial_integral_coefficients(6) = b*f4*g1
        radial_integral_coefficients(7) = b*f1*g4
        radial_integral_coefficients(8) = b*f3*g2
        radial_integral_coefficients(9) = b*f2*g4
        radial_integral_coefficients(10) = b*f3*g1
        radial_integral_coefficients(11) = b*f1*g3
        radial_integral_coefficients(12) = b*f4*g2
        return
    !
    !   3.0  IBR = 2  Degenerate case: only one non-zero R-integral
    !
    case (2)
        if (ia/=ic .or. ib==id) then
            if (ia==ic .or. ib/=id) call stop_invalid_orbital_combination
            !
            ik = ib
            ib = ia
            ia = ik
            ik = id
            id = ic
            ic = ik
            !
            kk = kb
            kb = ka
            ka = kk
            kk = kd
            kd = kc
            kc = kk
        end if
        !
        if (mod(tensor_rank,2) /= 1) return
        dk = tensor_rank*(tensor_rank + 1)
        h = cre(ka,tensor_rank,kc)*cre(kb,tensor_rank,kd)/dk
        radial_integral_coefficients(1) = h*dble(4*ka*(kb + kd))
        return
    !
    !   4.0  IBR = 3. Direct magnetic F-integrals
    !
    case (3)
        if (ia/=ic .or. ib/=id) call stop_invalid_orbital_combination
        if (mod(tensor_rank,2) /= 1) return
        dk = tensor_rank*(tensor_rank + 1)
        h = cre(ka,tensor_rank,ka)*cre(kb,tensor_rank,kb)/dk
        radial_integral_coefficients(1) = h*dble(16*ka*kb)
        return
    !
    !   5.0   IBR = 4. Exchange magnetic G- and H-integrals
    !
    case (4)
        if (ia/=id .or. ib/=ic) call stop_invalid_orbital_combination
        if (radial_integral_order - tensor_rank >= 0) then
            if (radial_integral_order - tensor_rank <= 0) then
                radial_integral_coefficients(1) = dble(ka + kb)*cre(ka,tensor_rank,kb)
                ip = abs(ka) - abs(kb) + tensor_rank + 1
                radial_integral_coefficients(1) = radial_integral_coefficients(1)*radial_integral_coefficients(1) &
                    /dble(tensor_rank*(tensor_rank + 1))
                if (mod(ip,2) /= 0) radial_integral_coefficients(1) = -radial_integral_coefficients(1)
                radial_integral_coefficients(3) = radial_integral_coefficients(1)
                radial_integral_coefficients(2) = radial_integral_coefficients(1) + radial_integral_coefficients(1)
                return
            end if
            !
            !   5.2  NU = K+1
            !
            dk = kb - ka
            gk = tensor_rank + 1
            fk = tensor_rank
            g1 = dk + gk
            g2 = dk - gk
            kk = tensor_rank + tensor_rank + 1
            h = cre(ka,tensor_rank,kb)**2
            if (ka*kb < 0) h = -h
            a = h*fk/gk/dble(kk*(kk + 2))
            radial_integral_coefficients(1) = -a*g1*g1
            radial_integral_coefficients(2) = -2.0_real64*a*g1*g2
            radial_integral_coefficients(3) = -a*g2*g2
            return
        end if
        !
        !   5.3  NU = K-1
        !
        dk = kb - ka
        fk = tensor_rank
        gk = tensor_rank + 1
        f1 = dk + fk
        f2 = dk - fk
        g1 = dk + gk
        g2 = dk - gk
        kk = tensor_rank + tensor_rank + 1
        h = cre(ka,tensor_rank,kb)**2
        if (ka*kb < 0) h = -h
        a = h*gk/fk/dble(kk*(kk - 2))
        radial_integral_coefficients(1) = -a*f2*f2
        radial_integral_coefficients(2) = -2.0_real64*a*f1*f2
        radial_integral_coefficients(3) = -a*f1*f1
        b = h/dble(kk*kk)
        b = b + b
        radial_integral_coefficients(4) = -b*f1*g2
        !     S(5) = S(4)
        radial_integral_coefficients(5) = -b*f2*g1
        radial_integral_coefficients(6) = -b*f1*g1
        radial_integral_coefficients(7) = -b*f2*g2
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
end subroutine cxk
