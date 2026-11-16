module npol_kernel
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use parameter_def, only: max_radial_points
    use npol_c, only: npol_config
    use npol_input, only: npol_validate_config
    implicit none (type, external)
    private
    public :: npol_check_grid, npol_build_grid, npol_compose
contains
    subroutine npol_check_grid(n, r, status, message)
        use, intrinsic :: iso_fortran_env, only: real64
        integer, intent(in) :: n
        real(real64), intent(in) :: r(:)
        integer, intent(out) :: status
        character(*), intent(out) :: message
        status = 1
        message = 'npol: N must be in 1:NNNP and fit the supplied arrays'
        if (n < 1 .or. n > max_radial_points .or. n > size(r)) return
        message = 'npol: radial grid must be finite and nonnegative'
        if (.not.all(ieee_is_finite(r(:n)))) return
        if (any(r(:n) < 0)) return
        status = 0
        message = ''
    end subroutine

    subroutine npol_build_grid(cfg, n, r, e1, e2, total, status, message)
        use, intrinsic :: iso_fortran_env, only: real64
        type(npol_config), intent(in) :: cfg
        integer, intent(in) :: n
        real(real64), intent(in) :: r(:)
        real(real64), intent(inout) :: e1(:), e2(:), total(:)
        integer, intent(out) :: status
        character(*), intent(out) :: message
        real(real64) :: a(max_radial_points), b(max_radial_points), v(max_radial_points)
        integer :: i
        call npol_check_grid(n, r, status, message)
        if (status /= 0) return
        status = 1
        message = 'npol: output arrays shorter than N'
        if (min(size(e1), size(e2), size(total)) < n) return
        call npol_validate_config(cfg, real(cfg%isotope_z,real64), status, message)
        if (status /= 0) return
        a(:n) = 0
        b(:n) = 0
        do i = 1, n
            if (cfg%enabled .and. cfg%use_e1) then
                call channel(r(i), cfg%k1, cfg%b1, 4, a(i), status)
                if (status /= 0) exit
            end if
            if (cfg%enabled .and. cfg%use_e2) then
                call channel(r(i), cfg%k2, cfg%b2, 6, b(i), status)
                if (status /= 0) exit
            end if
            if (abs(a(i)) > huge(1.0_real64)-abs(b(i))) then
                status = 1
                exit
            end if
            v(i) = a(i) + b(i)
        end do
        if (status /= 0) then
            message = 'npol: potential is not representable as a finite value'
            return
        end if
        e1(:n) = a(:n)
        e2(:n) = b(:n)
        total(:n) = v(:n)
    end subroutine

    ! Evaluate K/(r**p+b**p) without forming potentially overflowing powers.
    ! Binary mantissas/exponents preserve subnormal results when representable.
    subroutine channel(r, k, b, p, v, status)
        use, intrinsic :: iso_fortran_env, only: real64
        real(real64), intent(in) :: r, k, b
        integer, intent(in) :: p
        real(real64), intent(out) :: v
        integer, intent(out) :: status
        real(real64) :: s, m, denominator
        integer :: e
        status = 0
        v = 0
        if (k == 0) return
        s = max(r,b)
        denominator = (r/s)**p + (b/s)**p
        m = fraction(k)/(fraction(s)**p*denominator)
        e = exponent(k) - p*exponent(s) + exponent(m)
        if (e > maxexponent(v)) then
            status = 1
            return
        end if
        if (e < minexponent(v)-digits(v)) return
        v = -scale(fraction(m), e)
        if (.not.ieee_is_finite(v)) status = 1
    end subroutine

    subroutine npol_compose(n, r, coul, v, total, status, message)
        use, intrinsic :: iso_fortran_env, only: real64
        integer, intent(in) :: n
        real(real64), intent(in) :: r(:), coul(:), v(:)
        real(real64), intent(inout) :: total(:)
        integer, intent(out) :: status
        character(*), intent(out) :: message
        real(real64) :: candidate(max_radial_points), shift
        integer :: i
        call npol_check_grid(n, r, status, message)
        if (status /= 0) return
        status = 1
        message = 'npol: composition arrays shorter than N'
        if (min(size(coul), size(v), size(total)) < n) return
        message = 'npol: nonfinite Coulomb or NP input'
        if (.not.all(ieee_is_finite(coul(:n))) .or. .not.all(ieee_is_finite(v(:n)))) return
        candidate(:n) = coul(:n)
        message = 'npol: composed potential is not finite'
        do i = 1, n
            if (v(i) == 0 .or. r(i) == 0) cycle
            if (r(i) > 1) then
                if (abs(v(i)) > huge(shift)/r(i)) return
            end if
            shift = -r(i)*v(i)
            if (shift > 0) then
                if (coul(i) > huge(shift)-shift) return
            else
                if (coul(i) < -huge(shift)-shift) return
            end if
            candidate(i) = coul(i) + shift
        end do
        total(:n) = candidate(:n)
        status = 0
        message = ''
    end subroutine
end module npol_kernel
