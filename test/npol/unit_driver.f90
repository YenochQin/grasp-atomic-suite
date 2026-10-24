program npol_unit_driver
    use, intrinsic :: ieee_arithmetic
    use vast_kind_param, only: double
    use parameter_def, only: NNNP
    use npol_C
    use npol_input
    use npol_kernel
    use npol_initialization
    use npol_test_support
    implicit none
    type(npol_config) :: cfg
    real(double) :: r(NNNP+1), a(NNNP+1), b(NNNP+1), v(NNNP+1)
    real(double) :: coul(NNNP+1), zz(NNNP+1), saved(NNNP+1), nan, inf
    integer :: status, n, i, u
    character(32) :: mode
    character(4096) :: path
    character(NPOL_MESSAGE_LEN) :: message
    call get_command_argument(1, mode)
    call get_command_argument(2, path)
    cfg = fixture()
    nan = ieee_value(0.0_double, ieee_quiet_nan)
    inf = ieee_value(0.0_double, ieee_positive_inf)
    r = 0
    r(:3) = [0.0_double, 2.0_double, 4.0_double]
    a = 123
    b = 123
    v = 123
    coul = 7
    zz = coul
    select case (trim(mode))
    case ('input')
        call npol_read_config(trim(path), cfg, status, message)
        if (status == 0) call npol_validate_config(cfg, 2.0_double, status, message)
        if (status /= 0) then
            write(*,'(A)') trim(message)
            stop 1
        end if
    case ('init')
        call npol_init(trim(path), 2.0_double, 3, r, zz, status, message)
        if (status /= 0) then
            call require(all(zz == coul), 'failed init preserves Coulomb and tails')
            call require(.not.npol_ready .and. .not.npol_settings%enabled, 'failed init stays disabled')
            call require(all(vnp_total == 0), 'failed init has no partial NP array')
            write(*,'(A)') trim(message)
            stop 1
        end if
        call require(npol_ready .and. npol_n == 3, 'initialization committed')
        call require(all(zz_coul(:3) == coul(:3)), 'baseline captured')
        call require(all(zz(:3) == coul(:3)-r(:3)*vnp_total(:3)), 'initialization composition')
        call require(all(zz(4:) == coul(4:)), 'init preserves tail')
        if (.not.npol_settings%enabled) call require(all(zz == coul), 'disabled exact copy')
    case ('kernel')
        call npol_build_grid(cfg, 3, r, a, b, v, status, message)
        call require(status == 0, message)
        call near(a(1), -1.0_double, 'E1 origin')
        call near(b(1), -2.0_double, 'E2 origin')
        call near(a(2), -0.5_double, 'E1 r=b')
        call near(b(2), -1.0_double, 'E2 r=b')
        call near(a(3), -16.0_double/272, 'E1 r>b')
        call near(b(3), -128.0_double/4160, 'E2 r>b')
        call require(all(v(:3) == a(:3)+b(:3)), 'channel sum')
        cfg%use_e1 = .false.
        cfg%k1 = nan
        cfg%b1 = nan
        call npol_build_grid(cfg, 3, r, a, b, v, status, message)
        call require(status == 0 .and. all(a(:3) == 0), 'unused E1 short circuit')
        cfg = fixture()
        cfg%use_e2 = .false.
        cfg%k2 = nan
        cfg%b2 = nan
        call npol_build_grid(cfg, 3, r, a, b, v, status, message)
        call require(status == 0 .and. all(b(:3) == 0), 'unused E2 short circuit')
        cfg%k1 = 0
        cfg%b1 = tiny(1.0_double)
        call npol_build_grid(cfg, 3, r, a, b, v, status, message)
        call require(status == 0 .and. all(v(:3) == 0), 'explicit zero, no denominator underflow')
        cfg%enabled = .false.
        cfg%b1 = nan
        call npol_build_grid(cfg, 3, r, a, b, v, status, message)
        call require(status == 0 .and. all(v(:3) == 0), 'disabled ignores channel values')
        cfg = fixture()
        cfg%use_e2 = .false.
        cfg%b1 = 1.0e75_double
        cfg%k1 = 1.0e300_double
        r(:3) = 1.0e80_double
        call npol_build_grid(cfg, 3, r, a, b, v, status, message)
        call require(status == 0, 'large powers do not overflow')
        call near(v(1)/1.0e-20_double, -1.0_double, 'scaled quotient')
        cfg%b1 = 1.0e80_double
        cfg%k1 = 1
        r(1) = 0
        call npol_build_grid(cfg, 1, r, a, b, v, status, message)
        call require(status == 0 .and. v(1) < 0, 'representable subnormal potential')
        cfg%b1 = tiny(1.0_double)
        saved = v
        call npol_build_grid(cfg, 1, r, a, b, v, status, message)
        call require(status /= 0 .and. all(v == saved), 'overflow fails transactionally')
    case ('compose')
        call write_fixture('npol.inp')
        call npol_init('npol.inp', 2.0_double, 3, r, zz, status, message)
        call require(status == 0, message)
        saved = zz
        do i = 1, 2
            call npol_compose(3, r, zz_coul, vnp_total, zz, status, message)
            call require(status == 0 .and. all(zz == saved), 'recomposition is idempotent')
        end do
        call near(zz(2), 10.0_double, 'minus r*V sign')
        open(newunit=u, file='npol.inp', status='replace')
        write(u,'(A)') '&npol schema_version=1, enabled=.false. /'
        close(u)
        zz = coul
        call npol_init('npol.inp', 2.0_double, 3, r, zz, status, message)
        call require(status == 0 .and. all(zz == coul), 'enabled to disabled, fresh Coulomb')
        call require(all(vnp_total == 0), 'disabled clears previous NP')
        call npol_reset()
        call require(npol_n == 0 .and. .not.npol_ready, 'reset state')
        call require(all(zz_coul == 0), 'reset baseline')
        call write_fixture('npol.inp')
        r(:2) = [0.0_double, 1.0_double]
        zz = 11
        call npol_init('npol.inp', 2.0_double, 2, r, zz, status, message)
        call require(status == 0 .and. npol_n == 2, 'new grid length')
        call require(all(zz_coul(:2) == 11) .and. all(zz_coul(3:) == 0), 'new baseline, no stale tail')
        call near(zz(2), 11.0_double+16.0_double/17+128.0_double/65, 'new grid potential')
        open(newunit=u, file='npol.inp', status='replace')
        write(u,'(A)') '&npol schema_version=99 /'
        close(u)
        zz = coul
        call npol_init('npol.inp', 2.0_double, 3, r, zz, status, message)
        call require(status /= 0 .and. all(zz == coul), 'reinit failure preserves fresh Coulomb')
        call require(.not.npol_settings%enabled .and. .not.npol_ready, 'reinit failure clears published state')
    case ('bounds')
        do n = 1, NNNP, NNNP-1
            call npol_build_grid(cfg, n, r, a, b, v, status, message)
            call require(status == 0, message)
            call require(all(a(n+1:) == 123), 'kernel leaves inactive tail')
            call npol_compose(n, r, coul, v, zz, status, message)
            call require(status == 0, message)
            call require(all(zz(n+1:) == 7), 'compose leaves inactive tail')
        end do
        do i = 1, 3
            select case(i)
            case(1)
                n = 0
            case(2)
                n = -1
            case(3)
                n = NNNP+1
            end select
            call npol_build_grid(cfg, n, r, a, b, v, status, message)
            call require(status /= 0, 'invalid N rejected by kernel')
            call npol_compose(n, r, coul, v, zz, status, message)
            call require(status /= 0, 'invalid N rejected by compose')
        end do
        do i = 1, 3
            select case(i)
            case(1)
                r(2) = -1
            case(2)
                r(2) = nan
            case(3)
                r(2) = inf
            end select
            call npol_build_grid(cfg, 3, r, a, b, v, status, message)
            call require(status /= 0, 'invalid radius rejected')
        end do
        r(2) = 2
        call npol_build_grid(cfg, 3, r, a(:2), b, v, status, message)
        call require(status /= 0, 'short output rejected')
        call npol_compose(3, r(:2), coul, v, zz, status, message)
        call require(status /= 0, 'short input rejected')
        v(2) = -huge(1.0_double)
        saved = zz
        call npol_compose(3, r, coul, v, zz, status, message)
        call require(status /= 0 .and. all(zz == saved), 'product overflow is transactional')
        r(2) = 1
        coul(2) = huge(1.0_double)
        call npol_compose(3, r, coul, v, zz, status, message)
        call require(status /= 0 .and. all(zz == saved), 'sum overflow is transactional')
        v(2) = nan
        call npol_compose(3, r, coul, v, zz, status, message)
        call require(status /= 0, 'nonfinite potential rejected')
    case default
        error stop 'unknown test mode'
    end select
    print '(A)', 'PASS '//trim(mode)
end program npol_unit_driver
