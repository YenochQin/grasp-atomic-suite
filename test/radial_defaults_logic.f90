! Exercise the same initialization/override sequence used by all applications.
program radial_defaults_logic
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use vast_kind_param, only: DOUBLE
    use suite_parameters, only: NNNP, FINITE_N, POINT_N, DEFAULT_ACCY
    use radial_grid_defaults, only: SET_RADIAL_DEFAULTS, &
        UPDATE_RADIAL_ACCURACY, VALIDATE_RADIAL_GRID
    use grid_C, only: RNT, H, HP, N
    use def_C, only: ACCY
    implicit none
    character(len=32) :: mode
    real(DOUBLE) :: radius_at_one, expected_accuracy

    call get_command_argument(1, mode)
    if (trim(mode) == 'invalid_charge') then
        call SET_RADIAL_DEFAULTS(2, 0.0D0)
        stop 0
    endif

    call SET_RADIAL_DEFAULTS(2, 1.0D0)
    radius_at_one = RNT
    call SET_RADIAL_DEFAULTS(2, 92.0D0)
    if (N /= FINITE_N) error stop 'finite-nucleus default point count'
    if (abs(RNT*92.0D0-radius_at_one) > 1.0D-18) &
        error stop 'finite-nucleus radius must scale as 1/Z'

    select case (trim(mode))
    case ('invalid_n')
        N = NNNP + 1
    case ('invalid_h')
        H = ieee_value(H, ieee_quiet_nan)
    case ('invalid_accy')
        ACCY = -1.0D0
    case ('')
        ! Simulate changing H interactively, then accepting the default ACCY.
        H = 0.025D0
        call UPDATE_RADIAL_ACCURACY
        expected_accuracy = H**6
        if (DEFAULT_ACCY > 0.0D0) expected_accuracy = DEFAULT_ACCY
        if (ACCY /= expected_accuracy) error stop 'accuracy did not follow final H'
        ! An explicit user tolerance must survive validation unchanged.
        ACCY = 1.0D-12
        call VALIDATE_RADIAL_GRID
        if (ACCY /= 1.0D-12) error stop 'explicit accuracy was overwritten'
        call SET_RADIAL_DEFAULTS(0, 1.0D0)
        radius_at_one = RNT
        call SET_RADIAL_DEFAULTS(0, 92.0D0)
        if (N /= POINT_N) error stop 'point-nucleus default point count'
        if (abs(RNT*92.0D0-radius_at_one) > 1.0D-16) &
            error stop 'point-nucleus radius must scale as 1/Z'
    case default
        error stop 'unknown test mode'
    end select
    call VALIDATE_RADIAL_GRID
end program radial_defaults_logic
