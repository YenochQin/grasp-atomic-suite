! Shared initialization for orbital optimization and property applications.
! Runtime state stays in the legacy grid_C/def_C modules for compatibility.
module radial_grid_defaults
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use vast_kind_param, only: DOUBLE
    use suite_parameters, only: NNNP, FINITE_RNT_SCALE, FINITE_H, FINITE_N, &
        POINT_RNT_SCALE, POINT_H, POINT_N, DEFAULT_HP, DEFAULT_ACCY
    use grid_C, only: RNT, H, HP, N
    use def_C, only: ACCY
    implicit none
    private
    public :: SET_RADIAL_DEFAULTS, UPDATE_RADIAL_ACCURACY, VALIDATE_RADIAL_GRID

contains

    subroutine SET_RADIAL_DEFAULTS(nuclear_parameters, charge)
        integer, intent(in) :: nuclear_parameters
        real(DOUBLE), intent(in) :: charge

        if (.not.ieee_is_finite(charge) .or. charge <= 0.0D0) &
            error stop 'radial defaults: nuclear charge must be positive and finite'
        if (nuclear_parameters == 0) then
            RNT = POINT_RNT_SCALE/charge
            H = POINT_H
            N = POINT_N
        else
            RNT = FINITE_RNT_SCALE/charge
            H = FINITE_H
            N = FINITE_N
        endif
        HP = DEFAULT_HP
        call UPDATE_RADIAL_ACCURACY
        call VALIDATE_RADIAL_GRID
    end subroutine SET_RADIAL_DEFAULTS

    subroutine UPDATE_RADIAL_ACCURACY
        if (.not.ieee_is_finite(DEFAULT_ACCY) .or. DEFAULT_ACCY < 0.0D0) &
            error stop 'radial defaults: DEFAULT_ACCY must be zero or positive and finite'
        if (DEFAULT_ACCY > 0.0D0) then
            ACCY = DEFAULT_ACCY
        else
            ACCY = H**6
        endif
    end subroutine UPDATE_RADIAL_ACCURACY

    subroutine VALIDATE_RADIAL_GRID
        if (N < 13 .or. N > NNNP) &
            error stop 'radial grid: N must be between 13 and NNNP'
        if (.not.ieee_is_finite(RNT) .or. RNT <= 0.0D0) &
            error stop 'radial grid: RNT must be positive and finite'
        if (.not.ieee_is_finite(H) .or. H <= 0.0D0) &
            error stop 'radial grid: H must be positive and finite'
        if (.not.ieee_is_finite(HP) .or. HP < 0.0D0) &
            error stop 'radial grid: HP must be nonnegative and finite'
        if (.not.ieee_is_finite(ACCY) .or. ACCY <= 0.0D0) &
            error stop 'radial grid: ACCY must be positive and finite'
    end subroutine VALIDATE_RADIAL_GRID
end module radial_grid_defaults
