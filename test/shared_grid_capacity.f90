! Verify that separately compiled shared modules agree on the radial capacity.
! This catches legacy .mod files shadowing a newly compiled parameter_def.
program shared_grid_capacity
    use parameter_def, only: NNNP, NNN1, KEYORB, NNNW, NNNWM1, NNNWM2
    use suite_parameters, only: CONFIG_NNNP => NNNP, CONFIG_NNN1 => NNN1, &
        CONFIG_KEYORB => KEYORB, CONFIG_NNNW => NNNW, &
        CONFIG_NNNWM1 => NNNWM1, CONFIG_NNNWM2 => NNNWM2
    use grid_C, only: R, RP, RPOR
    use pote_C, only: XP, XQ, YP
    use tatb_C, only: TA, TB
    implicit none

    if (any([NNNP, NNN1, KEYORB, NNNW, NNNWM1, NNNWM2] /= &
            [CONFIG_NNNP, CONFIG_NNN1, CONFIG_KEYORB, CONFIG_NNNW, CONFIG_NNNWM1, CONFIG_NNNWM2])) &
        error stop 'parameter_def does not re-export the configured limits'
    if (NNN1 /= NNNP + 10) error stop 'inconsistent radial capacity parameters'
    if (NNNWM1 /= NNNW - 1 .or. NNNWM2 /= NNNW - 2) error stop 'inconsistent orbital capacities'
    if (size(R) /= NNN1 .or. size(RP) /= NNN1 .or. size(RPOR) /= NNN1) &
        error stop 'grid_C was compiled with a different radial capacity'
    if (size(XP) /= NNNP .or. size(XQ) /= NNNP .or. size(YP) /= NNNP) &
        error stop 'pote_C was compiled with a different radial capacity'
    if (size(TA) /= NNN1 .or. size(TB) /= NNN1) &
        error stop 'tatb_C was compiled with a different radial capacity'
    print *, 'Shared radial capacity:', NNNP, NNN1
end program shared_grid_capacity
