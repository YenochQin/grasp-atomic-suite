! Verify that separately compiled shared modules agree on the radial capacity.
! This catches legacy .mod files shadowing a newly compiled parameter_def.
program shared_grid_capacity
    use parameter_def, only: NNNP, NNN1, KEYORB, NNNW, NNNWM1, NNNWM2
    use suite_parameters, only: CONFIG_NNNP => max_radial_points, CONFIG_NNN1 => radial_workspace_capacity, &
       & CONFIG_KEYORB => orbital_index_radix, CONFIG_NNNW => max_orbitals, CONFIG_NNNWM1 => max_orbitals_minus_one, &
       & CONFIG_NNNWM2 => max_orbitals_minus_two
    use grid_c, only: R => radial_grid, RP => radial_grid_derivative, RPOR => radial_derivative_over_radius
    use pote_c, only: XP, XQ, YP
    use tatb_c, only: TA, TB
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
