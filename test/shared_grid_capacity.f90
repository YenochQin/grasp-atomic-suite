! Verify that separately compiled shared modules agree on the radial capacity.
! This catches legacy .mod files shadowing a newly compiled parameter_def.
program shared_grid_capacity
    use parameter_def, only: NNNP, NNN1
    use grid_C, only: R, RP, RPOR
    use pote_C, only: XP, XQ, YP
    use tatb_C, only: TA, TB
    implicit none

    if (NNN1 /= NNNP + 10) error stop 'inconsistent radial capacity parameters'
    if (size(R) /= NNN1 .or. size(RP) /= NNN1 .or. size(RPOR) /= NNN1) &
        error stop 'grid_C was compiled with a different radial capacity'
    if (size(XP) /= NNNP .or. size(XQ) /= NNNP .or. size(YP) /= NNNP) &
        error stop 'pote_C was compiled with a different radial capacity'
    if (size(TA) /= NNN1 .or. size(TB) /= NNN1) &
        error stop 'tatb_C was compiled with a different radial capacity'
    print *, 'Shared radial capacity:', NNNP, NNN1
end program shared_grid_capacity
