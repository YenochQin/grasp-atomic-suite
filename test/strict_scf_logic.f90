program strict_scf_logic
    use vast_kind_param, only: double
    use strict_scf, only: update_strict_scf
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    real(double), parameter :: accy = 1.5625d-8
    real(double) :: residual(3), energy, previous, scf_max, relative, nan
    logical :: fixed(3), orbital_ok, energy_ok, converged
    integer :: streak, nit

    nan = ieee_value(0.d0, ieee_quiet_nan)
    fixed = [.true., .false., .false.]
    energy = -15615.d0
    previous = energy
    streak = 0

    ! Regression: the original code stopped on this energy plateau even
    ! though weakly occupied correlation orbitals were still changing.
    residual = [nan, 2.d-5, 3.d-5]
    do nit = 1, 20
        call step(.true., nit)
        if (orbital_ok .or. converged .or. streak /= 0) stop 1
    end do
    if (.not.energy_ok) stop 2

    ! No first-cycle energy comparison and no inherited restart streak.
    residual = [nan, accy, 0.d0]
    streak = 2
    call step(.true., 1)
    if (.not.orbital_ok .or. energy_ok .or. converged .or. streak /= 0) stop 3
    call step(.true., 2)
    if (.not.energy_ok .or. converged .or. streak /= 1) stop 4

    ! Orbital-only convergence does not suffice either; a failure must reset.
    previous = energy + 1.d-4
    call step(.true., 3)
    if (.not.orbital_ok .or. energy_ok .or. converged .or. streak /= 0) stop 5
    previous = energy
    call step(.true., 4)
    if (converged .or. streak /= 1) stop 6
    call step(.true., 5)
    if (.not.converged .or. streak /= 2) stop 7

    ! Nonfinite varied orbitals / energies cannot produce false success.
    residual(2) = nan
    call step(.true., 6)
    if (orbital_ok .or. converged .or. streak /= 0) stop 8
    residual(2) = 0.d0
    previous = nan
    call step(.true., 7)
    if (energy_ok .or. converged) stop 9
    previous = energy
    energy = nan
    call step(.true., 8)
    if (energy_ok .or. converged) stop 10
    energy = 0.d0
    call step(.true., 9)
    if (energy_ok .or. converged) stop 11

    ! EAL must not read/require an undefined weighted energy.
    energy = nan
    previous = nan
    streak = 0
    call step(.false., 1)
    if (converged .or. streak /= 1) stop 12
    call step(.false., 2)
    if (.not.converged .or. streak /= 2) stop 13

    ! A frozen-only CI run has no orbital updates but still uses the EOL gate.
    fixed = .true.
    residual = nan
    energy = -15615.d0
    previous = energy
    do nit = 1, 3
        call step(.true., nit)
        if (.not.orbital_ok .or. scf_max /= 0.d0) stop 14
        if (converged .neqv. (nit == 3)) stop 15
    end do
    print *, 'Strict SCF stopping regression passed.'

contains
    subroutine step(eol, iteration)
        logical, intent(in) :: eol
        integer, intent(in) :: iteration
        call update_strict_scf(eol, iteration, accy, residual, fixed, &
            energy, previous, streak, scf_max, relative, &
            orbital_ok, energy_ok, converged)
    end subroutine step
end program strict_scf_logic
