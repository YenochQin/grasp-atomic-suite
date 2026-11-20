!***********************************************************************
!                                                                      *
subroutine update_trial_orbital_energy(orbital_index, first_oscillation_sign, effective_principal_n, &
    & node_count_difference, energy_correction, &
    & solve_failed, orbital_sign_reversal_count)
    ! Module dependencies
    !                                                                      *
    !   This  subroutine implements Part 2 of Algorithm 7.1 in  C Froese   *
    !   Fischer,  Comput  Phys  Rep, 3 (1986) 273-326. (The present code   *
    !   actually  implements  the version used in her program MCHF where   *
    !   differences occur.)                                                *
    !                                                                      *
    !   Call(s) to: [RSCF92]: OUTBND.                                      *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford    Last Update: 08 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use radial_dirac_workspace, only: max_radial_points, mtp0, p, p0, q, q0, tf, tg, xu, xv
    use orbital_csf_basis, only: max_orbitals, orbital_energy_parameters, orbital_principal_n
    use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
        nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya

    ! Procedure interfaces
    use energy_outside_bounds_interface, only: energy_outside_bounds
    implicit none (type, external)

    ! Arguments
    integer , intent(in) :: orbital_index
    integer , intent(in) :: effective_principal_n
    integer , intent(in) :: node_count_difference
    integer , intent(out) :: orbital_sign_reversal_count
    real(real64) , intent(in) :: first_oscillation_sign
    real(real64) , intent(inout) :: energy_correction
    logical , intent(out) :: solve_failed

    ! Local constants
    real(real64), parameter :: p02 = 2.0e-2_real64
    real(real64), parameter :: p05 = 5.0e-2_real64
    real(real64), parameter :: p001 = 1.0e-3_real64
    real(real64), parameter :: p00001 = 1.0e-5_real64
    real(real64), parameter :: d2p5 = 2.5_real64

    ! Local declarations
    integer :: i
    real(real64) :: eps, abdele, delebe, dn, dnprme, dnpn25, etry, delta
    !-----------------------------------------------
    !
    !
    !   Determine if the iterative process has succeeded
    !
    if (first_oscillation_sign>0.0_real64 .and. node_count_difference==0) then
        solve_failed = .false.
    else
        solve_failed = .true.
    end if
    !
    !   Inversion counter
    !
    orbital_sign_reversal_count = 0
    !
    !   If unsuccessful, obtain a new estimate of the eigenvalue
    !
    if (solve_failed) then
        !
        !   Define quantities used later
        !
        eps = orbital_energy_parameters(orbital_index)
        abdele = abs(energy_correction)
        delebe = abs(energy_correction/eps)
        dn = dble(orbital_principal_n(orbital_index))
        dnprme = dble(effective_principal_n)
        dnpn25 = (dnprme/dn)**d2p5
        !
        if (abs(node_count_difference)==1 .and. delebe>p02 .or. node_count_difference==0 .and. delebe>=p00001 .and. &
            abdele>=p001) then
            do
                etry = eps + energy_correction
                if (energy_outside_bounds(etry) .and. node_count_difference/=0) etry = eps*dnpn25
                if (energy_outside_bounds(etry)) etry = eps - energy_correction
                if (energy_outside_bounds(etry)) then
                    energy_correction = 0.5_real64*energy_correction
                else
                    exit
                end if
            end do
        else if (node_count_difference == 0) then
            etry = eps
            orbital_sign_reversal_count = 1
            p0 = -p0
            p(:mtp0) = -p(:mtp0)
            q(:mtp0) = -q(:mtp0)
            solve_failed = .false.
        else if (first_oscillation_sign < 0.0_real64) then
            etry = eps*dnpn25
            if (energy_outside_bounds(etry)) etry = 0.5_real64*(epsmax + eps)
            if (energy_outside_bounds(etry)) etry = 0.5_real64*(epsmin + eps)
        else if (node_count_difference < 0) then
            delta = 1.0_real64 - eps/epsmax
            epsmax = eps
            if (delta < p05) then
                emax = emax*dnpn25
            else
                emax = eps*dnpn25
            end if
            etry = emax
            if (energy_outside_bounds(etry)) etry = 0.5_real64*(epsmax + epsmin)
        else if (node_count_difference > 0) then
            delta = 1.0_real64 - epsmin/eps
            epsmin = eps
            if (delta < p05) then
                emin = emin*dnpn25
            else
                emin = eps*dnpn25
            end if
            etry = emin
            if (energy_outside_bounds(etry)) etry = 0.5_real64*(epsmax + epsmin)
        end if
        orbital_energy_parameters(orbital_index) = etry
    end if
    !
    return
end subroutine update_trial_orbital_energy
