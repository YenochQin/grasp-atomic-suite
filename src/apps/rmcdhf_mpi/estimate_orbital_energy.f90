!***********************************************************************
!                                                                      *
subroutine estimate_orbital_energy(orbital_index)
    ! Module dependencies
    !                                                                      *
    !   This  subprogram implements  Part 1 of Algorithm 7.1 of C Froese   *
    !   Fischer, Comput Phys Rep, 3 (1986) 320-321.                        *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford    Last update: 26 Sep 1993   *
    !                                                                      *
    !***********************************************************************
    use atomic_calculation_state, only: electron_count, nuclear_charge, speed_of_light_au
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use orbital_csf_basis, only: max_orbitals, orbital_kappa, orbital_principal_n
    use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
        nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Arguments
    integer , intent(in) :: orbital_index

    ! Local declarations
    integer :: npj, nakabs
    real(real64) :: alpha, csq, fnrel, fkabs, fkap2, zalpha, gamma, ebym
    !-----------------------------------------------
    !
    !   Initializations
    !
    alpha = 1.0_real64/speed_of_light_au
    csq = speed_of_light_au*speed_of_light_au
    npj = orbital_principal_n(orbital_index)
    nakabs = abs(orbital_kappa(orbital_index))
    fnrel = dble(npj - nakabs)
    fkabs = dble(nakabs)
    fkap2 = fkabs*fkabs
    !
    !   Set ZINF, the asymptotic charge seen by the electron
    !
    !     ZINF = DBLE (IONCTY+1)
    !
    !   Changed on 07/06/93 by WPW
    !
    zinf = nuclear_charge + dble((-electron_count) + 1)
    !
    !   Set the lower bound
    !
    zalpha = zinf*alpha
    if (zalpha < fkabs) then
        gamma = sqrt(fkap2 - zalpha*zalpha)
        ebym = 1.0_real64/sqrt(1.0_real64 + (zalpha/(gamma + fnrel + 0.5_real64))**2)
        epsmin = (1.0_real64 - ebym)*csq
    else
        epsmin = 0.25_real64*csq/dble(npj*npj)
    end if
    emin = epsmin
    !
    !   Set the upper bound
    !
    zalpha = nuclear_charge*alpha
    if (zalpha < fkabs) then
        gamma = sqrt(fkap2 - zalpha*zalpha)
        ebym = 1.0_real64/sqrt(1.0_real64 + (zalpha/(gamma + fnrel - 0.5_real64))**2)
        epsmax = (1.0_real64 - ebym)*csq
    else
        epsmax = csq + csq
    end if
    emax = epsmax
    !
    return
end subroutine estimate_orbital_energy
