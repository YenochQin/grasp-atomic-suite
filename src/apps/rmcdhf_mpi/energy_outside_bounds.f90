
!***********************************************************************
!                                                                      *
logical function energy_outside_bounds (etry)
    ! Module dependencies
    !                                                                      *
    !   This  subprogram determines whether the trial eigenvalue etry is   *
    !   within the bounds (EPSMIN,EPSMAX)                                  *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
        nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
    implicit none (type, external)

    ! Arguments
    real(real64) , intent(in) :: etry
    !
    !
    if (etry>epsmin .and. etry<epsmax) then
        energy_outside_bounds = .false.
    else
        energy_outside_bounds = .true.
    end if
    !
    return
end function energy_outside_bounds
