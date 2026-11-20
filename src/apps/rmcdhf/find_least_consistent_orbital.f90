!***********************************************************************
!                                                                      *
subroutine find_least_consistent_orbital(orbital_index)
    ! Module dependencies
    !                                                                      *
    !   This subroutine finds the least self-consistent orbital            *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford    Last update: 10 Dec 1992   *
    !                                                                      *
    ! J initialized to zero
    ! XHH 1997.02.14
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use fixed_orbital_flags, only: fixed_orbital_count, is_orbital_fixed, max_orbitals
    use orbital_csf_basis, only: max_orbitals, orbital_count
    use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
        nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
    implicit none (type, external)

    ! Arguments
    integer , intent(out) :: orbital_index

    ! Local declarations
    integer :: i
    real(real64) :: dlrgst
    !-----------------------------------------------
    !
    !
    orbital_index = 0
    dlrgst = 0.0_real64
    do i = 1, orbital_count
        if (is_orbital_fixed(i)) cycle
        if (scnsty(i) <= dlrgst) cycle
        dlrgst = scnsty(i)
        orbital_index = i
    end do
    !
    return
end subroutine find_least_consistent_orbital
