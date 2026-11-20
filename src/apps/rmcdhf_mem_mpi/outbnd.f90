
!***********************************************************************
!                                                                      *
      logical function outbnd (etry)
!                                                                      *
!   This  subprogram determines whether the trial eigenvalue etry is   *
!   within the bounds (EPSMIN,EPSMAX)                                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
          nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(real64) , intent(in) :: etry
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!
      if (etry>epsmin .and. etry<epsmax) then
         outbnd = .false.
      else
         outbnd = .true.
      end if
!
      return
      end function outbnd
