
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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use scf_c
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
      endif
!
      return
      end function outbnd
