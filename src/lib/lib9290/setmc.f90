!***********************************************************************
!                                                                      *
      subroutine setmc
!                                                                      *
!   This subprogram sets machine-dependent parameters.                 *
!                                                                      *
!   Call(s) to: [LAPACK]: DLAMCH.                                      *
!                                                                      *
!   Written by Farid A Parpia              Last updated: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:37   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use debug_c
      use def_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use suite_linear_algebra, only: dlamch
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(real64) :: dnum
!-----------------------------------------------
!
!
!   Set the machine-dependent parameters:
!
!   TENMAX - Maximum size of exponent of 10
!
      tenmax = dlamch('L')
!
!   EXPMAX - Maximum size of exponent of e
!
      dnum = dlamch('O')
      expmax = log(dnum)
!
!   EXPMIN - Minimum size of exponent of e
!
      dnum = dlamch('U')
      expmin = log(dnum)
!
!   PRECIS - Machine precision
!
      precis = dlamch('E')
!
!   Debug printout
!
      if (ldbpg(1)) write (99, 300) tenmax, expmax, expmin, precis
!
      return
!
  300 format(/,'From SUBROUTINE SETMC:'/,' TENMAX (maximum exponent of 10): ',&
         f5.0,/,' EXPMAX (maximum exponent of e): ',1p,1d19.12,/,&
         ' EXPMIN (minimum exponent of e): ',1d19.12,/,&
         ' PRECIS (machine precision): ',1d19.12)
      return
!
      end subroutine setmc
