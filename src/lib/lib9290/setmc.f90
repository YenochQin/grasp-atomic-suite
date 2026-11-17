module setmc_i
    implicit none (type, external)
    private
    public :: setmc
contains
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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use debug_c, only: cutoff, double, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use def_c, only: double, expmax, expmin, precis, tenmax
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
      tenmax = dlamch("L")
!
!   EXPMAX - Maximum size of exponent of e
!
      dnum = dlamch("O")
      expmax = log(dnum)
!
!   EXPMIN - Minimum size of exponent of e
!
      dnum = dlamch("U")
      expmin = log(dnum)
!
!   PRECIS - Machine precision
!
      precis = dlamch("E")
!
!   Debug printout
!
      if (ldbpg(1)) write (99, 300) tenmax, expmax, expmin, precis
!
      return
!
  300 format(/,"From SUBROUTINE SETMC:"/," TENMAX (maximum exponent of 10): ",&
         f5.0,/," EXPMAX (maximum exponent of e): ",1p,1d19.12,/,&
         " EXPMIN (minimum exponent of e): ",1d19.12,/,&
         " PRECIS (machine precision): ",1d19.12)
      return
!
      end subroutine setmc
end module setmc_i
