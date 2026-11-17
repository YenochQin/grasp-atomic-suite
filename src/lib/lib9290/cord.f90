module cord_i
    implicit none (type, external)
    private
    public :: cord
contains
!***********************************************************************
!                                                                      *
      subroutine cord(ja, jb, ja1, ipca, jb1)
!-----------------------------------------------
!                                                                      *
!   Computes the MCP coefficients for contributions involving closed   *
!   shells.  The  standard formulae are given in I P Grant, Advances   *
!   in Physics  19 (1970) 747, Eq. (8.33).  In this segment JA1, JB1   *
!   point to the JLIST array, IA1, IB1 to the full list of orbitals.   *
!                                                                      *
!   Call(s) to: CLRX, SPEAK.                                           *
!                                                                      *
!                                           Last update: 15 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use orb_c
      use m_c, only: double, jjc1, jjc2, jjq1, jjq2, jlist, klist, max_orbitals, ncore, npeel, nq1, nq2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use speak_i, only: speak
      use clrx_i, only: clrx
      implicit none (type, external)
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ja
      integer  :: jb
      integer, intent(in) :: ja1
      integer, intent(in) :: ipca
      integer, intent(in) :: jb1
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(real64), parameter :: eps = 1.0e-10_real64
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ia1,ib1,ns,kap1,j1,nqs1,numax,nu,kap2,j2,nqs2,numin
      real(real64) :: x,const,gam
!-----------------------------------------------
!
!
!   Set quantum numbers required.
!
      if (ipca == 2) then
         ia1 = klist(ja1)
      else
         ia1 = jlist(ja1)
      end if
      ib1 = klist(jb1)
!
!   Force IA1 to be greater than IB1
!
      if (ia1 > ib1) then
         ns = ia1
         ia1 = ib1
         ib1 = ns
      end if
!
      kap1 = orbital_kappa(ia1)
      j1 = iabs(kap1)
      nqs1 = nq1(ia1)
!
      if (ia1 == ib1) then
!
!   Case when IA1 .EQ. IB1
!
         x = dble(nqs1*(nqs1 - 1)/2)
         call speak (ja, jb, ia1, ib1, ia1, ib1, 0, x)
         numax = j1 + j1 - 2
         if (numax <= 0) return
         const = dble(nqs1*nqs1/2)
         do nu = 2, numax, 2
            gam = clrx(kap1,nu,kap1)
            x = -const*gam*gam
            if (abs(x) < eps) cycle
            call speak (ja, jb, ia1, ib1, ia1, ib1, nu, x)
         end do
!
!   Case when IA1 .NE. IB1
!
      else
!
         kap2 = orbital_kappa(ib1)
         j2 = abs(kap2)
         nqs2 = nq1(ib1)
         const = dble(nqs1*nqs2)
         call speak (ja, jb, ia1, ib1, ia1, ib1, 0, const)
         numin = abs(j1 - j2)
         numax = j1 + j2 - 1
         if (kap1*kap2 < 0) numin = numin + 1
         do nu = numin, numax, 2
            gam = clrx(kap1,nu,kap2)
            x = -const*gam*gam
            if (abs(x) < eps) cycle
            call speak (ja, jb, ia1, ib1, ib1, ia1, nu, x)
         end do
!
      end if
!
      return
      end subroutine cord
end module cord_i
