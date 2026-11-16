!***********************************************************************
!                                                                      *
      subroutine newe(j, sgn, nprime, mx, deleps, fail, inv)
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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use int_c
      use orb_c
      use scf_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use outbnd_i, only: outbnd
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: j
      integer, intent(in) :: nprime
      integer, intent(in) :: mx
      integer, intent(out) :: inv
      real(real64), intent(in) :: sgn
      real(real64), intent(inout) :: deleps
      logical, intent(out) :: fail
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(real64), parameter :: p02 = 2.0e-2_real64
      real(real64), parameter :: p05 = 5.0e-2_real64
      real(real64), parameter :: p001 = 1.0e-3_real64
      real(real64), parameter :: p00001 = 1.0e-5_real64
      real(real64), parameter :: d2p5 = 2.5_real64
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(real64) :: eps, abdele, delebe, dn, dnprme, dnpn25, etry, delta
!-----------------------------------------------
!
!
!   Determine if the iterative process has succeeded
!
      if (sgn>0.0_real64 .and. mx==0) then
         fail = .false.
      else
         fail = .true.
      end if
!
!   Inversion counter
!
      inv = 0
!
!   If unsuccessful, obtain a new estimate of the eigenvalue
!
      if (fail) then
!
!   Define quantities used later
!
         eps = orbital_energy_parameters(j)
         abdele = abs(deleps)
         delebe = abs(deleps/eps)
         dn = dble(orbital_principal_n(j))
         dnprme = dble(nprime)
         dnpn25 = (dnprme/dn)**d2p5
!
         if (abs(mx)==1 .and. delebe>p02 .or. mx==0 .and. delebe>=p00001 .and. &
            abdele>=p001) then
    1       continue
            etry = eps + deleps
            if (outbnd(etry) .and. mx/=0) etry = eps*dnpn25
            if (outbnd(etry)) etry = eps - deleps
            if (outbnd(etry)) then
               deleps = 0.5_real64*deleps
               goto 1
            end if
         else if (mx == 0) then
            etry = eps
            inv = 1
            p0 = -p0
            p(:mtp0) = -p(:mtp0)
            q(:mtp0) = -q(:mtp0)
            fail = .false.
         else if (sgn < 0.0_real64) then
            etry = eps*dnpn25
            if (outbnd(etry)) etry = 0.5_real64*(epsmax + eps)
            if (outbnd(etry)) etry = 0.5_real64*(epsmin + eps)
         else if (mx < 0) then
            delta = 1.0_real64 - eps/epsmax
            epsmax = eps
            if (delta < p05) then
               emax = emax*dnpn25
            else
               emax = eps*dnpn25
            end if
            etry = emax
            if (outbnd(etry)) etry = 0.5_real64*(epsmax + epsmin)
         else if (mx > 0) then
            delta = 1.0_real64 - epsmin/eps
            epsmin = eps
            if (delta < p05) then
               emin = emin*dnpn25
            else
               emin = eps*dnpn25
            end if
            etry = emin
            if (outbnd(etry)) etry = 0.5_real64*(epsmax + epsmin)
         end if
         orbital_energy_parameters(j) = etry
      end if
!
      return
      end subroutine newe
