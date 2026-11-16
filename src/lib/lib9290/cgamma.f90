module cgamma_i
    implicit none (type, external)
    private
    public :: cgamma
contains
!***********************************************************************
!                                                                      *
      subroutine cgamma(argr, argi, resr, resi)
!                                                                      *
!   This subroutine returns in RES the complex Gamma function of the   *
!   complex argument ARG.  The suffixes R and I respectively distin-   *
!   guish the real and imaginary parts OF both RES and ARG.            *
!                                                                      *
!   Only RESR is nonzero if ARGI is zero.                              *
!                                                                      *
!   The  ARCTAN function required must return angles (in radians) in   *
!   the range  [0,2*\pi).                                              *
!                                                                      *
!   Subprogram(S) required: ARCTAN.                                    *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use def_c
      use arctan_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(real64), intent(in) :: argi, argr
      real(real64), intent(out) :: resr, resi
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(real64) :: fn(7), fd(7)
      real(real64) :: hlntpi, twoi, diff, argum, clngi, facneg, argur, ovlfac, &
         clngr, fac, obasq, argui, argui2, ovlfr, ovlfi, termr, termi, argur2, &
         obasqr, obasqi, zfacr, zfaci
      logical,save :: first, negarg
!----------------------------------------------------------------------*
!                                                                      *
!   These are the Bernoulli numbers B02, B04, ..., B14, expressed as   *
!   rational numbers. From abramowitz and stegun, p. 810.              *
!
      data fn/ 1.0_real64, -1.0_real64, 1.0_real64, -1.0_real64, 5.0_real64, -691.0_real64, 7.0_real64/
      data fd/ 6.0_real64, 30.0_real64, 42.0_real64, 30.0_real64, 66.0_real64, 2730.0_real64, 6.0_real64/
!
!----------------------------------------------------------------------*
!
      data hlntpi/ 1.0_real64/
!
      data first/ .true./
!
!   On the first entry to this routine, set up the constants required
!   for the reflection formula (cf. Abramowitz and Stegun 6.1.17) and
!   Stirling's approximation (cf. Abramowitz and Stegun 6.1.40).
!
      if (first) then
!
         hlntpi = 0.5_real64*log(pi + pi)
!
         do i = 1, 7
            fn(i) = fn(i)/fd(i)
            twoi = dble(i + i)
            fn(i) = fn(i)/(twoi*(twoi - 1.0_real64))
         end do
!
         first = .false.
!
      endif
!
!   Cases where the argument is real
!
      if (argi == 0.0_real64) then
!
!   Cases where the argument is real and negative
!
         if (argr <= 0.0_real64) then
!
!   Stop with an error message if the argument is too near a pole
!
            diff = abs(dble(nint(argr)) - argr)
            if (diff <= precis + precis) then
               write (*, 300) argr, argi
               stop
            else
!
!   Otherwise use the reflection formula (Abramowitz and Stegun 6.1.17)
!   to ensure that the argument is suitable for Stirling's formula
!
               argum = pi/(-argr*sin(pi*argr))
               if (argum < 0.0_real64) then
                  argum = -argum
                  clngi = pi
               else
                  clngi = 0.0_real64
               endif
               facneg = log(argum)
               argur = -argr
               negarg = .true.
!
            endif
!
!   Cases where the argument is real and positive
!
         else
!
            clngi = 0.0_real64
            argur = argr
            negarg = .false.
!
         endif
!
!   Use Abramowitz and Stegun formula 6.1.15 to ensure that
!   the argument in Stirling's formula is greater than 10
!
         ovlfac = 1.0_real64
    2    continue
         if (argur < 10.0_real64) then
            ovlfac = ovlfac*argur
            argur = argur + 1.0_real64
            go to 2
         endif
!
!   Now use Stirling's formula to compute Log (Gamma (ARGUM))
!
         clngr = (argur - 0.5_real64)*log(argur) - argur + hlntpi
         fac = argur
         obasq = 1.0_real64/(argur*argur)
         do i = 1, 7
            fac = fac*obasq
            clngr = clngr + fn(i)*fac
         end do
!
!   Include the contributions from the recurrence and reflection
!   formulae
!
         clngr = clngr - log(ovlfac)
         if (negarg) clngr = facneg - clngr
!
      else
!
!   Cases where the argument is complex
!
         argur = argr
         argui = argi
         argui2 = argui*argui
!
!   Use the recurrence formula (Abramowitz and Stegun 6.1.15)
!   to ensure that the magnitude of the argument in Stirling's
!   formula is greater than 10
!
         ovlfr = 1.0_real64
         ovlfi = 0.0_real64
    4    continue
         argum = sqrt(argur*argur + argui2)
         if (argum < 10.0_real64) then
            termr = ovlfr*argur - ovlfi*argui
            termi = ovlfr*argui + ovlfi*argur
            ovlfr = termr
            ovlfi = termi
            argur = argur + 1.0_real64
            go to 4
         endif
!
!   Now use Stirling's formula to compute Log (Gamma (ARGUM))
!
         argur2 = argur*argur
         termr = 0.5_real64*log(argur2 + argui2)
         termi = arctan(argui,argur)
         clngr = (argur - 0.5_real64)*termr - argui*termi - argur + hlntpi
         clngi = (argur - 0.5_real64)*termi + argui*termr - argui
         fac = (argur2 + argui2)**(-2)
         obasqr = (argur2 - argui2)*fac
         obasqi = -2.0_real64*argur*argui*fac
         zfacr = argur
         zfaci = argui
         do i = 1, 7
            termr = zfacr*obasqr - zfaci*obasqi
            termi = zfacr*obasqi + zfaci*obasqr
            fac = fn(i)
            clngr = clngr + termr*fac
            clngi = clngi + termi*fac
            zfacr = termr
            zfaci = termi
         end do
!
!   Add in the relevant pieces from the recurrence formula
!
         clngr = clngr - 0.5_real64*log(ovlfr*ovlfr + ovlfi*ovlfi)
         clngi = clngi - arctan(ovlfi,ovlfr)
!
      endif
!
!   Now exponentiate the complex Log Gamma function to get
!   the complex Gamma function
!
      if (clngr<=expmax .and. clngr>=expmin) then
         fac = exp(clngr)
      else
         write (*, 301) clngr
         stop
      endif
      resr = fac*cos(clngi)
      resi = fac*sin(clngi)
!
      return
!
  300 format('CGAMMA: Argument (',1p,1d19.12,',',1d19.12,')',&
         ' too close to a pole.')
  301 format('CGAMMA: Argument to exponential function',' (',1p,1d19.12,&
         ') out of range.')
      return
!
      end subroutine cgamma
end module cgamma_i
