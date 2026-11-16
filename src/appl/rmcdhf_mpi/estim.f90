!***********************************************************************
!                                                                      *
      subroutine estim(j)
!                                                                      *
!   This  subprogram implements  Part 1 of Algorithm 7.1 of C Froese   *
!   Fischer, Comput Phys Rep, 3 (1986) 320-321.                        *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 26 Sep 1993   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use def_c
      use grid_c
      use orb_c
      use scf_c
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: npj, nakabs
      real(double) :: alpha, csq, fnrel, fkabs, fkap2, zalpha, gamma, ebym
!-----------------------------------------------
!
!   Initializations
!
      alpha = 1.0d00/c
      csq = c*c
      npj = np(j)
      nakabs = abs(nak(j))
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
      zinf = z + dble((-nelec) + 1)
!
!   Set the lower bound
!
      zalpha = zinf*alpha
      if (zalpha < fkabs) then
         gamma = sqrt(fkap2 - zalpha*zalpha)
         ebym = 1.0d00/sqrt(1.0d00 + (zalpha/(gamma + fnrel + 0.5d00))**2)
         epsmin = (1.0d00 - ebym)*csq
      else
         epsmin = 0.25d00*csq/dble(npj*npj)
      endif
      emin = epsmin
!
!   Set the upper bound
!
      zalpha = z*alpha
      if (zalpha < fkabs) then
         gamma = sqrt(fkap2 - zalpha*zalpha)
         ebym = 1.0d00/sqrt(1.0d00 + (zalpha/(gamma + fnrel - 0.5d00))**2)
         epsmax = (1.0d00 - ebym)*csq
      else
         epsmax = csq + csq
      endif
      emax = epsmax
!
      return
      end subroutine estim
