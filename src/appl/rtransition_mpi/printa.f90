!***********************************************************************
!                                                                      *
      subroutine printa(asfa, asfb, i, j, omega, factor, lines, lsame)
!                                                                      *
!   This  routine  prints the basic oscillator strength  information   *
!   for transitions between level I and level J.                       *
!                                                                      *
!                                         Last revision: 28 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:39:37   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use def_c
      use cuto_c
      use jlabl_c, labj=>jlbr, labp=> jlbp
      use osc_c
      use syma_c
      use prnt_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical, intent(in) :: lsame
      integer  :: i, j
      integer, intent(inout) :: lines
      real(double), intent(in) :: asfa, asfb
      real(double), intent(in) :: omega, factor
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ipar, jpar, isigna, isignb
      real(double), dimension(10) :: dblfac
      real(double) :: stfac, dll1, omc, faau, fosc, fbau, eng, ggfactor, &
         acsq, absq, ac, ab, bc, bb, oscc, oscb, sa, sb, ams, am, bm, oscm
      character(len=4) :: jlabi, jlabj, jpari, jparj
      character(len=2) :: f1,f2
!-----------------------------------------------
!
!
      data dblfac/ 3.0000000000d00, 1.5000000000d01, 1.0500000000d02, &
         9.4500000000d02, 1.0395000000d04, 1.3513500000d05, 2.0270250000d06, &
         3.4459425000d07, 6.5472907500d08, 1.3749310575d10/
!
!   Evaluate statistical factors and constants
!
      stfac = iatjpoii(i)
      dll1 = dble(lk + lk + 1)
      stfac = stfac/dll1
!
      omc = omega/vacuum_speed_of_light_au
!     OMC = OMEGA/C
      faau = 2.0d00*omc*stfac/dble(iatjpoff(j))
      fosc = vacuum_speed_of_light_au*stfac/omc
!     FOSC = C*STFAC/OMC
      fbau = pi*fosc/omega
      eng = omega*factor
      if (ltc(1)) eng = factor/omega
!
!   J/pi labels for levels
!
      jlabi = labj(iatjpoii(i))
      jlabj = labj(iatjpoff(j))
      ipar = (iasparii(i) + 3)/2
      jpar = (iasparff(j) + 3)/2
      jpari = labp(ipar)
      jparj = labp(jpar)
!GG
      ggfactor = vacuum_speed_of_light_au**(2*lk - 2)*dble(lk)*dblfac(lk)**2/((2.0d00*dble(lk) + &
         1.0d00)*(dble(lk) + 1.0d00)*abs(omega)**(2*lk - 1))
!GG-end
!
!   Calculate Einstein A and B coefficients and oscillator strengths
!
      if (kk == 0) then
!
!   Electric multipoles
!
!   In atomic units
!
         if (asfa < 0) then
            isigna = -1
         else
            isigna = 1
         endif
         if (asfb < 0) then
            isignb = -1
         else
            isignb = 1
         endif

         acsq = asfa**2
         absq = asfb**2
         ac = acsq*faau
         ab = absq*faau
         bc = acsq*fbau
         bb = absq*fbau
         oscc = acsq*fosc
         oscb = absq*fosc
         sa = oscc*ggfactor
         sb = oscb*ggfactor
!GG-end
!
!   Convert to SI units if option 5 not set
!
         if (.not.ltc(7)) then
            ac = ac*einstein_a_au_to_si
            ab = ab*einstein_a_au_to_si
            bc = bc*einstein_b_au_to_si
            bb = bb*einstein_b_au_to_si
         endif
!
!   Accumulate total of A coefficients
!
         totc(j) = totc(j) + ac
         totb(j) = totb(j) + ab
!
!   Print information if both AC and AB are greater than CUTOFF
!
         if (abs(ac)>=cutoff .and. abs(ab)>=cutoff) then
            if (eng .lt. -1.0d-9) then
               if (lsame) then
                  f1 = 'f '
                  f2 = 'f '
               else
                  f1 = 'f2'
                  f2 = 'f1'
               end if
!
               write (24, 300) f1,ivecff(j),jlabj,jparj,f2,ivecii(i),jlabi,jpari, &
               -eng, -ac, -oscc, -sa !, ASFA
               write (24, 301) -ab, -oscb, -sb !, ASFB
            else if (eng .gt. 1.d-9) then
               if (lsame) then
                  f1 = 'f '
                  f2 = 'f '
               else
                  f1 = 'f1'
                  f2 = 'f2'
               end if
               write (24,300) f1,ivecii(i),jlabi,jpari,f2,ivecff(j), &
                     jlabj,jparj,eng,ac*iatjpoff(j)/iatjpoii(i),     &
                     oscc,sa
               write (24,301) ab*iatjpoff(j)/iatjpoii(i),oscb,sb
            end if
            lines = lines + 3
         endif
!
      else
!
!   Magnetic multipoles
!
!   In atomic units
!
         if (asfa < 0) then
            isigna = -1
         else
            isigna = 1
         endif
         ams = asfa**2
         am = ams*faau
         bm = ams*fbau
         oscm = ams*fosc
         sa = oscm*ggfactor*vacuum_speed_of_light_au*vacuum_speed_of_light_au*4
!
!   Convert to SI units if option 5 not set
!
         if (.not.ltc(7)) then
            am = am*einstein_a_au_to_si
            bm = am*einstein_b_au_to_si
         endif
!
!   Accumulate total of A coefficients
!
         totc(j) = totc(j) + am
         totb(j) = totb(j) + am
!
!   Print information if AM is greater than CUTOFF
!
         if (abs(am) >= cutoff) then
            if (eng .le. -1.0d-9) then
               if (lsame) then
                  f1 = 'f '
                  f2 = 'f '
               else
                  f1 = 'f2'
                  f2 = 'f1'
               end if
               write (24,302) f1,ivecff(j),jlabj,jparj,f2,ivecii(i),&
                     jlabi,jpari,-eng,-am,-oscm,-sa
            else if (eng .gt. 1.d-09) then
               if (lsame) then
                  f1 = 'f '
                  f2 = 'f '
               else
                  f1 = 'f1'
                  f2 = 'f2'
               end if
               write (24,302) f1,ivecii(i),jlabi,jpari,f2,ivecff(j),&
                     jlabj,jparj,eng,am*iatjpoff(j)/iatjpoii(i),    &
                     oscm,sa
            end if
            lines = lines+2
         endif
!
      endif
!
      return
!
!cjb  format for highly charged ions         F11.2 -> F13.2
! 300 FORMAT(1X,A2,I3,1X,2A4,A2,I3,1X,2A4,0P,F11.2,' C',1P,  &
  300 format(1x,a2,i3,1x,2a4,a2,i3,1x,2a4,0p,f13.2,' C',1p,  &
           3d13.5)
  301 format(42x,' B',1p,3d13.5)
!cjb  format for highly charged ions         F11.2 -> F13.2
! 302 FORMAT(1X,A2,I3,1X,2A4,A2,I3,1X,2A4,0P,F11.2,' M',1P,  &
  302 format(1x,a2,i3,1x,2a4,a2,i3,1x,2a4,0p,f13.2,' M',1p,  &
           3d13.5)
      return
!
      end subroutine printa
