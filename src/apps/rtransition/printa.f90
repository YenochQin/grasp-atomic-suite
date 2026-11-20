!***********************************************************************
!                                                                      *
      subroutine printa(asfa, asfb, i, j, omega, factor, lines, lsame)
      use calculation_file_units, only: summary_output_unit
!                                                                      *
!   This  routine  prints the basic oscillator strength  information   *
!   for transitions between level I and level J.                       *
!                                                                      *
!                                         Last revision: 28 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use atomic_calculation_state, only: einstein_a_au_to_si, einstein_b_au_to_si, pi, vacuum_speed_of_light_au
      use transition_output_cutoff, only: cutoff
      use angular_symmetry_labels, only: right_justified_j_labels => jlbr, parity_labels => jlbp
      use transition_matrix_workspace, only: nsdim, hb1, hb2, hc1, hc2, hm1, hm2, jja, jjb, lk, kk, ntdim, &
          xsldr, totc, totb, isldr, isldr1, nint, ninteg, nptr, lab, nkp, kp, ltc
      use level_symmetry_metadata, only: iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
      use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical, intent(in) :: lsame
      integer :: i
      integer, intent(in) :: j
      integer, intent(inout) :: lines
      real(real64), intent(in) :: asfa, asfb
      real(real64), intent(in) :: omega, factor
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ipar, jpar, isigna, isignb
      real(real64) :: dblfac(10)
      real(real64) :: stfac, dll1, omc, faau, fosc, fbau, eng, ggfactor, &
         acsq, absq, ac, ab, bc, bb, oscc, oscb, sa, sb, ams, am, bm, oscm
      character(len=4) :: jlabi, jlabj, jpari, jparj
      character(len=2) :: f1,f2
!-----------------------------------------------
!
!
!  *****  DBLFAC(I) = (2I+1)!!
!
      data dblfac/ 3.0000000000_real64, 1.5000000000e1_real64, 1.0500000000e2_real64, &
         9.4500000000e2_real64, 1.0395000000e4_real64, 1.3513500000e5_real64, 2.0270250000e6_real64, &
         3.4459425000e7_real64, 6.5472907500e8_real64, 1.3749310575e10_real64/
!
!   Evaluate statistical factors and constants
!
      stfac = iatjpoii(i)
      dll1 = dble(lk + lk + 1)
      stfac = stfac/dll1
!
      omc = omega/vacuum_speed_of_light_au
!     OMC = OMEGA/C
      faau = 2.0_real64*omc*stfac/dble(iatjpoff(j))
      fosc = vacuum_speed_of_light_au*stfac/omc
!     FOSC = C*STFAC/OMC
      fbau = pi*fosc/omega
      eng = omega*factor
      if (ltc(1)) eng = factor/omega
!
!   J/pi labels for levels
!
      jlabi = right_justified_j_labels(iatjpoii(i))
      jlabj = right_justified_j_labels(iatjpoff(j))
      ipar = (iasparii(i) + 3)/2
      jpar = (iasparff(j) + 3)/2
      jpari = parity_labels(ipar)
      jparj = parity_labels(jpar)
!GG
      ggfactor = vacuum_speed_of_light_au**(2*lk - 2)*dble(lk)*dblfac(lk)**2/((2.0_real64*dble(lk) + &
         1.0_real64)*(dble(lk) + 1.0_real64)*abs(omega)**(2*lk - 1))
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
         end if
         if (asfb < 0) then
            isignb = -1
         else
            isignb = 1
         end if

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
!
!   Convert to SI units if option 5 not set
!
         if (.not.ltc(7)) then
            ac = ac*einstein_a_au_to_si
            ab = ab*einstein_a_au_to_si
            bc = bc*einstein_b_au_to_si
            bb = bb*einstein_b_au_to_si
         end if
!
!   Accumulate total of A coefficients
!
         totc(j) = totc(j) + ac
         totb(j) = totb(j) + ab
!
!   Print information if both AC and AB are greater than CUTOFF
!
         if (abs(ac)>=cutoff .and. abs(ab)>=cutoff) then
            if (-1.0e9_real64 < eng .and. eng < -1.0e-2_real64) then
               if (lsame) then
                  f1 = "f "
                  f2 = "f "
               else
                  f1 = "f2"
                  f2 = "f1"
               end if
!
            write (summary_output_unit, '(1x,a2,i3,1x,2a4,a2,i3,1x,2a4,0p,f13.2," C",1p,3d13.5)') &
               f1,ivecff(j),jlabj,jparj,f2,ivecii(i),jlabi,jpari, &
               -eng, -ac, -oscc, -sa !, ASFA
            write (summary_output_unit, '(42x," B",1p,3d13.5)') -ab, -oscb, -sb !, ASFB
            else if ( 1.0e9_real64 > eng .and. eng > 1.0e-2_real64) then
               if (lsame) then
                  f1 = "f "
                  f2 = "f "
               else
                  f1 = "f1"
                  f2 = "f2"
               end if
               write (summary_output_unit,'(1x,a2,i3,1x,2a4,a2,i3,1x,2a4,0p,f13.2," C",1p,3d13.5)') &
                  f1,ivecii(i),jlabi,jpari,f2,ivecff(j), &
                     jlabj,jparj,eng,ac*iatjpoff(j)/iatjpoii(i),     &
                     oscc,sa
               write (summary_output_unit,'(42x," B",1p,3d13.5)') ab*iatjpoff(j)/iatjpoii(i),oscb,sb
            end if
            lines = lines + 3
         end if
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
         end if
         ams = asfa**2
         am = ams*faau
         bm = ams*fbau
         oscm = ams*fosc
         sa = oscm*ggfactor*vacuum_speed_of_light_au*vacuum_speed_of_light_au*4
!GG-end
!
!   Convert to SI units if option 5 not set
!
         if (.not.ltc(7)) then
            am = am*einstein_a_au_to_si
            bm = am*einstein_b_au_to_si
         end if
!
!   Accumulate total of A coefficients
!
         totc(j) = totc(j) + am
         totb(j) = totb(j) + am
!
!   Print information if AM is greater than CUTOFF
!
         if (abs(am) >= cutoff) then
            if (-1.0e9_real64 < eng .and. eng < -1.0e-2_real64) then
               if (lsame) then
                  f1 = "f "
                  f2 = "f "
               else
                  f1 = "f2"
                  f2 = "f1"
               end if
               write (summary_output_unit,'(1x,a2,i3,1x,2a4,a2,i3,1x,2a4,0p,f13.2," M",1p,3d13.5)') &
                  f1,ivecff(j),jlabj,jparj,f2,ivecii(i),&
                     jlabi,jpari,-eng,-am,-oscm,-sa
!     :               ,ASFA (M-value)
!            WRITE (24,*) 'Relative sign',ISIGNA
            else if ( 1.0e9_real64 > eng .and. eng > 1.0e-2_real64) then
               if (lsame) then
                  f1 = "f "
                  f2 = "f "
               else
                  f1 = "f1"
                  f2 = "f2"
               end if
               write (summary_output_unit,'(1x,a2,i3,1x,2a4,a2,i3,1x,2a4,0p,f13.2," M",1p,3d13.5)') &
                  f1,ivecii(i),jlabi,jpari,f2,ivecff(j),&
                     jlabj,jparj,eng,am*iatjpoff(j)/iatjpoii(i),    &
                     oscm,sa
!     :              ,ASFA (M-value)
            end if
            lines = lines+2
         end if
!
      end if
!
      return
!
!cjb  format for highly charged ions         F11.2 -> F13.2
!cjb  format for highly charged ions         F11.2 -> F13.2
      return
!
      end subroutine printa
