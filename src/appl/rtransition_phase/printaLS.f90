!***********************************************************************
!                                                                      *
      subroutine printals (inum_ii,inum_ff,asfa,asfb,i,j,omega,factor)
!                                                                      *
!   This  routine  prints the basic oscillator strength  information   *
!   for transitions between level I and level J in LSJ coupling.       *
!                                                                      *
!   Written by G. Gaigalas,                                            *
!   NIST                                                  May 2011     *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use def_c
      use cuto_c
      use osc_c
      use syma_c
      use jj2lsjbio_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: i, j, inum_ii, inum_ff
      real(double), intent(in) :: asfa, asfb
      real(double), intent(in) :: omega, factor
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!      INTEGER :: IPAR, JPAR, ISIGNA, ISIGNB
      character (len=1) :: im
      real(double), dimension(10) :: dblfac
      real(double) :: d,dd,ab,angs,angsa,ac,am,ams,sa,sb,oscc,oscb,  &
                      sigma,faau,stfac,dll1,omc,fosc,eng,ggfactor,   &
                      acsq,absq,oscm
!
!  *****  DBLFAC(I) = (2I+1)!!
!
      data dblfac/ 3.0000000000d00, 1.5000000000d01, 1.0500000000d02, &
         9.4500000000d02, 1.0395000000d04, 1.3513500000d05, 2.0270250000d06, &
         3.4459425000d07, 6.5472907500d08, 1.3749310575d10/
!-----------------------------------------------
!
      d    = dabs(rlev_ener_2(inum_ff) - rlev_ener_1(inum_ii))
      if(dabs(d - dabs(omega)) > 0.00000001) then
       print*, '                                                      '
       print*, ' INCORRECT INPUT, PROGRAM STOP!!!                     '
       print*, '                                                      '
       print*, ' lsj.lbl files from the jj2lsj runs are inconsistent  '
       print*, ' with the files used in the bioscl calculation.       '
       print*, ' The bioscl calculation was for rci wave functions    '
       print*, ' but the jj2lsj run was for rscf wave functions or    '
       print*, ' vice versa.                                          '
       print*, '                                                      '
       print*, ' If the bioscl calculation was for RCI wave functions '
       print*, ' run jj2lsj to produce labels for RCI wave functions. '
       print*, ' Then rerun rtransition.                                   '
       print*, '                                                      '
       print*, ' If the bioscl calculation was for rscf wave functions'
       print*, ' run jj2lsj to produce labels for rscf wave functions.'
       print*, ' Then rerun rtransition.                                   '
       stop
      end if
      dd    = d*aucm
      angs  =  1.0d08 / dd
      angsa = angs
      if(angs .gt. 2000.d0) then
         sigma = (1.d8/angs)**2
         angsa = angs/(1.d0 + 8342.13d-8 + 206030./(130.d+8 - sigma) &
                       +15997./(38.9d+8 - sigma))
      end if
!
!   Evaluate statistical factors and constants
!
      stfac = iatjpoii(i)
      dll1 = dble (lk+lk+1)
      stfac = stfac/dll1
!
      omc = omega/cvac
      faau = 2.0d00*omc*stfac/dble (iatjpoff(j))
      fosc = cvac*stfac/omc
      eng = omega*factor
      if (ltc(1)) eng = factor/omega
!
      ggfactor = cvac**(2*lk-2)*dble(lk)*dblfac(lk)**2/    &
                 ((2.0d00 * dble(lk)+1.0d00)*            &
                 (dble(lk)+1.0d00)*abs(omega)**(2*lk-1))
!
!   Calculate Einstein A and B coefficients and oscillator strengths
!
      if (kk .eq. 0) then
!
!   Electric multipoles
!
!   In atomic units
!
         acsq = asfa**2
         absq = asfb**2
         ac = acsq*faau
         ab = absq*faau
         oscc = acsq*fosc
         oscb = absq*fosc
         sa=oscc*ggfactor
         sb=oscb*ggfactor
!
!   Convert to SI units if option 5 not set
!
         if (.not. ltc(7)) then
            ac = ac*fasi
            ab = ab*fasi
         endif
!
!   Print information if both AC and AB are greater than CUTOFF
!
         if((abs(ac).ge.cutoff) .and. (abs(ab).ge.cutoff)) then
            im = 'E'
            if (-1.d+9 < eng .and. eng .lt. -1.0d-2) then
               write (32,'(A1)')  ""
               write (32,'(A1)')  ""
               if(dabs(rlev_ener_1(inum_ii)) .lt. 9999.5) then
                  write(32,'(I4,F14.8,2X,A)')                      &
                  iatjpoii(i)-1,rlev_ener_1(inum_ii),              &
                  trim(string_csf1(inum_ii))
               else if(dabs(rlev_ener_1(inum_ii)) .lt. 99999.5) then
                   write(32,'(I4,F14.7,2X,A)')                     &
                   iatjpoii(i)-1,rlev_ener_1(inum_ii),             &
                   trim(string_csf1(inum_ii))
               else
                   write(32,'(I4,F14.6,2X,A)')                     &
                   iatjpoii(i)-1,rlev_ener_1(inum_ii),             &
                   trim(string_csf1(inum_ii))
               end if
               if(dabs(rlev_ener_2(inum_ff)) .lt. 9999.5) then
                  write(32,'(I4,F14.8,2X,A)')                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               else if(dabs(rlev_ener_2(inum_ff)) .lt. 99999.5) then
                  write(32,'(I4,F14.7,2X,A)')                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               else
                  write(32,'(I4,F14.6,2X,A)')                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               end if
               write (32,38) dd,angs,angsa,                        &
                   im,lk,-sb,-oscb,-ab,                            &
                   dabs(dabs(sb)-dabs(sa))/max(dabs(sb),dabs(sa))
               write (32,39)  -sa,-oscc,-ac
            else if ( 1.d+9 .gt. eng .and. eng .gt. 1.d-2) then
               write (32,'(A1)')  ""
               write (32,'(A1)')  ""
               if(dabs(rlev_ener_2(inum_ff)) .lt. 9999.5) then
                  write(32,'(I4,F14.8,2X,A)')                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               else if(dabs(rlev_ener_2(inum_ff)) .lt. 99999.5) then
                  write(32,'(I4,F14.7,2X,A)')                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               else
                  write(32,'(I4,F14.6,2X,A)')                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               end if
               if(dabs(rlev_ener_1(inum_ii)) .lt. 9999.5) then
                  write(32,'(I4,F14.8,2X,A)')                      &
                  iatjpoii(i)-1,rlev_ener_1(inum_ii),              &
                  trim(string_csf1(inum_ii))
               else if(dabs(rlev_ener_1(inum_ii)) .lt. 99999.5) then
                   write(32,'(I4,F14.7,2X,A)')                     &
                   iatjpoii(i)-1,rlev_ener_1(inum_ii),             &
                   trim(string_csf1(inum_ii))
               else
                   write(32,'(I4,F14.6,2X,A)')                     &
                   iatjpoii(i)-1,rlev_ener_1(inum_ii),             &
                   trim(string_csf1(inum_ii))
               end if
               write (32,38) dd,angs,angsa,                        &
                   im,lk,sb,oscb,ab*iatjpoff(j)/iatjpoii(i),       &
                   dabs(dabs(sb)-dabs(sa))/max(dabs(sb),dabs(sa))
               write (32,39)  sa,oscc,ac*iatjpoff(j)/iatjpoii(i)
            end if
         endif
      else
!
!   Magnetic multipoles
!
!   In atomic units
!
         ams = asfa**2
         am = ams*faau
         oscm = ams*fosc
         sa=oscm*ggfactor*cvac*cvac*4
!
!   Convert to SI units if option 5 not set
!
         if (.not. ltc(7)) am = am*fasi
!
!   Print information if AM is greater than CUTOFF
!
         if (abs(am) .ge. cutoff) then
            im = 'M'
            if (-1.d+9 < eng .and. eng .lt. -1.0d-2) then
               write (32,'(A1)')  ""
               write (32,'(A1)')  ""
               if(dabs(rlev_ener_1(inum_ii)) .lt. 9999.5) then
                  write(32,'(I4,F14.8,2X,A)')                      &
                  iatjpoii(i)-1,rlev_ener_1(inum_ii),              &
                  trim(string_csf1(inum_ii))
               else if(dabs(rlev_ener_1(inum_ii)) .lt. 99999.5) then
                   write(32,'(I4,F14.7,2X,A)')                     &
                   iatjpoii(i)-1,rlev_ener_1(inum_ii),             &
                   trim(string_csf1(inum_ii))
               else
                   write(32,'(I4,F14.6,2X,A)')                     &
                   iatjpoii(i)-1,rlev_ener_1(inum_ii),             &
                   trim(string_csf1(inum_ii))
               end if
               if(dabs(rlev_ener_2(inum_ff)) .lt. 9999.5) then
                  write(32,'(I4,F14.8,2X,A)')                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               else if(dabs(rlev_ener_2(inum_ff)) .lt. 99999.5) then
                  write(32,'(I4,F14.7,2X,A)')                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               else
                  write(32,'(I4,F14.6,2X,A)')                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               end if
               write (32,40) dd,angs,angsa,im,lk,-sa,-oscm,-am
            else if ( 1.d+9 .gt. eng .and. eng .gt. 1.d-2) then
               write (32,'(A1)')  ""
               write (32,'(A1)')  ""
               if(dabs(rlev_ener_2(inum_ff)) .lt. 9999.5) then
                  write(32,'(I4,F14.8,2X,A)')                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               else if(dabs(rlev_ener_2(inum_ff)) .lt. 99999.5) then
                  write(32,'(I4,F14.7,2X,A)')                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               else
                  write(32,'(I4,F14.6,2X,A)')                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               end if
               if(dabs(rlev_ener_1(inum_ii)) .lt. 9999.5) then
                  write(32,'(I4,F14.8,2X,A)')                      &
                  iatjpoii(i)-1,rlev_ener_1(inum_ii),              &
                  trim(string_csf1(inum_ii))
               else if(dabs(rlev_ener_1(inum_ii)) .lt. 99999.5) then
                   write(32,'(I4,F14.7,2X,A)')                     &
                   iatjpoii(i)-1,rlev_ener_1(inum_ii),             &
                   trim(string_csf1(inum_ii))
               else
                   write(32,'(I4,F14.6,2X,A)')                     &
                   iatjpoii(i)-1,rlev_ener_1(inum_ii),             &
                   trim(string_csf1(inum_ii))
               end if
               write (32,40) dd,angs,angsa,                        &
                   im,lk,sa,oscm,am*iatjpoff(j)/iatjpoii(i)
            end if
         endif
      endif
      return
   38 format(f11.2,' CM-1',2x,f11.2,' ANGS(VAC)',2x,f11.2,' ANGS(AIR)'/&
         1x,a1,i1,2x,'S = ',1pd12.5,3x,'GF = ',d12.5,3x,'AKI = ',d12.5,&
         3x,'dT = ',0pf8.5)
   39 format(9x,1pd12.5,8x,d12.5,9x,d12.5)
   40 format(f11.2,' CM-1',2x,f11.2,' ANGS(VAC)',2x,f11.2,' ANGS(AIR)'/&
          1x,a1,i1,2x,'S = ',1pd12.5,3x,'GF = ',d12.5,3x,'AKI = ',d12.5)
      end subroutine printals
