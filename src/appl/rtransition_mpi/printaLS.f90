!***********************************************************************
!                                                                      *
      subroutine printals (inum_ii,inum_ff,asfa,asfb,i,j,omega,factor)
      use, intrinsic :: iso_fortran_env, only: real32
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
      use, intrinsic :: iso_fortran_env, only: real64
      use atomic_calculation_state, only: einstein_a_au_to_si, hartree_to_wavenumber_cm, vacuum_speed_of_light_au
      use transition_output_cutoff, only: cutoff
      use transition_matrix_workspace, only: nsdim, hb1, hb2, hc1, hc2, hm1, hm2, jja, jjb, lk, kk, ntdim, &
          xsldr, totc, totb, isldr, isldr1, nint, ninteg, nptr, lab, nkp, kp, ltc
      use level_symmetry_metadata, only: iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
      use transition_lsj_labels, only: rlev_ener_1, rlev_ener_2, string_csf1, string_csf2, nvectoti, nvectotf, &
          iopen_status1, iopen_status2, lev_pos_1, lev_j_1, lev_par_1, lev_pos_2, lev_j_2, lev_par_2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: i, j, inum_ii, inum_ff
      real(real64), intent(in) :: asfa, asfb
      real(real64), intent(in) :: omega, factor
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!      INTEGER :: IPAR, JPAR, ISIGNA, ISIGNB
      character (len=1) :: im
      real(real64) :: dblfac(10)
      real(real64) :: d,dd,ab,angs,angsa,ac,am,ams,sa,sb,oscc,oscb,  &
                      sigma,faau,stfac,dll1,omc,fosc,eng,ggfactor,   &
                      acsq,absq,oscm
!
!  *****  DBLFAC(I) = (2I+1)!!
!
      data dblfac/ 3.0000000000_real64, 1.5000000000e1_real64, 1.0500000000e2_real64, &
         9.4500000000e2_real64, 1.0395000000e4_real64, 1.3513500000e5_real64, 2.0270250000e6_real64, &
         3.4459425000e7_real64, 6.5472907500e8_real64, 1.3749310575e10_real64/
!-----------------------------------------------
!
      d    = dabs(rlev_ener_2(inum_ff) - rlev_ener_1(inum_ii))
      if(dabs(d - dabs(omega)) > 0.00000001_real64) then
       print*, "                                                      "
       print*, " INCORRECT INPUT, PROGRAM STOP!!!                     "
       print*, "                                                      "
       print*, " lsj.lbl files from the jj2lsj runs are inconsistent  "
       print*, " with the files used in the bioscl calculation.       "
       print*, " The bioscl calculation was for rci wave functions    "
       print*, " but the jj2lsj run was for rscf wave functions or    "
       print*, " vice versa.                                          "
       print*, "                                                      "
       print*, " If the bioscl calculation was for RCI wave functions "
       print*, " run jj2lsj to produce labels for RCI wave functions. "
       print*, " Then rerun rtransition.                              "
       print*, "                                                      "
       print*, " If the bioscl calculation was for rscf wave functions"
       print*, " run jj2lsj to produce labels for rscf wave functions."
       print*, " Then rerun rtransition.                              "
       stop
      end if
      dd    = d*hartree_to_wavenumber_cm
      angs  =  1.0e8_real64 / dd
      angsa = angs
      if(angs > 2000.0_real64) then
         sigma = (1.0e8_real64/angs)**2
         angsa = angs/(1.0_real64 + 8342.13e-8_real64 + 206030.0_real32/(130.0e8_real64 - sigma) &
                       +15997.0_real32/(38.9e8_real64 - sigma))
      end if
!
!   Evaluate statistical factors and constants
!
      stfac = iatjpoii(i)
      dll1 = dble (lk+lk+1)
      stfac = stfac/dll1
!
      omc = omega/vacuum_speed_of_light_au
      faau = 2.0_real64*omc*stfac/dble (iatjpoff(j))
      fosc = vacuum_speed_of_light_au*stfac/omc
      eng = omega*factor
      if (ltc(1)) eng = factor/omega
!
      ggfactor = vacuum_speed_of_light_au**(2*lk-2)*dble(lk)*dblfac(lk)**2/    &
                 ((2.0_real64 * dble(lk)+1.0_real64)*            &
                 (dble(lk)+1.0_real64)*abs(omega)**(2*lk-1))
!
!   Calculate Einstein A and B coefficients and oscillator strengths
!
      if (kk == 0) then
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
            ac = ac*einstein_a_au_to_si
            ab = ab*einstein_a_au_to_si
         end if
!
!   Print information if both AC and AB are greater than CUTOFF
!
         if((abs(ac)>=cutoff) .and. (abs(ab)>=cutoff)) then
            im = "E"
            if (-1.0e9_real64 < eng .and. eng < -1.0e-2_real64) then
               write (32,"(A1)")  ""
               write (32,"(A1)")  ""
               if(dabs(rlev_ener_1(inum_ii)) < 9999.5_real32) then
                  write(32,"(I4,F14.8,2X,A)")                      &
                  iatjpoii(i)-1,rlev_ener_1(inum_ii),              &
                  trim(string_csf1(inum_ii))
               else if(dabs(rlev_ener_1(inum_ii)) < 99999.5_real32) then
                   write(32,"(I4,F14.7,2X,A)")                     &
                   iatjpoii(i)-1,rlev_ener_1(inum_ii),             &
                   trim(string_csf1(inum_ii))
               else
                   write(32,"(I4,F14.6,2X,A)")                     &
                   iatjpoii(i)-1,rlev_ener_1(inum_ii),             &
                   trim(string_csf1(inum_ii))
               end if
               if(dabs(rlev_ener_2(inum_ff)) < 9999.5_real32) then
                  write(32,"(I4,F14.8,2X,A)")                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               else if(dabs(rlev_ener_2(inum_ff)) < 99999.5_real32) then
                  write(32,"(I4,F14.7,2X,A)")                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               else
                  write(32,"(I4,F14.6,2X,A)")                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               end if
               write (32,'(f11.2," CM-1",2x,f11.2," ANGS(VAC)",2x,f11.2," ANGS(AIR)"/1x,a1,i1,2x,"S = ",1pd12.5,3x,&
                  &"GF = ",d12.5,3x,"AKI = ",d12.5,3x,"dT = ",0pf8.5)') dd,angs,angsa,                        &
                   im,lk,-sb,-oscb,-ab,                            &
                   dabs(dabs(sb)-dabs(sa))/max(dabs(sb),dabs(sa))
               write (32,'(9x,1pd12.5,8x,d12.5,9x,d12.5)')  -sa,-oscc,-ac
            else if ( 1.0e9_real64 > eng .and. eng > 1.0e-2_real64) then
               write (32,"(A1)")  ""
               write (32,"(A1)")  ""
               if(dabs(rlev_ener_2(inum_ff)) < 9999.5_real32) then
                  write(32,"(I4,F14.8,2X,A)")                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               else if(dabs(rlev_ener_2(inum_ff)) < 99999.5_real32) then
                  write(32,"(I4,F14.7,2X,A)")                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               else
                  write(32,"(I4,F14.6,2X,A)")                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               end if
               if(dabs(rlev_ener_1(inum_ii)) < 9999.5_real32) then
                  write(32,"(I4,F14.8,2X,A)")                      &
                  iatjpoii(i)-1,rlev_ener_1(inum_ii),              &
                  trim(string_csf1(inum_ii))
               else if(dabs(rlev_ener_1(inum_ii)) < 99999.5_real32) then
                   write(32,"(I4,F14.7,2X,A)")                     &
                   iatjpoii(i)-1,rlev_ener_1(inum_ii),             &
                   trim(string_csf1(inum_ii))
               else
                   write(32,"(I4,F14.6,2X,A)")                     &
                   iatjpoii(i)-1,rlev_ener_1(inum_ii),             &
                   trim(string_csf1(inum_ii))
               end if
               write (32,'(f11.2," CM-1",2x,f11.2," ANGS(VAC)",2x,f11.2," ANGS(AIR)"/1x,a1,i1,2x,"S = ",1pd12.5,3x,&
                  &"GF = ",d12.5,3x,"AKI = ",d12.5,3x,"dT = ",0pf8.5)') dd,angs,angsa,                        &
                   im,lk,sb,oscb,ab*iatjpoff(j)/iatjpoii(i),       &
                   dabs(dabs(sb)-dabs(sa))/max(dabs(sb),dabs(sa))
               write (32,'(9x,1pd12.5,8x,d12.5,9x,d12.5)')  sa,oscc,ac*iatjpoff(j)/iatjpoii(i)
            end if
         end if
      else
!
!   Magnetic multipoles
!
!   In atomic units
!
         ams = asfa**2
         am = ams*faau
         oscm = ams*fosc
         sa=oscm*ggfactor*vacuum_speed_of_light_au*vacuum_speed_of_light_au*4
!
!   Convert to SI units if option 5 not set
!
         if (.not. ltc(7)) am = am*einstein_a_au_to_si
!
!   Print information if AM is greater than CUTOFF
!
         if (abs(am) >= cutoff) then
            im = "M"
            if (-1.0e9_real64 < eng .and. eng < -1.0e-2_real64) then
               write (32,"(A1)")  ""
               write (32,"(A1)")  ""
               if(dabs(rlev_ener_1(inum_ii)) < 9999.5_real32) then
                  write(32,"(I4,F14.8,2X,A)")                      &
                  iatjpoii(i)-1,rlev_ener_1(inum_ii),              &
                  trim(string_csf1(inum_ii))
               else if(dabs(rlev_ener_1(inum_ii)) < 99999.5_real32) then
                   write(32,"(I4,F14.7,2X,A)")                     &
                   iatjpoii(i)-1,rlev_ener_1(inum_ii),             &
                   trim(string_csf1(inum_ii))
               else
                   write(32,"(I4,F14.6,2X,A)")                     &
                   iatjpoii(i)-1,rlev_ener_1(inum_ii),             &
                   trim(string_csf1(inum_ii))
               end if
               if(dabs(rlev_ener_2(inum_ff)) < 9999.5_real32) then
                  write(32,"(I4,F14.8,2X,A)")                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               else if(dabs(rlev_ener_2(inum_ff)) < 99999.5_real32) then
                  write(32,"(I4,F14.7,2X,A)")                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               else
                  write(32,"(I4,F14.6,2X,A)")                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               end if
               write (32,'(f11.2," CM-1",2x,f11.2," ANGS(VAC)",2x,f11.2," ANGS(AIR)"/1x,a1,i1,2x,"S = ",1pd12.5,3x,&
                  &"GF = ",d12.5,3x,"AKI = ",d12.5)') dd,angs,angsa,im,lk,-sa,-oscm,-am
            else if ( 1.0e9_real64 > eng .and. eng > 1.0e-2_real64) then
               write (32,"(A1)")  ""
               write (32,"(A1)")  ""
               if(dabs(rlev_ener_2(inum_ff)) < 9999.5_real32) then
                  write(32,"(I4,F14.8,2X,A)")                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               else if(dabs(rlev_ener_2(inum_ff)) < 99999.5_real32) then
                  write(32,"(I4,F14.7,2X,A)")                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               else
                  write(32,"(I4,F14.6,2X,A)")                      &
                  iatjpoff(j)-1,rlev_ener_2(inum_ff),              &
                  trim(string_csf2(inum_ff))
               end if
               if(dabs(rlev_ener_1(inum_ii)) < 9999.5_real32) then
                  write(32,"(I4,F14.8,2X,A)")                      &
                  iatjpoii(i)-1,rlev_ener_1(inum_ii),              &
                  trim(string_csf1(inum_ii))
               else if(dabs(rlev_ener_1(inum_ii)) < 99999.5_real32) then
                   write(32,"(I4,F14.7,2X,A)")                     &
                   iatjpoii(i)-1,rlev_ener_1(inum_ii),             &
                   trim(string_csf1(inum_ii))
               else
                   write(32,"(I4,F14.6,2X,A)")                     &
                   iatjpoii(i)-1,rlev_ener_1(inum_ii),             &
                   trim(string_csf1(inum_ii))
               end if
               write (32,'(f11.2," CM-1",2x,f11.2," ANGS(VAC)",2x,f11.2," ANGS(AIR)"/1x,a1,i1,2x,"S = ",1pd12.5,3x,&
                  &"GF = ",d12.5,3x,"AKI = ",d12.5)') dd,angs,angsa,                        &
                   im,lk,sa,oscm,am*iatjpoff(j)/iatjpoii(i)
            end if
         end if
      end if
      return
      end subroutine printals
