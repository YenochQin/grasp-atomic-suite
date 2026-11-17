!***********************************************************************
!                                                                      *
      subroutine csfm (asfa,asfb,lev1,lev2)
!                                                                      *
!   This routine calculates  the CSF Coulomb, Babuskin, and magnetic   *
!   matrix elements for  a transition  between  levels  separated by   *
!   energy OMEGA.                                                      *
!                                                                      *
!   Modified for different initial and final state orbitals            *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: orbital_index_radix
      use biorb_c
      use debug_c, only: cutoff, double, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use eigv_c
      use jlabl_c, labj=>jlbr, labp=>jlbp
      use orb_c
      use osc_c
      use prnt_c, only: byte, double, ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use syma_c, only: double, iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use spme_i, only: spme
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(real64), intent(out) :: asfa, asfb
      integer, intent(in) :: lev1, lev2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, parameter :: key=orbital_index_radix, nca=65536
      integer :: ja, jb, m, nclr, ncup, ncount, nfi, ipar, jpar, &
                 ia, ib, j
      real(real64) :: hcoul1, hbab1, hmag1, eps, couvx, coeff, contr
      character(len=4) :: jlabi, jlabj, jpari, jparj
!-----------------------------------------------------------------------
!
      eps = 1.0e-10_real64
!
      asfa = 0.0_real64
      asfb = 0.0_real64
      ncount = 0
      nfi = 0
!
!Rasa -- start
!
!   J/pi labels for levels
!
      jlabi = labj(iatjpoii(lev1))
      jlabj = labj(iatjpoff(lev2))
      ipar = (iasparii(lev1)+3)/2
      jpar = (iasparff(lev2)+3)/2
      jpari = labp(ipar)
      jparj = labp(jpar)

      if (ldbpr(18)) then
          write (*,302) ivecff(lev2), jlabj, jparj, ivecii(lev1),    &
       jlabi, jpari
          write (*,301)
      end if
! Rasa -- end
      do m = 1,ninteg
        ja = lab(m)
        jb = mod (ja,key)
        ja = ja/key
        nclr = nptr(m)+1
        ncup = nptr(m+1)
        if (nclr > ncup) goto 4
        if (mod(orbital_l(ja)+orbital_l(jb)+kk+lk,2) /= 0) goto 4
        call spme (ja,jb,hcoul1,hbab1,hmag1)
!
        do j = nclr,ncup
!
          ia = isldr(j)
          ib = isldr1(j)
!         IB = MOD (IA,NCA)
!         IA = IA/NCA
!
!  Check on consistency
!
          if (ia>ncfii) then
            goto 131
          else if (ib>ncfff) then
            goto 131
          end if
!
!  Observ that IA and IB refer to the merged list
!  whereas EVECII and EVECFF refers to the initial and
!  final state lists. Therefore we have IB -> IB-NCFII
!
          couvx =evecii(ia+(lev1-1)*ncfii)*                          &
                 evecff(ib+(lev2-1)*ncfff)
!     :           EVECFF(IB-NCFII+(LEV2-1)*NCFFF)
          coeff = xsldr(j)
          if (abs(couvx) > eps) then
            if (kk == 0) then
              asfa = asfa+hcoul1*coeff*couvx
              asfb = asfb+hbab1*coeff*couvx
!Rasa -- start
              if (ldbpr(18) .and. (dabs(hcoul1*coeff*couvx)          &
                .gt. cutoff)) then
                  write (*,300) ia, evecii(ia+(lev1-1)*ncfii),       &
                      ib-ncfii, evecff(ib-ncfii+(lev2-1)*ncfff),     &
                      orbital_principal_n(ja),orbital_symmetry_labels(ja),orbital_principal_n(jb), &
                         & orbital_symmetry_labels(jb),coeff,hcoul1,"C",  &
                      hcoul1*coeff*couvx
                  write(*,303) hbab1,"B",hbab1*coeff*couvx
              end if
!Rasa -- end
            else
              contr=hmag1*coeff*couvx
              asfa = asfa+hmag1*coeff*couvx
!Rasa -- start
              if (ldbpr(18) .and. (abs(contr) > cutoff)) then
                  write (*,300) ia, evecii(ia+(lev1-1)*ncfii),       &
                      ib-ncfii, evecff(ib-ncfii+(lev2-1)*ncfff),     &
                      orbital_principal_n(ja),orbital_symmetry_labels(ja),orbital_principal_n(jb), &
                         & orbital_symmetry_labels(jb),coeff,hmag1," ",   &
                      contr
              end if
!Rasa -- end
            end if
          end if
 131    continue
        end do
    4 end do
!
      return
!
  300 format (i5,2x,d12.5,i5,2x,d12.5, 2(1x,1i2,1a2),1p,d15.6,1x,d15.6, &
      1x,1a,1x,d15.6)
  301 format (4x,"i",3x,"Coeff.(i)",5x,"f",2x,"Coeff.(f)",2x,            &
      "Orb(i)",1x,"Orb(f)",2x,"MCT coeff.",4x,"Radial Factor",          &
      5x,"Contribution")
  302 format ("Upper level ",i3,a4,a4,10x," Lower level",i3,a4,a4)
          write (*,302) ivecff(lev2), jlabj, jparj, lev1, jlabi, jpari
  303 format (62x,1p,d15.6,1x,1a,1x,d15.6)
  304 format(65x,"Total:")
  305 format (62x,1p,16x,1a,1x,d15.6)
!
      end subroutine csfm
