!***********************************************************************
!                                                                      *
      subroutine spme(i, j, hcoul, hbab, hmag)
!                                                                      *
!   This routine calculates the reduced matrix elements for pair I,J   *
!   in either Coulomb/Babuskin gauge or for magnetic case.             *
!                                                                      *
!   These are defined in the  Brink and Satchler sense --- i.e. com-   *
!   patible with Pyper et al.'s MCT paper but not with Grant's radi-   *
!   ative transitions paper.                                           *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use bess_c, only: bessj, bessn, bj, double, max_radial_points, tc, td, wij
      use biorb_c
      use debug_c, only: cutoff, double, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use def_c
      use grid_c
      use orb_c
      use osc_c
      use tatb_c, only: double, mtp, radial_workspace_capacity, ta, tb
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use clrx_i, only: clrx
      use quad_i, only: quad
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: i,j
      real(real64), intent(out) :: hcoul, hbab, hmag
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nkji, nkjj, ipj, naki, nakj, ii, lp, lm
      real(real64) :: eps, hgauge, tji, tjj, fact, form, fl, flp, dfki, dfkj, &
         clp, clm, ciplp, cimlp, ciplm, cimlm, value, cjl, cip, tarom, esther, &
         tekel, peres, tamar, enoch, seth, shem, darius, cyrus, daniel, boaz, &
         esau, aaron, gimel
!-----------------------------------------------

!
      eps = 1.0e-10_real64
!
      if (ldbpr(12) .or. ldbpr(13) .or. ldbpr(14) .or. ldbpr(15)) write (99, &
         303) i, j
!
      hmag = 0.0_real64
      hcoul = 0.0_real64
      hgauge = 0.0_real64
      hbab = 0.0_real64
!
!  Evaluate factor multiplying Mbar(a,b)
!
      nkji = orbital_two_j(i)
      nkjj = orbital_two_j(j)
      tji = dble(nkji)
      tjj = dble(nkjj)
      fact = 1.0_real64
      ipj = (nkji + 1)/2 + nkjj
      if (mod(ipj,2) /= 0) fact = -fact
      naki = orbital_kappa(i)
      nakj = orbital_kappa(j)
      form = fact*sqrt(tjj + 1.0_real64)*clrx(nakj,lk,naki)
      fl = dble(lk)
      flp = fl + 1.0_real64
      dfki = dble(naki)
      dfkj = dble(nakj)
!
      mtp = min(mfii(nnii(i)),mfff(nnff(j)))
!
      if (kk == 0) then
!
         if (abs(form) > eps) then
!
!  To pick the right initial and final state radial functions convert
!  from orbital order of the merged list to the orbital orders of the
!  initial and final state lists, respectively.
!  N(I) -> NNII(I), N(J) -> NNFF(J)
!
            do ii = 1, radial_point_count
               tb(ii) = pfii(ii,nnii(i))*qfff(ii,nnff(j)) + qfii(ii,nnii(i))*&
                  pfff(ii,nnff(j))
               tc(ii) = pfii(ii,nnii(i))*qfff(ii,nnff(j)) - qfii(ii,nnii(i))*&
                  pfff(ii,nnff(j))
! zou
               tb(ii) = tb(ii)*speed_of_light_au/vacuum_speed_of_light_au
               tc(ii) = tc(ii)*speed_of_light_au/vacuum_speed_of_light_au
               td(ii) = pfii(ii,nnii(i))*pfff(ii,nnff(j)) + qfii(ii,nnii(i))*&
                  qfff(ii,nnff(j))*(speed_of_light_au/vacuum_speed_of_light_au)**2
!           TD(II) = PFII(II,NNII(I))*PFFF(II,NNFF(J)) +
!    :               QFII(II,NNII(I))*QFFF(II,NNFF(J))
! zou
            end do
            lp = lk + 1
            lm = lk - 1
!
!  Calculate Coulomb coefficients
!
            clp = sqrt(fl/flp)
            clm = -sqrt(flp/fl)
            ciplp = clp*(dfki - dfkj)
            cimlp = clp*flp
            ciplm = clm*(dfki - dfkj)
            cimlm = -clm*fl
!
!  Tabulate Coulomb integrand
!
            ta(1) = 0.0_real64
            ta(2:mtp) = radial_grid_derivative(2:mtp)*(bj(2:mtp,3)*(ciplp*tb(2:mtp)+cimlp*tc(2:mtp)&
               )+bj(2:mtp,1)*(ciplm*tb(2:mtp)+cimlm*tc(2:mtp)))
            call quad (value)
            hcoul = form*value
!
!  Calculate gauge dependent coefficients
!
            cjl = -(fl + flp)
            cip = dfki - dfkj
!
!  Tabulate gauge dependent integrand
!
            ta(1) = 0.0_real64
            ta(2:mtp) = radial_grid_derivative(2:mtp)*(bj(2:mtp,2)*(cjl*td(2:mtp))+bj(2:mtp,3)*(cip&
               *tb(2:mtp)+flp*tc(2:mtp))+bj(2:mtp,1)*(cip*tb(2:mtp)-fl*tc(2:mtp&
               )))
!
!  Print gauge dependent integrand if requested
!
            if (ldbpr(13)) write (99, 300) i, j, (ii,ta(ii),ii=1,radial_point_count)
            call quad (value)
            hgauge = form*value
            hbab = hcoul + sqrt(flp/fl)*hgauge
!
!  Print Coulomb and gauge dependent integrals if requested
!
            if (ldbpr(12)) write (99, 301) i, j, hcoul, hgauge, hbab
!
!  Calculate the contibutions from various terms if requested
!
            if (ldbpr(14) .or. ldbpr(15)) then
!
               ta(1) = 0.0_real64
               ta(2:mtp) = td(2:mtp)*bj(2:mtp,2)*radial_grid_derivative(2:mtp)
               if (ldbpr(15)) write (99, 304) (ii,ta(ii),ii=1,radial_point_count)
               call quad (tarom)
               if (ldbpr(14)) write (99, 305) tarom
               esther = -cjl
               if (ldbpr(14)) write (99, 306) esther
               tekel = esther*tarom
               if (ldbpr(14)) write (99, 307) tekel
!
               ta(1) = 0.0_real64
               ta(2:mtp) = -tb(2:mtp)*bj(2:mtp,3)*radial_grid_derivative(2:mtp)
               if (ldbpr(15)) write (99, 304) (ii,ta(ii),ii=1,radial_point_count)
               call quad (peres)
               if (ldbpr(14)) then
                  write (99, 308) peres
                  write (99, 306) cip
               end if
               tamar = cip*peres
               if (ldbpr(14)) write (99, 307) tamar
!
               ta(1) = 0.0_real64
               ta(2:mtp) = -tb(2:mtp)*bj(2:mtp,1)*radial_grid_derivative(2:mtp)
               if (ldbpr(15)) write (99, 304) (ii,ta(ii),ii=1,radial_point_count)
               call quad (enoch)
               if (ldbpr(14)) then
                  write (99, 309) enoch
                  write (99, 306) cip
               end if
               seth = cip*enoch
               if (ldbpr(14)) write (99, 307) seth
!
               ta(1) = 0.0_real64
               ta(2:mtp) = bj(2:mtp,1)*tc(2:mtp)*radial_grid_derivative(2:mtp)
               if (ldbpr(15)) write (99, 304) (ii,ta(ii),ii=1,radial_point_count)
               call quad (shem)
               if (ldbpr(14)) then
                  write (99, 310) shem
                  write (99, 306) fl
               end if
               darius = fl*shem
               if (ldbpr(14)) write (99, 307) darius
!
               ta(1) = 0.0_real64
               ta(2:mtp) = -bj(2:mtp,3)*tc(2:mtp)*radial_grid_derivative(2:mtp)
               if (ldbpr(15)) write (99, 304) (ii,ta(ii),ii=1,radial_point_count)
               call quad (cyrus)
               if (ldbpr(14)) write (99, 311) cyrus
               daniel = flp
               if (ldbpr(14)) write (99, 306) daniel
               boaz = daniel*cyrus
               if (ldbpr(14)) write (99, 307) boaz
!
               esau = tamar + seth + darius + boaz
               if (ldbpr(14)) write (99, 312) esau, tekel
               aaron = esau + tekel
               gimel = -value
               if (ldbpr(14)) then
                  write (99, 313) aaron, gimel
                  write (99, 314)
                  write (99, 315) peres, cip, tamar
                  write (99, 316) enoch, cip, seth
                  write (99, 317) shem, fl, darius
                  write (99, 318) cyrus, daniel, boaz
                  write (99, 319) esau
                  write (99, 320) tarom, esther, tekel
                  write (99, 321) gimel, aaron
               end if
!
            end if
!
         end if
!
         if (ldbpr(12)) write (99, 301) i, j, hcoul, hgauge, hbab
!
      else
!
!   Tabulate magnetic integrand
!
         if (abs(form) > eps) then
!
            ta(1) = 0.0_real64
            do ii = 2, mtp
!zou
               ta(ii) = (pfii(ii,nnii(i))*qfff(ii,nnff(j))+qfii(ii,nnii(i))*&
                  pfff(ii,nnff(j)))*bj(ii,2)*radial_grid_derivative(ii)*speed_of_light_au/vacuum_speed_of_light_au
!           TA(II) =  (PFII(II,NNII(I))*QFFF(II,NNFF(J)) +
!    :      QFII(II,NNII(I))*PFFF(II,NNFF(J)))*BJ(II,2)*RP(II)
!zou
            end do
            call quad (value)
            if (ldbpr(14)) write (99, 322) value
            hmag = -value*(fl + flp)*(dfki + dfkj)*form/sqrt(fl*flp)
!
         end if
!
         if (ldbpr(12)) write (99, 302) hmag
!
      end if
!
      return
!
  300 format(/,1x,20x,"Local form of gauge dependent integral for"," orbitals",&
         i3," and",i3,/,/,100(1x,7(i3,2x,1p,d11.3,2x),/))
  301 format(/,1x,"Orbital pair (",i2,",",i2,") integrals:","Coul. gauge ",1p,d&
         13.3,3x,"Gauge contribution",d13.3,3x,"Bab.gauge",d13.3)
  302 format(/,1x,"Magnetic single particle matrix element ",1p,d13.3)
  303 format(/,1x,30("+"),"   SPME called for orbitals ",i4,"   and",i4,3x,30(&
         "+"))
  304 format(/,100(1x,7(i3,2x,1p,d11.3,2x),/),/)
  305 format(/," Integral J(L)  =  ",1p,d28.12)
  306 format(" Factor multiplying this ",1p,d20.12)
  307 format(" Resultant contribution  ",1p,d20.12,/)
  308 format(" I + (L+1) integral  = ",1p,d24.12)
  309 format(" I + (L-1) integral  = ",1p,d24.12)
  310 format(" I - (L-1) integral  = ",1p,d24.12)
  311 format(" I - (L+1) integral  = ",1p,d24.12)
  312 format(" These last four add up to ",1p,d20.12,&
         "    compared to the first which is ",d20.12)
  313 format(" These all add up to ",1p,d20.12," which should equal ",d20.12)
  314 format(/,18x,"Integral",11x,"Coefficient",8x,"Contribution")
  315 format(" I + (L+1)",1p,3d20.10)
  316 format(" I + (L-1)",1p,3d20.10)
  317 format(" I - (L-1)",1p,3d20.10)
  318 format(" I - (L+1)",1p,3d20.10)
  319 format(54x,16("-"),/,50x,1p,d20.10)
  320 format(" J (L)",4x,1p,3d20.10)
  321 format(54x,16("-"),/," Correct figure",1p,d20.10,15x,d20.10)
  322 format(" I + (L)  integral",1p,d20.10)
      return
!
      end subroutine spme
