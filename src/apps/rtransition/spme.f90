!***********************************************************************
!                                                                      *
      subroutine spme(first_orbital, second_orbital, hcoul, hbab, hmag)
      use debug_output_settings, only: angular_debug_unit, ensure_debug_output
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
      use bessel_function_workspace, only: bessj, bessn, bj, max_radial_points, tc, td, wij
      use biorthogonal_orbital_sets, only: max_orbitals, nhff, nhii, eff, gamaff, eii, gamaii, ncfff, nwff, ncfii, nwii, npff, &
          nakff, npii, nakii, nklff, nkjff, nklii, nkjii, nnff, nnii
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use atomic_calculation_state, only: max_radial_points, speed_of_light_au, vacuum_speed_of_light_au
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use orbital_csf_basis, only: max_orbitals, orbital_kappa, orbital_two_j
      use transition_matrix_workspace, only: nsdim, hb1, hb2, hc1, hc2, hm1, hm2, jja, jjb, lk, kk, ntdim, &
          xsldr, totc, totb, isldr, isldr1, nint, ninteg, nptr, lab, nkp, kp, ltc
      use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use relativistic_three_j, only: clrx
      use radial_quadrature, only: quad
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: first_orbital,second_orbital
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
      if (ldbpr(12) .or. ldbpr(13) .or. ldbpr(14) .or. ldbpr(15)) then
         call ensure_debug_output()
         write (angular_debug_unit, &
         '(/,1x,30("+"),"   SPME called for orbitals ",i4,"   and",i4,3x,30("+"))') first_orbital, second_orbital
      end if
!
      hmag = 0.0_real64
      hcoul = 0.0_real64
      hgauge = 0.0_real64
      hbab = 0.0_real64
!
!  Evaluate factor multiplying Mbar(a,b)
!
      nkji = orbital_two_j(first_orbital)
      nkjj = orbital_two_j(second_orbital)
      tji = dble(nkji)
      tjj = dble(nkjj)
      fact = 1.0_real64
      ipj = (nkji + 1)/2 + nkjj
      if (mod(ipj,2) /= 0) fact = -fact
      naki = orbital_kappa(first_orbital)
      nakj = orbital_kappa(second_orbital)
      form = fact*sqrt(tjj + 1.0_real64)*clrx(nakj,lk,naki)
      fl = dble(lk)
      flp = fl + 1.0_real64
      dfki = dble(naki)
      dfkj = dble(nakj)
!
      mtp = min(mfii(nnii(first_orbital)),mfff(nnff(second_orbital)))
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
               tb(ii) = pfii(ii,nnii(first_orbital))*qfff(ii,nnff(second_orbital)) + qfii(ii,nnii(first_orbital))*&
                  pfff(ii,nnff(second_orbital))
               tc(ii) = pfii(ii,nnii(first_orbital))*qfff(ii,nnff(second_orbital)) - qfii(ii,nnii(first_orbital))*&
                  pfff(ii,nnff(second_orbital))
! zou
               tb(ii) = tb(ii)*speed_of_light_au/vacuum_speed_of_light_au
               tc(ii) = tc(ii)*speed_of_light_au/vacuum_speed_of_light_au
               td(ii) = pfii(ii,nnii(first_orbital))*pfff(ii,nnff(second_orbital)) + qfii(ii,nnii(first_orbital))*&
                  qfff(ii,nnff(second_orbital))*(speed_of_light_au/vacuum_speed_of_light_au)**2
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
            if (ldbpr(13)) then
               write (angular_debug_unit, '(/,1x,20x,"Local form of gauge dependent integral for"," orbitals",i3,&
               &" and",i3,/,/,100(1x,7(i3,2x,1p,d11.3,2x),/))') first_orbital, second_orbital, (ii,ta(ii),ii=1,radial_point_count)
            end if
            call quad (value)
            hgauge = form*value
            hbab = hcoul + sqrt(flp/fl)*hgauge
!
!  Print Coulomb and gauge dependent integrals if requested
!
            if (ldbpr(12)) then
               write (angular_debug_unit, '(/,1x,"Orbital pair (",i2,",",i2,") integrals:","Coul. gauge ",1p,d13.3,3x,&
               &"Gauge contribution",d13.3,3x,"Bab.gauge",d13.3)') first_orbital, second_orbital, hcoul, hgauge, hbab
            end if
!
!  Calculate the contibutions from various terms if requested
!
            if (ldbpr(14) .or. ldbpr(15)) then
!
               ta(1) = 0.0_real64
               ta(2:mtp) = td(2:mtp)*bj(2:mtp,2)*radial_grid_derivative(2:mtp)
               if (ldbpr(15)) write (angular_debug_unit, '(/,100(1x,7(i3,2x,1p,d11.3,2x),/),/)') (ii,ta(ii),ii=1,radial_point_count)
               call quad (tarom)
               if (ldbpr(14)) write (angular_debug_unit, '(/," Integral J(L)  =  ",1p,d28.12)') tarom
               esther = -cjl
               if (ldbpr(14)) write (angular_debug_unit, '(" Factor multiplying this ",1p,d20.12)') esther
               tekel = esther*tarom
               if (ldbpr(14)) write (angular_debug_unit, '(" Resultant contribution  ",1p,d20.12,/)') tekel
!
               ta(1) = 0.0_real64
               ta(2:mtp) = -tb(2:mtp)*bj(2:mtp,3)*radial_grid_derivative(2:mtp)
               if (ldbpr(15)) write (angular_debug_unit, '(/,100(1x,7(i3,2x,1p,d11.3,2x),/),/)') (ii,ta(ii),ii=1,radial_point_count)
               call quad (peres)
               if (ldbpr(14)) then
                  write (angular_debug_unit, '(" I + (L+1) integral  = ",1p,d24.12)') peres
                  write (angular_debug_unit, '(" Factor multiplying this ",1p,d20.12)') cip
               end if
               tamar = cip*peres
               if (ldbpr(14)) write (angular_debug_unit, '(" Resultant contribution  ",1p,d20.12,/)') tamar
!
               ta(1) = 0.0_real64
               ta(2:mtp) = -tb(2:mtp)*bj(2:mtp,1)*radial_grid_derivative(2:mtp)
               if (ldbpr(15)) write (angular_debug_unit, '(/,100(1x,7(i3,2x,1p,d11.3,2x),/),/)') (ii,ta(ii),ii=1,radial_point_count)
               call quad (enoch)
               if (ldbpr(14)) then
                  write (angular_debug_unit, '(" I + (L-1) integral  = ",1p,d24.12)') enoch
                  write (angular_debug_unit, '(" Factor multiplying this ",1p,d20.12)') cip
               end if
               seth = cip*enoch
               if (ldbpr(14)) write (angular_debug_unit, '(" Resultant contribution  ",1p,d20.12,/)') seth
!
               ta(1) = 0.0_real64
               ta(2:mtp) = bj(2:mtp,1)*tc(2:mtp)*radial_grid_derivative(2:mtp)
               if (ldbpr(15)) write (angular_debug_unit, '(/,100(1x,7(i3,2x,1p,d11.3,2x),/),/)') (ii,ta(ii),ii=1,radial_point_count)
               call quad (shem)
               if (ldbpr(14)) then
                  write (angular_debug_unit, '(" I - (L-1) integral  = ",1p,d24.12)') shem
                  write (angular_debug_unit, '(" Factor multiplying this ",1p,d20.12)') fl
               end if
               darius = fl*shem
               if (ldbpr(14)) write (angular_debug_unit, '(" Resultant contribution  ",1p,d20.12,/)') darius
!
               ta(1) = 0.0_real64
               ta(2:mtp) = -bj(2:mtp,3)*tc(2:mtp)*radial_grid_derivative(2:mtp)
               if (ldbpr(15)) write (angular_debug_unit, '(/,100(1x,7(i3,2x,1p,d11.3,2x),/),/)') (ii,ta(ii),ii=1,radial_point_count)
               call quad (cyrus)
               if (ldbpr(14)) write (angular_debug_unit, '(" I - (L+1) integral  = ",1p,d24.12)') cyrus
               daniel = flp
               if (ldbpr(14)) write (angular_debug_unit, '(" Factor multiplying this ",1p,d20.12)') daniel
               boaz = daniel*cyrus
               if (ldbpr(14)) write (angular_debug_unit, '(" Resultant contribution  ",1p,d20.12,/)') boaz
!
               esau = tamar + seth + darius + boaz
               if (ldbpr(14)) then
                  write (angular_debug_unit, '(" These last four add up to ",1p,d20.12,&
                  &"    compared to the first which is ",d20.12)') esau, tekel
               end if
               aaron = esau + tekel
               gimel = -value
               if (ldbpr(14)) then
                  write (angular_debug_unit, '(" These all add up to ",1p,d20.12," which should equal ",d20.12)') aaron, gimel
                  write (angular_debug_unit, '(/,18x,"Integral",11x,"Coefficient",8x,"Contribution")')
                  write (angular_debug_unit, '(" I + (L+1)",1p,3d20.10)') peres, cip, tamar
                  write (angular_debug_unit, '(" I + (L-1)",1p,3d20.10)') enoch, cip, seth
                  write (angular_debug_unit, '(" I - (L-1)",1p,3d20.10)') shem, fl, darius
                  write (angular_debug_unit, '(" I - (L+1)",1p,3d20.10)') cyrus, daniel, boaz
                  write (angular_debug_unit, '(54x,16("-"),/,50x,1p,d20.10)') esau
                  write (angular_debug_unit, '(" J (L)",4x,1p,3d20.10)') tarom, esther, tekel
                  write (angular_debug_unit, '(54x,16("-"),/," Correct figure",1p,d20.10,15x,d20.10)') gimel, aaron
               end if
!
            end if
!
         end if
!
         if (ldbpr(12)) then
            write (angular_debug_unit, '(/,1x,"Orbital pair (",i2,",",i2,") integrals:","Coul. gauge ",1p,d13.3,3x,&
            &"Gauge contribution",d13.3,3x,"Bab.gauge",d13.3)') first_orbital, second_orbital, hcoul, hgauge, hbab
         end if
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
               ta(ii) = (pfii(ii,nnii(first_orbital))*qfff(ii,nnff(second_orbital))+qfii(ii,nnii(first_orbital))*&
                  pfff(ii,nnff(second_orbital)))*bj(ii,2)*radial_grid_derivative(ii)*speed_of_light_au/vacuum_speed_of_light_au
!           TA(II) =  (PFII(II,NNII(I))*QFFF(II,NNFF(J)) +
!    :      QFII(II,NNII(I))*PFFF(II,NNFF(J)))*BJ(II,2)*RP(II)
!zou
            end do
            call quad (value)
            if (ldbpr(14)) write (angular_debug_unit, '(" I + (L)  integral",1p,d20.10)') value
            hmag = -value*(fl + flp)*(dfki + dfkj)*form/sqrt(fl*flp)
!
         end if
!
         if (ldbpr(12)) write (angular_debug_unit, '(/,1x,"Magnetic single particle matrix element ",1p,d13.3)') hmag
!
      end if
!
      return
!
      return
!
      end subroutine spme
