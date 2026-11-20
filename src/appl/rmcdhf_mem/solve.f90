!***********************************************************************
!                                                                      *
      subroutine solve(j, fail, inv, jp, nnp)
!                                                                      *
!   This subroutine  performs  step  2 in Algorithm 5.2 and 5.3 of C   *
!   Froese Fischer, Comput Phys Rep 3 (1986) 295. Some minor changes   *
!   have been made.                                                    *
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!      J     : (Input) The serial number of the orbital                *
!      JP    : (Output) The join point                                 *
!      FAIL  : (Output) If .true., the iterations did not yield an     *
!              acceptable solution (methods 1 and 2)                   *
!                                                                      *
!   Call(s) to: [RSCF92]: COUNT, ESTIM, IN, NEWE, OUT, cofpot,         *
!                         SETXUV, SETXV, SETXZ, START                  *
!               [LIB92]: DCBSRW, QUAD, SETPOT.                         *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 26 Sep 1993   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,    only: max_radial_points
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use atomic_calculation_state, only: speed_of_light_au, nsolv
      use radial_grid_state, only: radial_grid_derivative
      use radial_dirac_workspace, only: p, q, p0, q0, mtp0
      use orbital_inversion_flags, only: max_orbitals, noinvt
      use expected_orbital_nodes, only: max_orbitals, nnodep
      use orbital_csf_basis, only: max_orbitals, orbital_energy_parameters, orbital_kappa, orbital_l, orbital_principal_n, &
          orbital_symmetry_labels
      use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
          nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
      use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use estim_i, only: estim
      use eigen_i, only: eigen
      use dirac_coulomb_radial_orbital, only: dcbsrw
      use radial_dirac_coefficients, only: setpot
      use setxz_i, only: setxz
      use radial_origin_solution, only: start
      use out_i, only: out
      use in_i, only: in
      use setxuv_i, only: setxuv
      use radial_quadrature, only: quad
      use setxv_i, only: setxv
      use prwf_i, only: prwf
      use radial_node_count, only: count
      use newe_i, only: newe
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: j
      integer, intent(inout) :: inv
      integer, intent(out) :: jp
      integer, intent(out) :: nnp
      logical  :: fail
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: kount, nakj, nloops, icase, mtph, mtpi, i, mtpc, mtpv, mtpvc, &
         loc1, loc2, loc, mx, nprime
      real(real64) :: ph(max_radial_points), qh(max_radial_points), pv(max_radial_points), qv(max_radial_points)
      real(real64) :: twocsq, elast, tenej, deleps, sgn, qjpoh, qjpih, qjpoi, &
         qjpii, dnorm, alfa, p0h, p0v, q0v, qjpov, qjpiv, crnorm, dvnorm, aa, &
         bb, cc, discr, qq, root1, root2, pmx, apm, pati, abpi, ratio, test1, &
         test2, dele
      logical :: check
!-----------------------------------------------
!
!   Initialization
!
      check = .not.noinvt(j)
      fail = .false.
      kount = 0
      nakj = orbital_kappa(j)
!     TWOCSQ = 2.0D 00*137.036**2
      twocsq = 2.0_real64*speed_of_light_au*speed_of_light_au
!
!   Debug header
!
      if (ldbpr(22)) write (99, '(/,/," Debug printout active; orbital: ",1i2,1a2)') &
         orbital_principal_n(j), orbital_symmetry_labels(j)
!
      nloops = max(nsolv,3*orbital_principal_n(j))
!
      call estim (j)
      elast = orbital_energy_parameters(j)
      orbital_energy_parameters(j) = eigen(j)
!
!   Checks on lower bounds
!
      if (orbital_energy_parameters(j) < epsmin) then
         if (orbital_energy_parameters(j) > 0.0_real64) then
            epsmin = orbital_energy_parameters(j)
         else
            write (*, '(" E(",1i2,1a2,") = ",1p,d11.4,"; adjusted to EPSMIN")') &
               orbital_principal_n(j), orbital_symmetry_labels(j), orbital_energy_parameters(j)
            orbital_energy_parameters(j) = epsmin
            if (abs(emin - elast)<=1.0e-6_real64 .and. method(j)<=2) then
               call dcbsrw (orbital_principal_n(j), nakj, zinf, orbital_energy_parameters(j), p0, p, q, mtp0)
               write (*, '(" Returned hydrogenic function for ",1i2,1a2," with"," effective charge ",f7.3)') &
                  orbital_principal_n(j), orbital_symmetry_labels(j), zinf
               return
            end if
         end if
      end if
!
!   Check on upper bound
!
      if (method(j) <= 2) then
         if (orbital_energy_parameters(j) > epsmax) then
            tenej = 10.0_real64*orbital_energy_parameters(j)
            epsmax = min(tenej,twocsq)
            emax = epsmax
            if (orbital_energy_parameters(j) > twocsq) then
               write (*, '(" E(",1i2,1a2,") = ",1p,d11.4,"; adjusted to TWOCSQ")') &
                  orbital_principal_n(j), orbital_symmetry_labels(j), orbital_energy_parameters(j)
               orbital_energy_parameters(j) = twocsq
            end if
         end if
      end if
!
!   Iteration loop begins here
!
    1 continue
      kount = kount + 1
!
      p0 = orbital_origin_parameters(j)
!
      if (kount > 1) then
!
!   Check that bounds are ordered correctly
!
         if (epsmax <= epsmin) write (*, '(" Warning: difficulty with node-counting procedure"/,&
            &" lower bound on energy (",1p,d11.4,") exceeds upper"," bound (",1d11.4,"; E(",1i2,1a2,") = ",1d11.4)') &
            epsmin, epsmax, orbital_principal_n(j), orbital_symmetry_labels(j), &
            & orbital_energy_parameters(j&
            )
!
         if (kount>nloops .or. epsmax-epsmin<1.0_real64/dble(orbital_principal_n(j))**3) then
            write (*, '(" Method ",1i1," unable to solve for ",1i2,1a2," orbital")') &
               method(j), orbital_principal_n(j), orbital_symmetry_labels(j)
            write (*, '(" Iteration number: ",1i2,", limit: ",1i2,/," Present estimate of P0; ",1d21.14,/,&
               &" Present estimate of E(J): ",1d21.14,", DELEPS: ",1d21.14,/," Lower bound on energy: ",1d21.14,&
               &", upper bound: ",1d21.14,/," Join point: ",1i4,", Maximum tabulation point:",1i4,/,&
               &" Number of nodes counted: ",1i2,", Correct number: ",1i2,/," Sign of P at first oscillation: ",&
               &f3.0)') kount - 1, nloops, p0, orbital_energy_parameters(j), deleps, epsmin, epsmax&
               , jp, mtp, nnp, nnodep(j), sgn
            fail = .true.
            return
         end if
      end if
!
!   Set up arrays TF and TG; find join point
!
      call setpot (j, jp)
!
!   Set right-hand side to zero to form homogeneous equations;
!   integrate homogeneous equations outwards and inwards; store
!   small component at join point each time
!
      call setxz (j)
      icase = 1
      call start (j, icase, p0, ph, q0, qh)
      call out (j, jp, ph, qh)
      qjpoh = qh(jp)
      call in (j, jp, ph, qh, mtph)
      qjpih = qh(jp)
!
!   Set up right-hand side for inhomogeneous equations; integrate
!   inhomogeneous equations outwards and inwards; store small
!   component at join point each time
!
      call setxuv (j)
      icase = 2
      call start (j, icase, p0, p, q0, q)
      call out (j, jp, p, q)
      qjpoi = q(jp)
      call in (j, jp, p, q, mtpi)
      qjpii = q(jp)
!
!   Determine energy adjustment for methods 1 and 2
!
      if (method(j) <= 2) then
         ta(1) = 0.0_real64
         ta(2:mtpi) = (p(2:mtpi)**2+q(2:mtpi)**2)*radial_grid_derivative(2:mtpi)
         mtp = mtpi
         call quad (dnorm)
!        DELEPS = 137.036*P(JP)*(QJPII-QJPOI)/DNORM
         deleps = speed_of_light_au*p(jp)*(qjpii - qjpoi)/dnorm
      end if
!
!   Generate the continuous solution
!
      mtpc = max(mtph,mtpi)
      alfa = -(qjpii - qjpoi)/(qjpih - qjpoh)
      p0h = p0
      p0 = p0*(1.0_real64 + alfa)
      p(:mtpc) = p(:mtpc) + alfa*ph(:mtpc)
      q(:mtpc) = q(:mtpc) + alfa*qh(:mtpc)
!
      if (method(j)==2 .or. method(j)==4) then
!
!   Set up right-hand side for variational equations; integrate
!   variational equations outwards and inwards; store small
!   component at join point each time
!
         p0v = 0.0_real64
         call setxv (j)
         icase = 3
         call start (j, icase, p0v, pv, q0v, qv)
         call out (j, jp, pv, qv)
         qjpov = qv(jp)
         call in (j, jp, pv, qv, mtpv)
         qjpiv = qv(jp)
!
!   Generate continuous solutions
!
         mtpvc = max(mtpc,mtpv)
         alfa = -(qjpiv - qjpov)/(qjpih - qjpoh)
         pv(:mtpvc) = pv(:mtpvc) + alfa*ph(:mtpvc)
         qv(:mtpvc) = qv(:mtpvc) + alfa*qh(:mtpvc)
!
         ta(1) = 0.0_real64
         ta(2:mtpc) = radial_grid_derivative(2:mtpc)*(p(2:mtpc)**2+q(2:mtpc)**2)
         mtp = mtpc
         call quad (dnorm)
!
         mtp = min(mtpc,mtpvc)
         ta(1) = 0.0_real64
         ta(2:mtp) = radial_grid_derivative(2:mtp)*(p(2:mtp)*pv(2:mtp)+q(2:mtp)*qv(2:mtp))
         call quad (crnorm)
!
         ta(1) = 0.0_real64
         ta(2:mtpvc) = radial_grid_derivative(2:mtpvc)*(pv(2:mtpvc)**2+qv(2:mtpvc)**2)
         mtp = mtpvc
         call quad (dvnorm)
!
!   Determine deleps required to normalize new solution to
!   first order: modified form of solution to a quadratic
!   equation (see Press et al.)
!
         aa = dvnorm
         bb = crnorm + crnorm
         cc = dnorm - 1.0_real64
         discr = bb*bb - 4.0_real64*aa*cc
         if (discr > 0.0_real64) then
            qq = -0.5_real64*(bb + sign(1.0_real64,bb)*sqrt(discr))
            root1 = cc/qq
            root2 = qq/aa
            pmx = 0.0_real64
            apm = 0.0_real64
            do i = 2, jp
               pati = p(i)
               if (pati > pmx) then
                  pmx = pati
                  loc1 = i
               end if
               abpi = abs(pati)
               if (abpi <= apm) cycle
               apm = abpi
               loc2 = i
            end do
            if (pmx /= 0.0_real64) then
               ratio = apm/abs(pmx)
               if (ratio < 10.0_real64) then
                  loc = loc1
               else
                  loc = loc2
               end if
            else
               loc = loc2
            end if
            test1 = p(loc) + root1*pv(loc)
            test2 = p(loc) + root2*pv(loc)
            if (test1>0.0_real64 .and. test2<0.0_real64) then
               dele = root1
            else if (test1<0.0_real64 .and. test2>0.0_real64) then
               dele = root2
            else if (test1>0.0_real64 .and. test2>0.0_real64) then
               if (test1 < test2) then
                  dele = root1
               else
                  dele = root2
               end if
            else if (test1<0.0_real64 .and. test2<0.0_real64) then
               if (test1 > test2) then
                  dele = root1
               else
                  dele = root2
               end if
            end if
         else
            dele = -bb/(aa + aa)
         end if
!
!   Generate new solution
!
         mtp0 = max(mtpc,mtpvc)
         p0 = p0 + dele*alfa*p0h
         p(2:mtp0) = p(2:mtp0) + dele*pv(2:mtp0)
         q(2:mtp0) = q(2:mtp0) + dele*qv(2:mtp0)
      else
         mtp0 = mtpc
      end if
!
!   Debug printout
!
      if (ldbpr(23)) call prwf (j)
!
!   Count nodes in large component function; determine sign at first
!   oscillation, effective quantum number; note that node counting
!   is never enforced on the small component
!
      call count (p, mtp0, nnp, sgn)
!
!   DEBUG PRINTOUT
!
      if (ldbpr(22)) write (99, '(" Iteration number: ",1i2,", limit: ",1i2,/," Present estimate of P0; ",1d21.14,/&
         &," Present estimate of E(J): ",1d21.14,", DELEPS: ",1d21.14,/," Lower bound on energy: ",1d21.14,&
         &", upper bound: ",1d21.14,/," Join point: ",1i4,", Maximum tabulation point:",1i4,/,&
         &" Number of nodes counted: ",1i2,", Correct number: ",1i2,/," Sign of P at first oscillation: ",f3.0)') &
         kount, nloops, p0, orbital_energy_parameters(j), deleps, epsmin, &
         epsmax, jp, mtp, nnp, nnodep(j), sgn
!
!   Proceed according to method
!
      if (method(j) > 2) then
         if (check .and. sgn<0.0_real64) then
            inv = 1
            p0 = -p0
            p(2:mtp0) = -p(2:mtp0)
            q(2:mtp0) = -q(2:mtp0)
         end if
      else
         mx = nnp - nnodep(j)
         nprime = nnp + orbital_l(j) + 1
         call newe (j, sgn, nprime, mx, deleps, fail, inv)
         if (fail) goto 1
      end if
!
!   Solution found
!
      return
!
      return
!
      end subroutine solve
