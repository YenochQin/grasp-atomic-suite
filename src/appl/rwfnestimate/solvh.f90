

!***********************************************************************
!                                                                      *
      subroutine solvh(iorb, fail)
!                                                                      *
!   This routine solves the homogeneous Dirac radial equation.         *
!                                                                      *
!   Arguments:  IORB : (Input) Index of orbital                        *
!               FAIL : (Output) .TRUE. if solution not obtained        *
!                                                                      *
!   The direct potential is assumed tabulated in the COMMON array YP   *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!               [RSCF92]: COUNT, SBSTEP, SETPOT, START, TAIL.          *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use def_c
      use grid_c
      use orb_c
      use pote_c, only: yp
      use tatb_c, only: ta, mtp
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setpot_i
      use start_i
      use sbstep_i
      use tail_i
      use count_i
      use quad_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: iorb
      logical , intent(out) :: fail
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mxk, npiorb, nkiorb, nakabs, nliorb, nnp, nrel, kount, jp, &
         itype, nstrt, i, npc
      real(double) :: epslon, csq, alpha, fkabs, fkap2, zalpha, gamma, ebym, &
         emin, emax, dele, eest, q0, pfjpo, qfjpo, pfjpi, qfjpi, ratio, sgn, &
         dnorm, dmsmch, dnfac
!-----------------------------------------------
!
      data mxk/ 75/
!
!   A solution is deemed continuous when the relative mismatch in the
!   small component is within EPSLON
!
      epslon = accy*0.1e00
!
!   Establish the number of nodes in the large component
!
      npiorb = orbital_principal_n(iorb)
      nkiorb = orbital_kappa(iorb)
      nakabs = abs(nkiorb)
      if (nkiorb < 0) then
         nliorb = nakabs - 1
      else
         nliorb = nakabs
      endif
      nnp = npiorb - nliorb - 1
!
!   Establish the bounds on, and an estimate of, the eigenvalue
!
      csq = speed_of_light_au*speed_of_light_au
      alpha = 1.0d00/speed_of_light_au
      nrel = npiorb - nakabs
      fkabs = dble(nakabs)
      fkap2 = fkabs*fkabs
!
      zalpha = yp(radial_point_count)*alpha
      if (zalpha < fkabs) then
         gamma = sqrt(fkap2 - zalpha*zalpha)
         ebym = 1.0d00/sqrt(1.0d00 + (zalpha/(gamma + nrel + 0.5d00))**2)
         emin = (1.0d00 - ebym)*csq
      else
         emin = 0.25d00*csq/dble(npiorb*npiorb)
      endif
!
      zalpha = nuclear_charge*alpha
!
      if (zalpha < fkabs) then
         gamma = sqrt(fkap2 - zalpha*zalpha)
         ebym = 1.0d00/sqrt(1.0d00 + (zalpha/(gamma + nrel))**2)
         orbital_energy_parameters(iorb) = (1.0d00 - ebym)*csq
      else
         orbital_energy_parameters(iorb) = csq
      endif
!
      if (zalpha < fkabs) then
         gamma = sqrt(fkap2 - zalpha*zalpha)
         ebym = 1.0d00/sqrt(1.0d00 + (zalpha/(gamma + nrel - 0.5d00))**2)
         emax = (1.0d00 - ebym)*csq
      else
         emax = csq + csq
      endif
!
      dele = 0.0d00
!
!   Initialize
!
      fail = .false.
      kount = -1
!
!   Iteration loop begins here
!
    1 continue
      kount = kount + 1
      if (kount > mxk) then
         fail = .true.
         return
      endif
!
!   Generate estimate of eigenvalue for this iteration
!
      eest = orbital_energy_parameters(iorb) + dele
      if (eest>emin .and. eest<emax) then
         orbital_energy_parameters(iorb) = eest
      else
         orbital_energy_parameters(iorb) = 0.5d00*(emin + emax)
      endif
!
!   E is the positive binding energy. An invalid trial energy cannot
!   produce the decaying bound-state tail used by TAIL.
!
      if (.not.ieee_is_finite(orbital_energy_parameters(iorb)) .or. &
          orbital_energy_parameters(iorb)<=0.0d00 .or. orbital_energy_parameters(iorb)>=2.0d00*csq) then
         fail = .true.
         return
      endif
!
!   Set up arrays TF and TG; find join point
!
      call setpot (iorb, jp)
!
!   Initialize outward integration
!
      itype = 1
      call start (iorb, itype, orbital_origin_parameters(iorb), radial_large_components(:,iorb), q0, &
         & radial_small_components(:,iorb))
!
!   Continue outward integration to classical turning point
!
      nstrt = 6
      call sbstep (iorb, nstrt, jp, radial_large_components(:,iorb), radial_small_components(:,iorb))
      pfjpo = radial_large_components(jp,iorb)
      qfjpo = radial_small_components(jp,iorb)
!
!   Initialize inward integration
!
      call tail (iorb, radial_large_components(:,iorb), radial_small_components(:,iorb), jp, mtp)
!
!   Continue inward integration to classical turning point
!
      nstrt = mtp
      call sbstep (iorb, nstrt, jp, radial_large_components(:,iorb), radial_small_components(:,iorb))
      pfjpi = radial_large_components(jp,iorb)
      qfjpi = radial_small_components(jp,iorb)
!
!   Make large component continuous, determine mismatch in small
!   component
!
      ratio = pfjpo/pfjpi
      radial_large_components(jp:radial_point_count,iorb) = radial_large_components(jp:radial_point_count,iorb)*ratio
      radial_small_components(jp:radial_point_count,iorb) = radial_small_components(jp:radial_point_count,iorb)*ratio
!
!   Count nodes
!
      call count (radial_large_components(:,iorb), mtp, npc, sgn)
!
!   Correct the energy estimate if the number of nodes is wrong
!
      if (npc > nnp) then
         emin = orbital_energy_parameters(iorb)
         dele = 0.5d00*(emax - emin)
         go to 1
      else if (npc < nnp) then
         emax = orbital_energy_parameters(iorb)
         dele = -0.5d00*(emax - emin)
         go to 1
      endif
!
!   Correct number of nodes
!
!   Compute 'norm' of solution
!
      ta(1) = 0.0d00
      ta(2:radial_point_count) = (radial_large_components(2:radial_point_count,iorb)**2+ &
         & radial_small_components(2:radial_point_count,iorb)**2)*radial_grid_derivative(2:radial_point_count)
      mtp = radial_point_count
      call quad (dnorm)
      if (.not.ieee_is_finite(dnorm) .or. dnorm<=0.0d00) then
         fail = .true.
         return
      endif
!
!   Determine correction to eigenvalue from magic formula
!   correct slope at origin
!
      qfjpi = qfjpi*ratio
      dmsmch = qfjpi - qfjpo
      if (.not.ieee_is_finite(dmsmch) .or. &
          .not.ieee_is_finite(qfjpo)) then
         fail = .true.
         return
      endif
      if (abs(dmsmch) > epslon*abs(qfjpo)) then
         dele = speed_of_light_au*radial_large_components(jp,iorb)*dmsmch/dnorm
         if (.not.ieee_is_finite(dele)) then
            fail = .true.
            return
         endif
!
!   A large correction (notably for the first h/i orbitals) can make
!   the old multiplicative upper bound negative. Keep the current
!   trial inside the updated interval; the next iteration bisects it
!   if the proposed correction falls outside the bounds.
!
         if (dele < 0.0d00) then
            emax = max(orbital_energy_parameters(iorb), emax*(1.0d00 - 0.2e00*abs(dele/orbital_energy_parameters(iorb))))
         else
            emin = min(orbital_energy_parameters(iorb), emin*(1.0d00 + 0.2e00*abs(dele/orbital_energy_parameters(iorb))))
         endif
         orbital_origin_parameters(iorb) = orbital_origin_parameters(iorb)/sqrt(dnorm)
         go to 1
      endif
!
!   Normalize
!
      dnfac = 1.0d00/sqrt(dnorm)
      orbital_origin_parameters(iorb) = orbital_origin_parameters(iorb)*dnfac
      radial_large_components(:radial_point_count,iorb) = radial_large_components(:radial_point_count,iorb)*dnfac
      radial_small_components(:radial_point_count,iorb) = radial_small_components(:radial_point_count,iorb)*dnfac
!
!   Find maximum tabulation point
!
      i = radial_point_count + 1
    5 continue
      i = i - 1
      if (abs(radial_large_components(i,iorb)) < epslon) then
         radial_large_components(i,iorb) = 0.0d00
         radial_small_components(i,iorb) = 0.0d00
         go to 5
      else
         orbital_last_grid_index(iorb) = i
      endif
!
      return
      end subroutine solvh
