!***********************************************************************
!                                                                      *
      subroutine qed(jstate, slfint, ucf)
!                                                                      *
!   This  routine estimates corrections to the  energy levels due to   *
!   self-energy.                                                       *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, IQ, QUAD, RALLOC, SCREEN.      *
!               [RCI92]: FZALF, HOVLAP.                                *
!                                                                      *
!                                           Last update: 30 Oct 1992   *
!   Modified by Xinghong He                 Last update: 24 Jun 1997   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: max_orbitals, max_radial_points
      use atomic_calculation_state, only: max_radial_points, nuclear_charge, pi, speed_of_light_au
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      !USE grid_C
      !USE horb_C
      use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
      use orbital_csf_basis, only: orbital_count, csf_count, orbital_principal_n, orbital_kappa
      use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      !USE tatb_C
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use csf_quantum_queries, only: iq
      use ratden_i, only: ratden
      use fzalf_i, only: fzalf
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: jstate
      real(real64), intent(out) :: slfint(max_orbitals)
      real(real64), intent(out) :: ucf(:)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: maxiter, j, i, ii, npj, kappa, mfj
      real(real64) :: ptemp(max_radial_points), qtemp(max_radial_points)
      real(real64) :: ucfj, zeff, ratio, valu
      character :: npchar, nakchar*2
!-----------------------------------------------
!
! Pre-set tolerable number for iteration in finding effective
! nuclear charge.
!
      maxiter = 20
!
! Modified so that UCFJ describes the current eigenstate
!
      do j = 1, orbital_count
         ucfj = 0.0_real64
!         DO 3 I = 1,NVEC
         i = jstate
         do ii = 1, csf_count
            ucfj = ucfj + dble(iq(j,ii))*csf_mixing_coefficients(ii + (i - 1)*csf_count)**2
         end do
!    3    CONTINUE
!     print *, ucfj,'ucf'
         ucf(j) = ucfj
! zou    UCF(J) = UCFJ/DBLE (NCF)
      end do
!
      do j = 1, orbital_count
!
         npj = orbital_principal_n(j)
!
         if (npj <= 8) then
!
!   Only orbitals with principal quantum number 8 or less can
!   be treated by this section of code
!
            kappa = orbital_kappa(j)
!
!   Begin by transferring the function to a temporary array
!
            mfj = orbital_last_grid_index(j)
!
            ptemp(1) = 0.0_real64
            qtemp(1) = 0.0_real64
            do i = 2, mfj
               ptemp(i) = radial_large_components(i,j)
               qtemp(i) = radial_small_components(i,j)
            end do
!
            zeff = nuclear_charge
            ratio = ratden(ptemp,qtemp,mfj,npj,kappa,zeff)
            valu = ratio*fzalf(npj,kappa,zeff)/dble(npj**3)
!
            slfint(j) = valu*zeff**4/(pi*speed_of_light_au**3)
!        print *, 'No. orb.=',j,' Zeff = ',zeff
!    &   , 'Scale= ',ratio
!    &   , 'S.E. = ',slfint(j)*2*13.6058,slfint(j)/ratio*2*13.6058
!
         else
!
!   The self-energy for orbitals with principal quantum number
!   greater than 8 is set to zero
!
            slfint(j) = 0.0_real64
!
         end if
!
      end do
!
!   Deallocate storage for the `generalised occupation numbers'
!
!
      return
      end subroutine qed


!***********************************************************************
!                                                                      *
      real(real64) function ratden (p, q, mtpo, np, kappa, z)
!                                                                      *
!   This subprogram computes the overlap of the orbital tabulated in   *
!   the arrays  P  and  Q  with maximum tabulation point  MTPO  with   *
!   a hydrogenic orbital with parameters  NP  KAPPA  Z .               *
!                                                                      *
!   Call(s) to: [LIB92]: DCBSRW, QUAD.                                 *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!
!...Switches:
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: max_radial_points
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use hydrogenic_radial_components, only: ph, qh
      use radial_quadrature_workspace, only: mtp, ta
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dirac_coulomb_radial_orbital, only: dcbsrw
      use radial_quadrature, only: quad
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: mtpo
      integer  :: np
      integer  :: kappa
      real(real64)  :: z
      real(real64), intent(in) :: p(max_radial_points)
      real(real64), intent(in) :: q(max_radial_points)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mtph, i, k
      real(real64) :: eh, pzh, result, result1
!-----------------------------------------------
!
!
!
!   Set up the hydrogenic orbital
!
      call dcbsrw (np, kappa, z, eh, pzh, ph, qh, mtph)
!
!   Compute the overlap
!
      mtp = min(mtph,mtpo)
      do i = 2, mtp
         if (radial_grid_derivative(i) > 0.0219_real64) cycle
         k = i
      end do
      mtp = k
      ta(1) = 0.0_real64
      ta(2:mtp) = (p(2:mtp)*p(2:mtp)+q(2:mtp)*q(2:mtp))*radial_grid_derivative(2:mtp)
!         TA(I) = (P(I)*PH(I)+Q(I)*QH(I))*RP(I)
      call quad (result)
      ta(2:mtp) = (ph(2:mtp)*ph(2:mtp)+qh(2:mtp)*qh(2:mtp))*radial_grid_derivative(2:mtp)
      call quad (result1)
!
      ratden = result/result1
!
      return
      end function ratden
