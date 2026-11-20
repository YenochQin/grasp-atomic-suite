!***********************************************************************
!                                                                      *
      subroutine qed_slfen(slfint)
!                                                                      *
!   This  routine estimates the F(Z\alpha) function of self energy for *
!   each orbital.                                                      *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, IQ, QUAD, RALLOC, SHIELD.      *
!               [RCI92]: FZALF, HOVLAP.                                *
!                                                                      *
!   Modified from subroutine QED by Yu Zou, Last update: 13 Mar 2000   *
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
      use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
      use orbital_csf_basis, only: orbital_count, orbital_principal_n, orbital_kappa
      use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
      use qed_orbital_cutoffs, only: nqedcut, nqedmax
      use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_integer, mpi_rank
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ratden_i, only: ratden
      use fzalf_i, only: fzalf
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(real64), intent(out) :: slfint(max_orbitals)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: maxiter, orbital_index, npj, kappa, last_radial_index, i, npjmax
!      REAL(DOUBLE), DIMENSION(1) :: UCF
      real(real64) :: ptemp(max_radial_points), qtemp(max_radial_points)
      real(real64) :: zeff, ratio, valu
      character :: npchar, nakchar*2
!-----------------------------------------------
!
! Pre-set tolerable number for iteration in finding effective
! nuclear charge.
!
      maxiter = 20
!
      if (mpi_rank == 0) then
         if (nqedcut==1) then
            npjmax = nqedmax
         else
            npjmax = 8
         end if
      end if
      call mpi_bcast (npjmax,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
!
!
      do orbital_index = 1, orbital_count
!
         npj = orbital_principal_n(orbital_index)
!
         if (npj <= npjmax) then
!
!   Only orbitals with principal quantum number 8 or less can
!   be treated by this section of code
!
            kappa = orbital_kappa(orbital_index)
!
!   Begin by transferring the function to a temporary array
!
            last_radial_index = orbital_last_grid_index(orbital_index)
!
            ptemp(1) = 0.0_real64
            qtemp(1) = 0.0_real64
            do i = 2, last_radial_index
               ptemp(i) = radial_large_components(i,orbital_index)
               qtemp(i) = radial_small_components(i,orbital_index)
            end do
            zeff = nuclear_charge
            ratio = ratden(ptemp,qtemp,last_radial_index,npj,kappa,zeff)
            valu = ratio*fzalf(npj,kappa,zeff)/dble(npj**3)
            slfint(orbital_index) = valu*zeff**4/(pi*speed_of_light_au**3)
!
         else
!
!   The self-energy for orbitals with principal quantum number
!   greater than 8 is set to zero
!
            slfint(orbital_index) = 0.0_real64
!
         end if
!
      end do
!
!   Deallocate storage for the `generalised occupation numbers'
!
!
      return
      end subroutine qed_slfen
