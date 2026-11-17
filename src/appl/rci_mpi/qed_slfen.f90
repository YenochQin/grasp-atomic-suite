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
      use parameter_def,   only: max_orbitals, max_radial_points
      use def_c
      use eigv_c
      use npar_c, only: double, nuclear_model_parameter_count, nuclear_model_parameters
      use orb_c, only: orbital_count, orbital_principal_n, orbital_kappa
      use prnt_c, only: byte, double, ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use wave_c
      use qedcut_c, only: nqedcut, nqedmax
      use mpi_c
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
      integer :: maxiter, j, npj, kappa, mfj, i, npjmax
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
      do j = 1, orbital_count
!
         npj = orbital_principal_n(j)
!
         if (npj <= npjmax) then
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
            zeff = nuclear_charge
            ratio = ratden(ptemp,qtemp,mfj,npj,kappa,zeff)
            valu = ratio*fzalf(npj,kappa,zeff)/dble(npj**3)
            slfint(j) = valu*zeff**4/(pi*speed_of_light_au**3)
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
      end subroutine qed_slfen
