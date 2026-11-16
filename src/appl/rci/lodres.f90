!***********************************************************************
!                                                                      *
      subroutine lodres
!                                                                      *
!   Loads the data from the  .res  file. A number of checks are made   *
!   to ensure correctness and consistency.                             *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, GETYN, SETQIC.                 *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Oct 1992   *
!   Block version by Xinghong He          Last revision:  1 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: max_radial_points
      use memory_man
      use decide_c
      use def_c, only: nelecr, electron_count, nuclear_charge, emn, speed_of_light_au, accy
      use grid_c
      use npar_c
      use npot_c, only: negative_radius_times_potential, first_grid_index_beyond_fermi_c
      use orb_c
      use wave_c
      use wfac_c
      use where_c
      use hblock_c
      use iccu_c
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      use setqic_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ncfres, nwres, nblockres, i, np10, j
      logical :: yes
!-----------------------------------------------
!
!     POINTER (pncfblk, ncfblk(0:*))
!
!     POINTER (piccutblk, iccutblk(1))
!
!-----------------------------------------------------------------------
      write (6, *) 'Calling lodres ...'
!
!   Read the basic parameters of the electron cloud; check these
!   against those deduced from the  .csl file
!
      read (imcdf) nelecr, ncfres, nwres, nblockres

      if (nelecr/=electron_count .or. ncfres/=csf_count .or. nwres/=orbital_count .or. nblockres/=symmetry_block_count&
         ) stop 'lodres: NELEC/NCF/NW does not match'
!
!   Read the nuclear parameters
!
      read (imcdf) nuclear_charge, emn
      read (imcdf) nuclear_model_parameter_count, (nuclear_model_parameters(i),i=1,nuclear_model_parameter_count)
      read (imcdf) radial_point_count, (negative_radius_times_potential(i),i=1,radial_point_count), first_grid_index_beyond_fermi_c

      if (radial_point_count > max_radial_points) stop 'lodres: N greater than NNNP'
!
!   Read the physical effects specifications
!   iccutblk() is now an array of length nblock.
!
      read (imcdf) speed_of_light_au, lfordr, (iccutblk(i),i=1,symmetry_block_count), ltrans, wfact, lvp, &
         lnms, lsms
!
!   Read the remaining parameters controlling the radial grid and the
!   grid arrays
!
      np10 = radial_point_count + 10
      read (imcdf) radial_grid_scale, transformed_grid_step, asymptotic_grid_step, (radial_grid(i),i=1,np10), &
         & (radial_grid_derivative(i),i=1,np10), (radial_derivative_over_radius(i),i=1,&
         np10)
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      accy = transformed_grid_step**6
!
!   Set up the coefficients for the numerical procedures
!
      call setqic
!
!   Allocate storage for the radial wavefunction arrays
!
      call alloc (radial_large_components, max_radial_points,orbital_count, 'PF', 'LODMIX')
      call alloc (radial_small_components, max_radial_points,orbital_count, 'QF', 'LODMIX')
!
!   Read the orbital wavefunctions and the associated arrays
!
      do j = 1, orbital_count
         read (imcdf) orbital_energy_parameters(j), orbital_origin_exponents(j), orbital_origin_parameters(j), &
            & orbital_last_grid_index(j)
         read (imcdf) (radial_large_components(i,j),i=1,orbital_last_grid_index(j)), (radial_small_components(i,j),i=1, &
            & orbital_last_grid_index(j))
      end do
!
!   Determine if the self-energy contribution is to be estimated
!
      write (stderr_unit, *) 'Estimate contributions from the self-energy?'
      lse = getyn()
      if (lse) then
        write(734,'(a)') 'y            ! Estimate contributions from the self-energy?'
      else
        write(734,'(a)') 'n            ! Estimate contributions from the self-energy?'
      end if

      return
      end subroutine lodres
