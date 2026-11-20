!***********************************************************************
!                                                                      *
      subroutine lodres
      use calculation_file_units, only: calculation_log_unit
!                                                                      *
!   Loads the data from the  .res  file. A number of checks are made   *
!   to ensure correctness and consistency.                             *
!   This subroutine is called by every node and the only differences   *
!   from the non-mpi version are:                                      *
!     1) STOP and stopmpi (not necessary)                              *
!     2) LSE needs coomunication - can be moved out of this routine.   *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, GETYN, SETQIC.                         *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Oct 1992   *
!   Block version by Xinghong He          Last revision:  1 Jun 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use grasp_parameters,   only: max_radial_points
      use array_memory_management, only: alloc
      use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
      use atomic_calculation_state, only: nelecr, electron_count, nuclear_charge, emn, speed_of_light_au, accy
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
      use nuclear_potential_state, only: negative_radius_times_potential, first_grid_index_beyond_fermi_c
      use orbital_csf_basis, only: csf_count, max_orbitals, orbital_count, orbital_energy_parameters, &
          orbital_origin_exponents
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
      use breit_frequency_scale, only: wfact
      use where_c, only: restart_matrix_unit
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use csf_block_cutoffs, only: iccut, iccutblk
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_logical, mpi_rank
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use yes_no_input, only: getyn
      use radial_stencil_weights, only: setqic
      use mpi_utilities, only: stopmpi
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ncfres, nwres, nblockres, i, np10, orbital_index
      logical :: yes
!-----------------------------------------------
!
!     POINTER (pncfblk, ncfblk(0:*))
!
!     POINTER (piccutblk, iccutblk(1))
!
!-----------------------------------------------------------------------
!
!   Read the basic parameters of the electron cloud; check these
!   against those deduced from the  .csl file
!
      read (restart_matrix_unit) nelecr, ncfres, nwres, nblockres

      if (nelecr/=electron_count .or. ncfres/=csf_count .or. nwres/=orbital_count .or. nblockres/=symmetry_block_count&
         ) call stopmpi ("lodres: NELEC/NCF/NW does not match",mpi_rank)
!
!   Read the nuclear parameters
!
      read (restart_matrix_unit) nuclear_charge, emn
      read (restart_matrix_unit) nuclear_model_parameter_count, (nuclear_model_parameters(i),i=1,nuclear_model_parameter_count)
      read (restart_matrix_unit) radial_point_count, (negative_radius_times_potential(i),i=1,radial_point_count), &
         first_grid_index_beyond_fermi_c

      if (radial_point_count > max_radial_points) call stopmpi ("lodres: N greater than NNNP",mpi_rank)
!
!   Read the physical effects specifications
!   iccutblk() is now an array of length nblock.
!
      read (restart_matrix_unit) speed_of_light_au, lfordr, (iccutblk(i),i=1,symmetry_block_count), ltrans, wfact, lvp, &
         lnms, lsms
!
!   Read the remaining parameters controlling the radial grid and the
!   grid arrays
!
      np10 = radial_point_count + 10
      read (restart_matrix_unit) radial_grid_scale, transformed_grid_step, asymptotic_grid_step, (radial_grid(i),i=1,np10), &
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
      call alloc (radial_large_components, max_radial_points,orbital_count, "PF", "LODMIX")
      call alloc (radial_small_components, max_radial_points,orbital_count, "QF", "LODMIX")
!
!   Read the orbital wavefunctions and the associated arrays
!
      do orbital_index = 1, orbital_count
         read (restart_matrix_unit) orbital_energy_parameters(orbital_index), orbital_origin_exponents(orbital_index), &
            & orbital_origin_parameters(orbital_index), &
            & orbital_last_grid_index(orbital_index)
         read (restart_matrix_unit) (radial_large_components(i,orbital_index),i=1,orbital_last_grid_index(orbital_index)), &
            & (radial_small_components(i, &
            orbital_index),i=1, &
            & orbital_last_grid_index(orbital_index))
      end do
!
!   Determine if the self-energy contribution is to be estimated
!
      if (mpi_rank == 0) then
         write (stderr_unit, *) "Estimate contributions from the self-energy?"
         lse = getyn()
         if (lse) then
           write(calculation_log_unit,"(a)") "y            ! Estimate contributions from the self-energy?"
         else
           write(calculation_log_unit,"(a)") "n            ! Estimate contributions from the self-energy?"
         end if
      end if
      call mpi_bcast (lse, 1, mpi_logical, 0, mpi_comm_world, mpi_error_code)
      return
      end subroutine lodres
