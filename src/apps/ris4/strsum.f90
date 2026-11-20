!***********************************************************************
!                                                                      *
      subroutine strsum
      use calculation_file_units, only: summary_output_unit
!                                                                      *
!   Generates the first part of  grasp92.sum  (on stream 24).          *
!                                                                      *
!   Call(s) to: [LIB92] CALEN, convrt2.                                *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 09 Dec 1992   *
!   Modified by Xinghong He               Last revision: 22 Dec 1997   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
      use atomic_calculation_state, only: electron_count, max_radial_points, nuclear_charge, speed_of_light_au
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
      use nuclear_potential_state, only: first_grid_index_beyond_fermi_c, max_radial_points, negative_radius_times_potential
      use orbital_csf_basis, only: csf_count, max_orbitals, orbital_count, orbital_energy_parameters, &
          orbital_origin_exponents, orbital_principal_n, orbital_symmetry_labels
      use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
      use breit_frequency_scale, only: wfact
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use level_symmetry_metadata, only: iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
      use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use integer_format_with_context, only: convrt2
      use engout_i, only: engout
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer   :: lenth, nb, iccut, orbital_index, iend, ibeg
      character :: record*256, cdata*26, clevel*2
!-----------------------------------------------
!
!   Get the date and time of day; make this information the
!   header of the summary file
!
!   Write out the basic dimensions of the electron cloud
!
      write (summary_output_unit, *)
      call convrt2 (electron_count, record, lenth, "strsum.nelec")
      write (summary_output_unit, *) "There are "//record(1:lenth)//" electrons in the cloud"
      call convrt2 (csf_count, record, lenth, "strsum.ncf")
      write (summary_output_unit, *) " in "//record(1:lenth)//" relativistic CSFs"
      call convrt2 (orbital_count, record, lenth, "strsum.nw")
      write (summary_output_unit, *) " based on "//record(1:lenth)//" relativistic subshells."
!
!   If the CSFs are not treated uniformly, write out an
!   informative message
!
      if (lfordr) then
         call convrt2 (iccut, record, lenth, "strsum.icccut")
         write (summary_output_unit, *) " CSFs 1--"//record(1:lenth)//" constitute"//   &
              " the zero-order space;"
      end if
!
!   Write out the nuclear parameters
!
      write (summary_output_unit, *)
      write (summary_output_unit, '("The atomic number is ",1f14.10,";")') nuclear_charge
      if (nuclear_model_parameter_count == 2) then
         write (summary_output_unit, *) " Fermi nucleus:"
         write (summary_output_unit, '("  c =",1p,1d19.12," Bohr radii,"/"  a =",1d19.12," Bohr radii;")') &
            nuclear_model_parameters(1), nuclear_model_parameters(2)
         call convrt2 (first_grid_index_beyond_fermi_c, record, lenth, "strsum.nnuc")
         write (summary_output_unit, *) " there are "//record(1:lenth)//&
            " tabulation points in the nucleus."
      else
         write (summary_output_unit, *) " point nucleus."
      end if
!
!   Write out the physical effects specifications
!
      write (summary_output_unit, *)
      write (summary_output_unit, '("Speed of light = ",1pd19.12," atomic units.")') speed_of_light_au
!
!   Write out the parameters of the radial grid
!
      write (summary_output_unit, *)
      if (asymptotic_grid_step == 0.0_real64) then
         write (summary_output_unit, '("Radial grid: R(I) = RNT*(exp((I-1)*H)-1),"," I = 1, ..., N;"//" RNT  = ",1p,d19.12,&
            &" Bohr radii;"/" H    = ",d19.12," Bohr radii;"/" N    = ",1i4,";")') &
            radial_grid_scale, transformed_grid_step, radial_point_count
      else
         write (summary_output_unit, '("Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,"," I = 1, ..., N;"//" RNT  = ",1p,&
            &d19.12,&
            &" Bohr radii;"/" H    = ",d19.12," Bohr radii;"/" HP   = ",d19.12," Bohr radii;"/" N    = ",1i4,";")') &
            radial_grid_scale, transformed_grid_step, asymptotic_grid_step, radial_point_count
      end if
      write (summary_output_unit, '(" R(1) = ",1p,1d19.12," Bohr radii;"/" R(2) = ",1d19.12," Bohr radii;"/" R(N) = ",1d19.12,&
         &" Bohr radii.")') radial_grid(1), radial_grid(2), radial_grid(radial_point_count)
!
!   Write out the orbital properties
!
      write (summary_output_unit, *)
      write (summary_output_unit, *) "Subshell radial wavefunction summary:"
      write (summary_output_unit, *)
      write (summary_output_unit, '(" Subshell",11x,"e",20x,"p0",18x,"gamma",19x,"P(2)",18x,"Q(2)",10x,"MTP")')
      write (summary_output_unit, *)
      do orbital_index = 1, orbital_count
         write (summary_output_unit, '(3x,1i2,1a2,1x,1p,5(3x,1d19.12),3x,1i3)') &
            orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index), orbital_energy_parameters(orbital_index), &
            & orbital_origin_parameters(orbital_index), orbital_origin_exponents(orbital_index), radial_large_components(2, &
               & orbital_index), &
            & radial_small_components(2,orbital_index), &
            orbital_last_grid_index(orbital_index)
      end do
!
!   Write the list of eigenpair indices
!
      write (summary_output_unit,*)
      call engout (hamiltonian_energy_shift,shifted_eigenvalues,iatjpo,iaspar,ivec,nvec,3)
!
      return
!
!
      end subroutine strsum
