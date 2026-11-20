!***********************************************************************
!                                                                      *
      subroutine strsum
!                                                                      *
!   Generates the first part of  grasp92.sum  (on stream 24).          *
!                                                                      *
!   Call(s) to: [LIB92] CALEN, convrt2.                                *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 09 Dec 1992   *
!   Modified by Xinghong He               Last revision: 22 Dec 1997   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
      use atomic_calculation_state, only: electron_count, emn, max_radial_points, nuclear_charge, speed_of_light_au
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
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use csf_block_cutoffs,          only: iccutblk
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use integer_format_with_context, only: convrt2
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer   :: lenth, nb, iccut, i, ntmp
      character :: record*256, cdata*26, clevel*2
!-----------------------------------------------
!     POINTER (pncfblk, ncfblk(0:*))
!
!     POINTER (piccutblk, iccutblk(1))

!
!   Get the date and time of day; make this information the
!   header of the summary file
!
!   Write out the basic dimensions of the electron cloud
!
      write (24, *)
      call convrt2 (electron_count, record, lenth, "strsum.nelec")
      write (24, *) "There are "//record(1:lenth)//" electrons in the cloud"
      call convrt2 (csf_count, record, lenth, "strsum.ncf")
      write (24, *) " in "//record(1:lenth)//" relativistic CSFs"
      call convrt2 (orbital_count, record, lenth, "strsum.nw")
      write (24, *) " based on "//record(1:lenth)//" relativistic subshells."
!
!   If the CSFs are not treated uniformly, write out an
!   informative message
!
      if (lfordr) then
         write (24, *)
         do nb = 1, symmetry_block_count
            iccut = iccutblk(nb)
            call convrt2 (iccut, record, lenth, "strsum.icccut")
            write (24, *) " CSFs 1--"//record(1:lenth)//" constitute"//&
               " the zero-order space;  nb = ", nb, " ncf = ", csf_counts_by_block(nb)
         end do
      end if
!
!   Write out the nuclear parameters
!
      write (24, *)
      write (24, '("The atomic number is ",1f14.10,";")') nuclear_charge
      if (emn == 0.0_real64) then
         write (24, *) " the nucleus is stationary;"
      else
         write (24, '(" the mass of the nucleus is ",1p,d19.12," electron masses;")') emn
      end if
      if (nuclear_model_parameter_count == 2) then
         write (24, *) " Fermi nucleus:"
         write (24, '("  c =",1p,1d19.12," Bohr radii,"/,"  a =",1d19.12," Bohr radii;")') &
            nuclear_model_parameters(1), nuclear_model_parameters(2)
         call convrt2 (first_grid_index_beyond_fermi_c, record, lenth, "strsum.nnuc")
         write (24, *) " there are "//record(1:lenth)//&
            " tabulation points in the nucleus."
      else
         write (24, *) " point nucleus."
      end if
!
!   Write out the physical effects specifications
!
      write (24, *)
      write (24, '("Speed of light = ",1p,d19.12," atomic units.")') speed_of_light_au
!
      write (24, *)
      if (ltrans .or. lvp .or. lnms .or. lsms) then
         write (24, *) "To H (Dirac Coulomb) is added"
         if (ltrans) write (24, '("  H (Transverse) --- factor multiplying the"," photon frequency: ",1p,d15.8,&
            &";")') wfact
         if (lvp) write (24, *) " H (Vacuum Polarisation);"
         if (lnms) write (24, *) " H (Normal Mass Shift);"
         if (lsms) write (24, *) " H (Specific Mass Shift);"
         write (24, *) " the total will be diagonalised."
      else
         write (24, *) "H (Dirac Coulomb) will be diagonalised by itself."
      end if
!
      if (lse) then
         write (24, *) &
            "Diagonal contributions from H (Self Energy) will be estimated"
         write (24, *) " from a screened hydrogenic approximation."
      end if
!
!   Write out the parameters of the radial grid
!
      write (24, *)
      if (asymptotic_grid_step == 0.0_real64) then
         write (24, '("Radial grid: R(I) = RNT*(exp((I-1)*H)-1),"," I = 1, ..., N;"/,/," RNT  = ",1p,d19.12,&
            &" Bohr radii;"/," H    = ",d19.12," Bohr radii;"/," N    = ",1i4,";")') &
            radial_grid_scale, transformed_grid_step, radial_point_count
      else
         write (24, '("Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,"," I = 1, ..., N;"/,/," RNT  = ",1p,&
            &d19.12," Bohr radii;"/," H    = ",d19.12," Bohr radii;"/," HP   = ",d19.12," Bohr radii;"/," N    = ",&
            &1i4,";")') radial_grid_scale, transformed_grid_step, asymptotic_grid_step, radial_point_count
      end if
      write (24, '(" R(1) = ",1p,1d19.12," Bohr radii;"/," R(2) = ",1d19.12," Bohr radii;"/," R(N) = ",1d19.12,&
         &" Bohr radii.")') radial_grid(1), radial_grid(2), radial_grid(radial_point_count)
!
!   Write out the orbital properties
!
      write (24, *)
      write (24, *) "Subshell radial wavefunction summary:"
      write (24, *)
      write (24, '("Subshell",6x,"e",13x,"p0",5x,"gamma",5x,"P(2)",7x,"Q(2)",4x,"MTP")')
      write (24, *)
      do i = 1, orbital_count
         write (24, '(1x,i2,a2,1x,1p,d17.10,1p,d11.3,0p,f6.2,1p,2(d11.3),i5)') &
            orbital_principal_n(i), orbital_symmetry_labels(i), orbital_energy_parameters(i), &
            & orbital_origin_parameters(i), orbital_origin_exponents(i), radial_large_components(2,i), &
            & radial_small_components(2,i), &
            orbital_last_grid_index(i)
      end do
!
!   Write the list of eigenpair indices
!
      write (24, *)
!
!   Find total number of eigenstates and print corresponding info
!
      ntmp = 0
      do i = 1, symmetry_block_count
         ntmp = ntmp + eigenvector_counts_by_block(i)
      end do

      write (24,*) ntmp, " levels will be computed"

      return
!
      return
!
      end subroutine strsum
