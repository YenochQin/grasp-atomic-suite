module write_input_summary_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_input_summary
contains

    !***********************************************************************
    !                                                                      *
    subroutine write_input_summary
        ! Module dependencies
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
        use, intrinsic :: iso_fortran_env, only: real64
        use hamiltonian_correction_flags, only: use_first_order_csf_perturbation, include_normal_mass_shift, &
            estimate_self_energy, include_specific_mass_shift, include_transverse_interaction, include_vacuum_polarization
        use atomic_calculation_state, only: electron_count, nuclear_mass_in_electron_masses, nuclear_charge, speed_of_light_au
        use radial_grid_state, only: radial_point_count, radial_grid, radial_grid_scale, transformed_grid_step, &
            asymptotic_grid_step
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
        use nuclear_potential_state, only: first_grid_index_beyond_fermi_c
        use orbital_csf_basis, only: csf_count, orbital_count, orbital_energy_parameters, orbital_origin_exponents, &
            orbital_principal_n, orbital_symmetry_labels
        use radial_wavefunctions, only: orbital_last_grid_index, orbital_origin_parameters, radial_large_components, &
            radial_small_components
        use breit_frequency_scale, only: transverse_photon_frequency_scale
        use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, symmetry_block_count
        use csf_block_cutoffs,          only: csf_cutoff_by_block
        use integer_format_with_context, only: format_integer_with_context
        implicit none (type, external)

        ! Local declarations
        integer   :: lenth, symmetry_block_index, iccut, i, ntmp
        character :: record*256
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
        write (summary_output_unit, *)
        call format_integer_with_context (electron_count, record, lenth, "write_input_summary: electron_count")
        write (summary_output_unit, *) "There are "//record(1:lenth)//" electrons in the cloud"
        call format_integer_with_context (csf_count, record, lenth, "write_input_summary: csf_count")
        write (summary_output_unit, *) " in "//record(1:lenth)//" relativistic CSFs"
        call format_integer_with_context (orbital_count, record, lenth, "write_input_summary: orbital_count")
        write (summary_output_unit, *) " based on "//record(1:lenth)//" relativistic subshells."
        !
        !   If the CSFs are not treated uniformly, write out an
        !   informative message
        !
        if (use_first_order_csf_perturbation) then
            write (summary_output_unit, *)
            do symmetry_block_index = 1, symmetry_block_count
                iccut = csf_cutoff_by_block(symmetry_block_index)
                call format_integer_with_context (iccut, record, lenth, "write_input_summary: iccut")
                write (summary_output_unit, *) " CSFs 1--"//record(1:lenth)//" constitute"//&
                    " the zero-order space;  nb = ", symmetry_block_index, " ncf = ", csf_counts_by_block(symmetry_block_index)
            end do
        end if
        !
        !   Write out the nuclear parameters
        !
        write (summary_output_unit, *)
        write (summary_output_unit, '("The atomic number is ",1f14.10,";")') nuclear_charge
        if (nuclear_mass_in_electron_masses == 0.0_real64) then
            write (summary_output_unit, *) " the nucleus is stationary;"
        else
            write (summary_output_unit, '(" the mass of the nucleus is ",1p,d19.12," electron masses;")') &
                nuclear_mass_in_electron_masses
        end if
        if (nuclear_model_parameter_count == 2) then
            write (summary_output_unit, *) " Fermi nucleus:"
            write (summary_output_unit, '("  c =",1p,1d19.12," Bohr radii,"/,"  a =",1d19.12," Bohr radii;")') &
                nuclear_model_parameters(1), nuclear_model_parameters(2)
            call format_integer_with_context (first_grid_index_beyond_fermi_c, record, lenth, &
                "write_input_summary: first_grid_index_beyond_fermi")
            write (summary_output_unit, *) " there are "//record(1:lenth)//&
                " tabulation points in the nucleus."
        else
            write (summary_output_unit, *) " point nucleus."
        end if
        !
        !   Write out the physical effects specifications
        !
        write (summary_output_unit, *)
        write (summary_output_unit, '("Speed of light = ",1p,d19.12," atomic units.")') speed_of_light_au
        !
        write (summary_output_unit, *)
        if (include_transverse_interaction .or. include_vacuum_polarization .or. include_normal_mass_shift .or. &
            include_specific_mass_shift) then
            write (summary_output_unit, *) "To H (Dirac Coulomb) is added"
            if (include_transverse_interaction) write (summary_output_unit, '("  H (Transverse) --- factor multiplying the", &
                & " photon frequency: ",1p,d15.8,&
                &";")') transverse_photon_frequency_scale
            if (include_vacuum_polarization) write (summary_output_unit, *) " H (Vacuum Polarisation);"
            if (include_normal_mass_shift) write (summary_output_unit, *) " H (Normal Mass Shift);"
            if (include_specific_mass_shift) write (summary_output_unit, *) " H (Specific Mass Shift);"
            write (summary_output_unit, *) " the total will be diagonalised."
        else
            write (summary_output_unit, *) "H (Dirac Coulomb) will be diagonalised by itself."
        end if
        !
        if (estimate_self_energy) then
            write (summary_output_unit, *) &
                "Diagonal contributions from H (Self Energy) will be estimated"
            write (summary_output_unit, *) " from a screened hydrogenic approximation."
        end if
        !
        !   Write out the parameters of the radial grid
        !
        write (summary_output_unit, *)
        if (asymptotic_grid_step == 0.0_real64) then
            write (summary_output_unit, '("Radial grid: R(I) = RNT*(exp((I-1)*H)-1),"," I = 1, ..., N;"/,/," RNT  = ",1p,d19.12,&
                &" Bohr radii;"/," H    = ",d19.12," Bohr radii;"/," N    = ",1i4,";")') &
                radial_grid_scale, transformed_grid_step, radial_point_count
        else
            write (summary_output_unit, '("Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,"," I = 1, ..., N;"/,/," RNT  = ",1p,&
                &d19.12," Bohr radii;"/," H    = ",d19.12," Bohr radii;"/," HP   = ",d19.12," Bohr radii;"/," N    = ",&
                &1i4,";")') radial_grid_scale, transformed_grid_step, asymptotic_grid_step, radial_point_count
        end if
        write (summary_output_unit, '(" R(1) = ",1p,1d19.12," Bohr radii;"/," R(2) = ",1d19.12," Bohr radii;"/," R(N) = ",1d19.12,&
            &" Bohr radii.")') radial_grid(1), radial_grid(2), radial_grid(radial_point_count)
        !
        !   Write out the orbital properties
        !
        write (summary_output_unit, *)
        write (summary_output_unit, *) "Subshell radial wavefunction summary:"
        write (summary_output_unit, *)
        write (summary_output_unit, '("Subshell",6x,"e",13x,"p0",5x,"gamma",5x,"P(2)",7x,"Q(2)",4x,"MTP")')
        write (summary_output_unit, *)
        do i = 1, orbital_count
            write (summary_output_unit, '(1x,i2,a2,1x,1p,d17.10,1p,d11.3,0p,f6.2,1p,2(d11.3),i5)') &
                orbital_principal_n(i), orbital_symmetry_labels(i), orbital_energy_parameters(i), &
                & orbital_origin_parameters(i), orbital_origin_exponents(i), radial_large_components(2,i), &
                & radial_small_components(2,i), &
                orbital_last_grid_index(i)
        end do
        !
        !   Write the list of eigenpair indices
        !
        write (summary_output_unit, *)
        !
        !   Find total number of eigenstates and print corresponding info
        !
        ntmp = 0
        do i = 1, symmetry_block_count
            ntmp = ntmp + eigenvector_counts_by_block(i)
        end do

        write (summary_output_unit,*) ntmp, " levels will be computed"

        return
        !
        return
        !
    end subroutine write_input_summary
end module write_input_summary_interface
