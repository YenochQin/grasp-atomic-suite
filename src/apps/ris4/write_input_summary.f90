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
        !...Translated by Gediminas Gaigalas 11/18/19
        use, intrinsic :: iso_fortran_env, only: real64
        use hamiltonian_correction_flags, only: use_first_order_csf_perturbation
        use atomic_calculation_state, only: electron_count, nuclear_charge, speed_of_light_au
        use radial_grid_state, only: radial_point_count, radial_grid, radial_grid_scale, transformed_grid_step, &
            asymptotic_grid_step
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
        use nuclear_potential_state, only: first_grid_index_beyond_fermi_c
        use orbital_csf_basis, only: csf_count, orbital_count, orbital_energy_parameters, orbital_origin_exponents, &
            orbital_principal_n, orbital_symmetry_labels
        use selected_level_indices, only: active_level_indices, active_level_count
        use radial_wavefunctions, only: orbital_last_grid_index, orbital_origin_parameters, radial_large_components, &
            radial_small_components
        use level_symmetry_metadata, only: level_parity_signs, level_two_j_plus_one
        use eigenstate_expansions, only: hamiltonian_energy_shift, shifted_eigenvalues
        use integer_format_with_context, only: format_integer_with_context

        ! Procedure interfaces
        use write_level_energies_interface, only: write_level_energies
        implicit none (type, external)

        ! Local declarations
        integer   :: lenth, iccut, orbital_index
        character :: record*256
        !-----------------------------------------------
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
            call format_integer_with_context (iccut, record, lenth, "write_input_summary: iccut")
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
                orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index), orbital_energy_parameters( &
                & orbital_index), &
                & orbital_origin_parameters(orbital_index), orbital_origin_exponents(orbital_index), radial_large_components(2, &
                & orbital_index), &
                & radial_small_components(2,orbital_index), &
                orbital_last_grid_index(orbital_index)
        end do
        !
        !   Write the list of eigenpair indices
        !
        write (summary_output_unit,*)
        call write_level_energies (hamiltonian_energy_shift,shifted_eigenvalues,level_two_j_plus_one,level_parity_signs, &
            active_level_indices,active_level_count,3)
        !
        return
        !
        !
    end subroutine write_input_summary
end module write_input_summary_interface
