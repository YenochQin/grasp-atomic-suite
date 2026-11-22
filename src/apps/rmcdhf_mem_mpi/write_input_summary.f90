module write_input_summary_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_input_summary
contains

    ! Generates the first part of  rscf92.sum  (on stream 24).
    !
    ! Written by Farid A. Parpia. Last revision: 26 Sep 1993
    subroutine write_input_summary
        ! Module dependencies
        use calculation_file_units, only: summary_output_unit

        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state, only: electron_count, nuclear_mass_in_electron_masses, requested_level_indices, &
            ionic_charge, selected_level_count, nuclear_charge, speed_of_light_au, requested_level_weights
        use first_order_csf_cutoffs, only: first_order_csf_cutoff
        use radial_grid_state, only: radial_point_count, radial_grid, radial_grid_scale, transformed_grid_step, &
            asymptotic_grid_step
        use rmcdhf_mcp_options, only: mcp_first_order_perturbation
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
        use nuclear_potential_state, only: first_grid_index_beyond_fermi_c
        use orbital_csf_basis, only: csf_count, orbital_count
        use integer_format, only: format_integer
        implicit none (type, external)

        ! Local declarations
        integer :: text_length, iend, i, ibeg
        character :: record*256, number_text*4
        ! Both the nuclear charge and the number of electrons are
        ! known at this point; load IONCTY with the ionicity
        ionic_charge = nint(nuclear_charge) - electron_count
        ! Write a timestamp in the summary header.
        !
        ! Report the electron, orbital, and CSF counts.
        write (summary_output_unit, *)
        call format_integer (electron_count, record, text_length)
        write (summary_output_unit, *) "There are "//record(1:text_length)//" electrons in the cloud"
        call format_integer (csf_count, record, text_length)
        write (summary_output_unit, *) " in "//record(1:text_length)//" relativistic CSFs"
        call format_integer (orbital_count, record, text_length)
        write (summary_output_unit, *) " based on "//record(1:text_length)//" relativistic subshells."
        ! Report the first-order cutoff when the calculation restricts CSF interactions.
        if (mcp_first_order_perturbation) then
            write (summary_output_unit, *)
            call format_integer (first_order_csf_cutoff, record, text_length)
            write (summary_output_unit, *) " CSFs 1--"//record(1:text_length)//" constitute"//&
                " the zero-order space;"
        end if
        ! Report the nuclear model and isotope parameters.
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
            call format_integer (first_grid_index_beyond_fermi_c, record, text_length)
            write (summary_output_unit, *) " there are "//record(1:text_length)//&
                " tabulation points in the nucleus."
        else
            write (summary_output_unit, *) " point nucleus."
        end if
        ! Report the enabled Hamiltonian corrections.
        write (summary_output_unit, *)
        write (summary_output_unit, '("Speed of light = ",3p,d19.12," atomic units.")') speed_of_light_au
        ! Report the radial grid parameters.
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
        write (summary_output_unit, *)
        ! (E)AL calculation, returns here
        if (selected_level_count == 0) then
            write (summary_output_unit, *) "(E)AL calculation."
            return
        end if
        ! Info exclusively for EOL calculations
        if (selected_level_count == 1) then
            write (summary_output_unit, *) "OL calculation."
            call format_integer (requested_level_indices(1), record, text_length)
            write (summary_output_unit, *) "Level "//record(1:text_length)//" will be optimised."
        else
            write (summary_output_unit, *) "EOL calculation."
            call format_integer (selected_level_count, record, text_length)
            write (summary_output_unit, *) record(1:text_length)//" levels will be optimised;"
            record(1:20) = " their indices are: "
            iend = 20
            do i = 1, selected_level_count
                ibeg = iend + 1
                call format_integer (requested_level_indices(i), number_text, text_length)
                if (i /= selected_level_count) then
                    iend = ibeg + text_length + 1
                    record(ibeg:iend) = number_text(1:text_length)//", "
                else
                    iend = ibeg + text_length
                    record(ibeg:iend) = number_text(1:text_length)//"."
                end if
                if (iend < 120) cycle
                write (summary_output_unit, *) record(1:iend)
                record(1:2) = "  "
                iend = 2
            end do
            if (iend /= 2) write (summary_output_unit, *) record(1:iend)
            if (requested_level_weights(1) == (-1.0_real64)) then
                write (summary_output_unit, *) "Each is assigned its statistical weight;"
            else if (requested_level_weights(1) == (-2.0_real64)) then
                write (summary_output_unit, *) "All levels are weighted equally;"
            else
                write (summary_output_unit, *) " weighted as follows:"
                write (summary_output_unit, *) (requested_level_weights(i),i=1,selected_level_count)
            end if
        end if

        return
    end subroutine write_input_summary
end module write_input_summary_interface
