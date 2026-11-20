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
        !   Generates the first part of  erwf.sum  (on stream 24).             *
        !                                                                      *
        !   Call(s) to: [LIB92] CALEN, CONVRT.                                 *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 17 Dec 1992   *
        !                                                                      *
        !***********************************************************************
        !
        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state, only: electron_count, nuclear_charge, speed_of_light_au
        use radial_grid_state, only: radial_point_count, radial_grid, radial_grid_scale, transformed_grid_step, &
            asymptotic_grid_step
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
        use nuclear_potential_state, only: first_grid_index_beyond_fermi_c
        use orbital_csf_basis, only: csf_count, orbital_count
        use calendar_clock, only: read_calendar_clock
        use integer_format, only: format_integer
        implicit none (type, external)

        ! Local declarations
        integer :: text_length
        character(len=10) :: time_text
        character(len=8) :: date_text
        character(len=256) :: record

        !-----------------------------------------------
        !
        !
        !   Get the date and time of day; make this information the
        !   header of the summary file
        !
        call read_calendar_clock (time_text, date_text)
        write (summary_output_unit, *) "ERWF run at ", time_text, " on ", date_text, "."
        !
        !   Write out the basic dimensions of the electron cloud
        !
        write (summary_output_unit, *)
        call format_integer (electron_count, record, text_length)
        write (summary_output_unit, *) "There are "//record(1:text_length)//" electrons in the cloud"
        call format_integer (csf_count, record, text_length)
        write (summary_output_unit, *) " in "//record(1:text_length)//" relativistic CSFs"
        call format_integer (orbital_count, record, text_length)
        write (summary_output_unit, *) " based on "//record(1:text_length)//" relativistic subshells."
        !
        !   Write out the nuclear parameters
        !
        write (summary_output_unit, *)
        write (summary_output_unit, '("The atomic number is ",1f14.10,";")') nuclear_charge
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
        !
        !   Write out the physical effects specifications
        !
        write (summary_output_unit, *)
        write (summary_output_unit, '("Speed of light = ",1p,d19.12," atomic units.")') speed_of_light_au
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
        write (summary_output_unit, *)
        !
        return
        !
        return
        !
    end subroutine write_input_summary
end module write_input_summary_interface
