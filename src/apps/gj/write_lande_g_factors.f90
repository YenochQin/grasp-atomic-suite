! Write the electronic Landé factor, its correction, and their sum for each selected level.
module write_lande_g_factors_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_lande_g_factors
contains

    ! Writes diagonal Landé-factor data to <name>.gj or <name>.cigj.
    subroutine write_lande_g_factors(state_name, mixing_file_kind, diagonal_lande_g_elements, diagonal_lande_correction_elements)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state,           only: vacuum_speed_of_light_au
        use angular_symmetry_labels,         only: format_j_label, parity_labels
        use isotope_properties,        only: nuclear_spin, nuclear_magnetic_dipole_moment, nuclear_electric_quadrupole_moment
        use selected_level_indices,          only: active_level_count, active_level_indices
        use level_symmetry_metadata,          only: level_two_j_plus_one, level_parity_signs
        use file_open, only: open_calculation_file
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: lande_output_unit
        character(len=24), intent(in) :: state_name
        integer, intent(in) :: mixing_file_kind
        real(real64), intent(in) :: diagonal_lande_g_elements(active_level_count), &
                & diagonal_lande_correction_elements(active_level_count)
        integer :: level_index, open_status, two_j_plus_one, state_name_end
        real(real64) :: electronic_j, inverse_j_norm, lande_g, lande_g_correction
        character :: file_path*256, file_form*11, file_status*3

        state_name_end = index(state_name," ")
        if (mixing_file_kind == 0) then
            file_path = state_name(1:state_name_end-1)//".cigj"
        else
            file_path = state_name(1:state_name_end-1)//".gj"
        end if
        file_form = "FORMATTED"
        file_status = "NEW"
        call open_calculation_file(lande_output_unit, file_path, file_form, file_status, open_status)
        if (open_status /= 0) then
            write (output_unit, *) "Error when opening", file_path
            stop
        end if
        write (lande_output_unit, '("Nuclear spin                        ",1p,d22.15," au")') nuclear_spin
        write (lande_output_unit, '("Nuclear magnetic dipole moment      ",1p,d22.15," n.m.")') nuclear_magnetic_dipole_moment
        write (lande_output_unit, '("Nuclear electric quadrupole moment  ",1p,d22.15," barns")') nuclear_electric_quadrupole_moment
        write (lande_output_unit, &
            & '(/,/," Interaction constants:"/,/," Level1  J Parity ",10x,"g_J",14x,"delta g_J",11x,"total g_J"/&
            &)')
        do level_index = 1, active_level_count
            two_j_plus_one = level_two_j_plus_one(level_index)
            if (two_j_plus_one <= 1) cycle
            electronic_j = 0.5_real64*dble(two_j_plus_one - 1)
            inverse_j_norm = sqrt(1.0_real64/(electronic_j*(electronic_j + 1.0_real64)))
            lande_g = vacuum_speed_of_light_au*inverse_j_norm*diagonal_lande_g_elements(level_index)
            lande_g_correction = 0.001160_real64*inverse_j_norm*diagonal_lande_correction_elements(level_index)
            write (lande_output_unit, '(1x,1i3,5x,2a,1p,3d20.10)') active_level_indices(level_index), &
                & format_j_label(two_j_plus_one), &
                & parity_labels((level_parity_signs(level_index)+3)/2), lande_g, lande_g_correction, &
                lande_g + lande_g_correction
        end do
        close(lande_output_unit)
        return
    end subroutine write_lande_g_factors
end module write_lande_g_factors_interface
