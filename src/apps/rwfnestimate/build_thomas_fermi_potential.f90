module build_thomas_fermi_potential_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: build_thomas_fermi_potential
contains

    !***********************************************************************
    !                                                                      *
    subroutine build_thomas_fermi_potential
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        !                                                                      *
        !   Calculation of the universal Thomas-Fermi potential.               *
        !                                                                      *
        !   Call(s) to: [LIB92]: DRAW.                                         *
        !                                                                      *
        !                                         Last revision: 09 Dec 1992   *
        !                                                                      *
        !***********************************************************************
        !
        use, intrinsic :: iso_fortran_env, only: real64
        use debug_output_settings, only: rmcdhf_debug_flags
        use atomic_calculation_state, only: nuclear_charge, electron_count
        use radial_grid_state, only: radial_point_count, radial_grid
        use nuclear_potential_state, only: negative_radius_times_potential
        use scf_radial_potentials, only: orbital_y_potential
        use radial_ascii_plot, only: plot_radial_orbital_ascii
        implicit none (type, external)

        ! Local declarations
        integer :: radial_index, complete_output_rows, row_count, output_row, first_column_index, second_column_index, &
            & third_column_index
        real(real64) :: screening_electron_count
        real(real64) :: third, residual_charge, inverse_screening_length, sqrt_scaled_radius, screening_numerator, &
            & screening_denominator, screening_amplitude
        !-----------------------------------------------
        !
        !
        third = 1.0_real64/3.0_real64
        !
        residual_charge = nuclear_charge - dble(electron_count - 1)
        screening_electron_count = max(nuclear_charge - residual_charge,0.0_real64)
        inverse_screening_length = screening_electron_count**third/0.8853_real64
        do radial_index = 1, radial_point_count
            !
            !   Rational function approximation to the universal Thomas-Fermi
            !   function
            !
            sqrt_scaled_radius = sqrt(inverse_screening_length*radial_grid(radial_index))
            screening_numerator = sqrt_scaled_radius*(0.60112_real64*sqrt_scaled_radius + 1.81061_real64) + 1.0_real64
            screening_denominator = sqrt_scaled_radius*(sqrt_scaled_radius*(sqrt_scaled_radius*(sqrt_scaled_radius* &
                & (0.04793_real64*sqrt_scaled_radius + 0.21465_real64) + 0.77112_real64) + &
                1.39515_real64) + 1.81061_real64) + 1.0_real64
            screening_amplitude = screening_numerator/screening_denominator
            orbital_y_potential(radial_index) = (negative_radius_times_potential(radial_index)-residual_charge)* &
                & screening_amplitude*screening_amplitude + residual_charge
        end do
        !
        !   Debug printout
        !
        if (rmcdhf_debug_flags(26)) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/,/,/," Thomas-Fermi potential"/,/,3(" --------- r --------- ------ -r*V(r) ------"))')
            complete_output_rows = radial_point_count/3
            if (3*complete_output_rows == radial_point_count) then
                row_count = complete_output_rows
            else
                row_count = complete_output_rows + 1
            end if
            do output_row = 1, row_count
                first_column_index = output_row
                second_column_index = first_column_index + row_count
                third_column_index = second_column_index + row_count
                if (third_column_index <= radial_point_count) then
                    call ensure_debug_output()
                    write (angular_debug_unit, '(1p,6(1x,1d21.14))') &
                        radial_grid(first_column_index), orbital_y_potential(first_column_index), &
                            & radial_grid(second_column_index), orbital_y_potential(second_column_index), &
                            radial_grid(third_column_index), orbital_y_potential(third_column_index)
                else if (second_column_index <= radial_point_count) then
                    call ensure_debug_output()
                    write (angular_debug_unit, '(1p,6(1x,1d21.14))') radial_grid(first_column_index), &
                        & orbital_y_potential(first_column_index), &
                        radial_grid(second_column_index), orbital_y_potential(second_column_index)
                else
                    call ensure_debug_output()
                    write (angular_debug_unit, '(1p,6(1x,1d21.14))') radial_grid(first_column_index), &
                        & orbital_y_potential(first_column_index)
                end if
            end do
            call plot_radial_orbital_ascii (orbital_y_potential, 1.0_real64, orbital_y_potential, 1.0_real64, radial_point_count)
        end if
        !
        return
        !
        return
        !
    end subroutine build_thomas_fermi_potential
end module build_thomas_fermi_potential_interface
