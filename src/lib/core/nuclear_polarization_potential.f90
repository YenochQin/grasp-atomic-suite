! Evaluate dipole and quadrupole polarization terms and compose them with the Coulomb grid.
module nuclear_polarization_potential
    ! Module dependencies
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use grasp_parameters, only: max_radial_points
    use nuclear_polarization_state, only: npol_config
    use nuclear_polarization_input, only: npol_validate_config
    implicit none (type, external)

    ! Module interface
    private
    public :: npol_check_grid, npol_build_grid, npol_compose
contains
    ! Check the active radial range and reject nonfinite or negative radii.
    subroutine npol_check_grid(grid_point_count, radial_grid, status, message)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64

        ! Arguments
        integer, intent(in) :: grid_point_count
        real(real64), intent(in) :: radial_grid(:)
        integer, intent(out) :: status
        character(*), intent(out) :: message
        status = 1
        message = "npol: N must be in 1:NNNP and fit the supplied arrays"
        if (grid_point_count < 1 .or. grid_point_count > max_radial_points .or. grid_point_count > size(radial_grid)) return
        message = "npol: radial grid must be finite and nonnegative"
        if (.not.all(ieee_is_finite(radial_grid(:grid_point_count)))) return
        if (any(radial_grid(:grid_point_count) < 0)) return
        status = 0
        message = ""
    end subroutine npol_check_grid

    subroutine npol_build_grid(config, grid_point_count, radial_grid, dipole_potential, quadrupole_potential, &
        & polarization_potential, status, message)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64

        ! Arguments
        type(npol_config), intent(in) :: config
        integer, intent(in) :: grid_point_count
        real(real64), intent(in) :: radial_grid(:)
        real(real64), intent(inout) :: dipole_potential(:), quadrupole_potential(:), polarization_potential(:)
        integer, intent(out) :: status
        character(*), intent(out) :: message

        ! Local declarations
        real(real64) :: dipole_candidate(max_radial_points), quadrupole_candidate(max_radial_points), &
            & total_candidate(max_radial_points)
        integer :: radial_index
        call npol_check_grid(grid_point_count, radial_grid, status, message)
        if (status /= 0) return
        status = 1
        message = "npol: output arrays shorter than N"
        if (min(size(dipole_potential), size(quadrupole_potential), size(polarization_potential)) < grid_point_count) return
        call npol_validate_config(config, real(config%isotope_z,real64), status, message)
        if (status /= 0) return
        dipole_candidate(:grid_point_count) = 0
        quadrupole_candidate(:grid_point_count) = 0
        do radial_index = 1, grid_point_count
            if (config%enabled .and. config%use_e1) then
                call channel(radial_grid(radial_index), config%k1, config%b1, 4, dipole_candidate(radial_index), status)
                if (status /= 0) exit
            end if
            if (config%enabled .and. config%use_e2) then
                call channel(radial_grid(radial_index), config%k2, config%b2, 6, quadrupole_candidate(radial_index), status)
                if (status /= 0) exit
            end if
            if (abs(dipole_candidate(radial_index)) > huge(1.0_real64)-abs(quadrupole_candidate(radial_index))) then
                status = 1
                exit
            end if
            total_candidate(radial_index) = dipole_candidate(radial_index) + quadrupole_candidate(radial_index)
        end do
        if (status /= 0) then
            message = "npol: potential is not representable as a finite value"
            return
        end if
        dipole_potential(:grid_point_count) = dipole_candidate(:grid_point_count)
        quadrupole_potential(:grid_point_count) = quadrupole_candidate(:grid_point_count)
        polarization_potential(:grid_point_count) = total_candidate(:grid_point_count)
    end subroutine npol_build_grid

    ! Scale both radii and combine binary exponents to avoid overflowing intermediate powers.
    subroutine channel(radius, channel_strength, cutoff_radius, denominator_power, channel_potential, status)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64

        ! Arguments
        real(real64), intent(in) :: radius, channel_strength, cutoff_radius
        integer, intent(in) :: denominator_power
        real(real64), intent(out) :: channel_potential
        integer, intent(out) :: status

        ! Local declarations
        real(real64) :: radius_scale, potential_mantissa, denominator
        integer :: potential_exponent
        status = 0
        channel_potential = 0
        if (channel_strength == 0) return
        radius_scale = max(radius,cutoff_radius)
        denominator = (radius/radius_scale)**denominator_power + (cutoff_radius/radius_scale)**denominator_power
        potential_mantissa = fraction(channel_strength)/(fraction(radius_scale)**denominator_power*denominator)
        potential_exponent = exponent(channel_strength) - denominator_power*exponent(radius_scale) + &
            & exponent(potential_mantissa)
        if (potential_exponent > maxexponent(channel_potential)) then
            status = 1
            return
        end if
        if (potential_exponent < minexponent(channel_potential)-digits(channel_potential)) return
        channel_potential = -scale(fraction(potential_mantissa), potential_exponent)
        if (.not.ieee_is_finite(channel_potential)) status = 1
    end subroutine channel

    ! Add -r*V_NP to a temporary Coulomb grid and publish it only after validation.
    subroutine npol_compose(grid_point_count, radial_grid, coulomb_grid, polarization_potential, composed_grid, status, message)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64

        ! Arguments
        integer, intent(in) :: grid_point_count
        real(real64), intent(in) :: radial_grid(:), coulomb_grid(:), polarization_potential(:)
        real(real64), intent(inout) :: composed_grid(:)
        integer, intent(out) :: status
        character(*), intent(out) :: message

        ! Local declarations
        real(real64) :: candidate_grid(max_radial_points), weighted_potential_shift
        integer :: radial_index
        call npol_check_grid(grid_point_count, radial_grid, status, message)
        if (status /= 0) return
        status = 1
        message = "npol: composition arrays shorter than N"
        if (min(size(coulomb_grid), size(polarization_potential), size(composed_grid)) < grid_point_count) return
        message = "npol: nonfinite Coulomb or NP input"
        if (.not.all(ieee_is_finite(coulomb_grid(:grid_point_count))) .or. &
            & .not.all(ieee_is_finite(polarization_potential(:grid_point_count)))) return
        candidate_grid(:grid_point_count) = coulomb_grid(:grid_point_count)
        message = "npol: composed potential is not finite"
        do radial_index = 1, grid_point_count
            if (polarization_potential(radial_index) == 0 .or. radial_grid(radial_index) == 0) cycle
            if (radial_grid(radial_index) > 1) then
                if (abs(polarization_potential(radial_index)) > huge(weighted_potential_shift)/radial_grid(radial_index)) return
            end if
            weighted_potential_shift = -radial_grid(radial_index)*polarization_potential(radial_index)
            if (weighted_potential_shift > 0) then
                if (coulomb_grid(radial_index) > huge(weighted_potential_shift)-weighted_potential_shift) return
            else
                if (coulomb_grid(radial_index) < -huge(weighted_potential_shift)-weighted_potential_shift) return
            end if
            candidate_grid(radial_index) = coulomb_grid(radial_index) + weighted_potential_shift
        end do
        composed_grid(:grid_point_count) = candidate_grid(:grid_point_count)
        status = 0
        message = ""
    end subroutine npol_compose
end module nuclear_polarization_potential
