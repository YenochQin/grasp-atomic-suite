! Orthogonalize orbitals of equal kappa with Gram-Schmidt, normalize them, and trim their tails.
module radial_orbital_orthogonalization
    implicit none (type, external)

    ! Module interface
    private
    public :: orthogonalize_radial_orbitals
contains
    subroutine orthogonalize_radial_orbitals
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only:  max_orbitals
        use debug_output_settings, only: rmcdhf_debug_flags
        use atomic_calculation_state, only: working_tolerance
        use orbital_csf_basis, only: max_orbitals, orbital_count, orbital_kappa, orbital_principal_n, orbital_symmetry_labels
        use radial_wavefunctions, only: orbital_origin_parameters, orbital_last_grid_index, radial_large_components, &
            radial_small_components
        use weighted_radial_overlap, only: integrate_weighted_radial_overlap
        implicit none (type, external)

        ! Local declarations
        integer :: projected_orbital_indices(max_orbitals)
        integer :: orbital_index, orbital_kappa_value, projection_count, combined_last_grid_index, previous_orbital_index, &
            & last_radial_index, projection_index
        real(real64) :: projection_overlaps(max_orbitals)
        real(real64) :: tail_cutoff, orbital_overlap, norm_squared, normalization_factor
        logical :: projected
        tail_cutoff = 0.01_real64*working_tolerance
        do orbital_index = 2, orbital_count
            orbital_kappa_value = orbital_kappa(orbital_index)
            projection_count = 0

            combined_last_grid_index = orbital_last_grid_index(orbital_index)
            projected = .false.
            do previous_orbital_index = 1, orbital_index - 1
                if (orbital_kappa(previous_orbital_index) /= orbital_kappa_value) cycle
                projected = .true.
                orbital_overlap = integrate_weighted_radial_overlap(orbital_index,previous_orbital_index,0)
                projection_count = projection_count + 1
                projected_orbital_indices(projection_count) = previous_orbital_index
                projection_overlaps(projection_count) = orbital_overlap
                ! Project the leading origin coefficient with the same overlap as P and Q.
                orbital_origin_parameters(orbital_index) = orbital_origin_parameters(orbital_index) - &
                    & orbital_overlap*orbital_origin_parameters(previous_orbital_index)
                last_radial_index = max(orbital_last_grid_index(orbital_index),orbital_last_grid_index(previous_orbital_index))
                combined_last_grid_index = max(combined_last_grid_index,orbital_last_grid_index(previous_orbital_index))

                radial_large_components(:last_radial_index,orbital_index) = radial_large_components(:last_radial_index, &
                    & orbital_index) - orbital_overlap*radial_large_components(:last_radial_index,previous_orbital_index)
                radial_small_components(:last_radial_index,orbital_index) = radial_small_components(:last_radial_index, &
                    & orbital_index) - orbital_overlap*radial_small_components(:last_radial_index,previous_orbital_index)
            end do

            if (projected) then
                last_radial_index = combined_last_grid_index

                ! Subtracting orbitals can extend the support. Include the full combined range in the norm.
                orbital_last_grid_index(orbital_index) = last_radial_index
                norm_squared = integrate_weighted_radial_overlap(orbital_index,orbital_index,0)
                normalization_factor = 1.0_real64/sqrt(norm_squared)
                orbital_origin_parameters(orbital_index) = normalization_factor*orbital_origin_parameters(orbital_index)
                radial_large_components(2:last_radial_index, &
                    & orbital_index) = normalization_factor*radial_large_components(2:last_radial_index,orbital_index)
                radial_small_components(2:last_radial_index, &
                    & orbital_index) = normalization_factor*radial_small_components(2:last_radial_index,orbital_index)
                do while (last_radial_index >= 1)
                    if (.not.(abs(radial_large_components(last_radial_index,orbital_index)) < tail_cutoff)) then
                        orbital_last_grid_index(orbital_index) = last_radial_index
                        exit
                    end if
                    radial_large_components(last_radial_index,orbital_index) = 0.0_real64
                    radial_small_components(last_radial_index,orbital_index) = 0.0_real64
                    last_radial_index = last_radial_index - 1
                end do
            end if
            if (.not.(rmcdhf_debug_flags(3) .and. projection_count>0)) cycle

            call ensure_debug_output()
            write (angular_debug_unit, '(1p,5(2x,1d10.3," = <",1i2,1a2,"|",1i2,1a2,">"))') &
                (projection_overlaps(projection_index),orbital_principal_n(orbital_index), &
                & orbital_symmetry_labels(orbital_index),orbital_principal_n(projected_orbital_indices(projection_index)), &
                & orbital_symmetry_labels(projected_orbital_indices(projection_index)),projection_index=1,projection_count)
        end do
        return
        return
    end subroutine orthogonalize_radial_orbitals
end module radial_orbital_orthogonalization
