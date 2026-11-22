! Collect angular coefficients and radial-integral identifiers for Breit matrix assembly.
! Store coefficients with magnitude above EPS and their labels in *
! the angular buffer. Print a diagnostic row if IBUG1 is nonzero.
!
! Last update: 14 Oct 1992
subroutine collect_breit_angular_coefficient(bra_csf_index, ket_csf_index, radial_integral_order, bra_orbital_1, bra_orbital_2, &
    ket_orbital_1, ket_orbital_2, integral_type_code, coefficient)
    ! Module dependencies
    use debug_output_settings, only: angular_debug_unit, ensure_debug_output
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: orbital_index_radix
    use angular_coefficient_buffer_state, only: angular_buffer_capacity, angular_coefficient_count, &
        angular_coefficient_labels, angular_coefficients
    use debug_output_settings, only: angular_coefficient_debug_flag
    use orbital_csf_basis, only: orbital_principal_n, orbital_symmetry_labels
    use numeric_constants, only: angular_tolerance
    use angular_coefficient_buffer, only: manage_angular_coefficient_buffer
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: bra_csf_index, ket_csf_index, radial_integral_order
    integer, intent(in) :: bra_orbital_1, bra_orbital_2, ket_orbital_1, ket_orbital_2, integral_type_code
    real(real64), intent(in) :: coefficient

    ! Local constants
    !
    !   Print coefficient if requested
    if (angular_coefficient_debug_flag /= 0) then
        call ensure_debug_output()
        write (angular_debug_unit, '(2(1x,1i2),4(1x,i2,a2),1x,1i2,1x,1i2,1x,1p,d19.12)') &
            bra_csf_index, ket_csf_index, &
            orbital_principal_n(bra_orbital_1), orbital_symmetry_labels(bra_orbital_1), &
            orbital_principal_n(bra_orbital_2), orbital_symmetry_labels(bra_orbital_2), &
            orbital_principal_n(ket_orbital_1), orbital_symmetry_labels(ket_orbital_1), &
            orbital_principal_n(ket_orbital_2), orbital_symmetry_labels(ket_orbital_2), &
            radial_integral_order, integral_type_code, coefficient
    end if
    ! Increment counter
    if(dabs(coefficient) > angular_tolerance) then
        angular_coefficient_count = angular_coefficient_count + 1
        ! Ensure that arrays are of adequate size; reallocate if necessary
        if (angular_coefficient_count > angular_buffer_capacity) call manage_angular_coefficient_buffer (2)
        ! Store integral indices and coefficient in COMMON/BUFFER/
        angular_coefficient_labels(1,angular_coefficient_count) = bra_orbital_1
        angular_coefficient_labels(2,angular_coefficient_count) = bra_orbital_2
        angular_coefficient_labels(3,angular_coefficient_count) = ket_orbital_1
        angular_coefficient_labels(4,angular_coefficient_count) = ket_orbital_2
        angular_coefficient_labels(5,angular_coefficient_count) = radial_integral_order
        angular_coefficient_labels(6,angular_coefficient_count) = integral_type_code
        angular_coefficients(angular_coefficient_count)   = coefficient
    end if
    return
    return
end subroutine collect_breit_angular_coefficient
