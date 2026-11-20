program breit_procedure_hooks
    use, intrinsic :: iso_fortran_env, only: real64
    use angular_coefficient_buffer_state, only: angular_coefficient_count, angular_coefficient_labels, angular_coefficients
    use angular_coefficient_buffer, only: manage_angular_coefficient_buffer
    use breit_angular_coefficient_buffer, only: collect_breit_angular_coefficient, store_breit_angular_coefficient
    use debug_output_settings, only: angular_coefficient_debug_flag
    implicit none (type, external)

    angular_coefficient_count = 0
    call emit_test_coefficient(collect_breit_angular_coefficient, 1.0_real64)
    if (angular_coefficient_count /= -1) error stop "The external Breit collector override was bypassed"

    angular_coefficient_count = 0
    angular_coefficient_debug_flag = 0
    call manage_angular_coefficient_buffer(1)
    call emit_test_coefficient(store_breit_angular_coefficient, 3.0_real64)
    if (angular_coefficient_count /= 1) error stop "The direct store used the external collector"
    if (any(angular_coefficient_labels(:,1) /= [3, 4, 5, 6, 2, 7])) error stop "Unexpected Breit coefficient labels"
    if (angular_coefficients(1) /= 3.0_real64) error stop "Unexpected Breit coefficient value"
    call manage_angular_coefficient_buffer(3)
contains
    subroutine emit_test_coefficient(entry, value)
        procedure(collect_breit_angular_coefficient), pointer, intent(in) :: entry
        real(real64), intent(in) :: value

        call entry(bra_csf_index=1, ket_csf_index=1, radial_integral_order=2, bra_orbital_1=3, bra_orbital_2=4, &
            ket_orbital_1=5, ket_orbital_2=6, integral_type_code=7, coefficient=value)
    end subroutine emit_test_coefficient
end program breit_procedure_hooks

! Stand in for the replacement that RCI supplies when linking the static library.
subroutine collect_breit_angular_coefficient(bra_csf_index, ket_csf_index, radial_integral_order, bra_orbital_1, bra_orbital_2, &
    ket_orbital_1, ket_orbital_2, integral_type_code, coefficient)
    use, intrinsic :: iso_fortran_env, only: real64
    use angular_coefficient_buffer_state, only: angular_coefficient_count
    implicit none (type, external)

    integer, intent(in) :: bra_csf_index, ket_csf_index, radial_integral_order
    integer, intent(in) :: bra_orbital_1, bra_orbital_2, ket_orbital_1, ket_orbital_2, integral_type_code
    real(real64), intent(in) :: coefficient

    if (any([bra_csf_index, ket_csf_index, radial_integral_order, bra_orbital_1, bra_orbital_2, &
        ket_orbital_1, ket_orbital_2, integral_type_code] /= [1, 1, 2, 3, 4, 5, 6, 7])) error stop "Wrong Breit collector arguments"
    angular_coefficient_count = angular_coefficient_count - int(coefficient)
end subroutine collect_breit_angular_coefficient
