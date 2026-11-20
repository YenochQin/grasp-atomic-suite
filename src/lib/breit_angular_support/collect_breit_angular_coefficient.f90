! Keep this default archive entry separate so RCI can supply its own collector.
subroutine collect_breit_angular_coefficient(bra_csf_index, ket_csf_index, radial_integral_order, bra_orbital_1, bra_orbital_2, &
    ket_orbital_1, ket_orbital_2, integral_type_code, coefficient)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use breit_angular_coefficient_buffer, only: store_breit_angular_coefficient
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: bra_csf_index, ket_csf_index, radial_integral_order
    integer, intent(in) :: bra_orbital_1, bra_orbital_2, ket_orbital_1, ket_orbital_2, integral_type_code
    real(real64), intent(in) :: coefficient

    call store_breit_angular_coefficient(bra_csf_index, ket_csf_index, radial_integral_order, bra_orbital_1, bra_orbital_2, &
        ket_orbital_1, ket_orbital_2, integral_type_code, coefficient)
end subroutine collect_breit_angular_coefficient
