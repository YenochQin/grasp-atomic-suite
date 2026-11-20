! A separate archive entry lets RCI replace TALK at link time.
subroutine talk(bra_csf_index, ket_csf_index, radial_integral_order, bra_orbital_1, bra_orbital_2, &
    ket_orbital_1, ket_orbital_2, integral_type_code, coefficient)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use breit_angular_coefficient_buffer, only: talk_fallback
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: bra_csf_index, ket_csf_index, radial_integral_order
    integer, intent(in) :: bra_orbital_1, bra_orbital_2, ket_orbital_1, ket_orbital_2, integral_type_code
    real(real64), intent(in) :: coefficient

    call talk_fallback(bra_csf_index, ket_csf_index, radial_integral_order, bra_orbital_1, bra_orbital_2, &
        ket_orbital_1, ket_orbital_2, integral_type_code, coefficient)
end subroutine talk
