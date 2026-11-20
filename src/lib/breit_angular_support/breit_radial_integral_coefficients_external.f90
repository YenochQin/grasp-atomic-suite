! A separate archive entry lets RCI replace CXK at link time.
subroutine cxk (radial_integral_coefficients, orbital_indices, twice_orbital_kappa, radial_integral_order, &
    tensor_rank, radial_integral_case, direct_exchange_code)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use breit_radial_integral_coefficients, only: cxk_fallback
    implicit none (type, external)

    ! Arguments
    real(real64), intent(inout) :: radial_integral_coefficients(12)
    integer, intent(in) :: orbital_indices(4), twice_orbital_kappa(4)
    integer, intent(in) :: radial_integral_order, tensor_rank, radial_integral_case, direct_exchange_code

    call cxk_fallback(radial_integral_coefficients, orbital_indices, twice_orbital_kappa, radial_integral_order, &
        tensor_rank, radial_integral_case, direct_exchange_code)
end subroutine cxk
