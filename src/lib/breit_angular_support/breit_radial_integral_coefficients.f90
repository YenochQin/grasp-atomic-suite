module breit_radial_integral_coefficients
    implicit none (type, external)

    ! Module interface
    private
    public :: cxk, cxk_fallback

    ! Procedure interfaces
    ! Derive the external interface from the fallback. RCI supplies cxk at link time.
    procedure(cxk_fallback) :: cxk
contains

    subroutine cxk_fallback (radial_integral_coefficients, orbital_indices, twice_orbital_kappa, radial_integral_order, &
        tensor_rank, radial_integral_case, direct_exchange_code)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        ! The radial order nu labels the integral. The tensor rank k labels the angular coupling.
        integer, intent(in) :: radial_integral_order
        integer, intent(in) :: tensor_rank
        integer, intent(in) :: radial_integral_case
        integer, intent(in) :: direct_exchange_code
        integer, intent(in) :: orbital_indices(4)
        integer, intent(in) :: twice_orbital_kappa(4)
        real(real64), intent(inout) :: radial_integral_coefficients(12)

        ! Stop if the Breit coefficient path is reached without an application CXK.
        stop "CXK: Error "
    end subroutine cxk_fallback
end module breit_radial_integral_coefficients
