! Stores switches for writes from the angular, radial, and property routines.
module debug_output_settings
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: angular_debug_unit = -1
    integer :: angular_coefficient_debug_flag, ibug2, ibug3, ibug4, ibug5, ibug6
    logical :: ldbpa(5)
    logical :: ldbpg(5)
    logical :: rmcdhf_debug_flags(30)
    real(real64) :: matrix_element_screening_cutoff  ! Magnitude cutoff for matrix elements and their contributions.
contains

    subroutine ensure_debug_output()
        ! Workers and standalone callers can emit traces without a named debug file.
        if (angular_debug_unit == -1) then
            open(newunit=angular_debug_unit, file="fort.99", status="unknown", form="formatted")
        end if
    end subroutine ensure_debug_output
end module debug_output_settings
