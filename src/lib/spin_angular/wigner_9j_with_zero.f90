! Reduce a Wigner 9j symbol containing zero to a single Wigner 6j symbol.
module wigner_9j_with_zero
    implicit none (type, external)

    ! Module interface
    private
    public :: wigner_9j_j_zero
contains
    subroutine wigner_9j_j_zero(two_j_11,two_j_12,two_j_13,two_j_21,two_j_22,two_j_23,two_j_31,two_j_32,two_j_33,nine_j_value)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, one
        use wigner_6j, only: wigner_6j_value
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: two_j_11, two_j_12, two_j_13, two_j_21, two_j_22, two_j_23, two_j_31, two_j_32, two_j_33
        real(real64), intent(out) :: nine_j_value

        ! Local declarations
        integer      :: twice_phase_exponent
        real(real64) :: six_j_value, normalization_squared
        if (two_j_11 == 0) then
            call wigner_6j_value (two_j_22, two_j_32, two_j_12, two_j_33, two_j_23, two_j_21, 0, six_j_value)
            normalization_squared = dble((two_j_12 + 1)*(two_j_21 + 1))
            twice_phase_exponent = two_j_32 + two_j_12 + two_j_23 + two_j_21
        else if (two_j_12 == 0) then
            call wigner_6j_value (two_j_23, two_j_33, two_j_13, two_j_31, two_j_21, two_j_22, 0, six_j_value)
            normalization_squared = dble((two_j_13 + 1)*(two_j_22 + 1))
            twice_phase_exponent = two_j_33 + two_j_13 + two_j_21 + two_j_22
        else if (two_j_13 == 0) then
            call wigner_6j_value (two_j_21, two_j_31, two_j_11, two_j_32, two_j_22, two_j_23, 0, six_j_value)
            normalization_squared = dble((two_j_11 + 1)*(two_j_23 + 1))
            twice_phase_exponent = two_j_31 + two_j_11 + two_j_22 + two_j_23
        else if (two_j_21 == 0) then
            call wigner_6j_value (two_j_32, two_j_12, two_j_22, two_j_13, two_j_33, two_j_31, 0, six_j_value)
            normalization_squared = dble((two_j_22 + 1)*(two_j_31 + 1))
            twice_phase_exponent = two_j_12 + two_j_22 + two_j_33 + two_j_31
        else if (two_j_22 == 0) then
            call wigner_6j_value (two_j_33, two_j_13, two_j_23, two_j_11, two_j_31, two_j_32, 0, six_j_value)
            normalization_squared = dble((two_j_23 + 1)*(two_j_32 + 1))
            twice_phase_exponent = two_j_13 + two_j_23 + two_j_31 + two_j_32
        else if (two_j_23 == 0) then
            call wigner_6j_value (two_j_31, two_j_11, two_j_21, two_j_12, two_j_32, two_j_33, 0, six_j_value)
            normalization_squared = dble((two_j_21 + 1)*(two_j_33 + 1))
            twice_phase_exponent = two_j_11 + two_j_21 + two_j_32 + two_j_33
        else if (two_j_31 == 0) then
            call wigner_6j_value (two_j_12, two_j_13, two_j_11, two_j_23, two_j_22, two_j_32, 0, six_j_value)
            normalization_squared = dble((two_j_11 + 1)*(two_j_32 + 1))
            twice_phase_exponent = two_j_13 + two_j_11 + two_j_22 + two_j_32
        else if (two_j_32 == 0) then
            call wigner_6j_value (two_j_13, two_j_11, two_j_12, two_j_21, two_j_23, two_j_33, 0, six_j_value)
            normalization_squared = dble((two_j_12 + 1)*(two_j_33 + 1))
            twice_phase_exponent = two_j_11 + two_j_12 + two_j_23 + two_j_33
        else if (two_j_33 == 0) then
            call wigner_6j_value (two_j_11, two_j_12, two_j_13, two_j_22, two_j_21, two_j_31, 0, six_j_value)
            normalization_squared = dble((two_j_13 + 1)*(two_j_31 + 1))
            twice_phase_exponent = two_j_12 + two_j_13 + two_j_21 + two_j_31
        else
            six_j_value = zero
            normalization_squared = one
        end if
        nine_j_value = six_j_value/dsqrt(normalization_squared)
        if (mod(twice_phase_exponent,4) /= 0) nine_j_value = -nine_j_value
        return
    end subroutine wigner_9j_j_zero
end module wigner_9j_with_zero
