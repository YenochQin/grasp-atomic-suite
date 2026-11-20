! Test a Wigner 9j symbol's triangles or sum its three-6j products.
module wigner_9j
    implicit none (type, external)

    ! Module interface
    private
    public :: wigner_9j_value
contains
    subroutine wigner_9j_value(two_j_11,two_j_12,two_j_13,two_j_21,two_j_22,two_j_23,two_j_31,two_j_32,two_j_33, &
        & check_triangles_only, &
        evaluation_status,nine_j_value)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero
        use angular_triangle_and_parity, only: satisfies_angular_triangle_and_parity
        use wigner_9j_with_zero, only: wigner_9j_j_zero
        use wigner_6j, only: wigner_6j_value
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: two_j_11, two_j_12, two_j_13, two_j_21, two_j_22, two_j_23, two_j_31, two_j_32, two_j_33
        integer, intent(in)  :: check_triangles_only
        integer, intent(out) :: evaluation_status
        real(real64), intent(inout) :: nine_j_value

        ! Local declarations
        integer :: triangle_bound_1, triangle_bound_2, triangle_bound_3, triangle_bound_4, triangle_bound_5, &
            triangle_bound_6, min_two_intermediate_j, max_two_intermediate_j, two_intermediate_j
        real(real64) :: first_six_j, second_six_j, third_six_j, weighted_six_j_product

        ! Check the three rows and three columns when no value is requested.
        if (check_triangles_only == 1) then
            evaluation_status = 0
            if (satisfies_angular_triangle_and_parity(two_j_11,two_j_12,two_j_13) == 0) return
            if (satisfies_angular_triangle_and_parity(two_j_21,two_j_22,two_j_23) == 0) return
            if (satisfies_angular_triangle_and_parity(two_j_31,two_j_32,two_j_33) == 0) return
            if (satisfies_angular_triangle_and_parity(two_j_11,two_j_21,two_j_31) == 0) return
            if (satisfies_angular_triangle_and_parity(two_j_12,two_j_22,two_j_32) == 0) return
            if (satisfies_angular_triangle_and_parity(two_j_13,two_j_23,two_j_33) == 0) return
            evaluation_status = 1
            return
        end if
        if (two_j_11*two_j_12*two_j_13*two_j_21*two_j_22*two_j_23*two_j_31*two_j_32*two_j_33 == 0) then
            evaluation_status = 1
            call wigner_9j_j_zero (two_j_11, two_j_12, two_j_13, two_j_21, two_j_22, two_j_23, two_j_31, two_j_32, two_j_33, &
                & nine_j_value)
        else
            triangle_bound_1 = iabs(two_j_11 - two_j_33)
            triangle_bound_2 = iabs(two_j_23 - two_j_12)
            triangle_bound_3 = iabs(two_j_21 - two_j_32)
            triangle_bound_4 = iabs(two_j_12 - two_j_23)
            triangle_bound_5 = iabs(two_j_32 - two_j_21)
            triangle_bound_6 = iabs(two_j_11 - two_j_33)

            ! Intersect the three lower triangle bounds for the summed 2x.
            min_two_intermediate_j = max0(triangle_bound_1,triangle_bound_2,triangle_bound_3,triangle_bound_4, &
                triangle_bound_5,triangle_bound_6)
            triangle_bound_1 = two_j_11 + two_j_33
            triangle_bound_2 = two_j_23 + two_j_12
            triangle_bound_3 = two_j_12 + two_j_23
            triangle_bound_4 = two_j_32 + two_j_21
            triangle_bound_5 = two_j_11 + two_j_33
            triangle_bound_6 = two_j_21 + two_j_32
            max_two_intermediate_j = min0(triangle_bound_1,triangle_bound_2,triangle_bound_3,triangle_bound_4, &
                triangle_bound_5,triangle_bound_6)
            evaluation_status = 1
            nine_j_value = zero
            ! The triangle parity fixes 2x modulo two, so admissible x values differ by one.
            do two_intermediate_j = min_two_intermediate_j, max_two_intermediate_j, 2
                call wigner_6j_value (two_j_11, two_j_12, two_j_13, two_j_23, two_j_33, two_intermediate_j, 0, first_six_j)
                call wigner_6j_value (two_j_21, two_j_22, two_j_23, two_j_12, two_intermediate_j, two_j_32, 0, second_six_j)
                call wigner_6j_value (two_j_31, two_j_32, two_j_33, two_intermediate_j, two_j_11, two_j_21, 0, third_six_j)
                weighted_six_j_product = first_six_j*second_six_j*third_six_j*dble(two_intermediate_j + 1)
                if (mod(two_intermediate_j,2) /= 0) weighted_six_j_product = -weighted_six_j_product
                nine_j_value = weighted_six_j_product + nine_j_value
            end do
        end if
        return
    end subroutine wigner_9j_value
end module wigner_9j
