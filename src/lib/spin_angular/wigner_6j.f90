! Evaluate a Wigner 6j symbol using symmetry, special cases and a small cache.
module wigner_6j
    implicit none (type, external)

    ! Module interface
    private
    public :: wigner_6j_value
contains
    subroutine wigner_6j_value(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6,check_triangle_conditions,six_j_value)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, one
        use wigner_6j_triangle_selection, only: allows_wigner_6j_triangles
        use wigner_6j_with_j_one_half, only: wigner_6j_j_one_half
        use wigner_6j_with_j_one, only: wigner_6j_j_one
        use wigner_6j_with_j_three_halves, only: wigner_6j_j_three_halves
        use wigner_6j_with_j_two, only: wigner_6j_j_two
        use racah_coefficients, only: evaluate_racah_coefficient
        use wigner_6j_with_j_three, only: wigner_6j_j_three
        use wigner_6j_with_j_four, only: wigner_6j_j_four
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, two_j_6
        integer, intent(in) :: check_triangle_conditions
        real(real64), intent(out) :: six_j_value

        ! Local declarations
        integer :: twice_phase_exponent
        real(real64) :: small_six_j_cache(0:4, 0:4, 0:4, 0:4, 0:4, 0:4)
        real(real64) :: uncached_six_j_sentinel, normalization_squared
        logical :: store_in_six_j_cache
        save :: small_six_j_cache

        ! Initial values
        data small_six_j_cache/ 15625*1.0e-20_real64/
        data uncached_six_j_sentinel/ 1.0e-20_real64/
!$OMP THREADPRIVATE(small_six_j_cache)
        six_j_value = zero
        ! Callers that have already established the triangle conditions can skip the check.
        if (check_triangle_conditions /= 0) then
            if (allows_wigner_6j_triangles(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 0) return
        end if
        store_in_six_j_cache = .false.

        ! Cache symbols whose six doubled arguments are all between 0 and 4.
        if (max0(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) <= 4) then
            six_j_value = small_six_j_cache(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6)
            if (six_j_value == uncached_six_j_sentinel) then
                store_in_six_j_cache = .true.
            else
                return
            end if
        end if

        ! A zero angular momentum reduces the symbol to a phase and a normalization.
        if (two_j_1*two_j_2*two_j_3*two_j_4*two_j_5*two_j_6 == 0) then
            if (two_j_1 == 0) then
                normalization_squared = dble((two_j_5 + 1)*(two_j_3 + 1))
                twice_phase_exponent = two_j_4 + two_j_5 + two_j_3
            else if (two_j_2 == 0) then
                normalization_squared = dble((two_j_4 + 1)*(two_j_3 + 1))
                twice_phase_exponent = two_j_1 + two_j_5 + two_j_6
            else if (two_j_3 == 0) then
                normalization_squared = dble((two_j_1 + 1)*(two_j_4 + 1))
                twice_phase_exponent = two_j_1 + two_j_5 + two_j_6
            else if (two_j_4 == 0) then
                normalization_squared = dble((two_j_2 + 1)*(two_j_3 + 1))
                twice_phase_exponent = two_j_1 + two_j_2 + two_j_3
            else if (two_j_5 == 0) then
                normalization_squared = dble((two_j_1 + 1)*(two_j_3 + 1))
                twice_phase_exponent = two_j_1 + two_j_2 + two_j_3
            else
                normalization_squared = dble((two_j_1 + 1)*(two_j_2 + 1))
                twice_phase_exponent = two_j_1 + two_j_2 + two_j_3
            end if
            six_j_value = one/dsqrt(normalization_squared)
            ! The sign is (-1)**phase. The stored exponent is doubled, hence the modulus of four.
            if (mod(twice_phase_exponent,4) /= 0) six_j_value = -six_j_value

        ! Use 6j symmetries to move the small argument to the sixth slot.
        else if (min0(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 1) then
            if (two_j_1 == 1) then
                call wigner_6j_j_one_half (two_j_5, two_j_3, two_j_4, two_j_2, two_j_6, 0, six_j_value)
            else if (two_j_2 == 1) then
                call wigner_6j_j_one_half (two_j_1, two_j_6, two_j_5, two_j_4, two_j_3, 0, six_j_value)
            else if (two_j_3 == 1) then
                call wigner_6j_j_one_half (two_j_1, two_j_5, two_j_6, two_j_4, two_j_2, 0, six_j_value)
            else if (two_j_4 == 1) then
                call wigner_6j_j_one_half (two_j_2, two_j_3, two_j_1, two_j_5, two_j_6, 0, six_j_value)
            else if (two_j_5 == 1) then
                call wigner_6j_j_one_half (two_j_1, two_j_3, two_j_2, two_j_4, two_j_6, 0, six_j_value)
            else
                call wigner_6j_j_one_half (two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, 0, six_j_value)
            end if
        else if (min0(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 2) then
            if (two_j_1 == 2) then
                call wigner_6j_j_one (two_j_5, two_j_3, two_j_4, two_j_2, two_j_6, 0, six_j_value)
            else if (two_j_2 == 2) then
                call wigner_6j_j_one (two_j_1, two_j_6, two_j_5, two_j_4, two_j_3, 0, six_j_value)
            else if (two_j_3 == 2) then
                call wigner_6j_j_one (two_j_1, two_j_5, two_j_6, two_j_4, two_j_2, 0, six_j_value)
            else if (two_j_4 == 2) then
                call wigner_6j_j_one (two_j_2, two_j_3, two_j_1, two_j_5, two_j_6, 0, six_j_value)
            else if (two_j_5 == 2) then
                call wigner_6j_j_one (two_j_1, two_j_3, two_j_2, two_j_4, two_j_6, 0, six_j_value)
            else
                call wigner_6j_j_one (two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, 0, six_j_value)
            end if
        else if (min0(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 3) then
            if (two_j_1 == 3) then
                call wigner_6j_j_three_halves (two_j_5, two_j_3, two_j_4, two_j_2, two_j_6, 0, six_j_value)
            else if (two_j_2 == 3) then
                call wigner_6j_j_three_halves (two_j_1, two_j_6, two_j_5, two_j_4, two_j_3, 0, six_j_value)
            else if (two_j_3 == 3) then
                call wigner_6j_j_three_halves (two_j_1, two_j_5, two_j_6, two_j_4, two_j_2, 0, six_j_value)
            else if (two_j_4 == 3) then
                call wigner_6j_j_three_halves (two_j_2, two_j_3, two_j_1, two_j_5, two_j_6, 0, six_j_value)
            else if (two_j_5 == 3) then
                call wigner_6j_j_three_halves (two_j_1, two_j_3, two_j_2, two_j_4, two_j_6, 0, six_j_value)
            else
                call wigner_6j_j_three_halves (two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, 0, six_j_value)
            end if
        else if (min0(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 4) then
            if (two_j_1 == 4) then
                call wigner_6j_j_two (two_j_5, two_j_3, two_j_4, two_j_2, two_j_6, 0, six_j_value)
            else if (two_j_2 == 4) then
                call wigner_6j_j_two (two_j_1, two_j_6, two_j_5, two_j_4, two_j_3, 0, six_j_value)
            else if (two_j_3 == 4) then
                call wigner_6j_j_two (two_j_1, two_j_5, two_j_6, two_j_4, two_j_2, 0, six_j_value)
            else if (two_j_4 == 4) then
                call wigner_6j_j_two (two_j_2, two_j_3, two_j_1, two_j_5, two_j_6, 0, six_j_value)
            else if (two_j_5 == 4) then
                call wigner_6j_j_two (two_j_1, two_j_3, two_j_2, two_j_4, two_j_6, 0, six_j_value)
            else
                call wigner_6j_j_two (two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, 0, six_j_value)
            end if
        else if (min0(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 5) then
            call evaluate_racah_coefficient (two_j_1, two_j_2, two_j_5, two_j_4, two_j_3, two_j_6, six_j_value)
            if (mod(two_j_1 + two_j_2 + two_j_5 + two_j_4,4) /= 0) six_j_value = -six_j_value
        else if (min0(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 6) then
            if (two_j_1 == 6) then
                call wigner_6j_j_three (two_j_5, two_j_3, two_j_4, two_j_2, two_j_6, 0, six_j_value)
            else if (two_j_2 == 6) then
                call wigner_6j_j_three (two_j_1, two_j_6, two_j_5, two_j_4, two_j_3, 0, six_j_value)
            else if (two_j_3 == 6) then
                call wigner_6j_j_three (two_j_1, two_j_5, two_j_6, two_j_4, two_j_2, 0, six_j_value)
            else if (two_j_4 == 6) then
                call wigner_6j_j_three (two_j_2, two_j_3, two_j_1, two_j_5, two_j_6, 0, six_j_value)
            else if (two_j_5 == 6) then
                call wigner_6j_j_three (two_j_1, two_j_3, two_j_2, two_j_4, two_j_6, 0, six_j_value)
            else
                call wigner_6j_j_three (two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, 0, six_j_value)
            end if
        else if (min0(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 7) then
            call evaluate_racah_coefficient (two_j_1, two_j_2, two_j_5, two_j_4, two_j_3, two_j_6, six_j_value)
            if (mod(two_j_1 + two_j_2 + two_j_5 + two_j_4,4) /= 0) six_j_value = -six_j_value
        else if (min0(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 8) then
            if (two_j_1 == 8) then
                call wigner_6j_j_four (two_j_5, two_j_3, two_j_4, two_j_2, two_j_6, 0, six_j_value)
            else if (two_j_2 == 8) then
                call wigner_6j_j_four (two_j_1, two_j_6, two_j_5, two_j_4, two_j_3, 0, six_j_value)
            else if (two_j_3 == 8) then
                call wigner_6j_j_four (two_j_1, two_j_5, two_j_6, two_j_4, two_j_2, 0, six_j_value)
            else if (two_j_4 == 8) then
                call wigner_6j_j_four (two_j_2, two_j_3, two_j_1, two_j_5, two_j_6, 0, six_j_value)
            else if (two_j_5 == 8) then
                call wigner_6j_j_four (two_j_1, two_j_3, two_j_2, two_j_4, two_j_6, 0, six_j_value)
            else
                call wigner_6j_j_four (two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, 0, six_j_value)
            end if
        else if (min0(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 9) then
            call evaluate_racah_coefficient (two_j_1, two_j_2, two_j_5, two_j_4, two_j_3, two_j_6, six_j_value)
            if (mod(two_j_1 + two_j_2 + two_j_5 + two_j_4,4) /= 0) six_j_value = -six_j_value
        else
            call evaluate_racah_coefficient(two_j_1,two_j_2,two_j_5,two_j_4,two_j_3,two_j_6,six_j_value)
            if (mod(two_j_1 + two_j_2 + two_j_5 + two_j_4,4) /= 0) six_j_value = -six_j_value
        end if
        if (store_in_six_j_cache) small_six_j_cache(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) = six_j_value
        return
    end subroutine wigner_6j_value
end module wigner_6j
