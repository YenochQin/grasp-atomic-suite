! Evaluate a Wigner 6j symbol with sixth angular momentum 2.
module wigner_6j_with_j_two
    implicit none (type, external)

    ! Module interface
    private
    public :: wigner_6j_j_two
contains
    subroutine wigner_6j_j_two(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,check_triangle_conditions,six_j_value)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, half, one, two, three, four
        use wigner_6j_triangle_selection, only: allows_wigner_6j_triangles
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: two_j_1, two_j_2, two_j_3, two_j_4, two_j_5
        integer, intent(in)       :: check_triangle_conditions
        real(real64), intent(out) :: six_j_value

        ! Local declarations
        integer      :: triangle_sum
        real(real64) :: real_triangle_sum, two_j_3_real, two_j_1_real, two_j_2_real, phase_sign, first_numerator_term, &
            second_numerator_term, third_numerator_term
        six_j_value = zero
        if (check_triangle_conditions /= 0) then
            if (allows_wigner_6j_triangles(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,4) == 0) return
        end if
        triangle_sum = (two_j_1 + two_j_2 + two_j_3)/2
        real_triangle_sum = dble(triangle_sum)
        two_j_3_real = dble(two_j_3)
        two_j_1_real = dble(two_j_1)
        two_j_2_real = dble(two_j_2)
        phase_sign = one
        if (mod(triangle_sum,2) /= 0) phase_sign = -phase_sign
        if (two_j_1 - two_j_5 == 4) then
            if (two_j_2 - two_j_4 == 4) then
                six_j_value = phase_sign*dsqrt((real_triangle_sum - two)*(real_triangle_sum - one)*real_triangle_sum* &
                    (real_triangle_sum + one)/((two_j_1_real - three)*(two_j_1_real&
                    - two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)))
                six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_3_real - three)*(real_triangle_sum - &
                    two_j_3_real - two)*(real_triangle_sum - two_j_3_real - one)*(real_triangle_sum - &
                    two_j_3_real)/((two_j_2_real - three)*(two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real* &
                    (two_j_2_real + one)))
            else if (two_j_4 - two_j_2 == 4) then
                six_j_value = phase_sign*dsqrt((real_triangle_sum - two_j_2_real - three)*(real_triangle_sum - &
                    two_j_2_real - two)*(real_triangle_sum - two_j_2_real - one)*(real_triangle_sum&
                    - two_j_2_real)/((two_j_2_real + one)*(two_j_2_real + two)*(two_j_2_real + three)*(two_j_2_real + four)* &
                    (two_j_2_real + two + &
                    three)))
                six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_1_real + one)*(real_triangle_sum - &
                    two_j_1_real + two)*(real_triangle_sum - two_j_1_real + three)*(real_triangle_sum - &
                    two_j_1_real + four)/((two_j_1_real - three)*(two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real* &
                    (two_j_1_real + one)))
            else if (two_j_2 - two_j_4 == 2) then
                six_j_value = phase_sign*two*dsqrt((real_triangle_sum - one)*real_triangle_sum*(real_triangle_sum + &
                    one)/((two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real*&
                    (two_j_2_real + one)*(two_j_2_real + two)))
                six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_3_real - two)*(real_triangle_sum - &
                    two_j_3_real - one)*(real_triangle_sum - two_j_3_real)*(real_triangle_sum - two_j_2_real)* &
                    (real_triangle_sum&
                    - two_j_1_real + one)/((two_j_1_real - three)*(two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real* &
                    (two_j_1_real + one)))
            else if (two_j_4 - two_j_2 == 2) then
                six_j_value = phase_sign*two*dsqrt((real_triangle_sum + one)*(real_triangle_sum - two_j_3_real)* &
                    (real_triangle_sum - two_j_2_real - two)*(real_triangle_sum - two_j_2_real - one&
                    )*(real_triangle_sum - two_j_2_real)/(two_j_2_real*(two_j_2_real + one)*(two_j_2_real + two)* &
                    (two_j_2_real + three)*(two_j_2_real + four)))
                six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_1_real + one)*(real_triangle_sum - &
                    two_j_1_real + two)*(real_triangle_sum - two_j_1_real + three)/((two_j_1_real - &
                    three)*(two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)))
            else if (two_j_2 - two_j_4 == 0) then
                six_j_value = phase_sign*dsqrt(two*three*real_triangle_sum*(real_triangle_sum + one)* &
                    (real_triangle_sum - two_j_3_real - one)*(real_triangle_sum - two_j_3_real)/((two_j_2_real&
                    - one)*two_j_2_real*(two_j_2_real + one)*(two_j_2_real + two)*(two_j_2_real + three)))
                six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - &
                    two_j_2_real)*(real_triangle_sum - two_j_1_real + one)*(real_triangle_sum - two_j_1_real + two)&
                    /((two_j_1_real - three)*(two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)))
            end if
        else if (two_j_5 - two_j_1 == 4) then
            if (two_j_2 - two_j_4 == 4) then
                six_j_value = phase_sign*dsqrt((real_triangle_sum - two_j_1_real - three)*(real_triangle_sum - &
                    two_j_1_real - two)*(real_triangle_sum - two_j_1_real - one)*(real_triangle_sum&
                    - two_j_1_real)/((two_j_1_real + one)*(two_j_1_real + two)*(two_j_1_real + three)*(two_j_1_real + four)* &
                    (two_j_1_real + two + &
                    three)))
                six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_2_real + one)*(real_triangle_sum - &
                    two_j_2_real + two)*(real_triangle_sum - two_j_2_real + three)*(real_triangle_sum - &
                    two_j_2_real + four)/((two_j_2_real - three)*(two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real* &
                    (two_j_2_real + one)))
            else if (two_j_4 - two_j_2 == 4) then
                six_j_value = phase_sign*dsqrt((real_triangle_sum - two_j_3_real + four)*(real_triangle_sum - &
                    two_j_3_real + three)*(real_triangle_sum - two_j_3_real + two)*(real_triangle_sum&
                    - two_j_3_real + one)/((two_j_1_real + three + two)*(two_j_1_real + four)*(two_j_1_real + three)* &
                    (two_j_1_real + two)*&
                    (two_j_1_real + one)))
                six_j_value = six_j_value*dsqrt((real_triangle_sum + three + two)*(real_triangle_sum + four)* &
                    (real_triangle_sum + three)*(real_triangle_sum + two&
                    )/((two_j_2_real + three + two)*(two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real + two)* &
                    (two_j_2_real + one))&
                    )
            else if (two_j_2 - two_j_4 == 2) then
                six_j_value = -phase_sign*dsqrt((real_triangle_sum - two_j_3_real + one)*(real_triangle_sum + two)* &
                    (real_triangle_sum - two_j_1_real - two)*(real_triangle_sum - two_j_1_real - &
                    one)*(real_triangle_sum - two_j_1_real)/((two_j_1_real + two + three)*(two_j_1_real + four)* &
                    (two_j_1_real + three)*(two_j_1_real + two&
                    )*(two_j_1_real + one)))
                six_j_value = six_j_value*two*dsqrt((real_triangle_sum - two_j_2_real + three)*(real_triangle_sum - &
                    two_j_2_real + two)*(real_triangle_sum - two_j_2_real + one)/((&
                    two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)*(two_j_2_real + two)))
            else if (two_j_4 - two_j_2 == 2) then
                six_j_value = -phase_sign*dsqrt((real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_2_real + &
                    one)*(real_triangle_sum - two_j_3_real + three)*(real_triangle_sum - two_j_3_real + &
                    two)*(real_triangle_sum - two_j_3_real + one)/((two_j_1_real + three + two)*(two_j_1_real + four)* &
                    (two_j_1_real + three)*(two_j_1_real&
                    + two)*(two_j_1_real + one)))
                six_j_value = six_j_value*two*dsqrt((real_triangle_sum + four)*(real_triangle_sum + three)* &
                    (real_triangle_sum + two)/((two_j_2_real + four)*(&
                    two_j_2_real + three)*(two_j_2_real + two)*(two_j_2_real + one)*two_j_2_real))
            else if (two_j_2 - two_j_4 == 0) then
                six_j_value = phase_sign*dsqrt(three*two*(real_triangle_sum - two_j_1_real)*(real_triangle_sum - &
                    two_j_1_real - one)*(real_triangle_sum - two_j_2_real + two)*(real_triangle_sum&
                    - two_j_2_real + one)/((two_j_1_real + three + two)*(two_j_1_real + four)*(two_j_1_real + three)* &
                    (two_j_1_real + two)*&
                    (two_j_1_real + one)))
                six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_3_real + two)*(real_triangle_sum - &
                    two_j_3_real + one)*(real_triangle_sum + three)*(real_triangle_sum + two)&
                    /((two_j_2_real + three)*(two_j_2_real + two)*(two_j_2_real + one)*two_j_2_real*(two_j_2_real - one)))
            end if
        else if (two_j_1 - two_j_5 == 2) then
            if (two_j_2 - two_j_4 == 4) then
                six_j_value = phase_sign*two*dsqrt((real_triangle_sum - one)*real_triangle_sum*(real_triangle_sum + &
                    one)/((two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real*&
                    (two_j_1_real + one)*(two_j_1_real + two)))
                six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_3_real - two)*(real_triangle_sum - &
                    two_j_3_real - one)*(real_triangle_sum - two_j_3_real)*(real_triangle_sum - two_j_1_real)* &
                    (real_triangle_sum&
                    - two_j_2_real + one)/((two_j_2_real - three)*(two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real* &
                    (two_j_2_real + one)))
            else if (two_j_4 - two_j_2 == 4) then
                six_j_value = -phase_sign*dsqrt((real_triangle_sum - two_j_3_real + one)*(real_triangle_sum + two)* &
                    (real_triangle_sum - two_j_2_real - two)*(real_triangle_sum - two_j_2_real - &
                    one)*(real_triangle_sum - two_j_2_real)/((two_j_2_real + two + three)*(two_j_2_real + four)* &
                    (two_j_2_real + three)*(two_j_2_real + two&
                    )*(two_j_2_real + one)))
                six_j_value = six_j_value*two*dsqrt((real_triangle_sum - two_j_1_real + three)*(real_triangle_sum - &
                    two_j_1_real + two)*(real_triangle_sum - two_j_1_real + one)/((&
                    two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two)))
            else if (two_j_2 - two_j_4 == 2) then
                six_j_value = phase_sign*((two_j_3_real + two_j_1_real)*(two_j_3_real - two_j_1_real + two) - &
                    (two_j_2_real - two)*(two_j_2_real - two_j_1_real + two))/dsqrt((two_j_1_real&
                    - two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two))
                six_j_value = six_j_value*dsqrt(real_triangle_sum*(real_triangle_sum + one)*(real_triangle_sum - &
                    two_j_3_real - one)*(real_triangle_sum - two_j_3_real)/((two_j_2_real - two)*(two_j_2_real&
                    - one)*two_j_2_real*(two_j_2_real + one)*(two_j_2_real + two)))
            else if (two_j_4 - two_j_2 == 2) then
                six_j_value = phase_sign*((two_j_3_real + two_j_2_real + four)*(two_j_3_real - two_j_2_real - two) - &
                    (two_j_1_real - two)*(two_j_1_real + two_j_2_real + four))/&
                    dsqrt(two_j_2_real*(two_j_2_real + one)*(two_j_2_real + two)*(two_j_2_real + three)*(two_j_2_real + four))
                six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - &
                    two_j_2_real)*(real_triangle_sum - two_j_1_real + one)*(real_triangle_sum - two_j_1_real + two)&
                    /((two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two)))
            else if (two_j_2 - two_j_4 == 0) then
                six_j_value = phase_sign*half*((two_j_3_real + two_j_2_real + two)*(two_j_3_real - two_j_2_real) - &
                    two_j_1_real*two_j_1_real + four)/dsqrt((two_j_2_real - one)*&
                    two_j_2_real*(two_j_2_real + one)*(two_j_2_real + two)*(two_j_2_real + three))
                six_j_value = six_j_value*dsqrt(three*two*(real_triangle_sum + one)*(real_triangle_sum - two_j_3_real)* &
                    (real_triangle_sum - two_j_2_real)*(real_triangle_sum - two_j_1_real + one)&
                    /((two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two)))
            end if
        else if (two_j_5 - two_j_1 == 2) then
            if (two_j_2 - two_j_4 == 4) then
                six_j_value = phase_sign*two*dsqrt((real_triangle_sum + one)*(real_triangle_sum - two_j_3_real)* &
                    (real_triangle_sum - two_j_1_real - two)*(real_triangle_sum - two_j_1_real - one&
                    )*(real_triangle_sum - two_j_1_real)/(two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two)* &
                    (two_j_1_real + three)*(two_j_1_real + four)))
                six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_2_real + one)*(real_triangle_sum - &
                    two_j_2_real + two)*(real_triangle_sum - two_j_2_real + three)/((two_j_2_real - &
                    three)*(two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)))
            else if (two_j_4 - two_j_2 == 4) then
                six_j_value = -phase_sign*dsqrt((real_triangle_sum - two_j_2_real)*(real_triangle_sum - two_j_1_real + &
                    one)*(real_triangle_sum - two_j_3_real + three)*(real_triangle_sum - two_j_3_real + &
                    two)*(real_triangle_sum - two_j_3_real + one)/((two_j_2_real + three + two)*(two_j_2_real + four)* &
                    (two_j_2_real + three)*(two_j_2_real&
                    + two)*(two_j_2_real + one)))
                six_j_value = six_j_value*two*dsqrt((real_triangle_sum + four)*(real_triangle_sum + three)* &
                    (real_triangle_sum + two)/((two_j_1_real + four)*(&
                    two_j_1_real + three)*(two_j_1_real + two)*(two_j_1_real + one)*two_j_1_real))
            else if (two_j_2 - two_j_4 == 2) then
                six_j_value = phase_sign*((two_j_3_real + two_j_1_real + four)*(two_j_3_real - two_j_1_real - two) - &
                    (two_j_2_real - two)*(two_j_1_real + two_j_2_real + four))/&
                    dsqrt(two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two)*(two_j_1_real + three)*(two_j_1_real + four))
                six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - &
                    two_j_1_real)*(real_triangle_sum - two_j_2_real + one)*(real_triangle_sum - two_j_2_real + two)&
                    /((two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)*(two_j_2_real + two)))
            else if (two_j_4 - two_j_2 == 2) then
                six_j_value = phase_sign*(three*(real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_2_real) - &
                    (real_triangle_sum - two_j_3_real)*(real_triangle_sum + four))/dsqrt((two_j_1_real&
                    + four)*(two_j_1_real + three)*(two_j_1_real + two)*(two_j_1_real + one)*two_j_1_real)
                six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_3_real + two)*(real_triangle_sum - &
                    two_j_3_real + one)*(real_triangle_sum + three)*(real_triangle_sum + two)&
                    /((two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real + two)*(two_j_2_real + one)*two_j_2_real))
            else if (two_j_2 - two_j_4 == 0) then
                six_j_value = -phase_sign*((real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - two_j_2_real) - &
                    (real_triangle_sum - two_j_3_real)*(real_triangle_sum + three))/dsqrt((&
                    two_j_1_real + four)*(two_j_1_real + three)*(two_j_1_real + two)*(two_j_1_real + one)*two_j_1_real)
                six_j_value = six_j_value*dsqrt(three*two*(real_triangle_sum - two_j_1_real)*(real_triangle_sum - &
                    two_j_2_real + one)*(real_triangle_sum - two_j_3_real + one)*(real_triangle_sum&
                    + two)/((two_j_2_real + three)*(two_j_2_real + two)*(two_j_2_real + one)*two_j_2_real*(two_j_2_real - &
                    one)))
            end if
        else if (two_j_5 - two_j_1 == 0) then
            if (two_j_2 - two_j_4 == 4) then
                six_j_value = phase_sign*dsqrt(three*two*real_triangle_sum*(real_triangle_sum + one)* &
                    (real_triangle_sum - two_j_3_real - one)*(real_triangle_sum - two_j_3_real)/((two_j_1_real&
                    - one)*two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two)*(two_j_1_real + three)))
                six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - &
                    two_j_1_real)*(real_triangle_sum - two_j_2_real + one)*(real_triangle_sum - two_j_2_real + two)&
                    /((two_j_2_real - three)*(two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)))
            else if (two_j_4 - two_j_2 == 4) then
                six_j_value = phase_sign*dsqrt(three*two*(real_triangle_sum - two_j_2_real)*(real_triangle_sum - &
                    two_j_2_real - one)*(real_triangle_sum - two_j_1_real + two)*(real_triangle_sum&
                    - two_j_1_real + one)/((two_j_2_real + three + two)*(two_j_2_real + four)*(two_j_2_real + three)* &
                    (two_j_2_real + two)*&
                    (two_j_2_real + one)))
                six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_3_real + two)*(real_triangle_sum - &
                    two_j_3_real + one)*(real_triangle_sum + three)*(real_triangle_sum + two)&
                    /((two_j_1_real + three)*(two_j_1_real + two)*(two_j_1_real + one)*two_j_1_real*(two_j_1_real - one)))
            else if (two_j_2 - two_j_4 == 2) then
                six_j_value = phase_sign*half*((two_j_3_real + two_j_1_real + two)*(two_j_3_real - two_j_1_real) - &
                    two_j_2_real*two_j_2_real + four)/dsqrt((two_j_1_real - one)*&
                    two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two)*(two_j_1_real + three))
                six_j_value = six_j_value*dsqrt(three*two*(real_triangle_sum + one)*(real_triangle_sum - two_j_3_real)* &
                    (real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_2_real + one)&
                    /((two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)*(two_j_2_real + two)))
            else if (two_j_4 - two_j_2 == 2) then
                six_j_value = -phase_sign*((real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - two_j_1_real) - &
                    (real_triangle_sum - two_j_3_real)*(real_triangle_sum + three))/dsqrt((&
                    two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real + two)*(two_j_2_real + one)*two_j_2_real)
                six_j_value = six_j_value*dsqrt(three*two*(real_triangle_sum - two_j_2_real)*(real_triangle_sum - &
                    two_j_1_real + one)*(real_triangle_sum - two_j_3_real + one)*(real_triangle_sum&
                    + two)/((two_j_1_real + three)*(two_j_1_real + two)*(two_j_1_real + one)*two_j_1_real*(two_j_1_real - &
                    one)))
            else if (two_j_2 - two_j_4 == 0) then
                first_numerator_term = (real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_1_real - one)* &
                    (real_triangle_sum - two_j_2_real)*(real_triangle_sum - two_j_2_real - one)
                second_numerator_term = four*(real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_2_real)* &
                    (real_triangle_sum - two_j_3_real)*(real_triangle_sum + two)
                third_numerator_term = (real_triangle_sum - two_j_3_real)*(real_triangle_sum - two_j_3_real - one)* &
                    (real_triangle_sum + three)*(real_triangle_sum + two)
                six_j_value = phase_sign*(first_numerator_term - second_numerator_term + third_numerator_term)/ &
                    dsqrt((two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two)*(two_j_1_real + &
                    three)*(two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)*(two_j_2_real + two)*(two_j_2_real + three))
            end if
        end if
        return
    end subroutine wigner_6j_j_two
end module wigner_6j_with_j_two
