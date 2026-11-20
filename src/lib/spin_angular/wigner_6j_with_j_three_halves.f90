! Evaluate a Wigner 6j symbol with sixth angular momentum three halves.
module wigner_6j_with_j_three_halves
    implicit none (type, external)

    ! Module interface
    private
    public :: sixj35
contains
    subroutine sixj35(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,check_triangle_conditions,six_j_value)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, one, two, three, four
        use wigner_6j_triangle_selection, only: ixjtik
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: two_j_1, two_j_2, two_j_3, two_j_4, two_j_5
        integer, intent(in)       :: check_triangle_conditions
        real(real64), intent(out) :: six_j_value

        ! Local declarations
        integer :: triangle_sum
        real(real64) :: real_triangle_sum, two_j_3_real, two_j_1_real, two_j_2_real, phase_sign
        six_j_value = zero
        if (check_triangle_conditions /= 0) then
            if (ixjtik(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,3) == 0) return
        end if
        triangle_sum = (two_j_1 + two_j_2 + two_j_3)/2
        real_triangle_sum = dble(triangle_sum)
        two_j_3_real = dble(two_j_3)
        two_j_1_real = dble(two_j_1)
        two_j_2_real = dble(two_j_2)
        phase_sign = one
        if (mod(triangle_sum,2) /= 0) phase_sign = -phase_sign
        if (two_j_1 - two_j_5 == 3) then
            if (two_j_2 - two_j_4 == 3) then
                six_j_value = phase_sign*dsqrt((real_triangle_sum - one)*real_triangle_sum*(real_triangle_sum + one)* &
                    (real_triangle_sum - two_j_3_real - two)*(real_triangle_sum - two_j_3_real - &
                    one)*(real_triangle_sum - two_j_3_real)/((two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real* &
                    (two_j_1_real + one)*(two_j_2_real - two)*(two_j_2_real - &
                    one)*two_j_2_real*(two_j_2_real + one)))
            else if (two_j_4 - two_j_2 == 3) then
                six_j_value = phase_sign*dsqrt((real_triangle_sum - two_j_2_real - two)*(real_triangle_sum - &
                    two_j_2_real - one)*(real_triangle_sum - two_j_2_real)*(real_triangle_sum - two_j_1_real + one&
                    )*(real_triangle_sum - two_j_1_real + two)*(real_triangle_sum - two_j_1_real + three)/((two_j_2_real + &
                    1)*(two_j_2_real + two)*(two_j_2_real + three)&
                    *(two_j_2_real + four)*(two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)))
            else if (two_j_2 - two_j_4 == 1) then
                six_j_value = phase_sign*dsqrt(three*real_triangle_sum*(real_triangle_sum + one)*(real_triangle_sum - &
                    two_j_3_real - one)*(real_triangle_sum - two_j_3_real)*(real_triangle_sum - two_j_2_real)&
                    *(real_triangle_sum - two_j_1_real + one)/((two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)* &
                    (two_j_2_real + two)*(two_j_1_real - two)*(two_j_1_real - &
                    one)*two_j_1_real*(two_j_1_real + one)))
            else if (two_j_4 - two_j_2 == 1) then
                six_j_value = phase_sign*dsqrt(three*(real_triangle_sum + one)*(real_triangle_sum - two_j_3_real)* &
                    (real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - two_j_2_real)*(&
                    real_triangle_sum - two_j_1_real + one)*(real_triangle_sum - two_j_1_real + two)/(two_j_2_real* &
                    (two_j_2_real + one)*(two_j_2_real + two)*(two_j_2_real + three)*&
                    (two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)))
            end if
        else if (two_j_5 - two_j_1 == 3) then
            if (two_j_2 - two_j_4 == 3) then
                six_j_value = phase_sign*dsqrt((real_triangle_sum - two_j_1_real - two)*(real_triangle_sum - &
                    two_j_1_real - one)*(real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_2_real + one&
                    )*(real_triangle_sum - two_j_2_real + two)*(real_triangle_sum - two_j_2_real + three)/((two_j_1_real + &
                    one)*(two_j_1_real + two)*(two_j_1_real + &
                    three)*(two_j_1_real + four)*(two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)))
            else if (two_j_4 - two_j_2 == 3) then
                six_j_value = -phase_sign*dsqrt((real_triangle_sum + two)*(real_triangle_sum + three)* &
                    (real_triangle_sum + four)*(real_triangle_sum - two_j_3_real + one)*&
                    (real_triangle_sum - two_j_3_real + two)*(real_triangle_sum - two_j_3_real + three)/((two_j_1_real + &
                    one)*(two_j_1_real + two)*(two_j_1_real + three)&
                    *(two_j_1_real + four)*(two_j_2_real + one)*(two_j_2_real + two)*(two_j_2_real + three)*(two_j_2_real + &
                    four)))
            else if (two_j_2 - two_j_4 == 1) then
                six_j_value = -phase_sign*dsqrt(three*(real_triangle_sum + two)*(real_triangle_sum - two_j_3_real + one) &
                    *(real_triangle_sum - two_j_2_real + one)*(real_triangle_sum&
                    - two_j_2_real + two)*(real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - two_j_1_real)/ &
                    ((two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)*(two_j_2_real + &
                    two)*(two_j_1_real + one)*(two_j_1_real + two)*(two_j_1_real + three)*(two_j_1_real + four)))
            else if (two_j_4 - two_j_2 == 1) then
                six_j_value = phase_sign*dsqrt(three*(real_triangle_sum + two)*(real_triangle_sum + three)* &
                    (real_triangle_sum - two_j_3_real + one)*(real_triangle_sum - two_j_3_real&
                    + two)*(real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_2_real + one)/(two_j_2_real* &
                    (two_j_2_real + one)*(two_j_2_real + two)*(two_j_2_real + &
                    three)*(two_j_1_real + one)*(two_j_1_real + two)*(two_j_1_real + three)*(two_j_1_real + four)))
            end if
        else if (two_j_1 - two_j_5 == 1) then
            if (two_j_2 - two_j_4 == 3) then
                six_j_value = phase_sign*dsqrt((three*real_triangle_sum*(real_triangle_sum + one)*(real_triangle_sum - &
                    two_j_3_real - one)*(real_triangle_sum - two_j_3_real)*(real_triangle_sum - two_j_1_real&
                    )*(real_triangle_sum - two_j_2_real + one))/((two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)* &
                    (two_j_1_real + two)*(two_j_2_real - two)*(two_j_2_real&
                    - one)*two_j_2_real*(two_j_2_real + one)))
            else if (two_j_4 - two_j_2 == 3) then
                six_j_value = -phase_sign*dsqrt(three*(real_triangle_sum + two)*(real_triangle_sum - two_j_3_real + one) &
                    *(real_triangle_sum - two_j_1_real + one)*(real_triangle_sum&
                    - two_j_1_real + two)*(real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - two_j_2_real)/ &
                    ((two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)*(two_j_1_real + &
                    two)*(two_j_2_real + one)*(two_j_2_real + two)*(two_j_2_real + three)*(two_j_2_real + four)))
            else if (two_j_2 - two_j_4 == 1) then
                six_j_value = phase_sign*(two*(real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_2_real) - &
                    (real_triangle_sum + two)*(real_triangle_sum - two_j_3_real - one))*dsqrt(&
                    (real_triangle_sum + one)*(real_triangle_sum - two_j_3_real)/((two_j_1_real - one)*two_j_1_real* &
                    (two_j_1_real + one)*(two_j_1_real + two)*(two_j_2_real - one)*two_j_2_real&
                    *(two_j_2_real + one)*(two_j_2_real + two)))
            else if (two_j_4 - two_j_2 == 1) then
                six_j_value = phase_sign*((real_triangle_sum - two_j_1_real + two)*(real_triangle_sum - two_j_2_real + &
                    one) - two*(real_triangle_sum - two_j_3_real + one)*(real_triangle_sum + &
                    one))*dsqrt((real_triangle_sum - two_j_2_real)*(real_triangle_sum - two_j_1_real + one)/(two_j_2_real* &
                    (two_j_2_real + one)*(two_j_2_real + two)*(two_j_2_real + &
                    three)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two)))
            end if
        else if (two_j_5 - two_j_1 == 1) then
            if (two_j_2 - two_j_4 == 3) then
                six_j_value = phase_sign*dsqrt(three*(real_triangle_sum + one)*(real_triangle_sum - two_j_3_real)* &
                    (real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - two_j_1_real)*(&
                    real_triangle_sum - two_j_2_real + one)*(real_triangle_sum - two_j_2_real + two)/(two_j_1_real* &
                    (two_j_1_real + one)*(two_j_1_real + two)*(two_j_1_real + three)*&
                    (two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)))
            else if (two_j_4 - two_j_2 == 3) then
                six_j_value = phase_sign*dsqrt(three*(real_triangle_sum + two)*(real_triangle_sum + three)* &
                    (real_triangle_sum - two_j_3_real + one)*(real_triangle_sum - two_j_3_real&
                    + two)*(real_triangle_sum - two_j_1_real + one)*(real_triangle_sum - two_j_2_real)/(two_j_1_real* &
                    (two_j_1_real + one)*(two_j_1_real + two)*(two_j_1_real + &
                    three)*(two_j_2_real + one)*(two_j_2_real + two)*(two_j_2_real + three)*(two_j_2_real + four)))
            else if (two_j_2 - two_j_4 == 1) then
                six_j_value = phase_sign*((real_triangle_sum - two_j_2_real + two)*(real_triangle_sum - two_j_1_real + &
                    one) - two*(real_triangle_sum - two_j_3_real + one)*(real_triangle_sum + &
                    one))*dsqrt((real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_2_real + one)/(two_j_1_real* &
                    (two_j_1_real + one)*(two_j_1_real + two)*(two_j_1_real + &
                    three)*(two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)*(two_j_2_real + two)))
            else if (two_j_4 - two_j_2 == 1) then
                six_j_value = -phase_sign*(two*(real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_2_real) - &
                    (real_triangle_sum + three)*(real_triangle_sum - two_j_3_real))*dsqrt((real_triangle_sum&
                    + two)*(real_triangle_sum - two_j_3_real + one)/(two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two) &
                    *(two_j_1_real + three)*two_j_2_real*(two_j_2_real&
                    + one)*(two_j_2_real + two)*(two_j_2_real + three)))
            end if
        end if
        return
    end subroutine sixj35
end module wigner_6j_with_j_three_halves
