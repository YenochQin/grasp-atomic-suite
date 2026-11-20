! Evaluate a Wigner 6j symbol with sixth angular momentum 3.
module wigner_6j_with_j_three
    implicit none (type, external)
    private
    public :: sixj3
contains
      subroutine sixj3(two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, check_triangle_conditions, six_j_value)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants,          only: zero, one, two, three, four, seven
      use wigner_6j_triangle_selection, only: ixjtik
      implicit none (type, external)
      integer, intent(in) :: two_j_1, two_j_2, two_j_3, two_j_4, two_j_5
      integer, intent(in) :: check_triangle_conditions
      real(real64), intent(out) :: six_j_value
      integer      :: triangle_sum
      real(real64) :: real_triangle_sum, two_j_3_real, two_j_1_real, two_j_2_real, phase_sign
      six_j_value = zero
      if (check_triangle_conditions /= 0) then
         if (ixjtik(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,6) == 0) return
      end if
      triangle_sum = (two_j_1 + two_j_2 + two_j_3)/2
      real_triangle_sum = dble(triangle_sum)
      two_j_3_real = dble(two_j_3)
      two_j_1_real = dble(two_j_1)
      two_j_2_real = dble(two_j_2)
      phase_sign = one
      if (mod(triangle_sum,2) /= 0) phase_sign = -phase_sign
      if (two_j_1 - two_j_5 == 6) then
         if (two_j_2 - two_j_4 == 6) then
            six_j_value = phase_sign*dsqrt((real_triangle_sum + one)*real_triangle_sum*(real_triangle_sum - one)* &
                (real_triangle_sum - two)*(real_triangle_sum - three)*(real_triangle_sum&
                - four)/((two_j_1_real + one)*two_j_1_real*(two_j_1_real - one)*(two_j_1_real - two)*(two_j_1_real - &
                    three)*(two_j_1_real - four)&
               *(two_j_1_real - three - two)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_3_real - three - two)*(real_triangle_sum - &
                two_j_3_real - four)*(real_triangle_sum - two_j_3_real - &
               three)*(real_triangle_sum - two_j_3_real - two)*(real_triangle_sum - two_j_3_real - one)* &
                   (real_triangle_sum - two_j_3_real)/((two_j_2_real + one)*two_j_2_real*(two_j_2_real - &
               one)*(two_j_2_real - two)*(two_j_2_real - three)*(two_j_2_real - four)*(two_j_2_real - three - two)))
         else if (two_j_4 - two_j_2 == 6) then
            six_j_value = phase_sign*dsqrt((real_triangle_sum - two_j_2_real)*(real_triangle_sum - two_j_2_real - &
                one)*(real_triangle_sum - two_j_2_real - two)*(real_triangle_sum - two_j_2_real - &
               three)*(real_triangle_sum - two_j_2_real - four)*(real_triangle_sum - two_j_2_real - three - two)/ &
                   ((two_j_1_real + one)*two_j_1_real*(two_j_1_real - &
               one)*(two_j_1_real - two)*(two_j_1_real - three)*(two_j_1_real - four)*(two_j_1_real - three - two)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_1_real + three + three)*(real_triangle_sum - &
                two_j_1_real + three + two)*(real_triangle_sum&
                - two_j_1_real + four)*(real_triangle_sum - two_j_1_real + three)*(real_triangle_sum - two_j_1_real + &
                    two)*(real_triangle_sum - two_j_1_real + one)/((two_j_2_real&
                + four + three)*(two_j_2_real + three + three)*(two_j_2_real + three + two)*(two_j_2_real + four&
               )*(two_j_2_real + three)*(two_j_2_real + two)*(two_j_2_real + one)))
         else if (two_j_2 - two_j_4 == 4) then
            six_j_value = phase_sign*dsqrt(two*three*(real_triangle_sum - two_j_2_real)*(real_triangle_sum - &
                two_j_1_real + one)*(real_triangle_sum - three)*(real_triangle_sum&
                - two)*(real_triangle_sum - one)*real_triangle_sum*(real_triangle_sum + one)/((two_j_1_real - &
                    three - two)*(two_j_1_real - four)*(&
               two_j_1_real - three)*(two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_3_real - four)*(real_triangle_sum - &
                two_j_3_real - three)*(real_triangle_sum - two_j_3_real - two)*(real_triangle_sum&
                - two_j_3_real - one)*(real_triangle_sum - two_j_3_real)/((two_j_2_real - four)*(two_j_2_real - three)* &
                    (two_j_2_real - two)*(two_j_2_real - one)&
               *two_j_2_real*(two_j_2_real + one)*(two_j_2_real + two)))
         else if (two_j_4 - two_j_2 == 4) then
            six_j_value = phase_sign*dsqrt(two*three*(real_triangle_sum + one)*(real_triangle_sum - two_j_3_real)* &
                (real_triangle_sum - two_j_2_real - four)*(real_triangle_sum - &
               two_j_2_real - three)*(real_triangle_sum - two_j_2_real - two)*(real_triangle_sum - two_j_2_real - one) &
                   *(real_triangle_sum - two_j_2_real)/((two_j_2_real + three + &
               three)*(two_j_2_real + three + two)*(two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real + two)* &
                   (two_j_2_real + &
               one)*two_j_2_real))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_1_real + three + two)*(real_triangle_sum - &
                two_j_1_real + four)*(real_triangle_sum - two_j_1_real + &
               three)*(real_triangle_sum - two_j_1_real + two)*(real_triangle_sum - two_j_1_real + one)/ &
                   ((two_j_1_real - three - two)*(two_j_1_real - &
               four)*(two_j_1_real - three)*(two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)))
         else if (two_j_2 - two_j_4 == 2) then
            six_j_value = phase_sign*dsqrt(three*(three + two)*(real_triangle_sum - two_j_2_real)* &
                (real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - two_j_1_real&
                + two)*(real_triangle_sum - two_j_1_real + one)*(real_triangle_sum - two)*(real_triangle_sum - one)* &
                    real_triangle_sum*(real_triangle_sum + one)/((two_j_1_real&
                - three - two)*(two_j_1_real - four)*(two_j_1_real - three)*(two_j_1_real - two)*(two_j_1_real - one)* &
                    two_j_1_real*(two_j_1_real&
                + one)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_3_real - three)*(real_triangle_sum - &
                two_j_3_real - two)*(real_triangle_sum - two_j_3_real - one)*(real_triangle_sum - &
               two_j_3_real)/((two_j_2_real - three)*(two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real* &
                   (two_j_2_real + one)*(two_j_2_real + two)*(two_j_2_real + &
               three)))
         else if (two_j_4 - two_j_2 == 2) then
            six_j_value = phase_sign*dsqrt(three*(three + two)*(real_triangle_sum + one)*real_triangle_sum* &
                (real_triangle_sum - two_j_3_real - one)*(real_triangle_sum&
                - two_j_3_real)*(real_triangle_sum - two_j_2_real - three)*(real_triangle_sum - two_j_2_real - two)* &
                    (real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - two_j_2_real)/((&
               two_j_1_real - three - two)*(two_j_1_real - four)*(two_j_1_real - three)*(two_j_1_real - two)* &
                   (two_j_1_real - one)*two_j_1_real*(two_j_1_real&
                + one)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_1_real + four)*(real_triangle_sum - &
                two_j_1_real + three)*(real_triangle_sum - two_j_1_real + two)*(real_triangle_sum&
                - two_j_1_real + one)/((two_j_2_real + three + two)*(two_j_2_real + four)*(two_j_2_real + three)* &
                    (two_j_2_real + two)*&
               (two_j_2_real + one)*two_j_2_real*(two_j_2_real - one)))
         else if (two_j_4 - two_j_2 == 0) then
            six_j_value = phase_sign*two*dsqrt((three + two)*(real_triangle_sum + one)*real_triangle_sum* &
                (real_triangle_sum - one)*(real_triangle_sum - two_j_3_real&
                - two)*(real_triangle_sum - two_j_3_real - one)*(real_triangle_sum - two_j_3_real)/((two_j_1_real - &
                    three - two)*(two_j_1_real - four)*(two_j_1_real&
                - three)*(two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_2_real - two)*(real_triangle_sum - &
                two_j_2_real - one)*(real_triangle_sum - two_j_2_real)*(real_triangle_sum - two_j_1_real + &
               three)*(real_triangle_sum - two_j_1_real + two)*(real_triangle_sum - two_j_1_real + one)/ &
                   ((two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real&
                + two)*(two_j_2_real + one)*two_j_2_real*(two_j_2_real - one)*(two_j_2_real - two)))
         end if
      else if (two_j_5 - two_j_1 == 6) then
         if (two_j_2 - two_j_4 == 6) then
            six_j_value = phase_sign*dsqrt((real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_1_real - &
                one)*(real_triangle_sum - two_j_1_real - two)*(real_triangle_sum - two_j_1_real - &
               three)*(real_triangle_sum - two_j_1_real - four)*(real_triangle_sum - two_j_1_real - three - two)/ &
                   ((two_j_2_real + one)*two_j_2_real*(two_j_2_real - &
               one)*(two_j_2_real - two)*(two_j_2_real - three)*(two_j_2_real - four)*(two_j_2_real - three - two)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_2_real + three + three)*(real_triangle_sum - &
                two_j_2_real + three + two)*(real_triangle_sum&
                - two_j_2_real + four)*(real_triangle_sum - two_j_2_real + three)*(real_triangle_sum - two_j_2_real + &
                    two)*(real_triangle_sum - two_j_2_real + one)/((two_j_1_real&
                + four + three)*(two_j_1_real + three + three)*(two_j_1_real + three + two)*(two_j_1_real + four&
               )*(two_j_1_real + three)*(two_j_1_real + two)*(two_j_1_real + one)))
         else if (two_j_4 - two_j_2 == 6) then
            six_j_value = phase_sign*dsqrt((real_triangle_sum - two_j_3_real + three + three)*(real_triangle_sum - &
                two_j_3_real + three + two)*(real_triangle_sum&
                - two_j_3_real + four)*(real_triangle_sum - two_j_3_real + three)*(real_triangle_sum - two_j_3_real + &
                    two)*(real_triangle_sum - two_j_3_real + one)/((two_j_1_real&
                + seven)*(two_j_1_real + three + three)*(two_j_1_real + three + two)*(two_j_1_real + four)* &
                    (two_j_1_real + &
               three)*(two_j_1_real + two)*(two_j_1_real + one)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum + seven)*(real_triangle_sum + three + three)* &
                (real_triangle_sum + three + two)*&
               (real_triangle_sum + four)*(real_triangle_sum + three)*(real_triangle_sum + two)/((two_j_2_real + &
                   seven)*(two_j_2_real + three + &
               three)*(two_j_2_real + three + two)*(two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real + two)* &
                   (two_j_2_real + &
               one)))
         else if (two_j_2 - two_j_4 == 4) then
            six_j_value = -phase_sign*dsqrt(two*three*(real_triangle_sum - two_j_3_real + one)*(real_triangle_sum + &
                two)*(real_triangle_sum - two_j_1_real - four)&
               *(real_triangle_sum - two_j_1_real - three)*(real_triangle_sum - two_j_1_real - two)* &
                   (real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - two_j_1_real)/((two_j_1_real + &
               seven)*(two_j_1_real + three + three)*(two_j_1_real + three + two)*(two_j_1_real + four)*(two_j_1_real + &
               three)*(two_j_1_real + two)*(two_j_1_real + one)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_2_real + three + two)*(real_triangle_sum - &
                two_j_2_real + four)*(real_triangle_sum - two_j_2_real + &
               three)*(real_triangle_sum - two_j_2_real + two)*(real_triangle_sum - two_j_2_real + one)/ &
                   ((two_j_2_real - four)*(two_j_2_real - three)*(two_j_2_real&
                - two)*(two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)*(two_j_2_real + two)))
         else if (two_j_4 - two_j_2 == 4) then
            six_j_value = -phase_sign*dsqrt(two*three*(real_triangle_sum - two_j_1_real)*(real_triangle_sum - &
                two_j_2_real + one)*(real_triangle_sum - two_j_3_real + three&
                + two)*(real_triangle_sum - two_j_3_real + four)*(real_triangle_sum - two_j_3_real + three)* &
                    (real_triangle_sum - two_j_3_real + two)*(real_triangle_sum - two_j_3_real&
                + one)/((two_j_1_real + seven)*(two_j_1_real + three + three)*(two_j_1_real + three + two)* &
                    (two_j_1_real + &
               four)*(two_j_1_real + three)*(two_j_1_real + two)*(two_j_1_real + one)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum + three + three)*(real_triangle_sum + three + two)* &
                (real_triangle_sum + four)*(&
               real_triangle_sum + three)*(real_triangle_sum + two)/((two_j_2_real + three + three)*(two_j_2_real + &
                   three + two)*(two_j_2_real&
                + four)*(two_j_2_real + three)*(two_j_2_real + two)*(two_j_2_real + one)*two_j_2_real))
         else if (two_j_2 - two_j_4 == 2) then
            six_j_value = phase_sign*dsqrt(three*(three + two)*(real_triangle_sum - two_j_3_real + one)* &
                (real_triangle_sum - two_j_3_real + two)*(&
               real_triangle_sum + three)*(real_triangle_sum + two)*(real_triangle_sum - two_j_1_real - three)* &
                   (real_triangle_sum - two_j_1_real - two)*(real_triangle_sum - two_j_1_real&
                - one)*(real_triangle_sum - two_j_1_real)/((two_j_1_real + seven)*(two_j_1_real + three + three)* &
                    (two_j_1_real + three + &
               two)*(two_j_1_real + four)*(two_j_1_real + three)*(two_j_1_real + two)*(two_j_1_real + one)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_2_real + four)*(real_triangle_sum - &
                two_j_2_real + three)*(real_triangle_sum - two_j_2_real + two)*(real_triangle_sum&
                - two_j_2_real + one)/((two_j_2_real - three)*(two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real* &
                    (two_j_2_real + one)*(two_j_2_real + &
               two)*(two_j_2_real + three)))
         else if (two_j_4 - two_j_2 == 2) then
            six_j_value = phase_sign*dsqrt(three*(three + two)*(real_triangle_sum - two_j_1_real)* &
                (real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - two_j_2_real&
                + two)*(real_triangle_sum - two_j_2_real + one)*(real_triangle_sum - two_j_3_real + four)* &
                    (real_triangle_sum - two_j_3_real + three)*(real_triangle_sum - two_j_3_real&
                + two)*(real_triangle_sum - two_j_3_real + one)/((two_j_1_real + seven)*(two_j_1_real + three + three) &
                    *(two_j_1_real + &
               three + two)*(two_j_1_real + four)*(two_j_1_real + three)*(two_j_1_real + two)*(two_j_1_real + one)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum + three + two)*(real_triangle_sum + four)* &
                (real_triangle_sum + three)*(real_triangle_sum + two&
               )/((two_j_2_real + three + two)*(two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real + two)* &
                   (two_j_2_real + one)*&
               two_j_2_real*(two_j_2_real - one)))
         else if (two_j_4 - two_j_2 == 0) then
            six_j_value = -phase_sign*two*dsqrt((three + two)*(real_triangle_sum - two_j_1_real)* &
                (real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - two_j_1_real&
                - two)*(real_triangle_sum - two_j_2_real + three)*(real_triangle_sum - two_j_2_real + two)* &
                    (real_triangle_sum - two_j_2_real + one)/((two_j_1_real + &
               seven)*(two_j_1_real + three + three)*(two_j_1_real + three + two)*(two_j_1_real + four)*(two_j_1_real + &
               three)*(two_j_1_real + two)*(two_j_1_real + one)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_3_real + three)*(real_triangle_sum - &
                two_j_3_real + two)*(real_triangle_sum - two_j_3_real + one)*(real_triangle_sum + &
               four)*(real_triangle_sum + three)*(real_triangle_sum + two)/((two_j_2_real + four)*(two_j_2_real + &
                   three)*(two_j_2_real + two)*&
               (two_j_2_real + one)*two_j_2_real*(two_j_2_real - one)*(two_j_2_real - two)))
         end if
      else if (two_j_1 - two_j_5 == 4) then
         if (two_j_2 - two_j_4 == 6) then
            six_j_value = phase_sign*dsqrt(two*three*(real_triangle_sum - two_j_1_real)*(real_triangle_sum - &
                two_j_2_real + one)*(real_triangle_sum - three)*(real_triangle_sum&
                - two)*(real_triangle_sum - one)*real_triangle_sum*(real_triangle_sum + one)/((two_j_2_real - &
                    three - two)*(two_j_2_real - four)*(&
               two_j_2_real - three)*(two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_3_real - four)*(real_triangle_sum - &
                two_j_3_real - three)*(real_triangle_sum - two_j_3_real - two)*(real_triangle_sum&
                - two_j_3_real - one)*(real_triangle_sum - two_j_3_real)/((two_j_1_real - four)*(two_j_1_real - three)* &
                    (two_j_1_real - two)*(two_j_1_real - one)&
               *two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two)))
         else if (two_j_4 - two_j_2 == 6) then
            six_j_value = -phase_sign*dsqrt(two*three*(real_triangle_sum - two_j_3_real + one)*(real_triangle_sum + &
                two)*(real_triangle_sum - two_j_2_real - four)&
               *(real_triangle_sum - two_j_2_real - three)*(real_triangle_sum - two_j_2_real - two)* &
                   (real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - two_j_2_real)/((two_j_2_real + &
               seven)*(two_j_2_real + three + three)*(two_j_2_real + three + two)*(two_j_2_real + four)*(two_j_2_real + &
               three)*(two_j_2_real + two)*(two_j_2_real + one)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_1_real + three + two)*(real_triangle_sum - &
                two_j_1_real + four)*(real_triangle_sum - two_j_1_real + &
               three)*(real_triangle_sum - two_j_1_real + two)*(real_triangle_sum - two_j_1_real + one)/ &
                   ((two_j_1_real - four)*(two_j_1_real - three)*(two_j_1_real&
                - two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two)))
         else if (two_j_2 - two_j_4 == 4) then
            six_j_value = phase_sign*((two + three)*(real_triangle_sum - two_j_2_real)*(real_triangle_sum - &
                two_j_1_real) - (real_triangle_sum + two)*(real_triangle_sum - two_j_3_real - &
               four))*dsqrt((real_triangle_sum - two)*(real_triangle_sum - one)*real_triangle_sum* &
                   (real_triangle_sum + one)/((two_j_1_real - four)*(two_j_1_real&
                - three)*(two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_3_real - three)*(real_triangle_sum - &
                two_j_3_real - two)*(real_triangle_sum - two_j_3_real - one)*(real_triangle_sum - &
               two_j_3_real)/((two_j_2_real - four)*(two_j_2_real - three)*(two_j_2_real - two)*(two_j_2_real - one)* &
                   two_j_2_real*(two_j_2_real + one)*(two_j_2_real + &
               two)))
         else if (two_j_4 - two_j_2 == 4) then
            six_j_value = -phase_sign*((two + three)*(real_triangle_sum + one)*(real_triangle_sum - two_j_3_real + &
                one) - (real_triangle_sum - two_j_2_real + one)&
               *(real_triangle_sum - two_j_1_real + three + two))*dsqrt((real_triangle_sum - two_j_2_real - three)* &
                   (real_triangle_sum - two_j_2_real - two)*(&
               real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - two_j_2_real)/((two_j_1_real - four)* &
                   (two_j_1_real - three)*(two_j_1_real - two)*(two_j_1_real - &
               one)*two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_1_real + four)*(real_triangle_sum - &
                two_j_1_real + three)*(real_triangle_sum - two_j_1_real + two)*(real_triangle_sum&
                - two_j_1_real + one)/((two_j_2_real + three + three)*(two_j_2_real + three + two)*(two_j_2_real + &
                    four)*(two_j_2_real&
                + three)*(two_j_2_real + two)*(two_j_2_real + one)*two_j_2_real))
         else if (two_j_2 - two_j_4 == 2) then
            six_j_value = phase_sign*(two*(real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - two_j_1_real) &
                - (real_triangle_sum + two)*(real_triangle_sum - two_j_3_real - three)&
               )*dsqrt(two*(three + two)*(real_triangle_sum - two_j_2_real)*(real_triangle_sum - two_j_1_real + one)* &
                   (real_triangle_sum - one)*real_triangle_sum*&
               (real_triangle_sum + one)/((two_j_1_real - four)*(two_j_1_real - three)*(two_j_1_real - two)* &
                   (two_j_1_real - one)*two_j_1_real*(two_j_1_real + &
               one)*(two_j_1_real + two)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_3_real - two)*(real_triangle_sum - &
                two_j_3_real - one)*(real_triangle_sum - two_j_3_real)/((two_j_2_real - three)*(&
               two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)*(two_j_2_real + two)* &
                   (two_j_2_real + three)))
         else if (two_j_4 - two_j_2 == 2) then
            six_j_value = -phase_sign*(two*real_triangle_sum*(real_triangle_sum - two_j_3_real + one) - &
                (real_triangle_sum - two_j_2_real + one)*(real_triangle_sum - two_j_1_real + four))*&
               dsqrt(two*(three + two)*(real_triangle_sum + one)*(real_triangle_sum - two_j_3_real)* &
                   (real_triangle_sum - two_j_2_real - two)*(real_triangle_sum&
                - two_j_2_real - one)*(real_triangle_sum - two_j_2_real)/((two_j_1_real - four)*(two_j_1_real - three)* &
                    (two_j_1_real - two)*(two_j_1_real - one)&
               *two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_1_real + three)*(real_triangle_sum - &
                two_j_1_real + two)*(real_triangle_sum - two_j_1_real + one)/((two_j_2_real + &
               three + two)*(two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real + two)*(two_j_2_real + one)* &
                   two_j_2_real*(two_j_2_real - &
               one)))
         else if (two_j_4 - two_j_2 == 0) then
            six_j_value = -phase_sign*((real_triangle_sum - one)*(real_triangle_sum - two_j_3_real + one) - &
                (real_triangle_sum - two_j_2_real + one)*(real_triangle_sum - two_j_1_real + &
               three))*dsqrt(two*three*(three + two)*(real_triangle_sum + one)*real_triangle_sum* &
                   (real_triangle_sum - two_j_3_real - &
               one)*(real_triangle_sum - two_j_3_real)/((two_j_1_real - four)*(two_j_1_real - three)*(two_j_1_real - &
                   two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real&
                + one)*(two_j_1_real + two)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - &
                two_j_2_real)*(real_triangle_sum - two_j_1_real + two)*(real_triangle_sum - two_j_1_real + one)&
               /((two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real + two)*(two_j_2_real + one)*two_j_2_real* &
                   (two_j_2_real - one)*(two_j_2_real - &
               two)))
         end if
      else if (two_j_5 - two_j_1 == 4) then
         if (two_j_2 - two_j_4 == 6) then
            six_j_value = phase_sign*dsqrt(two*three*(real_triangle_sum + one)*(real_triangle_sum - two_j_3_real)* &
                (real_triangle_sum - two_j_1_real - four)*(real_triangle_sum - &
               two_j_1_real - three)*(real_triangle_sum - two_j_1_real - two)*(real_triangle_sum - two_j_1_real - one) &
                   *(real_triangle_sum - two_j_1_real)/((two_j_1_real + three + &
               three)*(two_j_1_real + three + two)*(two_j_1_real + four)*(two_j_1_real + three)*(two_j_1_real + two)* &
                   (two_j_1_real + &
               one)*two_j_1_real))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_2_real + three + two)*(real_triangle_sum - &
                two_j_2_real + four)*(real_triangle_sum - two_j_2_real + &
               three)*(real_triangle_sum - two_j_2_real + two)*(real_triangle_sum - two_j_2_real + one)/ &
                   ((two_j_2_real - three - two)*(two_j_2_real - &
               four)*(two_j_2_real - three)*(two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)))
         else if (two_j_4 - two_j_2 == 6) then
            six_j_value = -phase_sign*dsqrt(two*three*(real_triangle_sum - two_j_2_real)*(real_triangle_sum - &
                two_j_1_real + one)*(real_triangle_sum - two_j_3_real + three&
                + two)*(real_triangle_sum - two_j_3_real + four)*(real_triangle_sum - two_j_3_real + three)* &
                    (real_triangle_sum - two_j_3_real + two)*(real_triangle_sum - two_j_3_real&
                + one)/((two_j_2_real + seven)*(two_j_2_real + three + three)*(two_j_2_real + three + two)* &
                    (two_j_2_real + &
               four)*(two_j_2_real + three)*(two_j_2_real + two)*(two_j_2_real + one)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum + three + three)*(real_triangle_sum + three + two)* &
                (real_triangle_sum + four)*(&
               real_triangle_sum + three)*(real_triangle_sum + two)/((two_j_1_real + three + three)*(two_j_1_real + &
                   three + two)*(two_j_1_real&
                + four)*(two_j_1_real + three)*(two_j_1_real + two)*(two_j_1_real + one)*two_j_1_real))
         else if (two_j_2 - two_j_4 == 4) then
            six_j_value = -phase_sign*((two + three)*(real_triangle_sum + one)*(real_triangle_sum - two_j_3_real + &
                one) - (real_triangle_sum - two_j_1_real + one)&
               *(real_triangle_sum - two_j_2_real + three + two))*dsqrt((real_triangle_sum - two_j_1_real - three)* &
                   (real_triangle_sum - two_j_1_real - two)*(&
               real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - two_j_1_real)/((two_j_2_real - four)* &
                   (two_j_2_real - three)*(two_j_2_real - two)*(two_j_2_real - &
               one)*two_j_2_real*(two_j_2_real + one)*(two_j_2_real + two)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_2_real + four)*(real_triangle_sum - &
                two_j_2_real + three)*(real_triangle_sum - two_j_2_real + two)*(real_triangle_sum&
                - two_j_2_real + one)/((two_j_1_real + three + three)*(two_j_1_real + three + two)*(two_j_1_real + &
                    four)*(two_j_1_real&
                + three)*(two_j_1_real + two)*(two_j_1_real + one)*two_j_1_real))
         else if (two_j_4 - two_j_2 == 4) then
            six_j_value = phase_sign*((two + three)*(real_triangle_sum - two_j_1_real)*(real_triangle_sum - &
                two_j_2_real) - (real_triangle_sum - two_j_3_real)*(real_triangle_sum + three + &
               three))*dsqrt((real_triangle_sum - two_j_3_real + four)*(real_triangle_sum - two_j_3_real + three)* &
                   (real_triangle_sum - two_j_3_real + two)*(&
               real_triangle_sum - two_j_3_real + one)/((two_j_1_real + three + three)*(two_j_1_real + three + two)* &
                   (two_j_1_real + four)*&
               (two_j_1_real + three)*(two_j_1_real + two)*(two_j_1_real + one)*two_j_1_real))
            six_j_value = six_j_value*dsqrt((real_triangle_sum + three + two)*(real_triangle_sum + four)* &
                (real_triangle_sum + three)*(real_triangle_sum + two&
               )/((two_j_2_real + three + three)*(two_j_2_real + three + two)*(two_j_2_real + four)*(two_j_2_real + &
                   three)*&
               (two_j_2_real + two)*(two_j_2_real + one)*two_j_2_real))
         else if (two_j_2 - two_j_4 == 2) then
            six_j_value = phase_sign*(two*(real_triangle_sum - two_j_3_real + two)*(real_triangle_sum + one) - &
                (real_triangle_sum - two_j_1_real + one)*(real_triangle_sum - two_j_2_real + &
               four))*dsqrt(two*(three + two)*(real_triangle_sum - two_j_3_real + one)*(real_triangle_sum + two)* &
                   (real_triangle_sum - two_j_1_real&
                - two)*(real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - two_j_1_real)/((two_j_1_real + &
                    three + three)*(two_j_1_real + three&
                + two)*(two_j_1_real + four)*(two_j_1_real + three)*(two_j_1_real + two)*(two_j_1_real + one)* &
                    two_j_1_real))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_2_real + three)*(real_triangle_sum - &
                two_j_2_real + two)*(real_triangle_sum - two_j_2_real + one)/((two_j_2_real - &
               three)*(two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)*(two_j_2_real + two)* &
                   (two_j_2_real + three)))
         else if (two_j_4 - two_j_2 == 2) then
            six_j_value = -phase_sign*(two*(real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - &
                two_j_2_real) - (real_triangle_sum - two_j_3_real)*(real_triangle_sum + three + two&
               ))*dsqrt(two*(three + two)*(real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_2_real + one) &
                   *(real_triangle_sum - two_j_3_real + &
               three)*(real_triangle_sum - two_j_3_real + two)*(real_triangle_sum - two_j_3_real + one)/ &
                   ((two_j_1_real + three + three)*(two_j_1_real + &
               three + two)*(two_j_1_real + four)*(two_j_1_real + three)*(two_j_1_real + two)*(two_j_1_real + one)* &
                   two_j_1_real))
            six_j_value = six_j_value*dsqrt((real_triangle_sum + four)*(real_triangle_sum + three)* &
                (real_triangle_sum + two)/((two_j_2_real + three + two&
               )*(two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real + two)*(two_j_2_real + one)*two_j_2_real* &
                   (two_j_2_real - one)))
         else if (two_j_4 - two_j_2 == 0) then
            six_j_value = phase_sign*((real_triangle_sum - two_j_1_real - two)*(real_triangle_sum - two_j_2_real) - &
                (real_triangle_sum - two_j_3_real)*(real_triangle_sum + four))*dsqrt(two&
               *three*(three + two)*(real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_1_real - one)* &
                   (real_triangle_sum - two_j_2_real + two)*(real_triangle_sum&
                - two_j_2_real + one)/((two_j_1_real + three + three)*(two_j_1_real + three + two)*(two_j_1_real + &
                    four)*(two_j_1_real&
                + three)*(two_j_1_real + two)*(two_j_1_real + one)*two_j_1_real))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_3_real + two)*(real_triangle_sum - &
                two_j_3_real + one)*(real_triangle_sum + three)*(real_triangle_sum + two)&
               /((two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real + two)*(two_j_2_real + one)*two_j_2_real* &
                   (two_j_2_real - one)*(two_j_2_real - &
               two)))
         end if
      else if (two_j_1 - two_j_5 == 2) then
         if (two_j_2 - two_j_4 == 6) then
            six_j_value = phase_sign*dsqrt(three*(three + two)*(real_triangle_sum - two_j_1_real)* &
                (real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - two_j_2_real&
                + two)*(real_triangle_sum - two_j_2_real + one)*(real_triangle_sum - two)*(real_triangle_sum - one)* &
                    real_triangle_sum*(real_triangle_sum + one)/((two_j_2_real&
                - three - two)*(two_j_2_real - four)*(two_j_2_real - three)*(two_j_2_real - two)*(two_j_2_real - one)* &
                    two_j_2_real*(two_j_2_real&
                + one)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_3_real - three)*(real_triangle_sum - &
                two_j_3_real - two)*(real_triangle_sum - two_j_3_real - one)*(real_triangle_sum - &
               two_j_3_real)/((two_j_1_real - three)*(two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real* &
                   (two_j_1_real + one)*(two_j_1_real + two)*(two_j_1_real + &
               three)))
         else if (two_j_4 - two_j_2 == 6) then
            six_j_value = phase_sign*dsqrt(three*(three + two)*(real_triangle_sum - two_j_3_real + one)* &
                (real_triangle_sum - two_j_3_real + two)*(&
               real_triangle_sum + three)*(real_triangle_sum + two)*(real_triangle_sum - two_j_2_real - three)* &
                   (real_triangle_sum - two_j_2_real - two)*(real_triangle_sum - two_j_2_real&
                - one)*(real_triangle_sum - two_j_2_real)/((two_j_2_real + seven)*(two_j_2_real + three + three)* &
                    (two_j_2_real + three + &
               two)*(two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real + two)*(two_j_2_real + one)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_1_real + four)*(real_triangle_sum - &
                two_j_1_real + three)*(real_triangle_sum - two_j_1_real + two)*(real_triangle_sum&
                - two_j_1_real + one)/((two_j_1_real - three)*(two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real* &
                    (two_j_1_real + one)*(two_j_1_real + &
               two)*(two_j_1_real + three)))
         else if (two_j_2 - two_j_4 == 4) then
            six_j_value = phase_sign*(two*(real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - two_j_2_real) &
                - (real_triangle_sum + two)*(real_triangle_sum - two_j_3_real - three)&
               )*dsqrt(two*(three + two)*(real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_2_real + one)* &
                   (real_triangle_sum - one)*real_triangle_sum*&
               (real_triangle_sum + one)/((two_j_2_real - four)*(two_j_2_real - three)*(two_j_2_real - two)* &
                   (two_j_2_real - one)*two_j_2_real*(two_j_2_real + &
               one)*(two_j_2_real + two)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_3_real - two)*(real_triangle_sum - &
                two_j_3_real - one)*(real_triangle_sum - two_j_3_real)/((two_j_1_real - three)*(&
               two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two)* &
                   (two_j_1_real + three)))
         else if (two_j_4 - two_j_2 == 4) then
            six_j_value = phase_sign*(two*(real_triangle_sum - two_j_3_real + two)*(real_triangle_sum + one) - &
                (real_triangle_sum - two_j_2_real + one)*(real_triangle_sum - two_j_1_real + &
               four))*dsqrt(two*(three + two)*(real_triangle_sum - two_j_3_real + one)*(real_triangle_sum + two)* &
                   (real_triangle_sum - two_j_2_real&
                - two)*(real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - two_j_2_real)/((two_j_2_real + &
                    three + three)*(two_j_2_real + three&
                + two)*(two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real + two)*(two_j_2_real + one)* &
                    two_j_2_real))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_1_real + three)*(real_triangle_sum - &
                two_j_1_real + two)*(real_triangle_sum - two_j_1_real + one)/((two_j_1_real - &
               three)*(two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two)* &
                   (two_j_1_real + three)))
         else if (two_j_2 - two_j_4 == 2) then
            six_j_value = phase_sign*(two*three*(real_triangle_sum - two_j_2_real)*(real_triangle_sum - &
                two_j_2_real - one)*(real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - two_j_1_real)&
                - two*four*(real_triangle_sum - two_j_2_real)*(real_triangle_sum - two_j_1_real)*(real_triangle_sum + &
                    two)*(real_triangle_sum - two_j_3_real - two) + (real_triangle_sum + &
               three)*(real_triangle_sum + two)*(real_triangle_sum - two_j_3_real - three)*(real_triangle_sum - &
                   two_j_3_real - two))
            six_j_value = six_j_value*dsqrt(real_triangle_sum*(real_triangle_sum + one)*(real_triangle_sum - &
                two_j_3_real - one)*(real_triangle_sum - two_j_3_real)/((two_j_1_real - three)*(two_j_1_real&
                - two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two)*(two_j_1_real + &
                    three)*(two_j_2_real - three)*&
               (two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)*(two_j_2_real + two)* &
                   (two_j_2_real + three)))
         else if (two_j_4 - two_j_2 == 2) then
            six_j_value = phase_sign*(two*three*(real_triangle_sum + one)*real_triangle_sum*(real_triangle_sum - &
                two_j_3_real + one)*(real_triangle_sum - two_j_3_real + two) - &
               two*four*(real_triangle_sum + one)*(real_triangle_sum - two_j_3_real + one)*(real_triangle_sum - &
                   two_j_2_real + one)*(real_triangle_sum - two_j_1_real + &
               three) + (real_triangle_sum - two_j_2_real + one)*(real_triangle_sum - two_j_2_real + two)* &
                   (real_triangle_sum - two_j_1_real + four)*(real_triangle_sum - two_j_1_real&
                + three))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - &
                two_j_2_real)*(real_triangle_sum - two_j_1_real + two)*(real_triangle_sum - two_j_1_real + one)&
               /((two_j_1_real - three)*(two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)* &
                   (two_j_1_real + two)*(two_j_1_real + &
               three)*(two_j_2_real + three + two)*(two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real + two)* &
                   (two_j_2_real + &
               one)*two_j_2_real*(two_j_2_real - one)))
         else if (two_j_4 - two_j_2 == 0) then
            six_j_value = phase_sign*two*(real_triangle_sum*(real_triangle_sum - one)*(real_triangle_sum - &
                two_j_3_real + one)*(real_triangle_sum - two_j_3_real + two) - three*&
               real_triangle_sum*(real_triangle_sum - two_j_3_real + one)*(real_triangle_sum - two_j_2_real + one)* &
                   (real_triangle_sum - two_j_1_real + two) + (real_triangle_sum - two_j_2_real + one)&
               *(real_triangle_sum - two_j_2_real + two)*(real_triangle_sum - two_j_1_real + three)* &
                   (real_triangle_sum - two_j_1_real + two))
            six_j_value = six_j_value*dsqrt(three*(real_triangle_sum + one)*(real_triangle_sum - two_j_3_real)* &
                (real_triangle_sum - two_j_2_real)*(real_triangle_sum - two_j_1_real + one)/((two_j_1_real&
                - three)*(two_j_1_real - two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)*(two_j_1_real + &
                    two)*(two_j_1_real + three)*&
               (two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real + two)*(two_j_2_real + one)*two_j_2_real* &
                   (two_j_2_real - one)*(two_j_2_real - two)&
               ))
         end if
      else if (two_j_5 - two_j_1 == 2) then
         if (two_j_2 - two_j_4 == 6) then
            six_j_value = phase_sign*dsqrt(three*(three + two)*(real_triangle_sum + one)*real_triangle_sum* &
                (real_triangle_sum - two_j_3_real - one)*(real_triangle_sum&
                - two_j_3_real)*(real_triangle_sum - two_j_1_real - three)*(real_triangle_sum - two_j_1_real - two)* &
                    (real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - two_j_1_real)/((&
               two_j_2_real - three - two)*(two_j_2_real - four)*(two_j_2_real - three)*(two_j_2_real - two)* &
                   (two_j_2_real - one)*two_j_2_real*(two_j_2_real&
                + one)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_2_real + four)*(real_triangle_sum - &
                two_j_2_real + three)*(real_triangle_sum - two_j_2_real + two)*(real_triangle_sum&
                - two_j_2_real + one)/((two_j_1_real + three + two)*(two_j_1_real + four)*(two_j_1_real + three)* &
                    (two_j_1_real + two)*&
               (two_j_1_real + one)*two_j_1_real*(two_j_1_real - one)))
         else if (two_j_4 - two_j_2 == 6) then
            six_j_value = phase_sign*dsqrt(three*(three + two)*(real_triangle_sum - two_j_2_real)* &
                (real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - two_j_1_real&
                + two)*(real_triangle_sum - two_j_1_real + one)*(real_triangle_sum - two_j_3_real + four)* &
                    (real_triangle_sum - two_j_3_real + three)*(real_triangle_sum - two_j_3_real&
                + two)*(real_triangle_sum - two_j_3_real + one)/((two_j_2_real + seven)*(two_j_2_real + three + three) &
                    *(two_j_2_real + &
               three + two)*(two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real + two)*(two_j_2_real + one)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum + three + two)*(real_triangle_sum + four)* &
                (real_triangle_sum + three)*(real_triangle_sum + two&
               )/((two_j_1_real + three + two)*(two_j_1_real + four)*(two_j_1_real + three)*(two_j_1_real + two)* &
                   (two_j_1_real + one)*&
               two_j_1_real*(two_j_1_real - one)))
         else if (two_j_2 - two_j_4 == 4) then
            six_j_value = -phase_sign*(two*real_triangle_sum*(real_triangle_sum - two_j_3_real + one) - &
                (real_triangle_sum - two_j_1_real + one)*(real_triangle_sum - two_j_2_real + four))*&
               dsqrt(two*(three + two)*(real_triangle_sum + one)*(real_triangle_sum - two_j_3_real)* &
                   (real_triangle_sum - two_j_1_real - two)*(real_triangle_sum&
                - two_j_1_real - one)*(real_triangle_sum - two_j_1_real)/((two_j_2_real - four)*(two_j_2_real - three)* &
                    (two_j_2_real - two)*(two_j_2_real - one)&
               *two_j_2_real*(two_j_2_real + one)*(two_j_2_real + two)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_2_real + three)*(real_triangle_sum - &
                two_j_2_real + two)*(real_triangle_sum - two_j_2_real + one)/((two_j_1_real + &
               three + two)*(two_j_1_real + four)*(two_j_1_real + three)*(two_j_1_real + two)*(two_j_1_real + one)* &
                   two_j_1_real*(two_j_1_real - &
               one)))
         else if (two_j_4 - two_j_2 == 4) then
            six_j_value = -phase_sign*(two*(real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - &
                two_j_1_real) - (real_triangle_sum - two_j_3_real)*(real_triangle_sum + three + two&
               ))*dsqrt(two*(three + two)*(real_triangle_sum - two_j_2_real)*(real_triangle_sum - two_j_1_real + one) &
                   *(real_triangle_sum - two_j_3_real + &
               three)*(real_triangle_sum - two_j_3_real + two)*(real_triangle_sum - two_j_3_real + one)/ &
                   ((two_j_2_real + three + three)*(two_j_2_real + &
               three + two)*(two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real + two)*(two_j_2_real + one)* &
                   two_j_2_real))
            six_j_value = six_j_value*dsqrt((real_triangle_sum + four)*(real_triangle_sum + three)* &
                (real_triangle_sum + two)/((two_j_1_real + three + two&
               )*(two_j_1_real + four)*(two_j_1_real + three)*(two_j_1_real + two)*(two_j_1_real + one)*two_j_1_real* &
                   (two_j_1_real - one)))
         else if (two_j_4 - two_j_2 == 2) then
            six_j_value = phase_sign*(two*three*(real_triangle_sum - two_j_1_real)*(real_triangle_sum - &
                two_j_1_real - one)*(real_triangle_sum - two_j_2_real)*(real_triangle_sum - two_j_2_real - one)&
                - two*four*(real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_2_real)*(real_triangle_sum - &
                    two_j_3_real)*(real_triangle_sum + four) + (real_triangle_sum - two_j_3_real)*(&
               real_triangle_sum - two_j_3_real - one)*(real_triangle_sum + three + two)*(real_triangle_sum + four))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_3_real + two)*(real_triangle_sum - &
                two_j_3_real + one)*(real_triangle_sum + three)*(real_triangle_sum + two)&
               /((two_j_1_real + three + two)*(two_j_1_real + four)*(two_j_1_real + three)*(two_j_1_real + two)* &
                   (two_j_1_real + one)*two_j_1_real&
               *(two_j_1_real - one)*(two_j_2_real + three + two)*(two_j_2_real + four)*(two_j_2_real + three)* &
                   (two_j_2_real + two)*(two_j_2_real&
                + one)*two_j_2_real*(two_j_2_real - one)))
         else if (two_j_2 - two_j_4 == 2) then
            six_j_value = phase_sign*(two*three*(real_triangle_sum + one)*real_triangle_sum*(real_triangle_sum - &
                two_j_3_real + one)*(real_triangle_sum - two_j_3_real + two) - &
               two*four*(real_triangle_sum + one)*(real_triangle_sum - two_j_3_real + one)*(real_triangle_sum - &
                   two_j_1_real + one)*(real_triangle_sum - two_j_2_real + &
               three) + (real_triangle_sum - two_j_1_real + one)*(real_triangle_sum - two_j_1_real + two)* &
                   (real_triangle_sum - two_j_2_real + four)*(real_triangle_sum - two_j_2_real&
                + three))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - &
                two_j_1_real)*(real_triangle_sum - two_j_2_real + two)*(real_triangle_sum - two_j_2_real + one)&
               /((two_j_2_real - three)*(two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)* &
                   (two_j_2_real + two)*(two_j_2_real + &
               three)*(two_j_1_real + three + two)*(two_j_1_real + four)*(two_j_1_real + three)*(two_j_1_real + two)* &
                   (two_j_1_real + &
               one)*two_j_1_real*(two_j_1_real - one)))
         else if (two_j_4 - two_j_2 == 0) then
            six_j_value = -phase_sign*two*((real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - &
                two_j_1_real - two)*(real_triangle_sum - two_j_2_real)*(real_triangle_sum - two_j_2_real - one&
               ) - three*(real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - two_j_2_real)* &
                   (real_triangle_sum - two_j_3_real)*(real_triangle_sum + three) + (real_triangle_sum - &
               two_j_3_real)*(real_triangle_sum - two_j_3_real - one)*(real_triangle_sum + four)*(real_triangle_sum + &
                   three))
            six_j_value = six_j_value*dsqrt(three*(real_triangle_sum - two_j_1_real)*(real_triangle_sum - &
                two_j_2_real + one)*(real_triangle_sum - two_j_3_real + one)*(real_triangle_sum + &
               two)/((two_j_1_real + three + two)*(two_j_1_real + four)*(two_j_1_real + three)*(two_j_1_real + two)* &
                   (two_j_1_real + &
               one)*two_j_1_real*(two_j_1_real - one)*(two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real + two)* &
                   (two_j_2_real + one)*two_j_2_real*(two_j_2_real&
                - one)*(two_j_2_real - two)))
         end if
      else if (two_j_1 - two_j_5 == 0) then
         if (two_j_2 - two_j_4 == 6) then
            six_j_value = phase_sign*two*dsqrt((three + two)*(real_triangle_sum + one)*real_triangle_sum* &
                (real_triangle_sum - one)*(real_triangle_sum - two_j_3_real&
                - two)*(real_triangle_sum - two_j_3_real - one)*(real_triangle_sum - two_j_3_real)/((two_j_2_real - &
                    three - two)*(two_j_2_real - four)*(two_j_2_real&
                - three)*(two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_1_real - two)*(real_triangle_sum - &
                two_j_1_real - one)*(real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_2_real + &
               three)*(real_triangle_sum - two_j_2_real + two)*(real_triangle_sum - two_j_2_real + one)/ &
                   ((two_j_1_real + four)*(two_j_1_real + three)*(two_j_1_real&
                + two)*(two_j_1_real + one)*two_j_1_real*(two_j_1_real - one)*(two_j_1_real - two)))
         else if (two_j_4 - two_j_2 == 6) then
            six_j_value = -phase_sign*two*dsqrt((three + two)*(real_triangle_sum - two_j_2_real)* &
                (real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - two_j_2_real&
                - two)*(real_triangle_sum - two_j_1_real + three)*(real_triangle_sum - two_j_1_real + two)* &
                    (real_triangle_sum - two_j_1_real + one)/((two_j_2_real + &
               seven)*(two_j_2_real + three + three)*(two_j_2_real + three + two)*(two_j_2_real + four)*(two_j_2_real + &
               three)*(two_j_2_real + two)*(two_j_2_real + one)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_3_real + three)*(real_triangle_sum - &
                two_j_3_real + two)*(real_triangle_sum - two_j_3_real + one)*(real_triangle_sum + &
               four)*(real_triangle_sum + three)*(real_triangle_sum + two)/((two_j_1_real + four)*(two_j_1_real + &
                   three)*(two_j_1_real + two)*&
               (two_j_1_real + one)*two_j_1_real*(two_j_1_real - one)*(two_j_1_real - two)))
         else if (two_j_2 - two_j_4 == 4) then
            six_j_value = -phase_sign*((real_triangle_sum - one)*(real_triangle_sum - two_j_3_real + one) - &
                (real_triangle_sum - two_j_1_real + one)*(real_triangle_sum - two_j_2_real + &
               three))*dsqrt(two*three*(three + two)*(real_triangle_sum + one)*real_triangle_sum* &
                   (real_triangle_sum - two_j_3_real - &
               one)*(real_triangle_sum - two_j_3_real)/((two_j_2_real - four)*(two_j_2_real - three)*(two_j_2_real - &
                   two)*(two_j_2_real - one)*two_j_2_real*(two_j_2_real&
                + one)*(two_j_2_real + two)))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - &
                two_j_1_real)*(real_triangle_sum - two_j_2_real + two)*(real_triangle_sum - two_j_2_real + one)&
               /((two_j_1_real + four)*(two_j_1_real + three)*(two_j_1_real + two)*(two_j_1_real + one)*two_j_1_real* &
                   (two_j_1_real - one)*(two_j_1_real - &
               two)))
         else if (two_j_4 - two_j_2 == 4) then
            six_j_value = phase_sign*((real_triangle_sum - two_j_2_real - two)*(real_triangle_sum - two_j_1_real) - &
                (real_triangle_sum - two_j_3_real)*(real_triangle_sum + four))*dsqrt(two&
               *three*(three + two)*(real_triangle_sum - two_j_2_real)*(real_triangle_sum - two_j_2_real - one)* &
                   (real_triangle_sum - two_j_1_real + two)*(real_triangle_sum&
                - two_j_1_real + one)/((two_j_2_real + three + three)*(two_j_2_real + three + two)*(two_j_2_real + &
                    four)*(two_j_2_real&
                + three)*(two_j_2_real + two)*(two_j_2_real + one)*two_j_2_real))
            six_j_value = six_j_value*dsqrt((real_triangle_sum - two_j_3_real + two)*(real_triangle_sum - &
                two_j_3_real + one)*(real_triangle_sum + three)*(real_triangle_sum + two)&
               /((two_j_1_real + four)*(two_j_1_real + three)*(two_j_1_real + two)*(two_j_1_real + one)*two_j_1_real* &
                   (two_j_1_real - one)*(two_j_1_real - &
               two)))
         else if (two_j_2 - two_j_4 == 2) then
            six_j_value = phase_sign*two*(real_triangle_sum*(real_triangle_sum - one)*(real_triangle_sum - &
                two_j_3_real + one)*(real_triangle_sum - two_j_3_real + two) - three*&
               real_triangle_sum*(real_triangle_sum - two_j_3_real + one)*(real_triangle_sum - two_j_1_real + one)* &
                   (real_triangle_sum - two_j_2_real + two) + (real_triangle_sum - two_j_1_real + one)&
               *(real_triangle_sum - two_j_1_real + two)*(real_triangle_sum - two_j_2_real + three)* &
                   (real_triangle_sum - two_j_2_real + two))
            six_j_value = six_j_value*dsqrt(three*(real_triangle_sum + one)*(real_triangle_sum - two_j_3_real)* &
                (real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_2_real + one)/((two_j_2_real&
                - three)*(two_j_2_real - two)*(two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)*(two_j_2_real + &
                    two)*(two_j_2_real + three)*&
               (two_j_1_real + four)*(two_j_1_real + three)*(two_j_1_real + two)*(two_j_1_real + one)*two_j_1_real* &
                   (two_j_1_real - one)*(two_j_1_real - two)&
               ))
         else if (two_j_4 - two_j_2 == 2) then
            six_j_value = -phase_sign*two*((real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - &
                two_j_2_real - two)*(real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_1_real - one&
               ) - three*(real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - two_j_1_real)* &
                   (real_triangle_sum - two_j_3_real)*(real_triangle_sum + three) + (real_triangle_sum - &
               two_j_3_real)*(real_triangle_sum - two_j_3_real - one)*(real_triangle_sum + four)*(real_triangle_sum + &
                   three))
            six_j_value = six_j_value*dsqrt(three*(real_triangle_sum - two_j_2_real)*(real_triangle_sum - &
                two_j_1_real + one)*(real_triangle_sum - two_j_3_real + one)*(real_triangle_sum + &
               two)/((two_j_2_real + three + two)*(two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real + two)* &
                   (two_j_2_real + &
               one)*two_j_2_real*(two_j_2_real - one)*(two_j_1_real + four)*(two_j_1_real + three)*(two_j_1_real + two)* &
                   (two_j_1_real + one)*two_j_1_real*(two_j_1_real&
                - one)*(two_j_1_real - two)))
         else if (two_j_2 - two_j_4 == 0) then
            six_j_value = phase_sign*((real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_1_real - one)* &
                (real_triangle_sum - two_j_1_real - two)*(real_triangle_sum - two_j_2_real)*(real_triangle_sum - &
                two_j_2_real&
                - one)*(real_triangle_sum - two_j_2_real - two) - three*three*(real_triangle_sum - two_j_1_real)* &
                    (real_triangle_sum - two_j_1_real - one)*(real_triangle_sum&
                - two_j_2_real)*(real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - two_j_3_real)* &
                    (real_triangle_sum + two) + three*three*(real_triangle_sum - two_j_1_real)*&
               (real_triangle_sum - two_j_2_real)*(real_triangle_sum - two_j_3_real)*(real_triangle_sum - &
                   two_j_3_real - one)*(real_triangle_sum + three)*(real_triangle_sum + two) - (real_triangle_sum&
                - two_j_3_real)*(real_triangle_sum - two_j_3_real - one)*(real_triangle_sum - two_j_3_real - two)* &
                    (real_triangle_sum + four)*(real_triangle_sum + three)*(real_triangle_sum&
                + two))
            six_j_value = six_j_value/dsqrt((two_j_1_real + four)*(two_j_1_real + three)*(two_j_1_real + two)* &
                (two_j_1_real + one)*two_j_1_real*(two_j_1_real - one&
               )*(two_j_1_real - two)*(two_j_2_real + four)*(two_j_2_real + three)*(two_j_2_real + two)*(two_j_2_real + &
                   one)*two_j_2_real*(two_j_2_real - &
               one)*(two_j_2_real - two))
         end if
      end if
      return
      end subroutine sixj3
end module wigner_6j_with_j_three
