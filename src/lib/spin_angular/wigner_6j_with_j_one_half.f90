! Evaluate a Wigner 6j symbol with sixth angular momentum one half.
module wigner_6j_with_j_one_half
    implicit none (type, external)
    private
    public :: sixj5
contains
      subroutine sixj5(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,check_triangle_conditions,six_j_value)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants,          only: zero, one, two
      use wigner_6j_triangle_selection, only: ixjtik
      implicit none (type, external)
      integer, intent(in) :: two_j_1, two_j_2, two_j_3, two_j_4, two_j_5
      integer, intent(in)       :: check_triangle_conditions
      real(real64), intent(out) :: six_j_value
      integer      :: triangle_sum
      real(real64) :: real_triangle_sum, two_j_3_real, two_j_2_real, two_j_1_real, phase_sign
      six_j_value = zero
      if (check_triangle_conditions /= 0) then
         if (ixjtik(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,1) == 0) return
      end if
      triangle_sum = (two_j_1 + two_j_2 + two_j_3)/2
      real_triangle_sum = dble(triangle_sum)
      two_j_3_real = dble(two_j_3)
      two_j_2_real = dble(two_j_2)
      two_j_1_real = dble(two_j_1)
      phase_sign = one
      if (mod(triangle_sum,2) /= 0) phase_sign = -phase_sign
      if (two_j_2 < two_j_4) then
         if (two_j_1 < two_j_5) then
            six_j_value = -phase_sign*dsqrt((real_triangle_sum + two)*(real_triangle_sum - two_j_3_real + one)/ &
                ((two_j_2_real + one)*(two_j_2_real + two)*(two_j_1_real&
                + one)*(two_j_1_real + two)))
         else if (two_j_1 > two_j_5) then
            six_j_value = phase_sign*dsqrt((real_triangle_sum - two_j_1_real + one)*(real_triangle_sum - &
                two_j_2_real)/((two_j_2_real + one)*(two_j_2_real + two)*two_j_1_real*(two_j_1_real + &
               one)))
         end if
      else if (two_j_2 > two_j_4) then
         if (two_j_1 < two_j_5) then
            six_j_value = phase_sign*dsqrt((real_triangle_sum - two_j_1_real)*(real_triangle_sum - two_j_2_real + &
                one)/(two_j_2_real*(two_j_2_real + one)*(two_j_1_real + one)*(two_j_1_real + &
               two)))
         else if (two_j_1 > two_j_5) then
            six_j_value = phase_sign*dsqrt((real_triangle_sum + one)*(real_triangle_sum - two_j_3_real)/ &
                (two_j_2_real*(two_j_2_real + one)*two_j_1_real*(two_j_1_real + one)))
         end if
      end if
      return
      end subroutine sixj5
end module wigner_6j_with_j_one_half
