! Evaluate a Wigner 6j symbol with sixth angular momentum 1.
module wigner_6j_with_j_one
    implicit none (type, external)
    private
    public :: sixj1
contains
      subroutine sixj1(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,check_triangle_conditions,six_j_value)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants,          only: zero, half, one, two, three
      use wigner_6j_triangle_selection, only: ixjtik
      implicit none (type, external)
      integer, intent(in) :: two_j_1, two_j_2, two_j_3, two_j_4, two_j_5
      integer, intent(in)       :: check_triangle_conditions
      real(real64), intent(out) :: six_j_value
      integer      :: triangle_sum
      real(real64) :: real_triangle_sum, phase_sign, two_j_3_real, two_j_2_real, two_j_1_real
      six_j_value = zero
      if (check_triangle_conditions /= 0) then
         if (ixjtik(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,2) == 0) return
      end if
      triangle_sum = (two_j_1 + two_j_2 + two_j_3)/2
      real_triangle_sum = dble(triangle_sum)
      phase_sign = one
      if (mod(triangle_sum,2) /= 0) phase_sign = -phase_sign
      two_j_3_real = dble(two_j_3)
      two_j_2_real = dble(two_j_2)
      two_j_1_real = dble(two_j_1)
      if (two_j_1 < two_j_5) then
         if (two_j_2 < two_j_4) then
            six_j_value = phase_sign*dsqrt((real_triangle_sum + two)*(real_triangle_sum + three)* &
                (real_triangle_sum - two_j_3_real + one)*(real_triangle_sum - two_j_3_real + two&
               )/((two_j_2_real + one)*(two_j_2_real + two)*(two_j_2_real + three)*(two_j_1_real + one)*(two_j_1_real + &
                   two)*(two_j_1_real + &
               three)))
         else if (two_j_2 == two_j_4) then
            six_j_value = (-phase_sign)*dsqrt(two*(real_triangle_sum + two)*(real_triangle_sum - two_j_1_real)* &
                (real_triangle_sum - two_j_2_real + one)*(real_triangle_sum - two_j_3_real + &
               one)/(two_j_2_real*(two_j_2_real + one)*(two_j_2_real + two)*(two_j_1_real + one)*(two_j_1_real + two)* &
                   (two_j_1_real + three)))
         else
            six_j_value = phase_sign*dsqrt((real_triangle_sum - two_j_1_real - one)*(real_triangle_sum - &
                two_j_1_real)*(real_triangle_sum - two_j_2_real + one)*(real_triangle_sum - two_j_2_real + two&
               )/((two_j_2_real - one)*two_j_2_real*(two_j_2_real + one)*(two_j_1_real + one)*(two_j_1_real + two)* &
                   (two_j_1_real + three)))
         end if
      else if (two_j_1 == two_j_5) then
         if (two_j_2 < two_j_4) then
            six_j_value = (-phase_sign)*dsqrt((real_triangle_sum + two)*(real_triangle_sum - two_j_1_real + one)* &
                (real_triangle_sum - two_j_2_real)*(real_triangle_sum - two_j_3_real + one)&
               *two/((two_j_2_real + one)*(two_j_2_real + two)*(two_j_2_real + three)*two_j_1_real*(two_j_1_real + one)* &
                   (two_j_1_real + two)))
         else if (two_j_2 == two_j_4) then
            six_j_value = (-phase_sign)*((two_j_2_real*two_j_2_real + two_j_1_real*two_j_1_real - two_j_3_real* &
                two_j_3_real)*half + two_j_2_real + two_j_1_real - two_j_3_real)/dsqrt(two_j_2_real*(two_j_2_real + one)*&
               (two_j_2_real + two)*two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two))
         else
            six_j_value = phase_sign*dsqrt((real_triangle_sum + one)*(real_triangle_sum - two_j_1_real)* &
                (real_triangle_sum - two_j_2_real + one)*(real_triangle_sum - two_j_3_real)*two/((two_j_2_real&
                - one)*two_j_2_real*(two_j_2_real + one)*two_j_1_real*(two_j_1_real + one)*(two_j_1_real + two)))
         end if
      else
         if (two_j_2 < two_j_4) then
            six_j_value = phase_sign*dsqrt((real_triangle_sum - two_j_1_real + one)*(real_triangle_sum - &
                two_j_1_real + two)*(real_triangle_sum - two_j_2_real - one)*(real_triangle_sum - two_j_2_real&
               )/((two_j_2_real + one)*(two_j_2_real + two)*(two_j_2_real + three)*(two_j_1_real - one)*two_j_1_real* &
                   (two_j_1_real + one)))
         else if (two_j_2 == two_j_4) then
            six_j_value = phase_sign*dsqrt((real_triangle_sum + one)*(real_triangle_sum - two_j_1_real + one)* &
                (real_triangle_sum - two_j_2_real)*(real_triangle_sum - two_j_3_real)*two/(two_j_2_real*(&
               two_j_2_real + one)*(two_j_2_real + two)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)))
         else
            six_j_value = phase_sign*dsqrt(real_triangle_sum*(real_triangle_sum + one)*(real_triangle_sum - &
                two_j_3_real - one)*(real_triangle_sum - two_j_3_real)/((two_j_2_real - one)*two_j_2_real*(&
               two_j_2_real + one)*(two_j_1_real - one)*two_j_1_real*(two_j_1_real + one)))
         end if
      end if
      return
      end subroutine sixj1
end module wigner_6j_with_j_one
