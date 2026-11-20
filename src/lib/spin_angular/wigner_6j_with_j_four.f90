! Evaluate a Wigner 6j symbol with sixth angular momentum 4.
module wigner_6j_with_j_four
    implicit none (type, external)
    private
    public :: sixj4
contains
      subroutine sixj4(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,check_triangle_conditions,six_j_value)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants,          only: zero, half, one, two, three, angular_tolerance => eps
      use wigner_6j_triangle_selection, only: ixjtik
      use racah_coefficients, only: dracah
      use wigner_6j_with_j_two, only: sixj2
      use wigner_6j_with_j_three, only: sixj3
      implicit none (type, external)
      integer, intent(in) :: two_j_1, two_j_2, two_j_3, two_j_4, two_j_5
      integer, intent(in) :: check_triangle_conditions
      real(real64), intent(out) :: six_j_value
      real(real64) :: j6_for_rank_3, j1, j2, j3, j4, j5, rank_4_weight, rank_2_weight, rank_3_weight, rank_2_six_j, &
          rank_3_six_j
      six_j_value = zero
      if (check_triangle_conditions /= 0) then
         if (ixjtik(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,8) == 0) return
      end if

      ! Use the Racah sum if the 6j symbol with sixth argument 3 is inadmissible.
      if (ixjtik(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,6) == 0) then
         call dracah (two_j_1, two_j_2, two_j_5, two_j_4, two_j_3, 8, six_j_value)
         if (mod(two_j_1 + two_j_2 + two_j_5 + two_j_4,4) /= 0) six_j_value = -six_j_value
      else
         j6_for_rank_3 = three
         j1 = dble(two_j_1)*half
         j2 = dble(two_j_2)*half
         j3 = dble(two_j_3)*half
         j4 = dble(two_j_4)*half
         j5 = dble(two_j_5)*half
         rank_4_weight = j6_for_rank_3*dsqrt((j6_for_rank_3+j4+j2+two)*(j6_for_rank_3-j4+j2+one)*(j6_for_rank_3+j4- &
             j2+one)*((-         &
            j6_for_rank_3)+j4+j2)*(j6_for_rank_3+j1+j5+two)*(j6_for_rank_3-j1+j5+one)*(j6_for_rank_3+j1-j5+one)* &
                (             &
            (-j6_for_rank_3)+j1+j5))
         rank_2_weight = (j6_for_rank_3+one)*dsqrt((j6_for_rank_3+j4+j2+one)*(j6_for_rank_3-j4+j2)*(j6_for_rank_3+ &
             j4-j2)*((-j6_for_rank_3)         &
             + j4+j2+one)*(j6_for_rank_3+j1+j5+one)*(j6_for_rank_3-j1+j5)*(j6_for_rank_3+j1-j5)*((-j6_for_rank_3)+ &
                 j1           &
             + j5+one))
         rank_3_weight = (two*j6_for_rank_3+one)*(two*(j6_for_rank_3*(j6_for_rank_3+one)*j3*(j3+one)-j4*(j4+one)* &
             j1*(j1+   &
            one)-j2*(j2+one)*j5*(j5+one))+(j6_for_rank_3*(j6_for_rank_3+one)-j4*(j4+one)-j2*(j2       &
             +one))*(j6_for_rank_3*(j6_for_rank_3+one)-j1*(j1+one)-j5*(j5+one)))
         if (dabs(rank_2_weight) < angular_tolerance) then
            rank_2_six_j = zero
         else
            call sixj2 (two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, 0, rank_2_six_j)
         end if
         if (dabs(rank_3_weight) < angular_tolerance) then
            rank_3_six_j = zero
         else
            call sixj3 (two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, 0, rank_3_six_j)
         end if
         six_j_value = (rank_3_weight*rank_3_six_j - rank_2_weight*rank_2_six_j)/rank_4_weight
      end if
      return
      end subroutine sixj4
end module wigner_6j_with_j_four
