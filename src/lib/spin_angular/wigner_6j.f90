! Evaluate a Wigner 6j symbol using symmetry, special cases and a small cache.
module wigner_6j
    implicit none (type, external)
    private
    public :: sixj
contains
      subroutine sixj(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6,check_triangle_conditions,six_j_value)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants,          only: zero, one
      use wigner_6j_triangle_selection, only: ixjtik
      use wigner_6j_with_j_one_half, only: sixj5
      use wigner_6j_with_j_one, only: sixj1
      use wigner_6j_with_j_three_halves, only: sixj35
      use wigner_6j_with_j_two, only: sixj2
      use racah_coefficients, only: dracah
      use wigner_6j_with_j_three, only: sixj3
      use wigner_6j_with_j_four, only: sixj4
      implicit none (type, external)
      integer, intent(in) :: two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, two_j_6
      integer, intent(in) :: check_triangle_conditions
      real(real64), intent(out) :: six_j_value
      integer :: twice_phase_exponent
      real(real64) :: small_six_j_cache(0:4, 0:4, 0:4, 0:4, 0:4, 0:4)
      real(real64) :: uncached_six_j_sentinel, normalization_squared
      logical :: store_in_six_j_cache
      save :: small_six_j_cache
      data small_six_j_cache/ 15625*1.0e-20_real64/
      data uncached_six_j_sentinel/ 1.0e-20_real64/
!$OMP THREADPRIVATE(small_six_j_cache)
      six_j_value = zero
      if (check_triangle_conditions /= 0) then
         if (ixjtik(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 0) return
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
         if (mod(twice_phase_exponent,4) /= 0) six_j_value = -six_j_value

      ! Use 6j symmetries to move the small argument to the sixth slot.
      else if (min0(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 1) then
         if (two_j_1 == 1) then
            call sixj5 (two_j_5, two_j_3, two_j_4, two_j_2, two_j_6, 0, six_j_value)
         else if (two_j_2 == 1) then
            call sixj5 (two_j_1, two_j_6, two_j_5, two_j_4, two_j_3, 0, six_j_value)
         else if (two_j_3 == 1) then
            call sixj5 (two_j_1, two_j_5, two_j_6, two_j_4, two_j_2, 0, six_j_value)
         else if (two_j_4 == 1) then
            call sixj5 (two_j_2, two_j_3, two_j_1, two_j_5, two_j_6, 0, six_j_value)
         else if (two_j_5 == 1) then
            call sixj5 (two_j_1, two_j_3, two_j_2, two_j_4, two_j_6, 0, six_j_value)
         else
            call sixj5 (two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, 0, six_j_value)
         end if
      else if (min0(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 2) then
         if (two_j_1 == 2) then
            call sixj1 (two_j_5, two_j_3, two_j_4, two_j_2, two_j_6, 0, six_j_value)
         else if (two_j_2 == 2) then
            call sixj1 (two_j_1, two_j_6, two_j_5, two_j_4, two_j_3, 0, six_j_value)
         else if (two_j_3 == 2) then
            call sixj1 (two_j_1, two_j_5, two_j_6, two_j_4, two_j_2, 0, six_j_value)
         else if (two_j_4 == 2) then
            call sixj1 (two_j_2, two_j_3, two_j_1, two_j_5, two_j_6, 0, six_j_value)
         else if (two_j_5 == 2) then
            call sixj1 (two_j_1, two_j_3, two_j_2, two_j_4, two_j_6, 0, six_j_value)
         else
            call sixj1 (two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, 0, six_j_value)
         end if
      else if (min0(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 3) then
         if (two_j_1 == 3) then
            call sixj35 (two_j_5, two_j_3, two_j_4, two_j_2, two_j_6, 0, six_j_value)
         else if (two_j_2 == 3) then
            call sixj35 (two_j_1, two_j_6, two_j_5, two_j_4, two_j_3, 0, six_j_value)
         else if (two_j_3 == 3) then
            call sixj35 (two_j_1, two_j_5, two_j_6, two_j_4, two_j_2, 0, six_j_value)
         else if (two_j_4 == 3) then
            call sixj35 (two_j_2, two_j_3, two_j_1, two_j_5, two_j_6, 0, six_j_value)
         else if (two_j_5 == 3) then
            call sixj35 (two_j_1, two_j_3, two_j_2, two_j_4, two_j_6, 0, six_j_value)
         else
            call sixj35 (two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, 0, six_j_value)
         end if
      else if (min0(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 4) then
         if (two_j_1 == 4) then
            call sixj2 (two_j_5, two_j_3, two_j_4, two_j_2, two_j_6, 0, six_j_value)
         else if (two_j_2 == 4) then
            call sixj2 (two_j_1, two_j_6, two_j_5, two_j_4, two_j_3, 0, six_j_value)
         else if (two_j_3 == 4) then
            call sixj2 (two_j_1, two_j_5, two_j_6, two_j_4, two_j_2, 0, six_j_value)
         else if (two_j_4 == 4) then
            call sixj2 (two_j_2, two_j_3, two_j_1, two_j_5, two_j_6, 0, six_j_value)
         else if (two_j_5 == 4) then
            call sixj2 (two_j_1, two_j_3, two_j_2, two_j_4, two_j_6, 0, six_j_value)
         else
            call sixj2 (two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, 0, six_j_value)
         end if
      else if (min0(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 5) then
         call dracah (two_j_1, two_j_2, two_j_5, two_j_4, two_j_3, two_j_6, six_j_value)
         if (mod(two_j_1 + two_j_2 + two_j_5 + two_j_4,4) /= 0) six_j_value = -six_j_value
      else if (min0(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 6) then
         if (two_j_1 == 6) then
            call sixj3 (two_j_5, two_j_3, two_j_4, two_j_2, two_j_6, 0, six_j_value)
         else if (two_j_2 == 6) then
            call sixj3 (two_j_1, two_j_6, two_j_5, two_j_4, two_j_3, 0, six_j_value)
         else if (two_j_3 == 6) then
            call sixj3 (two_j_1, two_j_5, two_j_6, two_j_4, two_j_2, 0, six_j_value)
         else if (two_j_4 == 6) then
            call sixj3 (two_j_2, two_j_3, two_j_1, two_j_5, two_j_6, 0, six_j_value)
         else if (two_j_5 == 6) then
            call sixj3 (two_j_1, two_j_3, two_j_2, two_j_4, two_j_6, 0, six_j_value)
         else
            call sixj3 (two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, 0, six_j_value)
         end if
      else if (min0(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 7) then
         call dracah (two_j_1, two_j_2, two_j_5, two_j_4, two_j_3, two_j_6, six_j_value)
         if (mod(two_j_1 + two_j_2 + two_j_5 + two_j_4,4) /= 0) six_j_value = -six_j_value
      else if (min0(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 8) then
         if (two_j_1 == 8) then
            call sixj4 (two_j_5, two_j_3, two_j_4, two_j_2, two_j_6, 0, six_j_value)
         else if (two_j_2 == 8) then
            call sixj4 (two_j_1, two_j_6, two_j_5, two_j_4, two_j_3, 0, six_j_value)
         else if (two_j_3 == 8) then
            call sixj4 (two_j_1, two_j_5, two_j_6, two_j_4, two_j_2, 0, six_j_value)
         else if (two_j_4 == 8) then
            call sixj4 (two_j_2, two_j_3, two_j_1, two_j_5, two_j_6, 0, six_j_value)
         else if (two_j_5 == 8) then
            call sixj4 (two_j_1, two_j_3, two_j_2, two_j_4, two_j_6, 0, six_j_value)
         else
            call sixj4 (two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, 0, six_j_value)
         end if
      else if (min0(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) == 9) then
         call dracah (two_j_1, two_j_2, two_j_5, two_j_4, two_j_3, two_j_6, six_j_value)
         if (mod(two_j_1 + two_j_2 + two_j_5 + two_j_4,4) /= 0) six_j_value = -six_j_value
      else
       call dracah(two_j_1,two_j_2,two_j_5,two_j_4,two_j_3,two_j_6,six_j_value)
         if (mod(two_j_1 + two_j_2 + two_j_5 + two_j_4,4) /= 0) six_j_value = -six_j_value
      end if
      if (store_in_six_j_cache) small_six_j_cache(two_j_1,two_j_2,two_j_3,two_j_4,two_j_5,two_j_6) = six_j_value
      return
      end subroutine sixj
end module wigner_6j
