! Find the first allowed rank and rank count for the direct and exchange Coulomb terms.
module coulomb_tensor_rank_range
    implicit none (type, external)
    private
    public :: skrc
contains
      subroutine skrc(orbital_indices, orbital_kappas, orbital_two_j_plus_one, first_direct_rank, direct_rank_count, &
         & first_exchange_rank, exchange_rank_count)
      implicit none (type, external)
      integer, intent(out) :: first_direct_rank
      integer, intent(out) :: direct_rank_count
      integer, intent(out) :: first_exchange_rank
      integer, intent(out) :: exchange_rank_count
      integer, intent(in) :: orbital_indices(4)
      integer, intent(in) :: orbital_kappas(4)
      integer, intent(in) :: orbital_two_j_plus_one(4)
      integer :: direct_kappa_sign_a, direct_kappa_sign_b, twice_direct_lower_a, twice_direct_lower_b, &
         & twice_direct_upper_a, twice_direct_upper_b, exchange_kappa_sign_a, exchange_kappa_sign_b, twice_exchange_lower_a, &
         twice_exchange_lower_b, twice_exchange_upper_a, twice_exchange_upper_b
      direct_rank_count = 0
      exchange_rank_count = 0
      direct_kappa_sign_a = 1
      if (orbital_kappas(1)*orbital_kappas(3) < 0) direct_kappa_sign_a = -1
      direct_kappa_sign_b = 1
      if (orbital_kappas(2)*orbital_kappas(4) < 0) direct_kappa_sign_b = -1
      twice_direct_lower_a = abs(orbital_two_j_plus_one(1)-orbital_two_j_plus_one(3))
      if (direct_kappa_sign_a < 0) twice_direct_lower_a = twice_direct_lower_a + 2
      twice_direct_lower_b = abs(orbital_two_j_plus_one(2)-orbital_two_j_plus_one(4))
      if (direct_kappa_sign_b < 0) twice_direct_lower_b = twice_direct_lower_b + 2
      if (mod((twice_direct_lower_a - twice_direct_lower_b)/2,2) == 0) then
         twice_direct_upper_a = orbital_two_j_plus_one(1) + orbital_two_j_plus_one(3) - 2
         if (direct_kappa_sign_a > 0) twice_direct_upper_a = twice_direct_upper_a - 2
         twice_direct_upper_b = orbital_two_j_plus_one(2) + orbital_two_j_plus_one(4) - 2
         if (direct_kappa_sign_b > 0) twice_direct_upper_b = twice_direct_upper_b - 2
         first_direct_rank = max(twice_direct_lower_a,twice_direct_lower_b)/2
         direct_rank_count = min(twice_direct_upper_a,twice_direct_upper_b)/2
         direct_rank_count = (direct_rank_count - first_direct_rank)/2 + 1
      end if
      if (orbital_indices(1)==orbital_indices(2) .or. orbital_indices(3)==orbital_indices(4)) return
      exchange_kappa_sign_a = 1
      if (orbital_kappas(1)*orbital_kappas(4) < 0) exchange_kappa_sign_a = -1
      exchange_kappa_sign_b = 1
      if (orbital_kappas(2)*orbital_kappas(3) < 0) exchange_kappa_sign_b = -1
      twice_exchange_lower_a = abs(orbital_two_j_plus_one(1)-orbital_two_j_plus_one(4))
      if (exchange_kappa_sign_a < 0) twice_exchange_lower_a = twice_exchange_lower_a + 2
      twice_exchange_lower_b = abs(orbital_two_j_plus_one(2)-orbital_two_j_plus_one(3))
      if (exchange_kappa_sign_b < 0) twice_exchange_lower_b = twice_exchange_lower_b + 2
      if (mod((twice_exchange_lower_a - twice_exchange_lower_b)/2,2) /= 0) return
      twice_exchange_upper_a = orbital_two_j_plus_one(1) + orbital_two_j_plus_one(4) - 2
      if (exchange_kappa_sign_a > 0) twice_exchange_upper_a = twice_exchange_upper_a - 2
      twice_exchange_upper_b = orbital_two_j_plus_one(2) + orbital_two_j_plus_one(3) - 2
      if (exchange_kappa_sign_b > 0) twice_exchange_upper_b = twice_exchange_upper_b - 2
      first_exchange_rank = max(twice_exchange_lower_a,twice_exchange_lower_b)/2
      exchange_rank_count = min(twice_exchange_upper_a,twice_exchange_upper_b)/2
      exchange_rank_count = (exchange_rank_count - first_exchange_rank)/2 + 1
      return
      end subroutine skrc
end module coulomb_tensor_rank_range
