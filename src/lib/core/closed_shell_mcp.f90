! Generate angular coefficients for direct and exchange Coulomb terms involving closed subshells.
module closed_shell_mcp
    implicit none (type, external)
    private
    public :: cord
contains
      subroutine cord(bra_csf_index, ket_csf_index, bra_shell_position, orbital_list_selector, ket_shell_position)
      use, intrinsic :: iso_fortran_env, only: real64
      use orbital_csf_basis, only: max_orbitals, orbital_kappa
      use csf_pair_angular_state, only: jjc1, jjc2, jjq1, jjq2, peel_orbital_indices => jlist, core_orbital_indices => klist, &
         & max_orbitals, ncore, npeel, bra_occupations => nq1, nq2
      use mcp_coefficient_output, only: speak
      use relativistic_three_j, only: clrx
      implicit none (type, external)
      integer  :: bra_csf_index
      integer  :: ket_csf_index
      integer, intent(in) :: bra_shell_position
      integer, intent(in) :: orbital_list_selector
      integer, intent(in) :: ket_shell_position
      real(real64), parameter :: coefficient_cutoff = 1.0e-10_real64
      integer :: first_orbital,second_orbital,swap_orbital,first_kappa,first_abs_kappa,first_occupation,maximum_rank, &
         & tensor_rank,second_kappa,second_abs_kappa,second_occupation,minimum_rank
      real(real64) :: mcp_coefficient,occupation_factor,three_j_symbol
      if (orbital_list_selector == 2) then
         first_orbital = core_orbital_indices(bra_shell_position)
      else
         first_orbital = peel_orbital_indices(bra_shell_position)
      end if
      second_orbital = core_orbital_indices(ket_shell_position)
      if (first_orbital > second_orbital) then
         swap_orbital = first_orbital
         first_orbital = second_orbital
         second_orbital = swap_orbital
      end if
      first_kappa = orbital_kappa(first_orbital)
      first_abs_kappa = iabs(first_kappa)
      first_occupation = bra_occupations(first_orbital)
      ! Equivalent electrons use pair counting for the monopole and only even higher ranks.
      if (first_orbital == second_orbital) then
         mcp_coefficient = dble(first_occupation*(first_occupation - 1)/2)
         call speak (bra_csf_index, ket_csf_index, first_orbital, second_orbital, first_orbital, second_orbital, 0, &
            & mcp_coefficient)
         maximum_rank = first_abs_kappa + first_abs_kappa - 2
         if (maximum_rank <= 0) return
         occupation_factor = dble(first_occupation*first_occupation/2)
         do tensor_rank = 2, maximum_rank, 2
            three_j_symbol = clrx(first_kappa,tensor_rank,first_kappa)
            mcp_coefficient = -occupation_factor*three_j_symbol*three_j_symbol
            if (abs(mcp_coefficient) < coefficient_cutoff) cycle
            call speak (bra_csf_index, ket_csf_index, first_orbital, second_orbital, first_orbital, second_orbital, &
               & tensor_rank, mcp_coefficient)
         end do
      else
         second_kappa = orbital_kappa(second_orbital)
         second_abs_kappa = abs(second_kappa)
         second_occupation = bra_occupations(second_orbital)
         occupation_factor = dble(first_occupation*second_occupation)
         call speak (bra_csf_index, ket_csf_index, first_orbital, second_orbital, first_orbital, second_orbital, 0, &
            & occupation_factor)
         ! The exchange ranks start at the triangle bound with the parity required by the orbital pair.
         minimum_rank = abs(first_abs_kappa - second_abs_kappa)
         maximum_rank = first_abs_kappa + second_abs_kappa - 1
         if (first_kappa*second_kappa < 0) minimum_rank = minimum_rank + 1
         do tensor_rank = minimum_rank, maximum_rank, 2
            three_j_symbol = clrx(first_kappa,tensor_rank,second_kappa)
            mcp_coefficient = -occupation_factor*three_j_symbol*three_j_symbol
            if (abs(mcp_coefficient) < coefficient_cutoff) cycle
            call speak (bra_csf_index, ket_csf_index, first_orbital, second_orbital, second_orbital, first_orbital, &
               & tensor_rank, mcp_coefficient)
         end do
      end if
      return
      end subroutine cord
end module closed_shell_mcp
