! Extract the occupations and coupling chains of a CSF pair into the angular workspace.
module csf_pair_quantum_numbers
    implicit none (type, external)
    private
    public :: setqna
contains
      subroutine setqna(bra_csf_index, ket_csf_index)
      use csf_pair_angular_state, only: bra_intermediate_couplings => jjc1, ket_intermediate_couplings => jjc2, &
         & bra_subshell_quantum_codes => jjq1, ket_subshell_quantum_codes => jjq2, peel_orbital_indices => jlist, &
         & core_orbital_indices => klist, max_orbitals, core_orbital_count => ncore, peel_orbital_count => npeel, &
         & bra_occupations => nq1, ket_occupations => nq2
      use orbital_csf_basis,           only: csf_count, orbital_count
      use csf_quantum_queries, only: iq, jqs, jcup
      use subshell_occupation_class, only: ichop
      implicit none (type, external)
      integer, intent(in) :: bra_csf_index
      integer, intent(in) :: ket_csf_index
      integer :: orbital_or_peel_index, quantum_number_selector, next_coupling_index, open_shell_count, &
         & first_peel_orbital, second_peel_orbital, peel_orbital
      do orbital_or_peel_index = 1, orbital_count
         bra_occupations(orbital_or_peel_index) = iq(orbital_or_peel_index,bra_csf_index)
         ket_occupations(orbital_or_peel_index) = iq(orbital_or_peel_index,ket_csf_index)
         do quantum_number_selector = 1, 3
            bra_subshell_quantum_codes(quantum_number_selector,orbital_or_peel_index) = jqs(quantum_number_selector, &
               & orbital_or_peel_index,bra_csf_index)
            ket_subshell_quantum_codes(quantum_number_selector,orbital_or_peel_index) = jqs(quantum_number_selector, &
               & orbital_or_peel_index,ket_csf_index)
         end do
      end do
      peel_orbital_count = 0
      core_orbital_count = 0
      do orbital_or_peel_index = 1, orbital_count
         if (ichop(orbital_or_peel_index,bra_csf_index)==(-1) .and. ichop(orbital_or_peel_index,ket_csf_index)==(-1)) cycle
         if (ichop(orbital_or_peel_index,bra_csf_index)==1 .and. ichop(orbital_or_peel_index,ket_csf_index)==1) then
            core_orbital_count = core_orbital_count + 1
            core_orbital_indices(core_orbital_count) = orbital_or_peel_index
         else
            peel_orbital_count = peel_orbital_count + 1
            peel_orbital_indices(peel_orbital_count) = orbital_or_peel_index
         end if
      end do
      if (peel_orbital_count <= 1) return
      next_coupling_index = 1
      open_shell_count = 0
      first_peel_orbital = peel_orbital_indices(1)
      second_peel_orbital = peel_orbital_indices(2)
      if (ichop(first_peel_orbital,bra_csf_index) /= 0) then
         bra_intermediate_couplings(1) = jqs(3,second_peel_orbital,bra_csf_index)
         if (ichop(second_peel_orbital,bra_csf_index) == 0) open_shell_count = 1
      else
         open_shell_count = 1
         if (ichop(second_peel_orbital,bra_csf_index) == 0) then
            bra_intermediate_couplings(1) = jcup(next_coupling_index,bra_csf_index)
            next_coupling_index = next_coupling_index + 1
         else
            bra_intermediate_couplings(1) = jqs(3,first_peel_orbital,bra_csf_index)
         end if
      end if
      do orbital_or_peel_index = 3, peel_orbital_count
         peel_orbital = peel_orbital_indices(orbital_or_peel_index)
         if (ichop(peel_orbital,bra_csf_index) /= 0) then
            bra_intermediate_couplings(orbital_or_peel_index-1) = bra_intermediate_couplings(orbital_or_peel_index-2)
         else
            if (open_shell_count /= 0) then
               bra_intermediate_couplings(orbital_or_peel_index-1) = jcup(next_coupling_index,bra_csf_index)
               next_coupling_index = next_coupling_index + 1
            else
               bra_intermediate_couplings(orbital_or_peel_index-1) = jqs(3,peel_orbital,bra_csf_index)
            end if
            open_shell_count = open_shell_count + 1
         end if
      end do
      next_coupling_index = 1
      open_shell_count = 0
      first_peel_orbital = peel_orbital_indices(1)
      second_peel_orbital = peel_orbital_indices(2)
      if (ichop(first_peel_orbital,ket_csf_index) /= 0) then
         ket_intermediate_couplings(1) = jqs(3,second_peel_orbital,ket_csf_index)
         if (ichop(second_peel_orbital,ket_csf_index) == 0) open_shell_count = 1
      else
         open_shell_count = 1
         if (ichop(second_peel_orbital,ket_csf_index) == 0) then
            ket_intermediate_couplings(1) = jcup(next_coupling_index,ket_csf_index)
            next_coupling_index = next_coupling_index + 1
         else
            ket_intermediate_couplings(1) = jqs(3,first_peel_orbital,ket_csf_index)
         end if
      end if
      do orbital_or_peel_index = 3, peel_orbital_count
         peel_orbital = peel_orbital_indices(orbital_or_peel_index)
         if (ichop(peel_orbital,ket_csf_index) /= 0) then
            ket_intermediate_couplings(orbital_or_peel_index-1) = ket_intermediate_couplings(orbital_or_peel_index-2)
         else
            if (open_shell_count /= 0) then
               ket_intermediate_couplings(orbital_or_peel_index-1) = jcup(next_coupling_index,ket_csf_index)
               next_coupling_index = next_coupling_index + 1
            else
               ket_intermediate_couplings(orbital_or_peel_index-1) = jqs(3,peel_orbital,ket_csf_index)
            end if
            open_shell_count = open_shell_count + 1
         end if
      end do
      return
      end subroutine setqna
end module csf_pair_quantum_numbers
