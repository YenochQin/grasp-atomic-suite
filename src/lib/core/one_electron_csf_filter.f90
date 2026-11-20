! Check whether CSF occupations can be connected by a one-electron operator.
module one_electron_csf_filter
    implicit none (type, external)
    private
    public :: ichkq1
contains
      integer function ichkq1 (bra_csf_index, ket_csf_index)
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use orbital_csf_basis, only: orbital_count
      use csf_quantum_queries, only: iq
      implicit none (type, external)
      integer, intent(in) :: bra_csf_index
      integer, intent(in) :: ket_csf_index
      integer :: differing_subshell_count, orbital_index, bra_occupation, ket_occupation

      ichkq1 = 0
      differing_subshell_count = 0
      do orbital_index = 1, orbital_count
         bra_occupation = iq(orbital_index,bra_csf_index)
         ket_occupation = iq(orbital_index,ket_csf_index)
         if (bra_occupation == ket_occupation) cycle
         differing_subshell_count = differing_subshell_count + 1
         if (differing_subshell_count > 2) return
         if (iabs(bra_occupation - ket_occupation) <= 1) cycle
         return
      end do
      if (differing_subshell_count==2 .or. differing_subshell_count==0) ichkq1 = 1
      return
      end function ichkq1
end module one_electron_csf_filter
