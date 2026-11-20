! Decode the quasispin, seniority and J of a tabulated subshell term.
module jj_subshell_term_decoding
    implicit none (type, external)
    private
    public :: rumtjj
contains
      subroutine rumtjj(term_index,two_orbital_j,two_q,two_seniority,two_term_j)
      use subshell_term_multiplicities, only: packed_term_codes_low_j => mt
      use high_j_subshell_term_data, only: packed_term_codes_high_j => mt11, packed_term_codes_two_j_9 => mt9
      use packed_integer_digit, only: jthn
      implicit none (type, external)
      integer, intent(in)  :: two_orbital_j, term_index
      integer, intent(out) :: two_q, two_seniority, two_term_j
      integer :: packed_term_code, unoffset_term_index

      ! Select the term table from 2j and decode the offset used for the few-electron j = 9/2 terms.
      if(two_orbital_j < 9) then
        packed_term_code=packed_term_codes_low_j(term_index)
      else if(two_orbital_j == 9) then
        if(term_index > 300) then
          unoffset_term_index=term_index-300
          packed_term_code=packed_term_codes_two_j_9(unoffset_term_index)
        else
          print*, "ERROR in RUMTJJ"
          stop
        end if
      else
        packed_term_code=packed_term_codes_high_j(term_index)
      end if
      two_q=jthn(packed_term_code,3,100)
      two_seniority=jthn(packed_term_code,2,100)
      two_term_j=jthn(packed_term_code,1,100)
      return
      end subroutine rumtjj
end module jj_subshell_term_decoding
