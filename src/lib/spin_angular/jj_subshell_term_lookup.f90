! Find the tabulated term that matches the given subshell quantum numbers.
module jj_subshell_term_lookup
    implicit none (type, external)
    private
    public :: nmtejj
contains
      integer function nmtejj(two_q,two_term_j,two_orbital_j,ket_occupation,bra_occupation) result(term_index)
      use subshell_term_multiplicities, only: packed_term_codes_low_j => mt
      use high_j_subshell_term_data, only: packed_term_codes_high_j => mt11, packed_term_codes_two_j_9 => mt9
      implicit none (type, external)
      integer, intent(in) :: two_q, two_term_j, two_orbital_j, ket_occupation, bra_occupation
      integer                :: two_seniority, packed_term_code, high_j_range_index, first_term_index, &
          last_term_index, table_index
      integer  :: low_j_first_terms(9), low_j_last_terms(9)
      integer :: high_j_first_terms(27), high_j_last_terms(27)
      data low_j_first_terms/1,0,3,0,6,0,12,0,26/
      data low_j_last_terms/2,0,5,0,11,0,25,0,63/
      data high_j_first_terms/1,0,8,0,16,0,25,0,35,0,46,0,58,0,71,0,85,0,100,0, &
      116,0,133,0,151,0,170/
      data high_j_last_terms/7,0,15,0,24,0,34,0,45,0,57,0,70,0,84,0,99,0,115,0,&
      132,0,150,0,169,0,189/
      term_index=0
      if(two_orbital_j > 37)return
      two_seniority=two_orbital_j+1-2*two_q

      ! Use three base-100 fields for 2Q, twice seniority and 2J.
      packed_term_code=(two_q*100+two_seniority)*100+two_term_j
      if(two_orbital_j < 9) then
         first_term_index=low_j_first_terms(two_orbital_j)
         if(first_term_index == 0)return
         last_term_index=low_j_last_terms(two_orbital_j)
         if(last_term_index == 0)return
         table_index=first_term_index
    2    if (packed_term_code == packed_term_codes_low_j(table_index)) goto 1
    3    table_index=table_index+1
         if(table_index > last_term_index)return
         goto 2
    1    term_index=table_index
      else if(two_orbital_j == 9) then
        if(max0(ket_occupation,bra_occupation) < 3) then
          first_term_index=1
          last_term_index=6
          table_index=first_term_index
    6     if (packed_term_code == packed_term_codes_two_j_9(table_index)) goto 5
    4     table_index=table_index+1
          if(table_index > last_term_index)return
          goto 6
    5     term_index=table_index+300
        else
          print*, "ERROR in FUNCTION NMTEJJ"
          stop
        end if
      else
        high_j_range_index=two_orbital_j-10
        first_term_index=high_j_first_terms(high_j_range_index)
        last_term_index=high_j_last_terms(high_j_range_index)
        table_index=first_term_index
   22   if (packed_term_code == packed_term_codes_high_j(table_index)) goto 21
   23   table_index=table_index+1
        if(table_index > last_term_index)return
        goto 22
   21   term_index=table_index
      end if
      return
      end function nmtejj
end module jj_subshell_term_lookup
