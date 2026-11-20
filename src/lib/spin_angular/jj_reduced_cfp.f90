! Look up a completely reduced creation or annihilation matrix element.
module jj_reduced_cfp
    implicit none (type, external)
    private
    public :: rmeajj
contains
      subroutine rmeajj(two_orbital_j,bra_term_index,two_bra_q,two_bra_term_j,ket_term_index,two_ket_q, &
          two_ket_term_j,reduced_cfp)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants,          only: zero
      use jj_reduced_coefficients, only: last_opposite_parity_terms_low_j => imgnjj, last_same_parity_terms_low_j => imgtjj, &
          first_opposite_parity_terms_low_j => impnjj, first_same_parity_terms_low_j => imptjj
      use jj_reduced_cfp_two_electron_terms, only: rmeajj11
      use jj_reduced_cfp_j_nine_halves, only: rmeajj9
      implicit none (type, external)
      integer,      intent(in)  :: two_orbital_j, bra_term_index, two_bra_q, two_bra_term_j, ket_term_index, &
          two_ket_q, two_ket_term_j
      real(real64), intent(out) :: reduced_cfp
      integer :: cfp_row, cfp_column, twice_phase_exponent,  lower_term_index, upper_term_index
      integer :: cfp_numerators_two_j_3(1, 2), cfp_denominators_two_j_3(1, 2)
      integer :: cfp_numerators_two_j_5(3, 3), cfp_denominators_two_j_5(3, 3)
      integer :: cfp_numerators_two_j_7(6, 8), cfp_denominators_two_j_7(6, 8)
      data cfp_numerators_two_j_3/12,20/
      data cfp_denominators_two_j_3/1,1/
      data cfp_numerators_two_j_5/-24,2*0,-30,120,-90,-54,-48,330/
      data cfp_denominators_two_j_5/4*1,2*7,1,2*7/
      data cfp_numerators_two_j_7/2*0,-40,3*0,-54,-33,-40,-65,30,0,88,-1,0,-1755,-40,  &
      0,198,-36,-72,4500,-234,-360,-78,156,0,108,-30,96,66,49,0,27,  &
      130,136,0,195,-104,-245,-624,1224,3*0,2720,-2448,-7752/
      data cfp_denominators_two_j_7/6*1,7,2*1,7,2*1,7,2*1,77,11,1,7,11,1,77,2*11,35,5,1, &
      91,1,13,2*5,2*1,2*7,1,11,1,3*11,3*1,143,77,91/
      reduced_cfp=zero
      if(two_orbital_j > 37) return
      if(two_orbital_j <= 9) then
        if(bra_term_index < 300) then
          if(first_same_parity_terms_low_j(bra_term_index) /= first_opposite_parity_terms_low_j(ket_term_index)) &
              return
          if(bra_term_index < ket_term_index) then
            lower_term_index=bra_term_index
            upper_term_index=ket_term_index
          else
            lower_term_index=ket_term_index
            upper_term_index=bra_term_index
          end if
        end if
      end if

      ! Recover the reduced element from tabulated signed squares and denominators.
      if(two_orbital_j == 1) then
        reduced_cfp=-2
      else if(two_orbital_j == 3)then
        cfp_row=lower_term_index-2
        cfp_column=upper_term_index-3
        reduced_cfp=-dsqrt(dble(cfp_numerators_two_j_3(cfp_row,cfp_column))/dble(cfp_denominators_two_j_3(cfp_row, &
            cfp_column)))
      else if(two_orbital_j == 5)then
        cfp_row=lower_term_index-5
        cfp_column=upper_term_index-8
        if(cfp_numerators_two_j_5(cfp_row,cfp_column) >= 0) then
          reduced_cfp=dsqrt(dble(cfp_numerators_two_j_5(cfp_row,cfp_column))/dble(cfp_denominators_two_j_5(cfp_row, &
              cfp_column)))
        else
          reduced_cfp=-dsqrt(-dble(cfp_numerators_two_j_5(cfp_row,cfp_column))/ &
              dble(cfp_denominators_two_j_5(cfp_row,cfp_column)))
        end if
      else if(two_orbital_j == 7)then
        cfp_row=lower_term_index-11
        cfp_column=upper_term_index-17
        if(cfp_numerators_two_j_7(cfp_row,cfp_column) >= 0) then
          reduced_cfp=dsqrt(dble(cfp_numerators_two_j_7(cfp_row,cfp_column))/dble(cfp_denominators_two_j_7(cfp_row, &
              cfp_column)))
        else
          reduced_cfp=-dsqrt(-dble(cfp_numerators_two_j_7(cfp_row,cfp_column))/ &
              dble(cfp_denominators_two_j_7(cfp_row,cfp_column)))
        end if
      else if(two_orbital_j == 9) then
        if(bra_term_index > 300) then
          if(ket_term_index < 300) then
            write(6,"(A)") " ERROR IN  RMEAJJ "
            stop
          end if
          if(two_orbital_j == two_ket_term_j) then
            call rmeajj11(bra_term_index,ket_term_index,two_orbital_j,reduced_cfp)
          else
            call rmeajj11(ket_term_index,bra_term_index,two_orbital_j,reduced_cfp)
            twice_phase_exponent=two_bra_q-two_ket_q+two_bra_term_j-two_ket_term_j+two_orbital_j-1
            if((twice_phase_exponent/4)*4 /= twice_phase_exponent)reduced_cfp=-reduced_cfp
          end if
        else
          call rmeajj9(bra_term_index,two_bra_q,two_bra_term_j,ket_term_index,two_ket_q,two_ket_term_j,reduced_cfp)
          if(bra_term_index > ket_term_index) then
            if(mod(two_bra_q+two_bra_term_j-two_ket_q-two_ket_term_j+two_orbital_j-1,4) /= 0) reduced_cfp=-reduced_cfp
          end if
          write(0,"(A)") " KLAIDA RMEAJJ SUB. "
          stop
        end if
      else
        if(two_orbital_j == two_ket_term_j) then
          call rmeajj11(bra_term_index,ket_term_index,two_orbital_j,reduced_cfp)
        else
          call rmeajj11(ket_term_index,bra_term_index,two_orbital_j,reduced_cfp)
          twice_phase_exponent=two_bra_q-two_ket_q+two_bra_term_j-two_ket_term_j+two_orbital_j-1
          if((twice_phase_exponent/4)*4 /= twice_phase_exponent)reduced_cfp=-reduced_cfp
        end if
      end if
      if(two_orbital_j < 9) then
        if(bra_term_index > ket_term_index) then
          if(mod(two_bra_q+two_bra_term_j-two_ket_q-two_ket_term_j+two_orbital_j-1,4) /= 0) reduced_cfp=-reduced_cfp
        end if
      end if
      return
      end subroutine rmeajj
end module jj_reduced_cfp
