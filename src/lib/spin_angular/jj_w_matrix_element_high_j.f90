! Sum two single-tensor elements over high-j intermediate subshell terms.
module jj_w_matrix_element_high_j
    implicit none (type, external)
    private
    public :: w1jjg
contains
      subroutine w1jjg(w_angular_rank,q_projection_1,q_projection_2,bra_shell_codes,bra_shell_values, &
          ket_shell_codes,ket_shell_values,w_matrix_element)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants,          only: zero, tenth, half, two, angular_tolerance => eps
      use jj_reduced_coefficients, only: last_opposite_parity_terms_low_j => imgnjj, last_same_parity_terms_low_j => imgtjj, &
          first_opposite_parity_terms_low_j => impnjj, first_same_parity_terms_low_j => imptjj
      use jj_reduced_coefficients_9, only: last_opposite_parity_terms_two_j_9 => imgnjj9, &
          last_same_parity_terms_two_j_9 => imgtjj9, first_opposite_parity_terms_two_j_9 => impnjj9, &
          first_same_parity_terms_two_j_9 => imptjj9
      use jj_reduced_coefficients_11, only: last_opposite_parity_terms_high_j => imgnjj11, &
          last_same_parity_terms_high_j => imgtjj11, first_opposite_parity_terms_high_j => impnjj11, &
          first_same_parity_terms_high_j => imptjj11
      use subshell_tensor_diagnostics, only: mes
      use jj_subshell_term_decoding, only: rumtjj
      use wigner_6j_triangle_selection, only: ixjtik
      use angular_triangle_and_parity, only: ittk
      use jj_a_matrix_element_few_electrons, only: a1jj
      use wigner_6j, only: sixj
      use quasispin_projection_selection, only: izas1
      implicit none (type, external)
      integer,      intent(in)               :: w_angular_rank
      integer,      intent(in) :: bra_shell_codes(7), ket_shell_codes(7)
      real(real64), intent(in)               :: q_projection_1, q_projection_2
      real(real64), intent(in) :: bra_shell_values(3), ket_shell_values(3)
      real(real64), intent(out)              :: w_matrix_element
      integer :: two_w_rank,occupation_change_2,two_abs_intermediate_projection,twice_phase_exponent, &
          intermediate_term_index,bra_term_lookup_index,ket_term_lookup_index,first_intermediate_term, &
          first_bra_range_term,last_intermediate_term,last_bra_range_term,two_intermediate_seniority
      integer :: intermediate_shell_codes(7)
      real(real64)               :: sum_over_terms, bra_a_matrix_element, six_j_value, ket_a_matrix_element
      real(real64) :: intermediate_shell_values(3)
      w_matrix_element=zero
      if(izas1(ket_shell_codes(7),ket_shell_values(3),bra_shell_codes(7),bra_shell_values(3)) == 0)return
      sum_over_terms=zero
      two_w_rank=w_angular_rank*2
      occupation_change_2=q_projection_2*two+q_projection_2*tenth
      if(ket_shell_codes(3) > 37) return
      if(ittk(ket_shell_codes(6),bra_shell_codes(6),two_w_rank) == 0)return
      bra_term_lookup_index=bra_shell_codes(1)
      ket_term_lookup_index=ket_shell_codes(1)
      if(ket_shell_codes(3) == 9) then
        if(ket_shell_codes(4) > 2) call mes(1)
        if(bra_shell_codes(4) > 2) call mes(1)
        bra_term_lookup_index=bra_term_lookup_index-300
        ket_term_lookup_index=ket_term_lookup_index-300
        first_bra_range_term=first_opposite_parity_terms_two_j_9(bra_term_lookup_index)
        first_intermediate_term=first_opposite_parity_terms_two_j_9(ket_term_lookup_index)
        if(first_bra_range_term /= first_intermediate_term)return
        last_bra_range_term=last_opposite_parity_terms_two_j_9(bra_term_lookup_index)
        last_intermediate_term=last_opposite_parity_terms_two_j_9(ket_term_lookup_index)
       else
        if(ket_shell_codes(4) > 2) call mes(1)
        if(bra_shell_codes(4) > 2) call mes(1)
        first_bra_range_term=first_opposite_parity_terms_high_j(bra_term_lookup_index)
        first_intermediate_term=first_opposite_parity_terms_high_j(ket_term_lookup_index)
        if(first_bra_range_term /= first_intermediate_term)return
        last_bra_range_term=last_opposite_parity_terms_high_j(bra_term_lookup_index)
        last_intermediate_term=last_opposite_parity_terms_high_j(ket_term_lookup_index)
      end if
      if(last_bra_range_term /= last_intermediate_term)return
      intermediate_shell_codes(2)=ket_shell_codes(2)
      intermediate_shell_codes(3)=ket_shell_codes(3)
      intermediate_shell_codes(4)=ket_shell_codes(4)+occupation_change_2
      intermediate_shell_values(3)=ket_shell_values(3)+half*dble(occupation_change_2)
      two_abs_intermediate_projection=two*dabs(intermediate_shell_values(3))+tenth

      ! Sum over admissible intermediate terms at the required occupation and projection.
      do intermediate_term_index=first_intermediate_term,last_intermediate_term
        call rumtjj(intermediate_term_index,intermediate_shell_codes(3),intermediate_shell_codes(7), &
            two_intermediate_seniority,intermediate_shell_codes(6))
        if(two_abs_intermediate_projection <= intermediate_shell_codes(7)) then
          if(ixjtik(bra_shell_codes(3),bra_shell_codes(3),two_w_rank,ket_shell_codes(6),bra_shell_codes(6), &
              intermediate_shell_codes(6)) /= 0) then
            intermediate_shell_codes(1)=intermediate_term_index
            intermediate_shell_values(2)=dble(intermediate_shell_codes(6))/two
            intermediate_shell_values(1)=dble(intermediate_shell_codes(7))/two
            call a1jj(bra_shell_codes,bra_shell_values,intermediate_shell_codes,intermediate_shell_values, &
                q_projection_1,bra_a_matrix_element)
            if(dabs(bra_a_matrix_element) > angular_tolerance) then
              call a1jj(intermediate_shell_codes,intermediate_shell_values,ket_shell_codes,ket_shell_values, &
                  q_projection_2,ket_a_matrix_element)
              if(dabs(ket_a_matrix_element) > angular_tolerance) then
                bra_a_matrix_element=bra_a_matrix_element*ket_a_matrix_element
                call sixj(bra_shell_codes(3),bra_shell_codes(3),two_w_rank,ket_shell_codes(6),bra_shell_codes(6), &
                    intermediate_shell_codes(6),0,six_j_value)
                sum_over_terms=sum_over_terms+bra_a_matrix_element*six_j_value
              end if
            end if
          end if
        end if
      end do
      w_matrix_element=sum_over_terms*dsqrt(dble((two_w_rank+1)))
      twice_phase_exponent=two_w_rank+bra_shell_codes(6)+ket_shell_codes(6)
      if(((twice_phase_exponent/4)*4) /= twice_phase_exponent)w_matrix_element=-w_matrix_element
      return
      end subroutine w1jjg
end module jj_w_matrix_element_high_j
