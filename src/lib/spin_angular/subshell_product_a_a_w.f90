! Multiply two single-tensor elements and a pair-tensor element on three subshells.
module subshell_product_a_a_w
    implicit none (type, external)
    private
    public :: gg1233
contains
      subroutine gg1233(bra_shell_codes_1,bra_shell_codes_2,bra_shell_codes_3,bra_shell_values_1, &
          bra_shell_values_2,bra_shell_values_3,ket_shell_codes_1,ket_shell_codes_2,ket_shell_codes_3, &
          ket_shell_values_1,  &
      ket_shell_values_2,ket_shell_values_3,w_angular_rank,q_projection_1,q_projection_2,q_projection_3, &
          q_projection_4,subshell_matrix_element)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants,          only: zero, tenth, angular_tolerance => eps
      use subshell_tensor_diagnostics, only: mes
      use quasispin_cg_rank_half, only: c0t5s
      use jj_reduced_cfp, only: rmeajj
      use jj_w_matrix_element, only: wj1
      implicit none (type, external)
      integer,      intent(in)               :: w_angular_rank
      integer,      intent(in) :: bra_shell_codes_1(7),bra_shell_codes_2(7),bra_shell_codes_3(7), &
          ket_shell_codes_1(7),ket_shell_codes_2(7),ket_shell_codes_3(7)
      real(real64), intent(in)               :: q_projection_1,q_projection_2,q_projection_3,q_projection_4
      real(real64), intent(in) :: bra_shell_values_1(3),bra_shell_values_2(3),bra_shell_values_3(3), &
          ket_shell_values_1(3),ket_shell_values_2(3),ket_shell_values_3(3)
      real(real64), intent(out)              :: subshell_matrix_element
      integer      :: occupation_change_1,occupation_change_2,occupation_change_3,occupation_change_4, &
          last_shell_occupation_change
      real(real64) :: cg_product, second_reduced_cfp, next_cg_factor, first_reduced_cfp, last_w_matrix_element
      subshell_matrix_element=zero
      if(bra_shell_codes_1(3) > 9) then
        if(bra_shell_codes_1(4) > 2) call mes(32)
        if(ket_shell_codes_1(4) > 2) call mes(32)
      end if
      if(bra_shell_codes_2(3) > 9) then
        if(bra_shell_codes_2(4) > 2) call mes(32)
        if(ket_shell_codes_2(4) > 2) call mes(32)
      end if
      if(bra_shell_codes_3(3) > 9) then
        if(bra_shell_codes_3(4) > 2) call mes(32)
        if(ket_shell_codes_3(4) > 2) call mes(32)
      end if
      occupation_change_1=q_projection_1+q_projection_1+tenth*q_projection_1
      if(bra_shell_codes_1(4) /= (ket_shell_codes_1(4)+occupation_change_1))return
      occupation_change_2=q_projection_2+q_projection_2+tenth*q_projection_2
      if(bra_shell_codes_2(4) /= (ket_shell_codes_2(4)+occupation_change_2))return
      occupation_change_3=q_projection_3+q_projection_3+tenth*q_projection_3
      occupation_change_4=q_projection_4+q_projection_4+tenth*q_projection_4
      last_shell_occupation_change=occupation_change_3+occupation_change_4
      if(bra_shell_codes_3(4) /= (ket_shell_codes_3(4)+last_shell_occupation_change))return
      call c0t5s(ket_shell_values_1(1),ket_shell_values_1(3),q_projection_1,bra_shell_values_1(1), &
          bra_shell_values_1(3),cg_product)
      if(dabs(cg_product) < angular_tolerance)return
      call c0t5s(ket_shell_values_2(1),ket_shell_values_2(3),q_projection_2,bra_shell_values_2(1), &
          bra_shell_values_2(3),next_cg_factor)
      if(dabs(next_cg_factor) < angular_tolerance)return
      cg_product=cg_product*next_cg_factor
      call rmeajj(bra_shell_codes_1(3),bra_shell_codes_1(1),bra_shell_codes_1(7),bra_shell_codes_1(6), &
          ket_shell_codes_1(1),ket_shell_codes_1(7),ket_shell_codes_1(6),first_reduced_cfp)
      if(dabs(first_reduced_cfp) < angular_tolerance)return
      call rmeajj(bra_shell_codes_2(3),bra_shell_codes_2(1),bra_shell_codes_2(7),bra_shell_codes_2(6), &
          ket_shell_codes_2(1),ket_shell_codes_2(7),ket_shell_codes_2(6),second_reduced_cfp)
      if(dabs(second_reduced_cfp) < angular_tolerance)return
      call wj1(bra_shell_codes_3,bra_shell_values_3,ket_shell_codes_3,ket_shell_values_3,w_angular_rank, &
          q_projection_3,q_projection_4,last_w_matrix_element)
      if(dabs(last_w_matrix_element) < angular_tolerance)return
      subshell_matrix_element=cg_product*last_w_matrix_element*first_reduced_cfp*second_reduced_cfp/ &
          dsqrt(dble((bra_shell_codes_1(7)+1)*(bra_shell_codes_2(7)+1)))
      return
      end subroutine gg1233
end module subshell_product_a_a_w
