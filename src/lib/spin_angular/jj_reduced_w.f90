! Select a table of completely reduced pair-tensor matrix elements by orbital j.
module jj_reduced_w
    implicit none (type, external)

    ! Module interface
    private
    public :: evaluate_jj_reduced_w
contains
    subroutine evaluate_jj_reduced_w(two_orbital_j,bra_term_index,ket_term_index,quasispin_rank,angular_rank,reduced_w_element)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: error_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use jj_reduced_w_j_one_half, only: evaluate_reduced_w_j_one_half
        use jj_reduced_w_j_three_halves, only: evaluate_reduced_w_j_three_halves
        use jj_reduced_w_j_five_halves, only: evaluate_reduced_w_j_five_halves
        use jj_reduced_w_j_seven_halves, only: evaluate_reduced_w_j_seven_halves
        implicit none (type, external)

        ! Arguments
        integer,      intent(in)  :: two_orbital_j, bra_term_index, ket_term_index, quasispin_rank, angular_rank
        real(real64), intent(out) :: reduced_w_element
        if(two_orbital_j == 1) then
            call evaluate_reduced_w_j_one_half(bra_term_index,ket_term_index,quasispin_rank,angular_rank,reduced_w_element)
        else if(two_orbital_j == 3) then
            call evaluate_reduced_w_j_three_halves(bra_term_index,ket_term_index,quasispin_rank,angular_rank,reduced_w_element)
        else if(two_orbital_j == 5) then
            call evaluate_reduced_w_j_five_halves(bra_term_index,ket_term_index,quasispin_rank,angular_rank,reduced_w_element)
        else if(two_orbital_j == 7) then
            call evaluate_reduced_w_j_seven_halves(bra_term_index,ket_term_index,quasispin_rank,angular_rank,reduced_w_element)
        else
            write(error_unit,"(A,I5)") " KLAIDA SUB. RWJJ J=",two_orbital_j
            stop
        end if
        return
    end subroutine evaluate_jj_reduced_w
end module jj_reduced_w
