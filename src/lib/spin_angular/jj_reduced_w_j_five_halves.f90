! Look up pair-tensor elements for j = 5/2 by quasispin and angular rank.
module jj_reduced_w_j_five_halves
    implicit none (type, external)

    ! Module interface
    private
    public :: evaluate_reduced_w_j_five_halves
contains
    subroutine evaluate_reduced_w_j_five_halves(bra_term_index,ket_term_index,quasispin_rank,angular_rank,reduced_w_element)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: error_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero
        use jj_reduced_coefficients, only: first_same_parity_terms_low_j
        implicit none (type, external)

        ! Arguments
        integer,      intent(in)  :: bra_term_index, ket_term_index, quasispin_rank, angular_rank
        real(real64), intent(out) :: reduced_w_element

        ! Local declarations
        integer                 :: bra_table_index, ket_table_index
        integer   :: w_numerators_0_0(6), w_numerators_0_1(6)
        integer :: w_numerators_1_2_group_1(3, 3),w_numerators_1_2_group_2(3, 3),w_numerators_0_3_group_1(3, 3), &
            w_numerators_0_3_group_2(3, 3),w_numerators_1_4_group_1(3, 3),w_numerators_1_4_group_2(3, &
            3),w_numerators_0_5_group_1(3, 3),w_numerators_0_5_group_2(3, 3)

        ! Initial values
        data w_numerators_0_0/54,12,30,12,30,54/
        data w_numerators_0_1/126,12,198,0,48,288/
        data w_numerators_1_2_group_1/420,360,270,360,2*0,-270,2*0/
        data w_numerators_1_2_group_2/0,1960,0,-1960,-1000,2430,0,2430,1980/
        data w_numerators_0_3_group_1/-882,3*0,384,-400,0,400,286/
        data w_numerators_0_3_group_2/4*0,162,-300,0,-300,22/
        data w_numerators_1_4_group_1/3780,-720,-4950,-720,2*0,4950,2*0/
        data w_numerators_1_4_group_2/2*0,3528,0,2430,1980,-3528,1980,-7722/
        data w_numerators_0_5_group_1/-1386,4*0,-440,0,440,-1430/
        data w_numerators_0_5_group_2/5*0,330,0,330,572/
        reduced_w_element=zero
        if(first_same_parity_terms_low_j(bra_term_index) /= first_same_parity_terms_low_j(ket_term_index)) return
        if(bra_term_index < 6 .or. bra_term_index > 11) return
        if(quasispin_rank == 0 .and. angular_rank == 0) then
            if(bra_term_index /= ket_term_index) return
            reduced_w_element=-dsqrt(dble(w_numerators_0_0(bra_term_index-5)))
        else if(quasispin_rank == 1 .and. angular_rank == 0) then
            if(bra_term_index /= ket_term_index) return
            if(bra_term_index == 6) then
                reduced_w_element=-dsqrt(dble(48))
            else if(bra_term_index == 9) then
                reduced_w_element=-dsqrt(dble(20))
            else if(bra_term_index == 10) then
                reduced_w_element=-dsqrt(dble(10))
            else if(bra_term_index == 11) then
                reduced_w_element=-dsqrt(dble(18))
            end if
        else if(quasispin_rank == 0 .and. angular_rank == 1) then
            if(bra_term_index /= ket_term_index) return
            reduced_w_element=-dsqrt(dble(w_numerators_0_1(bra_term_index-5))/dble(7))
        else if(quasispin_rank == 1 .and. angular_rank == 2) then
            if(bra_term_index < 9) then
                bra_table_index=bra_term_index-5
                ket_table_index=ket_term_index-5
                if(w_numerators_1_2_group_1(bra_table_index,ket_table_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_1_2_group_1(bra_table_index,ket_table_index))/dble(7))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_1_2_group_1(bra_table_index,ket_table_index))/dble(7))
                end if
            else
                bra_table_index=bra_term_index-8
                ket_table_index=ket_term_index-8
                if(w_numerators_1_2_group_2(bra_table_index,ket_table_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_1_2_group_2(bra_table_index,ket_table_index))/dble(49))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_1_2_group_2(bra_table_index,ket_table_index))/dble(49))
                end if
            end if
        else if(quasispin_rank == 0 .and. angular_rank == 3) then
            if(bra_term_index < 9) then
                bra_table_index=bra_term_index-5
                ket_table_index=ket_term_index-5
                if(w_numerators_0_3_group_1(bra_table_index,ket_table_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_0_3_group_1(bra_table_index,ket_table_index))/dble(21))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_0_3_group_1(bra_table_index,ket_table_index))/dble(21))
                end if
            else
                bra_table_index=bra_term_index-8
                ket_table_index=ket_term_index-8
                if(w_numerators_0_3_group_2(bra_table_index,ket_table_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_0_3_group_2(bra_table_index,ket_table_index))/dble(7))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_0_3_group_2(bra_table_index,ket_table_index))/dble(7))
                end if
            end if
        else if(quasispin_rank == 1 .and. angular_rank == 4) then
            if(bra_term_index < 9) then
                bra_table_index=bra_term_index-5
                ket_table_index=ket_term_index-5
                if(w_numerators_1_4_group_1(bra_table_index,ket_table_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_1_4_group_1(bra_table_index,ket_table_index))/dble(35))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_1_4_group_1(bra_table_index,ket_table_index))/dble(35))
                end if
            else
                bra_table_index=bra_term_index-8
                ket_table_index=ket_term_index-8
                if(w_numerators_1_4_group_2(bra_table_index,ket_table_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_1_4_group_2(bra_table_index,ket_table_index))/dble(49))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_1_4_group_2(bra_table_index,ket_table_index))/dble(49))
                end if
            end if
        else if(quasispin_rank == 0 .and. angular_rank == 5) then
            if(bra_term_index < 9) then
                bra_table_index=bra_term_index-5
                ket_table_index=ket_term_index-5
                if(w_numerators_0_5_group_1(bra_table_index,ket_table_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_0_5_group_1(bra_table_index,ket_table_index))/dble(21))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_0_5_group_1(bra_table_index,ket_table_index))/dble(21))
                end if
            else
                bra_table_index=bra_term_index-8
                ket_table_index=ket_term_index-8
                if(w_numerators_0_5_group_2(bra_table_index,ket_table_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_0_5_group_2(bra_table_index,ket_table_index))/dble(7))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_0_5_group_2(bra_table_index,ket_table_index))/dble(7))
                end if
            end if
        else
            write(error_unit,"(A,4I5)") " J1 J2 = ",bra_term_index,ket_term_index
            write(error_unit,"(A)") " ERROR IN SUB. RMEW5JJ "
            stop
        end if
        return
    end subroutine evaluate_reduced_w_j_five_halves
end module jj_reduced_w_j_five_halves
