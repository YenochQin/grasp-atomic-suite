! Look up higher-rank j = 7/2 pair-tensor elements and apply the term-order phase.
module jj_reduced_w_j_seven_halves_higher_ranks
    implicit none (type, external)

    ! Module interface
    private
    public :: rmew7bjj
contains
    subroutine rmew7bjj(bra_term_index,ket_term_index,quasispin_rank,angular_rank,reduced_w_element)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: error_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero
        use jj_reduced_coefficients, only: last_opposite_parity_terms_low_j => imgnjj, last_same_parity_terms_low_j => imgtjj, &
            first_opposite_parity_terms_low_j => impnjj, first_same_parity_terms_low_j => imptjj
        use jj_subshell_term_decoding, only: rumtjj
        implicit none (type, external)

        ! Arguments
        integer,      intent(in)  :: bra_term_index, ket_term_index, quasispin_rank, angular_rank
        real(real64), intent(out) :: reduced_w_element

        ! Local declarations
        integer :: two_orbital_j, two_bra_q, two_bra_seniority, two_bra_term_j, two_ket_q, two_ket_seniority, &
            two_ket_term_j, lower_table_index, upper_table_index, twice_phase_exponent, packed_pair_index, term_group
        integer  :: first_group_row_offsets(6)
        integer  :: second_group_row_offsets(8)
        integer :: w_numerators_1_2_group_1(21),w_denominators_1_2_group_1(21),w_numerators_0_3_group_1(21), &
            w_denominators_0_3_group_1(21),w_numerators_1_4_group_1(21), &
            w_denominators_1_4_group_1(21),w_numerators_0_5_group_1(21),w_denominators_0_5_group_1(21), &
            w_numerators_1_6_group_1(21),w_denominators_1_6_group_1(21),w_numerators_0_7_group_1(21), &
            w_denominators_0_7_group_1(21)
        integer :: w_numerators_1_2_group_2(36),w_denominators_1_2_group_2(36),w_numerators_0_3_group_2(36), &
            w_denominators_0_3_group_2(36),w_numerators_1_4_group_2(36), &
            w_denominators_1_4_group_2(36),w_numerators_0_5_group_2(36),w_denominators_0_5_group_2(36), &
            w_numerators_1_6_group_2(36),w_denominators_1_6_group_2(36),w_numerators_0_7_group_2(36), &
            w_denominators_0_7_group_2(36)

        ! Initial values
        data first_group_row_offsets/0,5,9,12,14,15/
        data second_group_row_offsets/0,7,13,18,22,25,27,28/
        data w_numerators_1_2_group_1/-252,1056,144,3*0,507,-88,-325,2*0,200,-520,80,  &
            0,-6125,-560,0,390,-3840,2040/
        data w_denominators_1_2_group_1/10,70,7,3*1,10,1,14,2*1,3,21,2*1,462,11,1,539,   &
            2*49/
        data w_numerators_1_2_group_2/0,-50,6*0,-1280,990,2640,1950,4*0,-480,4*0,-360, &
            1872,42,390,3*0,48,2*0,234,0,1040,340,0/
        data w_denominators_1_2_group_2/8*1,4*49,4*1,49,4*1,539,49,1,11,6*1,7,1,11,7,1/
        data w_numerators_0_3_group_1/-1188,-196,0,-234,2*0,189,0,-1911,1470,0,-56,3*0,&
            394805,5250,53760,-78,17408,12920/
        data w_denominators_0_3_group_1/70,10,1,7,2*1,110,1,242,121,5*1,22022,121,1573,  &
            2*847,1001/
        data w_numerators_0_3_group_2/8*0,110,0,-240,4*0,-1920,0,52,-224,2*0,32490,2*0,&
            -7644,0,12,224,2*0,-52,0,-2040,-364,0,1292/
        data w_denominators_0_3_group_2/8*1,7,1,7,4*1,77,1,7,11,2*1,847,2*1,121,1,70,10, &
            2*1,70,1,77,121,1,77/
        data w_numerators_1_4_group_1/0,54,-528,546,-378,0,-4335,-96,11271,3822,0,120, &
            12000,-624,-960,-30345,-210,-228480,580476,146880,-627912/
        data w_denominators_1_4_group_1/1,2*7,2*11,1,110,11,1694,121,2*1,77,2*11,3146,   &
            121,1573,2*5929,7007/
        data w_numerators_1_4_group_2/3*0,-90,4*0,2640,480,-360,20592,-42,390,2*0,     &
            468180,2*0,21840,0,-359424,-6750,0,917280,-10710,2*0,36,2*0,&
            -858,0,-5304,-69768,0/
        data w_denominators_1_4_group_2/8*1,2*49,2*539,1,11,2*1,5929,2*1,847,1,5929,539, &
            1,5929,121,2*1,11,2*1,7,1,121,847,1/
        data w_numerators_0_5_group_1/3*0,144,14,0,975,0,-50421,336,-7000,-88,3*0,     &
            -6845,-3360,14280,-1836,-103360,-28424/
        data w_denominators_0_5_group_1/3*1,7,2*1,22,1,2002,11,143,4*1,26026,143,1859,77,&
            1001,143143/
        data w_numerators_0_5_group_2/10*0,390,2*0,-70,3*0,32,14,2*0,-576,2*0,-210,0,  &
            63888,-154,0,-8568,-176,0,1938,1088,0,-28424/
        data w_denominators_0_5_group_2/10*1,7,6*1,7,3*1,77,2*1,11,1,1183,13,1,169,7,1,  &
            91,11,1,1183/
        data w_numerators_1_6_group_1/3*0,576,390,-144,0,520,-49,-12480,-408,520,      &
            -1960,-1664,3264,552250,-43520,-38760,-2652,15504,38760/
        data w_denominators_1_6_group_1/3*1,11,77,7,1,11,2*121,11,3,33,2*11,4719,847,    &
            11011,2*121,143/
        data w_numerators_1_6_group_2/6*0,-130,3*0,390,-48,-234,1040,340,0,3120,2*0,   &
            -1170,0,18720,-36,858,-5304,-69768,2*0,1020,4*0,-33592,     &
            31654,0/
        data w_denominators_1_6_group_2/10*1,11,1,7,11,7,1,121,2*1,121,1,121,11,7,121,   &
            847,2*1,11,4*1,2*121,1/
        data w_numerators_0_7_group_1/4*0,162,272,2*0,11025,4624,-1632,-120,3*0,       &
            1224510,306000,12558240,-6460,77520,-297160/
        data w_denominators_0_7_group_1/4*1,2*7,2*1,1573,121,143,4*1,20449,11011,143143, &
            121,1573,1859/
        data w_numerators_0_7_group_2/13*0,60,4*0,1600,0,2040,4410,2*0,18360,0,-18816, &
            -11016,0,11628,34,0,7752,9690,0,222870/
        data w_denominators_0_7_group_2/18*1,77,1,77,121,2*1,121,1,845,455,1,1183,5,1,   &
            143,121,1,1859/
        reduced_w_element = zero
        two_orbital_j = 7
        call rumtjj(bra_term_index,two_orbital_j,two_bra_q,two_bra_seniority,two_bra_term_j)
        call rumtjj(ket_term_index,two_orbital_j,two_ket_q,two_ket_seniority,two_ket_term_j)
        if(bra_term_index > ket_term_index) then
            lower_table_index=ket_term_index
            upper_table_index=bra_term_index
            twice_phase_exponent=two_ket_term_j-two_bra_term_j+two_ket_q-two_bra_q
        else
            lower_table_index=bra_term_index
            upper_table_index=ket_term_index
            twice_phase_exponent=4
        end if
        if(bra_term_index > 17) then
            lower_table_index=lower_table_index-17
            upper_table_index=upper_table_index-17
            packed_pair_index=second_group_row_offsets(lower_table_index)+upper_table_index
            term_group=2
        else
            lower_table_index=lower_table_index-11
            upper_table_index=upper_table_index-11
            term_group=1
            packed_pair_index=first_group_row_offsets(lower_table_index)+upper_table_index
        end if
        if(quasispin_rank == 1 .and. angular_rank == 2) then
            if(term_group == 1) then
                if(w_numerators_1_2_group_1(packed_pair_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_1_2_group_1(packed_pair_index))/ &
                        dble(w_denominators_1_2_group_1(packed_pair_index)))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_1_2_group_1(packed_pair_index))/ &
                        dble(w_denominators_1_2_group_1(packed_pair_index)))
                end if
            else
                if(w_numerators_1_2_group_2(packed_pair_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_1_2_group_2(packed_pair_index))/ &
                        dble(w_denominators_1_2_group_2(packed_pair_index)))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_1_2_group_2(packed_pair_index))/ &
                        dble(w_denominators_1_2_group_2(packed_pair_index)))
                end if
            end if
        else if(quasispin_rank == 0 .and. angular_rank == 3) then
            if(term_group == 1) then
                if(w_numerators_0_3_group_1(packed_pair_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_0_3_group_1(packed_pair_index))/ &
                        dble(w_denominators_0_3_group_1(packed_pair_index)))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_0_3_group_1(packed_pair_index))/ &
                        dble(w_denominators_0_3_group_1(packed_pair_index)))
                end if
            else
                if(w_numerators_0_3_group_2(packed_pair_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_0_3_group_2(packed_pair_index))/ &
                        dble(w_denominators_0_3_group_2(packed_pair_index)))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_0_3_group_2(packed_pair_index))/ &
                        dble(w_denominators_0_3_group_2(packed_pair_index)))
                end if
            end if
        else if(quasispin_rank == 1 .and. angular_rank == 4) then
            if(term_group == 1) then
                if(w_numerators_1_4_group_1(packed_pair_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_1_4_group_1(packed_pair_index))/ &
                        dble(w_denominators_1_4_group_1(packed_pair_index)))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_1_4_group_1(packed_pair_index))/ &
                        dble(w_denominators_1_4_group_1(packed_pair_index)))
                end if
            else
                if(w_numerators_1_4_group_2(packed_pair_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_1_4_group_2(packed_pair_index))/ &
                        dble(w_denominators_1_4_group_2(packed_pair_index)))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_1_4_group_2(packed_pair_index))/ &
                        dble(w_denominators_1_4_group_2(packed_pair_index)))
                end if
            end if
        else if(quasispin_rank == 0 .and. angular_rank == 5) then
            if(term_group == 1) then
                if(w_numerators_0_5_group_1(packed_pair_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_0_5_group_1(packed_pair_index))/ &
                        dble(w_denominators_0_5_group_1(packed_pair_index)))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_0_5_group_1(packed_pair_index))/ &
                        dble(w_denominators_0_5_group_1(packed_pair_index)))
                end if
            else
                if(w_numerators_0_5_group_2(packed_pair_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_0_5_group_2(packed_pair_index))/ &
                        dble(w_denominators_0_5_group_2(packed_pair_index)))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_0_5_group_2(packed_pair_index))/ &
                        dble(w_denominators_0_5_group_2(packed_pair_index)))
                end if
            end if
        else if(quasispin_rank == 1 .and. angular_rank == 6) then
            if(term_group == 1) then
                if(w_numerators_1_6_group_1(packed_pair_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_1_6_group_1(packed_pair_index))/ &
                        dble(w_denominators_1_6_group_1(packed_pair_index)))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_1_6_group_1(packed_pair_index))/ &
                        dble(w_denominators_1_6_group_1(packed_pair_index)))
                end if
            else
                if(w_numerators_1_6_group_2(packed_pair_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_1_6_group_2(packed_pair_index))/ &
                        dble(w_denominators_1_6_group_2(packed_pair_index)))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_1_6_group_2(packed_pair_index))/ &
                        dble(w_denominators_1_6_group_2(packed_pair_index)))
                end if
            end if
        else if(quasispin_rank == 0 .and. angular_rank == 7) then
            if(term_group == 1) then
                if(w_numerators_0_7_group_1(packed_pair_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_0_7_group_1(packed_pair_index))/ &
                        dble(w_denominators_0_7_group_1(packed_pair_index)))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_0_7_group_1(packed_pair_index))/ &
                        dble(w_denominators_0_7_group_1(packed_pair_index)))
                end if
            else
                if(w_numerators_0_7_group_2(packed_pair_index) >= 0) then
                    reduced_w_element=dsqrt(dble(w_numerators_0_7_group_2(packed_pair_index))/ &
                        dble(w_denominators_0_7_group_2(packed_pair_index)))
                else
                    reduced_w_element=-dsqrt(-dble(w_numerators_0_7_group_2(packed_pair_index))/ &
                        dble(w_denominators_0_7_group_2(packed_pair_index)))
                end if
            end if
        else
            write(error_unit,"(A,4I5)") " J1 J2 = ",bra_term_index,ket_term_index
            write(error_unit,"(A)") " ERROR IN SUB. RMEW7JJ "
            stop
        end if
        if(mod(twice_phase_exponent,4) /= 0)reduced_w_element=-reduced_w_element
        return
    end subroutine rmew7bjj
end module jj_reduced_w_j_seven_halves_higher_ranks
