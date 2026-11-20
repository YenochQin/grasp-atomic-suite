! Look up the completely reduced j = 9/2 creation or annihilation matrix element.
module jj_reduced_cfp_j_nine_halves
    implicit none (type, external)

    ! Module interface
    private
    public :: rmeajj9
contains
    subroutine rmeajj9(bra_term_index,two_bra_q,two_bra_term_j,ket_term_index,two_ket_q,two_ket_term_j,reduced_cfp)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: error_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero
        implicit none (type, external)

        ! Arguments
        integer,      intent(in)  :: bra_term_index, two_bra_q, two_bra_term_j, ket_term_index, two_ket_q, &
            two_ket_term_j
        real(real64), intent(out) :: reduced_cfp

        ! Local declarations
        integer :: cfp_row, cfp_column
        integer :: cfp_numerators_col_3(20),cfp_numerators_col_4(20),cfp_numerators_col_5(20), &
            cfp_numerators_col_6(20),cfp_numerators_col_7(20),cfp_denominators_col_7(20),cfp_numerators_col_8(20),&
            cfp_denominators_col_8(20),cfp_numerators_col_9(20),cfp_numerators_col_10(20),cfp_numerators_col_11(20), &
            cfp_denominators_col_11(20),cfp_numerators_col_12(20),cfp_denominators_col_12(20), &
            cfp_numerators_col_13(20),cfp_denominators_col_13(20),cfp_numerators_col_14(20),   &
            cfp_numerators_col_15(20),cfp_denominators_col_15(20),cfp_numerators_col_16(20),cfp_denominators_col_16(20), &
            cfp_numerators_col_17(20),cfp_denominators_col_17(20),cfp_numerators_col_18(20),cfp_denominators_col_18(20)

        ! Initial values
        data cfp_numerators_col_3/2*0,660,0,1664,0,1650,130,0,816,0,1680,8*0/
        data cfp_numerators_col_4/2*0,-55770,-20592,6292,-13260,0,12740,-19110,     &
            79968,-47880,-25410,-19278,7*0/
        data cfp_numerators_col_5/0,-9009,-4056,-260,585,16575,0,25200,4200,-7140,  &
            -475,504,-14280,13566,-4250,5*0/
        data cfp_numerators_col_6/0,-4992,2028,0,-1920,0,-12870,7350,0,8160,        &
            0,1512,0,-3648,0,-9000,4*0/
        data cfp_numerators_col_7/-2184,-63,-59904,-302460,5265,848691,0,-145152,   &
            217728,1049580,287337,5184,261120,691866,-750,0,-76608,3*0/
        data cfp_denominators_col_7/253,23,19481,19481,161,253253,1,36179,            &
            36179,2*36179,299,3*36179,1,3289,3*1/
        data cfp_numerators_col_8/1224,2652,188598,-31824,204,-12996,0,25500,-38250,&
            -768,81396,3213,-60543,-3066144,727776,207,-41553,3*0/
        data cfp_denominators_col_8/1265,115,13915,2783,23,13915,1,3*2783,13915,115,  &
            2783,236555,47311,17,21505,3*1/
        data cfp_numerators_col_9/3380,-2340,5460,-1400,3150,-3570,0,9000,1500,2550,&
            -5320,-10710,1050,1140,6300,8550,-330,5750,2*0/
        data cfp_numerators_col_10/0,2160,6240,0,-32,0,-21450,-9610,0,2688,         &
            0,7140,0,20160,0,6156,2*0,-10164,0/
        data cfp_numerators_col_11/0,132,-52728,196,-50,-24990,0,-21160,31740,26250,&
            12920,-357,9583,-344988,5700,171,-15,-4830,-84,0/
        data cfp_denominators_col_11/1,5,6655,1331,11,1331,1,4*1331,55,1331,6655,     &
            1331,11,2*121,11,1/
        data cfp_numerators_col_12/2*0,-209950,77520,12920,-25688,0,-4522,2261,     &
            -48640,-285144,931,273885,-112908,-2138580,137781,6654375,  &
            -59616,284089,0/
        data cfp_denominators_col_12/2*1,2*9317,231,9317,1,3993,1331,22627,9317,44,   &
            90508,429913,158389,3740,156332,39083,17765,1/
        data cfp_numerators_col_13/2*0,1530,13056,29376,720,0,-1890,-315,101124,    &
            -4560,-13965,-35131,13500,-685900,-1197,-28875,-5060,759,0/
        data cfp_denominators_col_13/2*1,77,3*1001,1,2*143,2431,1001,572,9724,2431,   &
            17017,2*884,221,17,1/
        data cfp_numerators_col_14/4*0,22848,0,-121550,4590,0,-32832,0,45220,0,     &
            -31680,0,82764,2*0,144716,0/
        data cfp_numerators_col_15/5*0,2128,0,-1938,2907,-4860,-17136,-1309,-8505,  &
            15876,420,1287,-6075,-132,-253,-650/
        data cfp_denominators_col_15/5*1,143,1,3*143,2717,52,572,247,13,20,988,19,95, &
            19/
        data cfp_numerators_col_16/7*0,570,95,4104,504,-1463,-39501,-60516,-1596,   &
            3933,621,-840,-805,390/
        data cfp_denominators_col_16/7*1,2*13,221,65,260,884,1105,221,68,3740,2*17,   &
            11/
        data cfp_numerators_col_17/9*0,-5796,17664,-5313,1771,-16632,-88,693,       &
            -94269,192500,30030,24570/
        data cfp_denominators_col_17/9*1,221,1235,65,221,1615,2*17,1615,7429,323,437/
        data cfp_numerators_col_18/13*0,-15000,280,-1170,-48750,-15600,3510,-33930/
        data cfp_denominators_col_18/13*1,323,2*17,3553,437,19,253/
        reduced_cfp=zero
        if(bra_term_index < ket_term_index) then
            cfp_row=bra_term_index-25
            cfp_column=ket_term_index-45
        else
            cfp_row=ket_term_index-25
            cfp_column=bra_term_index-45
        end if
        if(cfp_column == 1) then
            if(cfp_row == 7) reduced_cfp=-dsqrt(dble(60))
        else if(cfp_column == 2) then
            if(cfp_row == 8) reduced_cfp=-dsqrt(dble(12))
            if(cfp_row == 9) reduced_cfp=-dsqrt(dble(8))
        else if(cfp_column == 3)then
            reduced_cfp=-dsqrt(dble(cfp_numerators_col_3(cfp_row))/dble(33))
        else if(cfp_column == 4)then
            if(cfp_numerators_col_4(cfp_row) >= 0) then
                reduced_cfp=dsqrt(dble(cfp_numerators_col_4(cfp_row))/dble(3003))
            else
                reduced_cfp=-dsqrt(-dble(cfp_numerators_col_4(cfp_row))/dble(3003))
            end if
        else if(cfp_column == 5)then
            if(cfp_numerators_col_5(cfp_row) >= 0) then
                reduced_cfp=dsqrt(dble(cfp_numerators_col_5(cfp_row))/dble(715))
            else
                reduced_cfp=-dsqrt(-dble(cfp_numerators_col_5(cfp_row))/dble(715))
            end if
        else if(cfp_column == 6)then
            if(cfp_numerators_col_6(cfp_row) >= 0) then
                reduced_cfp=dsqrt(dble(cfp_numerators_col_6(cfp_row))/dble(143))
            else
                reduced_cfp=-dsqrt(-dble(cfp_numerators_col_6(cfp_row))/dble(143))
            end if
        else if(cfp_column == 7)then
            if(cfp_numerators_col_7(cfp_row) >= 0) then
                reduced_cfp=dsqrt(dble(cfp_numerators_col_7(cfp_row))/dble(cfp_denominators_col_7(cfp_row)))
            else
                reduced_cfp=-dsqrt(-dble(cfp_numerators_col_7(cfp_row))/dble(cfp_denominators_col_7(cfp_row)))
            end if
        else if(cfp_column == 8)then
            if(cfp_numerators_col_8(cfp_row) >= 0) then
                reduced_cfp=dsqrt(dble(cfp_numerators_col_8(cfp_row))/dble(cfp_denominators_col_8(cfp_row)))
            else
                reduced_cfp=-dsqrt(-dble(cfp_numerators_col_8(cfp_row))/dble(cfp_denominators_col_8(cfp_row)))
            end if
        else if(cfp_column == 9)then
            if(cfp_numerators_col_9(cfp_row) >= 0) then
                reduced_cfp=dsqrt(dble(cfp_numerators_col_9(cfp_row))/dble(325))
            else
                reduced_cfp=-dsqrt(-dble(cfp_numerators_col_9(cfp_row))/dble(325))
            end if
        else if(cfp_column == 10)then
            if(cfp_numerators_col_10(cfp_row) >= 0) then
                reduced_cfp=dsqrt(dble(cfp_numerators_col_10(cfp_row))/dble(165))
            else
                reduced_cfp=-dsqrt(-dble(cfp_numerators_col_10(cfp_row))/dble(165))
            end if
        else if(cfp_column == 11)then
            if(cfp_numerators_col_11(cfp_row) >= 0) then
                reduced_cfp=dsqrt(dble(cfp_numerators_col_11(cfp_row))/dble(cfp_denominators_col_11(cfp_row)))
            else
                reduced_cfp=-dsqrt(-dble(cfp_numerators_col_11(cfp_row))/dble(cfp_denominators_col_11(cfp_row)))
            end if
        else if(cfp_column == 12)then
            if(cfp_numerators_col_12(cfp_row) >= 0) then
                reduced_cfp=dsqrt(dble(cfp_numerators_col_12(cfp_row))/dble(cfp_denominators_col_12(cfp_row)))
            else
                reduced_cfp=-dsqrt(-dble(cfp_numerators_col_12(cfp_row))/dble(cfp_denominators_col_12(cfp_row)))
            end if
        else if(cfp_column == 13)then
            if(cfp_numerators_col_13(cfp_row) >= 0) then
                reduced_cfp=dsqrt(dble(cfp_numerators_col_13(cfp_row))/dble(cfp_denominators_col_13(cfp_row)))
            else
                reduced_cfp=-dsqrt(-dble(cfp_numerators_col_13(cfp_row))/dble(cfp_denominators_col_13(cfp_row)))
            end if
        else if(cfp_column == 14)then
            if(cfp_numerators_col_14(cfp_row) >= 0) then
                reduced_cfp=dsqrt(dble(cfp_numerators_col_14(cfp_row))/dble(715))
            else
                reduced_cfp=-dsqrt(-dble(cfp_numerators_col_14(cfp_row))/dble(715))
            end if
        else if(cfp_column == 15)then
            if(cfp_numerators_col_15(cfp_row) >= 0) then
                reduced_cfp=dsqrt(dble(cfp_numerators_col_15(cfp_row))/dble(cfp_denominators_col_15(cfp_row)))
            else
                reduced_cfp=-dsqrt(-dble(cfp_numerators_col_15(cfp_row))/dble(cfp_denominators_col_15(cfp_row)))
            end if
        else if(cfp_column == 16)then
            if(cfp_numerators_col_16(cfp_row) >= 0) then
                reduced_cfp=dsqrt(dble(cfp_numerators_col_16(cfp_row))/dble(cfp_denominators_col_16(cfp_row)))
            else
                reduced_cfp=-dsqrt(-dble(cfp_numerators_col_16(cfp_row))/dble(cfp_denominators_col_16(cfp_row)))
            end if
        else if(cfp_column == 17)then
            if(cfp_numerators_col_17(cfp_row) >= 0) then
                reduced_cfp=dsqrt(dble(cfp_numerators_col_17(cfp_row))/dble(cfp_denominators_col_17(cfp_row)))
            else
                reduced_cfp=-dsqrt(-dble(cfp_numerators_col_17(cfp_row))/dble(cfp_denominators_col_17(cfp_row)))
            end if
        else if(cfp_column == 18)then
            if(cfp_numerators_col_18(cfp_row) >= 0) then
                reduced_cfp=dsqrt(dble(cfp_numerators_col_18(cfp_row))/dble(cfp_denominators_col_18(cfp_row)))
            else
                reduced_cfp=-dsqrt(-dble(cfp_numerators_col_18(cfp_row))/dble(cfp_denominators_col_18(cfp_row)))
            end if
        else
            write(error_unit,"(A,4I5)") " IT ITS JT JTS= ",bra_term_index,ket_term_index,cfp_row,cfp_column
            write(error_unit,"(A)") " ERROR IN SUB. RMEAJJ9 "
            stop
        end if
        return
    end subroutine rmeajj9
end module jj_reduced_cfp_j_nine_halves
