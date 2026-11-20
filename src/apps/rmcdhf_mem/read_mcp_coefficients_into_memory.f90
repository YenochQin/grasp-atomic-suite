!*******************************************************************
!                                                                  *
subroutine read_mcp_coefficients_into_memory(nfile,total_contribution_count,ncontr)
    ! Module dependencies
    !                                                                  *
    !   Written by  G. Gaigalas               Vilnius, September 2021  *
    !                                                                  *
    !*******************************************************************
    !
    use, intrinsic :: iso_fortran_env, only: int64
    use rmcdhf_mcp_memory, only: coeff, coeff_31, coeff_32, coeff_33, coeff_34, coeff_35, coeff_36, coeff_37, coeff_38, &
        & coeff_39, coeff_40, coeff_41, coeff_42, coeff_43, coeff_44, coeff_45, coeff_46, coeff_47, coeff_48, &
        & coeff_49, coeff_50, iclmn, iclmn_31, iclmn_32, iclmn_33, iclmn_34, iclmn_35, iclmn_36, iclmn_37, &
        & iclmn_38, iclmn_39, iclmn_40, iclmn_41, iclmn_42, iclmn_43, iclmn_44, iclmn_45, iclmn_46, iclmn_47, &
        & iclmn_48, iclmn_49, iclmn_50, indx, indx_31, indx_32, indx_33, indx_34, indx_35, indx_36, indx_37, &
        & indx_38, indx_39, indx_40, indx_41, indx_42, indx_43, indx_44, indx_45, indx_46, indx_47, indx_48, &
        & indx_49, indx_50
    implicit none (type, external)

    ! Arguments
    integer(int64), intent(in) :: total_contribution_count
    integer, intent(in) :: nfile
    integer :: ncontr

    ! Local declarations
    integer :: i
    !-----------------------------------------------

    do i = 1, ncontr
        !*** Read column index, sparse matrix index, and
        !*** coefficient for all contributions from this integral
        if(nfile == 31) then
            iclmn(i) = iclmn_31(total_contribution_count+i)
            indx(i) = indx_31(total_contribution_count+i)
            coeff(i) = coeff_31(total_contribution_count+i)
        else if(nfile == 32) then
            iclmn(i) = iclmn_32(total_contribution_count+i)
            indx(i) = indx_32(total_contribution_count+i)
            coeff(i) = coeff_32(total_contribution_count+i)
        else if(nfile == 33) then
            iclmn(i) = iclmn_33(total_contribution_count+i)
            indx(i) = indx_33(total_contribution_count+i)
            coeff(i) = coeff_33(total_contribution_count+i)
        else if(nfile == 34) then
            iclmn(i) = iclmn_34(total_contribution_count+i)
            indx(i) = indx_34(total_contribution_count+i)
            coeff(i) = coeff_34(total_contribution_count+i)
        else if(nfile == 35) then
            iclmn(i) = iclmn_35(total_contribution_count+i)
            indx(i) = indx_35(total_contribution_count+i)
            coeff(i) = coeff_35(total_contribution_count+i)
        else if(nfile == 36) then
            iclmn(i) = iclmn_36(total_contribution_count+i)
            indx(i) = indx_36(total_contribution_count+i)
            coeff(i) = coeff_36(total_contribution_count+i)
        else if(nfile == 37) then
            iclmn(i) = iclmn_37(total_contribution_count+i)
            indx(i) = indx_37(total_contribution_count+i)
            coeff(i) = coeff_37(total_contribution_count+i)
        else if(nfile == 38) then
            iclmn(i) = iclmn_38(total_contribution_count+i)
            indx(i) = indx_38(total_contribution_count+i)
            coeff(i) = coeff_38(total_contribution_count+i)
        else if(nfile == 39) then
            iclmn(i) = iclmn_39(total_contribution_count+i)
            indx(i) = indx_39(total_contribution_count+i)
            coeff(i) = coeff_39(total_contribution_count+i)
        else if(nfile == 40) then
            iclmn(i) = iclmn_40(total_contribution_count+i)
            indx(i) = indx_40(total_contribution_count+i)
            coeff(i) = coeff_40(total_contribution_count+i)
        else if(nfile == 41) then
            iclmn(i) = iclmn_41(total_contribution_count+i)
            indx(i) = indx_41(total_contribution_count+i)
            coeff(i) = coeff_41(total_contribution_count+i)
        else if(nfile == 42) then
            iclmn(i) = iclmn_42(total_contribution_count+i)
            indx(i) = indx_42(total_contribution_count+i)
            coeff(i) = coeff_42(total_contribution_count+i)
        else if(nfile == 43) then
            iclmn(i) = iclmn_43(total_contribution_count+i)
            indx(i) = indx_43(total_contribution_count+i)
            coeff(i) = coeff_43(total_contribution_count+i)
        else if(nfile == 44) then
            iclmn(i) = iclmn_44(total_contribution_count+i)
            indx(i) = indx_44(total_contribution_count+i)
            coeff(i) = coeff_44(total_contribution_count+i)
        else if(nfile == 45) then
            iclmn(i) = iclmn_45(total_contribution_count+i)
            indx(i) = indx_45(total_contribution_count+i)
            coeff(i) = coeff_45(total_contribution_count+i)
        else if(nfile == 46) then
            iclmn(i) = iclmn_46(total_contribution_count+i)
            indx(i) = indx_46(total_contribution_count+i)
            coeff(i) = coeff_46(total_contribution_count+i)
        else if(nfile == 47) then
            iclmn(i) = iclmn_47(total_contribution_count+i)
            indx(i) = indx_47(total_contribution_count+i)
            coeff(i) = coeff_47(total_contribution_count+i)
        else if(nfile == 48) then
            iclmn(i) = iclmn_48(total_contribution_count+i)
            indx(i) = indx_48(total_contribution_count+i)
            coeff(i) = coeff_48(total_contribution_count+i)
        else if(nfile == 49) then
            iclmn(i) = iclmn_49(total_contribution_count+i)
            indx(i) = indx_49(total_contribution_count+i)
            coeff(i) = coeff_49(total_contribution_count+i)
        else if(nfile == 50) then
            iclmn(i) = iclmn_50(total_contribution_count+i)
            indx(i) = indx_50(total_contribution_count+i)
            coeff(i) = coeff_50(total_contribution_count+i)
        else
            print*, "Error in read3_mem"
            stop
        end if
    end do

    return
end subroutine read_mcp_coefficients_into_memory
