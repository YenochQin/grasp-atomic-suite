module read_mcp_coefficients_into_memory_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_mcp_coefficients_into_memory
contains

    !*******************************************************************
    !                                                                  *
    subroutine read_mcp_coefficients_into_memory(mcp_channel,total_contribution_count,contribution_count)
        ! Module dependencies
        !                                                                  *
        !   Written by  G. Gaigalas               Vilnius, September 2021  *
        !                                                                  *
        !*******************************************************************
        !
        use, intrinsic :: iso_fortran_env, only: int64
        use rmcdhf_mcp_memory, only: coeff, angular_contribution_coefficients_mcp_31, angular_contribution_coefficients_mcp_32, &
            & angular_contribution_coefficients_mcp_33, angular_contribution_coefficients_mcp_34, &
                & angular_contribution_coefficients_mcp_35, angular_contribution_coefficients_mcp_36, &
                    & angular_contribution_coefficients_mcp_37, angular_contribution_coefficients_mcp_38, &
            & angular_contribution_coefficients_mcp_39, angular_contribution_coefficients_mcp_40, &
                & angular_contribution_coefficients_mcp_41, angular_contribution_coefficients_mcp_42, &
                    & angular_contribution_coefficients_mcp_43, angular_contribution_coefficients_mcp_44, &
                        & angular_contribution_coefficients_mcp_45, angular_contribution_coefficients_mcp_46, &
                            & angular_contribution_coefficients_mcp_47, angular_contribution_coefficients_mcp_48, &
            & angular_contribution_coefficients_mcp_49, angular_contribution_coefficients_mcp_50, iclmn, &
                & contribution_csf_column_indices_mcp_31, contribution_csf_column_indices_mcp_32, &
                    & contribution_csf_column_indices_mcp_33, contribution_csf_column_indices_mcp_34, &
                        & contribution_csf_column_indices_mcp_35, contribution_csf_column_indices_mcp_36, &
                            & contribution_csf_column_indices_mcp_37, &
            & contribution_csf_column_indices_mcp_38, contribution_csf_column_indices_mcp_39, &
                & contribution_csf_column_indices_mcp_40, contribution_csf_column_indices_mcp_41, &
                    & contribution_csf_column_indices_mcp_42, contribution_csf_column_indices_mcp_43, &
                        & contribution_csf_column_indices_mcp_44, contribution_csf_column_indices_mcp_45, &
                            & contribution_csf_column_indices_mcp_46, contribution_csf_column_indices_mcp_47, &
            & contribution_csf_column_indices_mcp_48, contribution_csf_column_indices_mcp_49, &
                & contribution_csf_column_indices_mcp_50, indx, contribution_sparse_matrix_indices_mcp_31, &
                    & contribution_sparse_matrix_indices_mcp_32, contribution_sparse_matrix_indices_mcp_33, &
                        & contribution_sparse_matrix_indices_mcp_34, contribution_sparse_matrix_indices_mcp_35, &
                            & contribution_sparse_matrix_indices_mcp_36, contribution_sparse_matrix_indices_mcp_37, &
            & contribution_sparse_matrix_indices_mcp_38, contribution_sparse_matrix_indices_mcp_39, &
                & contribution_sparse_matrix_indices_mcp_40, contribution_sparse_matrix_indices_mcp_41, &
                    & contribution_sparse_matrix_indices_mcp_42, contribution_sparse_matrix_indices_mcp_43, &
                        & contribution_sparse_matrix_indices_mcp_44, contribution_sparse_matrix_indices_mcp_45, &
                            & contribution_sparse_matrix_indices_mcp_46, contribution_sparse_matrix_indices_mcp_47, &
                                & contribution_sparse_matrix_indices_mcp_48, &
            & contribution_sparse_matrix_indices_mcp_49, contribution_sparse_matrix_indices_mcp_50
        implicit none (type, external)

        ! Arguments
        integer(int64), intent(in) :: total_contribution_count
        integer, intent(in) :: mcp_channel
        integer, intent(in) :: contribution_count

        ! Local declarations
        integer :: i
        !-----------------------------------------------

        do i = 1, contribution_count
            !*** Read column index, sparse matrix index, and
            !*** coefficient for all contributions from this integral
            if(mcp_channel == 31) then
                iclmn(i) = contribution_csf_column_indices_mcp_31(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_31(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_31(total_contribution_count+i)
            else if(mcp_channel == 32) then
                iclmn(i) = contribution_csf_column_indices_mcp_32(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_32(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_32(total_contribution_count+i)
            else if(mcp_channel == 33) then
                iclmn(i) = contribution_csf_column_indices_mcp_33(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_33(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_33(total_contribution_count+i)
            else if(mcp_channel == 34) then
                iclmn(i) = contribution_csf_column_indices_mcp_34(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_34(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_34(total_contribution_count+i)
            else if(mcp_channel == 35) then
                iclmn(i) = contribution_csf_column_indices_mcp_35(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_35(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_35(total_contribution_count+i)
            else if(mcp_channel == 36) then
                iclmn(i) = contribution_csf_column_indices_mcp_36(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_36(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_36(total_contribution_count+i)
            else if(mcp_channel == 37) then
                iclmn(i) = contribution_csf_column_indices_mcp_37(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_37(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_37(total_contribution_count+i)
            else if(mcp_channel == 38) then
                iclmn(i) = contribution_csf_column_indices_mcp_38(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_38(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_38(total_contribution_count+i)
            else if(mcp_channel == 39) then
                iclmn(i) = contribution_csf_column_indices_mcp_39(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_39(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_39(total_contribution_count+i)
            else if(mcp_channel == 40) then
                iclmn(i) = contribution_csf_column_indices_mcp_40(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_40(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_40(total_contribution_count+i)
            else if(mcp_channel == 41) then
                iclmn(i) = contribution_csf_column_indices_mcp_41(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_41(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_41(total_contribution_count+i)
            else if(mcp_channel == 42) then
                iclmn(i) = contribution_csf_column_indices_mcp_42(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_42(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_42(total_contribution_count+i)
            else if(mcp_channel == 43) then
                iclmn(i) = contribution_csf_column_indices_mcp_43(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_43(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_43(total_contribution_count+i)
            else if(mcp_channel == 44) then
                iclmn(i) = contribution_csf_column_indices_mcp_44(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_44(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_44(total_contribution_count+i)
            else if(mcp_channel == 45) then
                iclmn(i) = contribution_csf_column_indices_mcp_45(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_45(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_45(total_contribution_count+i)
            else if(mcp_channel == 46) then
                iclmn(i) = contribution_csf_column_indices_mcp_46(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_46(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_46(total_contribution_count+i)
            else if(mcp_channel == 47) then
                iclmn(i) = contribution_csf_column_indices_mcp_47(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_47(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_47(total_contribution_count+i)
            else if(mcp_channel == 48) then
                iclmn(i) = contribution_csf_column_indices_mcp_48(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_48(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_48(total_contribution_count+i)
            else if(mcp_channel == 49) then
                iclmn(i) = contribution_csf_column_indices_mcp_49(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_49(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_49(total_contribution_count+i)
            else if(mcp_channel == 50) then
                iclmn(i) = contribution_csf_column_indices_mcp_50(total_contribution_count+i)
                indx(i) = contribution_sparse_matrix_indices_mcp_50(total_contribution_count+i)
                coeff(i) = angular_contribution_coefficients_mcp_50(total_contribution_count+i)
            else
                print*, "Error in read_mcp_coefficients_into_memory"
                stop
            end if
        end do

        return
    end subroutine read_mcp_coefficients_into_memory
end module read_mcp_coefficients_into_memory_interface
