module load_mcp_records_into_memory_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: load_mcp_records_into_memory
contains

    !*******************************************************************
    !                                                                  *
    subroutine load_mcp_records_into_memory(nfile, total_secondary_contribution_count, lab, ncontr)
        ! Module dependencies
        !                                                                  *
        !   Written by  G. Gaigalas               Vilnius, September 2021  *
        !                                                                  *
        !*******************************************************************
        !
        use, intrinsic :: iso_fortran_env, only: int64
        use rmcdhf_mcp_memory, only: lab_31, lab_32, lab_33, lab_34, lab_35, lab_36, lab_37, lab_38, lab_39, lab_40, &
            & lab_41, lab_42, lab_43, lab_44, lab_45, lab_46, lab_47, lab_48, lab_49, lab_50, &
                & integral_contribution_counts_by_label_mcp_31, integral_contribution_counts_by_label_mcp_32, &
            & integral_contribution_counts_by_label_mcp_33, integral_contribution_counts_by_label_mcp_34, &
                & integral_contribution_counts_by_label_mcp_35, integral_contribution_counts_by_label_mcp_36, &
                    & integral_contribution_counts_by_label_mcp_37, integral_contribution_counts_by_label_mcp_38, &
                        & integral_contribution_counts_by_label_mcp_39, integral_contribution_counts_by_label_mcp_40, &
                            & integral_contribution_counts_by_label_mcp_41, &
            & integral_contribution_counts_by_label_mcp_42, integral_contribution_counts_by_label_mcp_43, &
                & integral_contribution_counts_by_label_mcp_44, integral_contribution_counts_by_label_mcp_45, &
                    & integral_contribution_counts_by_label_mcp_46, integral_contribution_counts_by_label_mcp_47, &
                        & integral_contribution_counts_by_label_mcp_48, integral_contribution_counts_by_label_mcp_49, &
                            & integral_contribution_counts_by_label_mcp_50
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: nfile
        integer, intent(out) :: lab, ncontr
        integer(int64), intent(in) :: total_secondary_contribution_count
        !-----------------------------------------------

        if(nfile == 31) then
            lab    = lab_31(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_31(total_secondary_contribution_count)
        else if(nfile == 32) then
            lab    = lab_32(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_32(total_secondary_contribution_count)
        else if(nfile == 33) then
            lab    = lab_33(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_33(total_secondary_contribution_count)
        else if(nfile == 34) then
            lab    = lab_34(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_34(total_secondary_contribution_count)
        else if(nfile == 35) then
            lab    = lab_35(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_35(total_secondary_contribution_count)
        else if(nfile == 36) then
            lab    = lab_36(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_36(total_secondary_contribution_count)
        else if(nfile == 37) then
            lab    = lab_37(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_37(total_secondary_contribution_count)
        else if(nfile == 38) then
            lab    = lab_38(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_38(total_secondary_contribution_count)
        else if(nfile == 39) then
            lab    = lab_39(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_39(total_secondary_contribution_count)
        else if(nfile == 40) then
            lab    = lab_40(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_40(total_secondary_contribution_count)
        else if(nfile == 41) then
            lab    = lab_41(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_41(total_secondary_contribution_count)
        else if(nfile == 42) then
            lab    = lab_42(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_42(total_secondary_contribution_count)
        else if(nfile == 43) then
            lab    = lab_43(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_43(total_secondary_contribution_count)
        else if(nfile == 44) then
            lab    = lab_44(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_44(total_secondary_contribution_count)
        else if(nfile == 45) then
            lab    = lab_45(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_45(total_secondary_contribution_count)
        else if(nfile == 46) then
            lab    = lab_46(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_46(total_secondary_contribution_count)
        else if(nfile == 47) then
            lab    = lab_47(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_47(total_secondary_contribution_count)
        else if(nfile == 48) then
            lab    = lab_48(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_48(total_secondary_contribution_count)
        else if(nfile == 49) then
            lab    = lab_49(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_49(total_secondary_contribution_count)
        else if(nfile == 50) then
            lab    = lab_50(total_secondary_contribution_count)
            ncontr = integral_contribution_counts_by_label_mcp_50(total_secondary_contribution_count)
        else
            print*, "Error in read1_memi  NFILE=",nfile
            stop
        end if

        return
    end subroutine load_mcp_records_into_memory
end module load_mcp_records_into_memory_interface
