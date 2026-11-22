! Load the MCP structure, integral labels, and contributions into memory for reuse by SCF iterations.
module cache_mcp_files_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: cache_mcp_files
contains

    ! Written by  G. Gaigalas               Vilnius, September 2021
    subroutine cache_mcp_files(block_count,mcp_channel)
        ! Module dependencies
        use mcp_file_metadata, only: mcp_input_units
        use, intrinsic :: iso_fortran_env, only: real64, int64
        use array_memory_management, only: alloc
        use mpi_io_errors, only: fail_input_io_mpi
        use rmcdhf_mcp_memory, only: angular_contribution_coefficients_mcp_31, angular_contribution_coefficients_mcp_32, &
            & angular_contribution_coefficients_mcp_33, angular_contribution_coefficients_mcp_34, &
                & angular_contribution_coefficients_mcp_35, angular_contribution_coefficients_mcp_36, &
                    & angular_contribution_coefficients_mcp_37, angular_contribution_coefficients_mcp_38, &
            & angular_contribution_coefficients_mcp_39, angular_contribution_coefficients_mcp_40, &
                & angular_contribution_coefficients_mcp_41, angular_contribution_coefficients_mcp_42, &
                    & angular_contribution_coefficients_mcp_43, angular_contribution_coefficients_mcp_44, &
                        & angular_contribution_coefficients_mcp_45, angular_contribution_coefficients_mcp_46, &
                            & angular_contribution_coefficients_mcp_47, angular_contribution_coefficients_mcp_48, &
            & angular_contribution_coefficients_mcp_49, angular_contribution_coefficients_mcp_50, &
                & contribution_csf_column_indices_mcp_31, contribution_csf_column_indices_mcp_32, &
                    & contribution_csf_column_indices_mcp_33, contribution_csf_column_indices_mcp_34, &
                        & contribution_csf_column_indices_mcp_35, contribution_csf_column_indices_mcp_36, &
                            & contribution_csf_column_indices_mcp_37, contribution_csf_column_indices_mcp_38, &
            & contribution_csf_column_indices_mcp_39, contribution_csf_column_indices_mcp_40, &
                & contribution_csf_column_indices_mcp_41, contribution_csf_column_indices_mcp_42, &
                    & contribution_csf_column_indices_mcp_43, contribution_csf_column_indices_mcp_44, &
                        & contribution_csf_column_indices_mcp_45, contribution_csf_column_indices_mcp_46, &
                            & contribution_csf_column_indices_mcp_47, contribution_csf_column_indices_mcp_48, &
            & contribution_csf_column_indices_mcp_49, contribution_csf_column_indices_mcp_50, iendc_30, imax1_30, imax2_30, &
                & imin1_30, imin2_30, contribution_sparse_matrix_indices_mcp_31, contribution_sparse_matrix_indices_mcp_32, &
                    & contribution_sparse_matrix_indices_mcp_33, contribution_sparse_matrix_indices_mcp_34, &
            & contribution_sparse_matrix_indices_mcp_35, contribution_sparse_matrix_indices_mcp_36, &
                & contribution_sparse_matrix_indices_mcp_37, contribution_sparse_matrix_indices_mcp_38, &
                    & contribution_sparse_matrix_indices_mcp_39, contribution_sparse_matrix_indices_mcp_40, &
                        & contribution_sparse_matrix_indices_mcp_41, contribution_sparse_matrix_indices_mcp_42, &
                            & contribution_sparse_matrix_indices_mcp_43, contribution_sparse_matrix_indices_mcp_44, &
                                & contribution_sparse_matrix_indices_mcp_45, &
            & contribution_sparse_matrix_indices_mcp_46, contribution_sparse_matrix_indices_mcp_47, &
                & contribution_sparse_matrix_indices_mcp_48, contribution_sparse_matrix_indices_mcp_49, &
                    & contribution_sparse_matrix_indices_mcp_50, irow_30, lab_31, lab_32, lab_33, lab_34, lab_35, lab_36, &
            & lab_37, lab_38, lab_39, lab_40, lab_41, lab_42, lab_43, lab_44, lab_45, lab_46, lab_47, lab_48, lab_49, &
            & lab_50, integral_contribution_counts_by_label_mcp_31, integral_contribution_counts_by_label_mcp_32, &
                & integral_contribution_counts_by_label_mcp_33, integral_contribution_counts_by_label_mcp_34, &
                    & integral_contribution_counts_by_label_mcp_35, integral_contribution_counts_by_label_mcp_36, &
                        & integral_contribution_counts_by_label_mcp_37, integral_contribution_counts_by_label_mcp_38, &
            & integral_contribution_counts_by_label_mcp_39, integral_contribution_counts_by_label_mcp_40, &
                & integral_contribution_counts_by_label_mcp_41, integral_contribution_counts_by_label_mcp_42, &
                    & integral_contribution_counts_by_label_mcp_43, integral_contribution_counts_by_label_mcp_44, &
                        & integral_contribution_counts_by_label_mcp_45, integral_contribution_counts_by_label_mcp_46, &
                            & integral_contribution_counts_by_label_mcp_47, &
            & integral_contribution_counts_by_label_mcp_48, integral_contribution_counts_by_label_mcp_49, &
                & integral_contribution_counts_by_label_mcp_50
        use standard_io_units, only: stderr_unit
        use mpi_state, only: mpi_process_count, mpi_rank
        use sparse_hamiltonian_storage,          only: stored_matrix_element_count

        ! Procedure interfaces
        use query_total_system_memory_interface, only: query_total_system_memory
        use query_available_system_memory_interface, only: query_available_system_memory
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: block_count, mcp_channel

        ! Local declarations
        character :: mcp_channel_label*3
        real(real64) :: discarded_angular_coefficient
        real(real64) :: total_system_memory_gib, available_system_memory_gib
        integer :: i, symmetry_block_index, file_block_index, file_csf_count, file_coefficient_count, io_status, block_csf_count
        integer(int64) :: total_contribution_count, total_secondary_contribution_count
        integer :: discarded_csf_column_index, probed_integral_label, probed_integral_contribution_count, &
            & discarded_sparse_matrix_index
        integer :: file_matrix_element_count, discarded_column_end, discarded_row_index

        total_contribution_count = 0
        total_secondary_contribution_count = 0
        do symmetry_block_index = 1, block_count
            if(mcp_channel == 30) then
                ! Read in IROW from file mcp.30 ***
                read (mcp_input_units(30), iostat=io_status) mcp_channel_label, file_block_index, block_csf_count
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(30), &
                    "MCP cache_mcp_files block or channel header")
                if (file_block_index /= symmetry_block_index) then
                    write (stderr_unit, *) "cache_mcp_files", " branch 1: block index mismatch"
                    call fail_input_io_mpi(0, mcp_input_units(mcp_channel), "MCP cache_mcp_files invalid block or channel")
                end if
                read (mcp_input_units(30), iostat=io_status) file_matrix_element_count
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(30), &
                    "MCP cache_mcp_files sparse element count")
                stored_matrix_element_count = int(file_matrix_element_count,kind=int64)
                if(symmetry_block_index == 1) then
                    call alloc (imin1_30,block_count,"IMIN1_30","cache_mcp_files")
                    call alloc (imin2_30,block_count,"IMIN2_30","cache_mcp_files")
                    call alloc (imax1_30,block_count,"IMAX1_30","cache_mcp_files")
                    call alloc (imax2_30,block_count,"IMAX2_30","cache_mcp_files")
                    imin1_30(symmetry_block_index) = 0
                    imax1_30(symmetry_block_index) = block_csf_count
                    imin2_30(symmetry_block_index) = 0
                    imax2_30(symmetry_block_index) = stored_matrix_element_count
                else
                    imin1_30(symmetry_block_index) = imax1_30(symmetry_block_index-1) + 0
                    imax1_30(symmetry_block_index) = imax1_30(symmetry_block_index-1) + block_csf_count
                    imin2_30(symmetry_block_index) = imax2_30(symmetry_block_index-1) + 0
                    imax2_30(symmetry_block_index) = imax2_30(symmetry_block_index-1) + stored_matrix_element_count
                end if
                read (mcp_input_units(30), iostat=io_status) (discarded_column_end,i=mpi_rank + &
                    & 1,block_csf_count,mpi_process_count), &
                    (discarded_row_index,i=1,stored_matrix_element_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(30), "MCP cache_mcp_files sparse indices")
            else
                read (mcp_input_units(mcp_channel), iostat=io_status) mcp_channel_label, file_block_index, file_csf_count, &
                    & file_coefficient_count
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files block or channel header")
                if (file_block_index /= symmetry_block_index) then
                    write (stderr_unit, *) "cache_mcp_files", " branch 2: block index mismatch"
                    call fail_input_io_mpi(0, mcp_input_units(mcp_channel), "MCP cache_mcp_files invalid block or channel")
                end if

                total_secondary_contribution_count = total_secondary_contribution_count + 1
                read (mcp_input_units(mcp_channel), iostat=io_status) probed_integral_label, probed_integral_contribution_count
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files count-pass label/count")
                do while (probed_integral_label /= 0 .or. probed_integral_contribution_count /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) (discarded_csf_column_index, &
                        & discarded_sparse_matrix_index,discarded_angular_coefficient,             &
                        i=1,probed_integral_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files count-pass coefficients")
                    total_contribution_count = total_contribution_count + probed_integral_contribution_count
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) probed_integral_label, probed_integral_contribution_count
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files count-pass label/count")
                end do
            end if
        end do

        if(mcp_channel == 30) then
            ! Total physical memory
            total_system_memory_gib = query_total_system_memory()  ! GB
            available_system_memory_gib = query_available_system_memory()  ! GB
            if(mpi_rank == 0) write(*,"(A53)")                                &
                "-----------------------------------------------------"
            if(mpi_rank == 0) write(*,"(A53)")                                &
                "Spin-angular coefficient are putting into the memory:"
            if(mpi_rank == 0) write(*,"(A53)")                                &
                "-----------------------------------------------------"
            if(mpi_rank == 0) write(*,"(A1)") ""
            if(mpi_rank == 0) write(*,"(A30,F10.2,A3)")                       &
                "Total memory on computer:",total_system_memory_gib," Gb"
            if(mpi_rank == 0) write(*,"(A30,F10.2,A3)")                       &
                "Free  memory on computer:",available_system_memory_gib," Gb"
            if(mpi_rank == 0) write(*,"(A1)") ""
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.30:"
            call alloc (iendc_30,imax1_30(block_count),"IENDC_30","cache_mcp_files")
            call alloc (irow_30, imax2_30(block_count),"IROW_30", "cache_mcp_files")
        else if(mcp_channel == 31) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.31:"
            call alloc (lab_31,   total_secondary_contribution_count,"LAB_31",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_31,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_31","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_31, total_contribution_count, &
                "contribution_csf_column_indices_mcp_31", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_31,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_31",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_31, total_contribution_count, &
                "angular_contribution_coefficients_mcp_31", "cache_mcp_files")
        else if(mcp_channel == 32) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.32:"
            call alloc (lab_32,   total_secondary_contribution_count,"LAB_32",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_32,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_32","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_32, total_contribution_count, &
                "contribution_csf_column_indices_mcp_32", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_32,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_32",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_32, total_contribution_count, &
                "angular_contribution_coefficients_mcp_32", "cache_mcp_files")
        else if(mcp_channel == 33) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.33:"
            call alloc (lab_33,   total_secondary_contribution_count,"LAB_33",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_33,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_33","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_33, total_contribution_count, &
                "contribution_csf_column_indices_mcp_33", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_33,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_33",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_33, total_contribution_count, &
                "angular_contribution_coefficients_mcp_33", "cache_mcp_files")
        else if(mcp_channel == 34) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.34:"
            call alloc (lab_34,   total_secondary_contribution_count,"LAB_34",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_34,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_34","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_34, total_contribution_count, &
                "contribution_csf_column_indices_mcp_34", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_34,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_34",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_34, total_contribution_count, &
                "angular_contribution_coefficients_mcp_34", "cache_mcp_files")
        else if(mcp_channel == 35) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.35:"
            call alloc (lab_35,   total_secondary_contribution_count,"LAB_35",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_35,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_35","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_35, total_contribution_count, &
                "contribution_csf_column_indices_mcp_35", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_35,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_35",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_35, total_contribution_count, &
                "angular_contribution_coefficients_mcp_35", "cache_mcp_files")
        else if(mcp_channel == 36) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.36:"
            call alloc (lab_36,   total_secondary_contribution_count,"LAB_36",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_36,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_36","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_36, total_contribution_count, &
                "contribution_csf_column_indices_mcp_36", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_36,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_36",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_36, total_contribution_count, &
                "angular_contribution_coefficients_mcp_36", "cache_mcp_files")
        else if(mcp_channel == 37) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.37:"
            call alloc (lab_37,   total_secondary_contribution_count,"LAB_37",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_37,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_37","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_37, total_contribution_count, &
                "contribution_csf_column_indices_mcp_37", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_37,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_37",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_37, total_contribution_count, &
                "angular_contribution_coefficients_mcp_37", "cache_mcp_files")
        else if(mcp_channel == 38) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.38:"
            call alloc (lab_38,   total_secondary_contribution_count,"LAB_38",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_38,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_38","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_38, total_contribution_count, &
                "contribution_csf_column_indices_mcp_38", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_38,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_38",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_38, total_contribution_count, &
                "angular_contribution_coefficients_mcp_38", "cache_mcp_files")
        else if(mcp_channel == 39) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.39:"
            call alloc (lab_39,   total_secondary_contribution_count,"LAB_39",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_39,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_39","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_39, total_contribution_count, &
                "contribution_csf_column_indices_mcp_39", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_39,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_39",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_39, total_contribution_count, &
                "angular_contribution_coefficients_mcp_39", "cache_mcp_files")
        else if(mcp_channel == 40) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.40:"
            call alloc (lab_40,   total_secondary_contribution_count,"LAB_40",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_40,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_40","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_40, total_contribution_count, &
                "contribution_csf_column_indices_mcp_40", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_40,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_40",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_40, total_contribution_count, &
                "angular_contribution_coefficients_mcp_40", "cache_mcp_files")
        else if(mcp_channel == 41) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.41:"
            call alloc (lab_41,   total_secondary_contribution_count,"LAB_41",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_41,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_41","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_41, total_contribution_count, &
                "contribution_csf_column_indices_mcp_41", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_41,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_41",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_41, total_contribution_count, &
                "angular_contribution_coefficients_mcp_41", "cache_mcp_files")
        else if(mcp_channel == 42) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.42:"
            call alloc (lab_42,   total_secondary_contribution_count,"LAB_42",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_42,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_42","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_42, total_contribution_count, &
                "contribution_csf_column_indices_mcp_42", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_42,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_42",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_42, total_contribution_count, &
                "angular_contribution_coefficients_mcp_42", "cache_mcp_files")
        else if(mcp_channel == 43) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.43:"
            call alloc (lab_43,   total_secondary_contribution_count,"LAB_43",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_43,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_43","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_43, total_contribution_count, &
                "contribution_csf_column_indices_mcp_43", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_43,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_43",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_43, total_contribution_count, &
                "angular_contribution_coefficients_mcp_43", "cache_mcp_files")
        else if(mcp_channel == 44) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.44:"
            call alloc (lab_44,   total_secondary_contribution_count,"LAB_44",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_44,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_44","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_44, total_contribution_count, &
                "contribution_csf_column_indices_mcp_44", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_44,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_44",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_44, total_contribution_count, &
                "angular_contribution_coefficients_mcp_44", "cache_mcp_files")
        else if(mcp_channel == 45) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.45:"
            call alloc (lab_45,   total_secondary_contribution_count,"LAB_45",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_45,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_45","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_45, total_contribution_count, &
                "contribution_csf_column_indices_mcp_45", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_45,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_45",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_45, total_contribution_count, &
                "angular_contribution_coefficients_mcp_45", "cache_mcp_files")
        else if(mcp_channel == 46) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.46:"
            call alloc (lab_46,   total_secondary_contribution_count,"LAB_46",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_46,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_46","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_46, total_contribution_count, &
                "contribution_csf_column_indices_mcp_46", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_46,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_46",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_46, total_contribution_count, &
                "angular_contribution_coefficients_mcp_46", "cache_mcp_files")
        else if(mcp_channel == 47) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.47:"
            call alloc (lab_47,   total_secondary_contribution_count,"LAB_47",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_47,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_47","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_47, total_contribution_count, &
                "contribution_csf_column_indices_mcp_47", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_47,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_47",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_47, total_contribution_count, &
                "angular_contribution_coefficients_mcp_47", "cache_mcp_files")
        else if(mcp_channel == 48) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.48:"
            call alloc (lab_48,   total_secondary_contribution_count,"LAB_48",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_48,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_48","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_48, total_contribution_count, &
                "contribution_csf_column_indices_mcp_48", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_48,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_48",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_48, total_contribution_count, &
                "angular_contribution_coefficients_mcp_48", "cache_mcp_files")
        else if(mcp_channel == 49) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.49:"
            call alloc (lab_49,   total_secondary_contribution_count,"LAB_49",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_49,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_49","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_49, total_contribution_count, &
                "contribution_csf_column_indices_mcp_49", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_49,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_49",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_49, total_contribution_count, &
                "angular_contribution_coefficients_mcp_49", "cache_mcp_files")
        else if(mcp_channel == 50) then
            if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.50:"
            call alloc (lab_50,   total_secondary_contribution_count,"LAB_50",   "cache_mcp_files")
            call alloc (integral_contribution_counts_by_label_mcp_50,total_secondary_contribution_count, &
                "integral_contribution_counts_by_label_mcp_50","cache_mcp_files")
            call alloc (contribution_csf_column_indices_mcp_50, total_contribution_count, &
                "contribution_csf_column_indices_mcp_50", "cache_mcp_files")
            call alloc (contribution_sparse_matrix_indices_mcp_50,  total_contribution_count, &
                "contribution_sparse_matrix_indices_mcp_50",  "cache_mcp_files")
            call alloc (angular_contribution_coefficients_mcp_50, total_contribution_count, &
                "angular_contribution_coefficients_mcp_50", "cache_mcp_files")
        else
            print*, "Error in cache_mcp_files"
            call fail_input_io_mpi(0, mcp_input_units(mcp_channel), "MCP cache_mcp_files invalid block or channel")
        end if
        rewind (mcp_input_units(mcp_channel), iostat=io_status)
        if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), "MCP cache_mcp_files rewind")
        if(mcp_channel == 30) then
            read (mcp_input_units(mcp_channel), iostat=io_status)
            if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                "MCP cache_mcp_files skip header or payload")
            read (mcp_input_units(mcp_channel), iostat=io_status)
            if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                "MCP cache_mcp_files skip header or payload")
            read (mcp_input_units(mcp_channel), iostat=io_status)
            if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                "MCP cache_mcp_files skip header or payload")
        end if
        read (mcp_input_units(mcp_channel), iostat=io_status)
        if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
            "MCP cache_mcp_files skip header or payload")
        read (mcp_input_units(mcp_channel), iostat=io_status)
        if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
            "MCP cache_mcp_files skip header or payload")
        read (mcp_input_units(mcp_channel), iostat=io_status)
        if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
            "MCP cache_mcp_files skip header or payload")
        total_contribution_count = 0
        total_secondary_contribution_count = 0
        do symmetry_block_index = 1, block_count
            if(mcp_channel == 30) then
                ! Read in IROW from file mcp.30 ***
                read (mcp_input_units(30), iostat=io_status) mcp_channel_label, file_block_index, block_csf_count
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(30), &
                    "MCP cache_mcp_files block or channel header")
                if (file_block_index /= symmetry_block_index) then
                    write (stderr_unit, *) "cache_mcp_files", " branch 3: block index mismatch"
                    call fail_input_io_mpi(0, mcp_input_units(mcp_channel), "MCP cache_mcp_files invalid block or channel")
                end if
            else
                read (mcp_input_units(mcp_channel), iostat=io_status) mcp_channel_label, file_block_index, file_csf_count, &
                    & file_coefficient_count
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files block or channel header")
                if (file_block_index /= symmetry_block_index) then
                    write (stderr_unit, *) "cache_mcp_files", " branch 4: block index mismatch"
                    call fail_input_io_mpi(0, mcp_input_units(mcp_channel), "MCP cache_mcp_files invalid block or channel")
                end if
                total_secondary_contribution_count = total_secondary_contribution_count + 1
            end if

            if(mcp_channel == 30) then
                read (mcp_input_units(30), iostat=io_status) file_matrix_element_count
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(30), &
                    "MCP cache_mcp_files sparse element count")
                stored_matrix_element_count = int(file_matrix_element_count,kind=int64)
                read (mcp_input_units(30), iostat=io_status) (iendc_30(i),i=                                  &
                    imin1_30(symmetry_block_index)+mpi_rank+1,imax1_30(symmetry_block_index),mpi_process_count), &
                    (irow_30(i),i=imin2_30(symmetry_block_index)+1,imax2_30(symmetry_block_index))
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(30), "MCP cache_mcp_files sparse indices")
            else if(mcp_channel == 31) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_31(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_31(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_31(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_31(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_31(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_31(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_31(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_31(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_31(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_31(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_31(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else if(mcp_channel == 32) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_32(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_32(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_32(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_32(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_32(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_32(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_32(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_32(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_32(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_32(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_32(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else if(mcp_channel == 33) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_33(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_33(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_33(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_33(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_33(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_33(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_33(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_33(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_33(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_33(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_33(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else if(mcp_channel == 34) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_34(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_34(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_34(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_34(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_34(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_34(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_34(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_34(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_34(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_34(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_34(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else if(mcp_channel == 35) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_35(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_35(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_35(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_35(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_35(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_35(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_35(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_35(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_35(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_35(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_35(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else if(mcp_channel == 36) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_36(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_36(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_36(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_36(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_36(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_36(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_36(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_36(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_36(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_36(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_36(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else if(mcp_channel == 37) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_37(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_37(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_37(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_37(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_37(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_37(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_37(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_37(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_37(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_37(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_37(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else if(mcp_channel == 38) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_38(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_38(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_38(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_38(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_38(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_38(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_38(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_38(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_38(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_38(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_38(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else if(mcp_channel == 39) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_39(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_39(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_39(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_39(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_39(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_39(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_39(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_39(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_39(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_39(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_39(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else if(mcp_channel == 40) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_40(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_40(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_40(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_40(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_40(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_40(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_40(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_40(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_40(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_40(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_40(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else if(mcp_channel == 41) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_41(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_41(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_41(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_41(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_41(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_41(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_41(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_41(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_41(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_41(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_41(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else if(mcp_channel == 42) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_42(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_42(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_42(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_42(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_42(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_42(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_42(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_42(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_42(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_42(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_42(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else if(mcp_channel == 43) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_43(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_43(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_43(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_43(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_43(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_43(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_43(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_43(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_43(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_43(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_43(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else if(mcp_channel == 44) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_44(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_44(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_44(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_44(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_44(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_44(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_44(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_44(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_44(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_44(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_44(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else if(mcp_channel == 45) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_45(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_45(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_45(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_45(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_45(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_45(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_45(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_45(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_45(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_45(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_45(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else if(mcp_channel == 46) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_46(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_46(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_46(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_46(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_46(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_46(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_46(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_46(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_46(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_46(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_46(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else if(mcp_channel == 47) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_47(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_47(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_47(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_47(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_47(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_47(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_47(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_47(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_47(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_47(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_47(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else if(mcp_channel == 48) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_48(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_48(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_48(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_48(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_48(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_48(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_48(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_48(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_48(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_48(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_48(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else if(mcp_channel == 49) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_49(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_49(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_49(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_49(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_49(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_49(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_49(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_49(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_49(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_49(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_49(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else if(mcp_channel == 50) then
                read (mcp_input_units(mcp_channel), iostat=io_status) lab_50(total_secondary_contribution_count),              &
                    integral_contribution_counts_by_label_mcp_50(total_secondary_contribution_count)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                    & "MCP cache_mcp_files load-pass label/count")
                do while (lab_50(total_secondary_contribution_count) /= 0 .or.                    &
                    integral_contribution_counts_by_label_mcp_50(total_secondary_contribution_count) /= 0)
                    ! Read the column, sparse matrix index, and coefficient for each contribution to this integral.
                    read (mcp_input_units(mcp_channel), iostat=io_status) &
                        & (contribution_csf_column_indices_mcp_50(total_contribution_count+i),                   &
                        contribution_sparse_matrix_indices_mcp_50(total_contribution_count+i),                   &
                        angular_contribution_coefficients_mcp_50(total_contribution_count+i),                   &
                        i=1,integral_contribution_counts_by_label_mcp_50(total_secondary_contribution_count))
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files coefficient payload")
                    total_contribution_count  = total_contribution_count + &
                        & integral_contribution_counts_by_label_mcp_50(total_secondary_contribution_count)
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    read (mcp_input_units(mcp_channel), iostat=io_status) lab_50(total_secondary_contribution_count),           &
                        integral_contribution_counts_by_label_mcp_50(total_secondary_contribution_count)
                    if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(mcp_channel), &
                        & "MCP cache_mcp_files load-pass label/count")
                end do
            else
                print*, "Error in cache_mcp_files"
                call fail_input_io_mpi(0, mcp_input_units(mcp_channel), "MCP cache_mcp_files invalid block or channel")
            end if

        end do
        available_system_memory_gib = query_available_system_memory()  ! GB
        if(mpi_rank == 0)                                                    &
            write(*,"(A30,F10.2,A3)") "Free memory on computer",available_system_memory_gib," Gb"
        close(mcp_input_units(mcp_channel))

        return
    end subroutine cache_mcp_files
end module cache_mcp_files_interface
