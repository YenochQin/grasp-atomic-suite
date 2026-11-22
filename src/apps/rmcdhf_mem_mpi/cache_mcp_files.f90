! Populate the in-memory MCDHF solver's arrays from one indexed MCP shard.
module cache_mcp_files_interface
    use iso_fortran_env, only: int64, real64
    use mpi_indexed_mcp, only: mcp_group, mcp_segment_size, contribution_bytes, &
        read_indexed_sparse, read_indexed_group, read_indexed_contributions, close_indexed_mcp
    use mpi_mcp_calculation_input, only: indexed_input
    use mpi_state, only: mpi_rank, mpi_process_count
    use mpi_io_errors, only: fail_input_io_mpi
    use rmcdhf_mcp_memory, only: imin1_30, imax1_30, imin2_30, imax2_30, iendc_30, irow_30
    use rmcdhf_mcp_memory, only: lab_31, integral_contribution_counts_by_label_mcp_31, &
        contribution_csf_column_indices_mcp_31, contribution_sparse_matrix_indices_mcp_31, &
        angular_contribution_coefficients_mcp_31
    use rmcdhf_mcp_memory, only: lab_32, integral_contribution_counts_by_label_mcp_32, &
        contribution_csf_column_indices_mcp_32, contribution_sparse_matrix_indices_mcp_32, &
        angular_contribution_coefficients_mcp_32
    use rmcdhf_mcp_memory, only: lab_33, integral_contribution_counts_by_label_mcp_33, &
        contribution_csf_column_indices_mcp_33, contribution_sparse_matrix_indices_mcp_33, &
        angular_contribution_coefficients_mcp_33
    use rmcdhf_mcp_memory, only: lab_34, integral_contribution_counts_by_label_mcp_34, &
        contribution_csf_column_indices_mcp_34, contribution_sparse_matrix_indices_mcp_34, &
        angular_contribution_coefficients_mcp_34
    use rmcdhf_mcp_memory, only: lab_35, integral_contribution_counts_by_label_mcp_35, &
        contribution_csf_column_indices_mcp_35, contribution_sparse_matrix_indices_mcp_35, &
        angular_contribution_coefficients_mcp_35
    use rmcdhf_mcp_memory, only: lab_36, integral_contribution_counts_by_label_mcp_36, &
        contribution_csf_column_indices_mcp_36, contribution_sparse_matrix_indices_mcp_36, &
        angular_contribution_coefficients_mcp_36
    use rmcdhf_mcp_memory, only: lab_37, integral_contribution_counts_by_label_mcp_37, &
        contribution_csf_column_indices_mcp_37, contribution_sparse_matrix_indices_mcp_37, &
        angular_contribution_coefficients_mcp_37
    use rmcdhf_mcp_memory, only: lab_38, integral_contribution_counts_by_label_mcp_38, &
        contribution_csf_column_indices_mcp_38, contribution_sparse_matrix_indices_mcp_38, &
        angular_contribution_coefficients_mcp_38
    use rmcdhf_mcp_memory, only: lab_39, integral_contribution_counts_by_label_mcp_39, &
        contribution_csf_column_indices_mcp_39, contribution_sparse_matrix_indices_mcp_39, &
        angular_contribution_coefficients_mcp_39
    use rmcdhf_mcp_memory, only: lab_40, integral_contribution_counts_by_label_mcp_40, &
        contribution_csf_column_indices_mcp_40, contribution_sparse_matrix_indices_mcp_40, &
        angular_contribution_coefficients_mcp_40
    use rmcdhf_mcp_memory, only: lab_41, integral_contribution_counts_by_label_mcp_41, &
        contribution_csf_column_indices_mcp_41, contribution_sparse_matrix_indices_mcp_41, &
        angular_contribution_coefficients_mcp_41
    use rmcdhf_mcp_memory, only: lab_42, integral_contribution_counts_by_label_mcp_42, &
        contribution_csf_column_indices_mcp_42, contribution_sparse_matrix_indices_mcp_42, &
        angular_contribution_coefficients_mcp_42
    use rmcdhf_mcp_memory, only: lab_43, integral_contribution_counts_by_label_mcp_43, &
        contribution_csf_column_indices_mcp_43, contribution_sparse_matrix_indices_mcp_43, &
        angular_contribution_coefficients_mcp_43
    use rmcdhf_mcp_memory, only: lab_44, integral_contribution_counts_by_label_mcp_44, &
        contribution_csf_column_indices_mcp_44, contribution_sparse_matrix_indices_mcp_44, &
        angular_contribution_coefficients_mcp_44
    use rmcdhf_mcp_memory, only: lab_45, integral_contribution_counts_by_label_mcp_45, &
        contribution_csf_column_indices_mcp_45, contribution_sparse_matrix_indices_mcp_45, &
        angular_contribution_coefficients_mcp_45
    use rmcdhf_mcp_memory, only: lab_46, integral_contribution_counts_by_label_mcp_46, &
        contribution_csf_column_indices_mcp_46, contribution_sparse_matrix_indices_mcp_46, &
        angular_contribution_coefficients_mcp_46
    use rmcdhf_mcp_memory, only: lab_47, integral_contribution_counts_by_label_mcp_47, &
        contribution_csf_column_indices_mcp_47, contribution_sparse_matrix_indices_mcp_47, &
        angular_contribution_coefficients_mcp_47
    use rmcdhf_mcp_memory, only: lab_48, integral_contribution_counts_by_label_mcp_48, &
        contribution_csf_column_indices_mcp_48, contribution_sparse_matrix_indices_mcp_48, &
        angular_contribution_coefficients_mcp_48
    use rmcdhf_mcp_memory, only: lab_49, integral_contribution_counts_by_label_mcp_49, &
        contribution_csf_column_indices_mcp_49, contribution_sparse_matrix_indices_mcp_49, &
        angular_contribution_coefficients_mcp_49
    use rmcdhf_mcp_memory, only: lab_50, integral_contribution_counts_by_label_mcp_50, &
        contribution_csf_column_indices_mcp_50, contribution_sparse_matrix_indices_mcp_50, &
        angular_contribution_coefficients_mcp_50
    implicit none (type, external)
    private
    public :: cache_mcp_files
contains

    subroutine cache_mcp_files()
        integer :: block, channel, blocks, status, column_start, column_end
        integer(int64) :: rows, columns, row_end, column_offset

        blocks = indexed_input%block_count
        if (indexed_input%max_rank > 18) &
            call fail_input_io_mpi(0,indexed_input%unit,"MCP in-memory solver supports up to V(18)")
        rows = sum(indexed_input%row_counts)
        columns = sum(indexed_input%block_sizes)
        if (max(rows,columns) >= huge(0)) &
            call fail_input_io_mpi(0,indexed_input%unit,"MCP in-memory sparse cache exceeds 32-bit offsets")
        allocate(imin1_30(blocks),imax1_30(blocks),imin2_30(blocks),imax2_30(blocks),iendc_30(columns),irow_30(rows),stat=status)
        if (status /= 0) call fail_input_io_mpi(status,indexed_input%unit,"MCP allocate in-memory sparse cache")
        row_end = 0
        column_offset = 0
        iendc_30 = 0
        do block = 1, blocks
            imin1_30(block) = int(column_offset)
            imin2_30(block) = int(row_end)
            column_offset = column_offset + indexed_input%block_sizes(block)
            row_end = row_end + indexed_input%row_counts(block)
            imax1_30(block) = int(column_offset)
            imax2_30(block) = int(row_end)
            column_end = imax1_30(block)
            column_start = column_end + 1
            if (indexed_input%column_counts(block) > 0) column_start = imin1_30(block) + mpi_rank + 1
            call read_indexed_sparse(indexed_input,block,iendc_30(column_start:column_end:mpi_process_count), &
                irow_30(imin2_30(block)+1:imax2_30(block)))
        end do
        do channel = 1, indexed_input%max_rank+2
            select case (channel)
            case (1)
                call load_channel(channel,lab_31,integral_contribution_counts_by_label_mcp_31, &
                    contribution_csf_column_indices_mcp_31,contribution_sparse_matrix_indices_mcp_31, &
                    angular_contribution_coefficients_mcp_31)
            case (2)
                call load_channel(channel,lab_32,integral_contribution_counts_by_label_mcp_32, &
                    contribution_csf_column_indices_mcp_32,contribution_sparse_matrix_indices_mcp_32, &
                    angular_contribution_coefficients_mcp_32)
            case (3)
                call load_channel(channel,lab_33,integral_contribution_counts_by_label_mcp_33, &
                    contribution_csf_column_indices_mcp_33,contribution_sparse_matrix_indices_mcp_33, &
                    angular_contribution_coefficients_mcp_33)
            case (4)
                call load_channel(channel,lab_34,integral_contribution_counts_by_label_mcp_34, &
                    contribution_csf_column_indices_mcp_34,contribution_sparse_matrix_indices_mcp_34, &
                    angular_contribution_coefficients_mcp_34)
            case (5)
                call load_channel(channel,lab_35,integral_contribution_counts_by_label_mcp_35, &
                    contribution_csf_column_indices_mcp_35,contribution_sparse_matrix_indices_mcp_35, &
                    angular_contribution_coefficients_mcp_35)
            case (6)
                call load_channel(channel,lab_36,integral_contribution_counts_by_label_mcp_36, &
                    contribution_csf_column_indices_mcp_36,contribution_sparse_matrix_indices_mcp_36, &
                    angular_contribution_coefficients_mcp_36)
            case (7)
                call load_channel(channel,lab_37,integral_contribution_counts_by_label_mcp_37, &
                    contribution_csf_column_indices_mcp_37,contribution_sparse_matrix_indices_mcp_37, &
                    angular_contribution_coefficients_mcp_37)
            case (8)
                call load_channel(channel,lab_38,integral_contribution_counts_by_label_mcp_38, &
                    contribution_csf_column_indices_mcp_38,contribution_sparse_matrix_indices_mcp_38, &
                    angular_contribution_coefficients_mcp_38)
            case (9)
                call load_channel(channel,lab_39,integral_contribution_counts_by_label_mcp_39, &
                    contribution_csf_column_indices_mcp_39,contribution_sparse_matrix_indices_mcp_39, &
                    angular_contribution_coefficients_mcp_39)
            case (10)
                call load_channel(channel,lab_40,integral_contribution_counts_by_label_mcp_40, &
                    contribution_csf_column_indices_mcp_40,contribution_sparse_matrix_indices_mcp_40, &
                    angular_contribution_coefficients_mcp_40)
            case (11)
                call load_channel(channel,lab_41,integral_contribution_counts_by_label_mcp_41, &
                    contribution_csf_column_indices_mcp_41,contribution_sparse_matrix_indices_mcp_41, &
                    angular_contribution_coefficients_mcp_41)
            case (12)
                call load_channel(channel,lab_42,integral_contribution_counts_by_label_mcp_42, &
                    contribution_csf_column_indices_mcp_42,contribution_sparse_matrix_indices_mcp_42, &
                    angular_contribution_coefficients_mcp_42)
            case (13)
                call load_channel(channel,lab_43,integral_contribution_counts_by_label_mcp_43, &
                    contribution_csf_column_indices_mcp_43,contribution_sparse_matrix_indices_mcp_43, &
                    angular_contribution_coefficients_mcp_43)
            case (14)
                call load_channel(channel,lab_44,integral_contribution_counts_by_label_mcp_44, &
                    contribution_csf_column_indices_mcp_44,contribution_sparse_matrix_indices_mcp_44, &
                    angular_contribution_coefficients_mcp_44)
            case (15)
                call load_channel(channel,lab_45,integral_contribution_counts_by_label_mcp_45, &
                    contribution_csf_column_indices_mcp_45,contribution_sparse_matrix_indices_mcp_45, &
                    angular_contribution_coefficients_mcp_45)
            case (16)
                call load_channel(channel,lab_46,integral_contribution_counts_by_label_mcp_46, &
                    contribution_csf_column_indices_mcp_46,contribution_sparse_matrix_indices_mcp_46, &
                    angular_contribution_coefficients_mcp_46)
            case (17)
                call load_channel(channel,lab_47,integral_contribution_counts_by_label_mcp_47, &
                    contribution_csf_column_indices_mcp_47,contribution_sparse_matrix_indices_mcp_47, &
                    angular_contribution_coefficients_mcp_47)
            case (18)
                call load_channel(channel,lab_48,integral_contribution_counts_by_label_mcp_48, &
                    contribution_csf_column_indices_mcp_48,contribution_sparse_matrix_indices_mcp_48, &
                    angular_contribution_coefficients_mcp_48)
            case (19)
                call load_channel(channel,lab_49,integral_contribution_counts_by_label_mcp_49, &
                    contribution_csf_column_indices_mcp_49,contribution_sparse_matrix_indices_mcp_49, &
                    angular_contribution_coefficients_mcp_49)
            case (20)
                call load_channel(channel,lab_50,integral_contribution_counts_by_label_mcp_50, &
                    contribution_csf_column_indices_mcp_50,contribution_sparse_matrix_indices_mcp_50, &
                    angular_contribution_coefficients_mcp_50)
            end select
        end do
        call close_indexed_mcp(indexed_input)
    end subroutine cache_mcp_files

    subroutine load_channel(channel, labels, counts, columns, positions, coefficients)
        integer, intent(in) :: channel
        integer, pointer, intent(out) :: labels(:), counts(:), columns(:), positions(:)
        real(real64), pointer, intent(out) :: coefficients(:)
        integer :: block, ngrp, ncoeff, status
        integer(int64) :: index, last_group, remaining, position, previous, group_index, stored, first, start, last
        type(mcp_group) :: group

        remaining = sum(indexed_input%group_counts(channel,:)) + indexed_input%block_count
        stored = sum(indexed_input%channel_counts(channel,:))
        ! The current in-memory kernels use default integers to advance through each channel.
        if (max(remaining,stored) >= huge(0)) &
            call fail_input_io_mpi(0,indexed_input%unit,"MCP in-memory channel exceeds 32-bit offsets")
        ngrp = int(remaining)
        ncoeff = int(stored)
        allocate(labels(ngrp),counts(ngrp),columns(ncoeff),positions(ncoeff),coefficients(ncoeff),stat=status)
        if (status /= 0) call fail_input_io_mpi(status,indexed_input%unit,"MCP allocate in-memory coefficient channel")
        group_index = 0
        stored = 0
        do block = 1, indexed_input%block_count
            previous = 0
            remaining = indexed_input%channel_counts(channel,block)
            last_group = indexed_input%first_groups(channel,block)+indexed_input%group_counts(channel,block)-1
            position = indexed_input%sparse_positions(block) + 4*(indexed_input%column_counts(block) + &
                indexed_input%row_counts(block)) + contribution_bytes*sum(indexed_input%channel_counts(:channel-1,block))
            do index = indexed_input%first_groups(channel,block), last_group
                call read_indexed_group(indexed_input,index,position,remaining,previous,group,channel==1)
                group_index = group_index + 1
                labels(group_index) = int(group%label)
                counts(group_index) = int(group%count)
                do first = 1, group%count, mcp_segment_size
                    start = stored + first
                    last = min(stored+group%count,start+mcp_segment_size-1)
                    call read_indexed_contributions(indexed_input,group,first,columns(start:last),positions(start:last), &
                        coefficients(start:last))
                    if (any(columns(start:last) <= mpi_rank) .or. any(columns(start:last) > indexed_input%block_sizes(block)) &
                        .or. any(positions(start:last) < 1) .or. any(positions(start:last) > indexed_input%row_counts(block))) &
                        call fail_input_io_mpi(0,indexed_input%unit,"MCP in-memory contribution index outside rank block")
                    if (any(mod(columns(start:last)-1,mpi_process_count) /= mpi_rank)) &
                        call fail_input_io_mpi(0,indexed_input%unit,"MCP in-memory contribution column owned by another rank")
                end do
                stored = stored + group%count
                remaining = remaining - group%count
                position = position + group%bytes
                previous = group%label
            end do
            if (remaining /= 0) call fail_input_io_mpi(0,indexed_input%unit,"MCP in-memory channel contribution total mismatch")
            group_index = group_index + 1
            labels(group_index) = 0
            counts(group_index) = 0
        end do
    end subroutine load_channel
end module cache_mcp_files_interface
