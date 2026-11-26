! Cache group contractions for the current EOL mixing coefficients and level weights.
module mcp_density_contractions
    use iso_fortran_env, only: int64, real64
    use grasp_parameters, only: orbital_index_radix
    use mcp_file_metadata, only: max_mcp_file_multipole_rank
    use hamiltonian_symmetry_blocks, only: symmetry_block_count
    use orbital_csf_basis, only: orbital_count
    use orbital_lagrange_multipliers, only: lagrange_multiplier_count, lagrange_orbital_pair_labels
    use fixed_orbital_flags, only: is_orbital_fixed
    use mpi_indexed_mcp, only: decode_mcp_label, mcp_segment_size, contribution_bytes
    use mpi_io_errors, only: fail_input_io_mpi
    use csf_density_matrix_element_interface, only: csf_density_matrix_element
    use mpi_mcp_calculation_input, only: indexed_input
    use rmcdhf_mcp_io, only: mcp_epoch, mcp_cursor, mcp_buffer, read_sparse, start_channel, next_group, skip_group, read_segment, &
        profile_channel, profile_cache, mcp_now, profile_stage, ensure_buffer, directory_buffer_bytes
    implicit none (type, external)
    private

    public :: begin_density_epoch, release_contractions, prepare_contractions, contraction_bounds, &
        contracted_labels, contracted_values, record_solver_workspace

    integer, allocatable :: contracted_labels(:)
    real(real64), allocatable :: contracted_values(:)
    integer(int64), allocatable :: group_bounds(:,:,:)
    integer(int64) :: solver_reserve = 0
    logical :: cache_ready = .false.
contains

    subroutine release_contractions()
        if (allocated(contracted_labels)) deallocate(contracted_labels)
        if (allocated(contracted_values)) deallocate(contracted_values)
        if (allocated(group_bounds)) deallocate(group_bounds)
        cache_ready = .false.
    end subroutine release_contractions

    subroutine begin_density_epoch()
        call release_contractions()
        mcp_epoch = mcp_epoch + 1
    end subroutine begin_density_epoch

    subroutine record_solver_workspace(csf_count, row_count, highest_root, mixing_size)
        integer, intent(in) :: csf_count, highest_root, mixing_size
        integer(int64), intent(in) :: row_count
        integer(int64) :: n, root, limit, submatrix, workspace
        n = csf_count
        root = highest_root
        limit = min(n,2*root+80)
        submatrix = min(n,4000_int64)
        ! Sparse matrix, mixing vectors, Davidson work arrays and packed LAPACK initialization.
        workspace = 12*row_count + 8*(n+1) + 16_int64*mixing_size
        if (n > 1) workspace = workspace + 8*(2*n*limit+2*limit**2+11*limit+root) + 24*n + 4*(7*limit+root) + &
            4*submatrix*(submatrix+1) + 8*submatrix*min(root,submatrix) + 96*submatrix
        solver_reserve = max(solver_reserve,workspace)
    end subroutine record_solver_workspace

    subroutine contraction_bounds(block, channel, first, last)
        integer, intent(in) :: block, channel
        integer(int64), intent(out) :: first, last
        first = group_bounds(1,channel,block)
        last = group_bounds(2,channel,block)
    end subroutine contraction_bounds

    real(real64) function contract_group(cursor, buffer, rows, eol) result(value)
        type(mcp_cursor), intent(inout) :: cursor
        type(mcp_buffer), intent(inout) :: buffer
        integer, intent(in) :: rows(:)
        logical, intent(in) :: eol
        integer :: length, i, row, column, block
        integer(int64) :: first
        real(real64) :: contribution
        block = cursor%block
        value = 0.0_real64
        do first = 1, cursor%group%count, mcp_segment_size
            length = int(min(int(mcp_segment_size,int64),cursor%group%count-first+1))
            call read_segment(cursor,buffer,first,length,size(rows,kind=int64))
            do i = 1, length
                row = rows(buffer%positions(i))
                column = buffer%columns(i)
                if (row < 1 .or. row > column) &
                    call fail_input_io_mpi(0,cursor%unit,"MCP contraction sparse row outside column")
                contribution = csf_density_matrix_element(eol,row,column,block)*buffer%values(i)
                if (row /= column) contribution = contribution + contribution
                value = value + contribution
            end do
        end do
    end function contract_group

    subroutine prepare_contractions(eol)
        logical, intent(in) :: eol
        integer :: channel, block, status, pair, label, orbitals(4), rows, columns, buffer_capacity
        integer(int64) :: structure_bytes, table_bytes, peak, total, bytes_used, row_capacity, column_capacity
        integer, allocatable :: sparse_rows(:), column_ends(:)
        logical :: needed(orbital_count), sparse_loaded
        real(real64) :: started
        type(mcp_cursor) :: cursor
        type(mcp_buffer) :: buffer

        if (cache_ready) return
        started = mcp_now()
        row_capacity = maxval(indexed_input%row_counts)
        column_capacity = maxval(indexed_input%column_counts)
        buffer_capacity = int(min(int(mcp_segment_size,int64),maxval(indexed_input%channel_counts)))
        allocate(sparse_rows(row_capacity),column_ends(column_capacity),contracted_labels(indexed_input%group_count), &
            contracted_values(indexed_input%group_count),group_bounds(2,31:32+max_mcp_file_multipole_rank,symmetry_block_count), &
            stat=status)
        if (status == 0) call ensure_buffer(buffer,buffer_capacity,status)
        if (status /= 0) call fail_input_io_mpi(status,indexed_input%unit,"MCP allocate density contraction workspace")
        structure_bytes = 4*(row_capacity+column_capacity)
        table_bytes = 12*indexed_input%group_count + 16_int64*(max_mcp_file_multipole_rank+2)*symmetry_block_count
        peak = max(solver_reserve,structure_bytes + table_bytes + contribution_bytes*buffer_capacity + directory_buffer_bytes)

        ! Ordinary and extra updates can use any free orbital. Multipliers also use both frozen/varied sides.
        needed = .not. is_orbital_fixed(:orbital_count)
        do pair = 1, lagrange_multiplier_count
            label = lagrange_orbital_pair_labels(pair)
            needed(mod(label,orbital_index_radix)) = .true.
            needed(label/orbital_index_radix) = .true.
        end do

        total = 0
        bytes_used = 0
        do block = 1, symmetry_block_count
            rows = int(indexed_input%row_counts(block))
            columns = int(indexed_input%column_counts(block))
            sparse_loaded = .false.
            do channel = 31, 32+max_mcp_file_multipole_rank
                group_bounds(1,channel,block) = total+1
                call start_channel(cursor,block,channel)
                do while (next_group(cursor))
                    label = int(cursor%group%label)
                    call decode_mcp_label(label,channel==31,orbitals)
                    if (.not. any(needed(orbitals))) then
                        call skip_group(cursor)
                        cycle
                    end if
                    if (.not. sparse_loaded) then
                        call read_sparse(block,column_ends(:columns),sparse_rows(:rows),"contraction")
                        bytes_used = max(bytes_used,4_int64*(int(rows,int64)+columns))
                        sparse_loaded = .true.
                    end if
                    total = total+1
                    contracted_labels(total) = label
                    ! The sum for an unselected block is still evaluated. This preserves signed zeros and IEEE behavior.
                    contracted_values(total) = contract_group(cursor,buffer,sparse_rows(:rows),eol)
                end do
                group_bounds(2,channel,block) = total
                call profile_channel("contraction",0,cursor)
            end do
        end do
        call profile_cache(structure_bytes,table_bytes,solver_reserve,peak,"allocated peak bounds")
        call profile_cache(bytes_used,12*total+16_int64*(max_mcp_file_multipole_rank+2)*symmetry_block_count, &
            solver_reserve,peak,"used")
        call profile_stage("contraction_scan",0,0,started)
        cache_ready = .true.
    end subroutine prepare_contractions
end module mcp_density_contractions
