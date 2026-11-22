! Read the same sparse blocks and integral groups for the matrix and orbital equations.
module rmcdhf_mcp_io
    use iso_fortran_env, only: int64, real64
    use mpi_indexed_mcp, only: mcp_group, read_indexed_sparse, read_indexed_directory, &
        read_indexed_contributions, close_indexed_mcp, contribution_bytes, mcp_group_directory_bytes, decode_mcp_label, &
        validate_indexed_group
    use mpi_io_errors, only: fail_input_io_mpi
    use mpi_state, only: mpi_rank, mpi_process_count, mpi_wtime
    use mcp_file_metadata, only: max_mcp_file_multipole_rank
    use hamiltonian_symmetry_blocks, only: csf_counts_by_block
    use mpi_mcp_calculation_input, only: indexed_input
    implicit none (type, external)
    private

    public :: read_mcp_options, start_mcp_logs, close_mcp_input
    public :: mcp_cursor, mcp_buffer, read_sparse, start_channel, &
        next_group, skip_group, read_segment, profile_channel, profile_sparse, profile_indexed_matrix_skip, profile_stage, &
        profile_cache, log_bits, log_radial_grid, mcp_epoch, verification_enabled, logging_enabled, mcp_now, ensure_buffer

    logical :: logging_enabled = .false., verification_enabled = .false.
    integer(int64) :: mcp_epoch = 0
    integer :: trace_unit = -1, bits_unit = -1
    character(len=*), parameter :: io_context = "MCP assemble_orbital_potential_coefficients"

    integer, parameter :: directory_batch_size = 256
    integer(int64), parameter, public :: directory_buffer_bytes = directory_batch_size*mcp_group_directory_bytes
    type :: mcp_cursor
        integer :: unit, block, channel
        integer(int64) :: index, groups_left, position, remaining, previous_label
        integer(int64) :: groups_read, groups_skipped, contributions, skipped_contributions, bytes_read, bytes_skipped, max_group
        type(mcp_group) :: group, directory_batch(directory_batch_size)
        integer :: batch_index, batch_count
    end type
    type :: mcp_buffer
        integer, allocatable :: columns(:), positions(:)
        real(real64), allocatable :: values(:)
    end type
contains

    subroutine read_mcp_options()
        integer :: i, status
        character(len=64) :: argument
        i = 1
        do while (i <= command_argument_count())
            call get_command_argument(i,argument,status=status)
            if (status /= 0) call fail_input_io_mpi(status,-1,"MCP command-line option")
            select case (argument)
            case ("--mcp-trace")
                logging_enabled = .true.
            case ("--mcp-verify")
                verification_enabled = .true.
            case default
                call fail_input_io_mpi(0,-1,"MCP unknown option "//trim(argument))
            end select
            i = i + 1
        end do
    end subroutine read_mcp_options

    subroutine start_mcp_logs(prefix)
        character(*), intent(in) :: prefix
        integer :: status
        if (logging_enabled) then
            open(newunit=trace_unit,file=trim(prefix)//".io.csv",status="replace",action="write",iostat=status)
            if (status /= 0) call fail_input_io_mpi(status,trace_unit,"MCP open I/O log")
            write(trace_unit,"(A)") "phase,source,epoch,rank,orbital,block,channel,sparse_blocks,groups_read,groups_skipped,"// &
                "contributions,skipped_contributions,bytes_read,bytes_skipped,max_group,wall_seconds,structure_bytes,"// &
                "table_bytes,solver_bytes,peak_estimate_bytes,reason"
        end if
        if (verification_enabled) then
            open(newunit=bits_unit,file=trim(prefix)//".bits",status="replace",action="write",iostat=status)
            if (status /= 0) call fail_input_io_mpi(status,bits_unit,"MCP open numerical comparison log")
        end if
    end subroutine start_mcp_logs

    subroutine close_mcp_input()
        call close_indexed_mcp(indexed_input)
        if (logging_enabled) close(trace_unit)
        if (verification_enabled) close(bits_unit)
    end subroutine close_mcp_input

    real(real64) function mcp_now()
        mcp_now = 0.0_real64
        if (logging_enabled) mcp_now = mpi_wtime()
    end function mcp_now

    subroutine profile(phase, source, orbital, block, channel, sparsity, groups, skipped, contributions, skipped_contributions, &
        bytes, skipped_bytes, largest, elapsed, structure, table, solver, peak, reason)
        character(*), intent(in) :: phase, source, reason
        integer, intent(in) :: orbital, block, channel, sparsity
        integer(int64), intent(in) :: groups, skipped, contributions, skipped_contributions, bytes, skipped_bytes, largest, &
            structure, table, solver, peak
        real(real64), intent(in) :: elapsed
        if (.not. logging_enabled) return
        write(trace_unit,'(2(A,","),13(I0,","),ES16.8,4(",",I0),",",A)') phase,source,mcp_epoch,mpi_rank,orbital,block,channel, &
            sparsity,groups,skipped,contributions,skipped_contributions,bytes,skipped_bytes,largest,elapsed, &
            structure,table,solver,peak,reason
    end subroutine profile

    subroutine profile_stage(phase, orbital, block, started)
        character(*), intent(in) :: phase
        integer, intent(in) :: orbital, block
        real(real64), intent(in) :: started
        if (.not. logging_enabled) return
        call profile(phase,"",orbital,block,0,0,0_int64,0_int64,0_int64,0_int64,0_int64,0_int64,0_int64, &
            mpi_wtime()-started,0_int64,0_int64,0_int64,0_int64,"")
    end subroutine profile_stage

    subroutine profile_cache(structure, table, solver, peak, reason)
        integer(int64), intent(in) :: structure, table, solver, peak
        character(*), intent(in) :: reason
        if (.not. logging_enabled) return
        call profile("cache","",0,0,0,0,0_int64,0_int64,0_int64,0_int64,0_int64,0_int64,0_int64, &
            0.0_real64,structure,table,solver,peak,reason)
    end subroutine profile_cache

    subroutine profile_channel(source, orbital, cursor)
        character(*), intent(in) :: source
        integer, intent(in) :: orbital
        type(mcp_cursor), intent(in) :: cursor
        if (.not. logging_enabled) return
        call profile("channel",source,orbital,cursor%block,cursor%channel,0,cursor%groups_read,cursor%groups_skipped, &
            cursor%contributions,cursor%skipped_contributions,cursor%bytes_read,cursor%bytes_skipped,cursor%max_group, &
            0.0_real64,0_int64,0_int64,0_int64,0_int64,"")
    end subroutine profile_channel

    subroutine profile_indexed_matrix_skip(block)
        integer, intent(in) :: block
        integer :: channel
        integer(int64) :: groups, contributions, bytes
        if (.not. logging_enabled) return
        bytes = 4*(indexed_input%column_counts(block)+indexed_input%row_counts(block))
        call profile("sparse","matrix_skip",0,block,30,0,0_int64,0_int64,0_int64,0_int64,0_int64,bytes,0_int64, &
            0.0_real64,0_int64,0_int64,0_int64,0_int64,"")
        do channel = 1, indexed_input%max_rank+2
            groups = indexed_input%group_counts(channel,block)
            contributions = indexed_input%channel_counts(channel,block)
            bytes = contribution_bytes*contributions + mcp_group_directory_bytes*groups
            call profile("channel","matrix_skip",0,block,channel+30,0,0_int64,groups,0_int64,contributions,0_int64,bytes,0_int64, &
                0.0_real64,0_int64,0_int64,0_int64,0_int64,"")
        end do
    end subroutine profile_indexed_matrix_skip

    subroutine profile_sparse(source, orbital, block, bytes, columns, rows)
        character(*), intent(in) :: source
        integer, intent(in) :: orbital, block
        integer(int64), intent(in) :: bytes, columns, rows
        if (.not. logging_enabled) return
        call profile("sparse",source,orbital,block,30,1,0_int64,0_int64,rows,0_int64,bytes,0_int64,0_int64, &
            0.0_real64,4*(columns+rows),0_int64,0_int64,0_int64,"")
    end subroutine profile_sparse

    subroutine log_bits(phase, orbital, block, channel, label, value)
        character(*), intent(in) :: phase
        integer, intent(in) :: orbital, block, channel
        integer(int64), intent(in) :: label
        real(real64), intent(in) :: value
        if (verification_enabled) write(bits_unit,'(A,5(1X,I0),1X,Z16.16)') phase,mcp_epoch,orbital,block,channel,label, &
            transfer(value,0_int64)
    end subroutine log_bits

    subroutine log_radial_grid(phase, orbital, points)
        use scf_radial_potentials, only: orbital_y_potential, large_exchange_source, small_exchange_source
        character(*), intent(in) :: phase
        integer, intent(in) :: orbital, points
        integer :: i
        if (.not. verification_enabled) return
        do i = 1, points
            call log_bits(phase,orbital,0,1,int(i,int64),orbital_y_potential(i))
            call log_bits(phase,orbital,0,2,int(i,int64),large_exchange_source(i))
            call log_bits(phase,orbital,0,3,int(i,int64),small_exchange_source(i))
        end do
    end subroutine log_radial_grid

    subroutine read_sparse(block, ends, rows, source)
        integer, intent(in) :: block
        integer, intent(out) :: ends(:), rows(:)
        character(*), intent(in) :: source
        integer(int64) :: columns, positions
        columns = size(ends,kind=int64)
        positions = size(rows,kind=int64)
        call read_indexed_sparse(indexed_input,block,ends,rows)
        call profile_sparse(source,0,block,4*(columns+positions),columns,positions)
    end subroutine read_sparse

    subroutine start_channel(cursor, block, channel)
        type(mcp_cursor), intent(out) :: cursor
        integer, intent(in) :: block, channel
        cursor%unit = indexed_input%unit
        cursor%block = block
        cursor%channel = channel
        cursor%previous_label = 0
        cursor%groups_read = 0
        cursor%groups_skipped = 0
        cursor%contributions = 0
        cursor%skipped_contributions = 0
        cursor%bytes_read = 0
        cursor%bytes_skipped = 0
        cursor%max_group = 0
        cursor%batch_count = 0
        cursor%batch_index = 1
        cursor%index = indexed_input%first_groups(channel-30,block)
        cursor%groups_left = indexed_input%group_counts(channel-30,block)
        cursor%remaining = indexed_input%channel_counts(channel-30,block)
        cursor%position = indexed_input%sparse_positions(block) + 4*(indexed_input%column_counts(block) + &
            indexed_input%row_counts(block)) + contribution_bytes*sum(indexed_input%channel_counts(:channel-31,block))
    end subroutine start_channel

    logical function next_group(cursor) result(found)
        type(mcp_cursor), intent(inout) :: cursor
        found = cursor%groups_left > 0
        if (.not. found) then
            if (cursor%remaining /= 0) &
                call fail_input_io_mpi(0,cursor%unit,io_context//" contribution total disagrees with channel directory")
            return
        end if
        if (cursor%batch_index > cursor%batch_count) then
            cursor%batch_count = int(min(int(directory_batch_size,int64),cursor%groups_left))
            call read_indexed_directory(indexed_input,cursor%index,cursor%directory_batch(:cursor%batch_count))
            cursor%batch_index = 1
            cursor%index = cursor%index + cursor%batch_count
            cursor%bytes_read = cursor%bytes_read + mcp_group_directory_bytes*cursor%batch_count
        end if
        cursor%group = cursor%directory_batch(cursor%batch_index)
        cursor%batch_index = cursor%batch_index + 1
        call validate_indexed_group(indexed_input,cursor%position,cursor%remaining,cursor%previous_label, &
            cursor%group,cursor%channel==31)
        cursor%groups_left = cursor%groups_left - 1
        cursor%position = cursor%position + cursor%group%bytes
        cursor%remaining = cursor%remaining - cursor%group%count
        cursor%previous_label = cursor%group%label
        cursor%max_group = max(cursor%max_group,cursor%group%count)
    end function next_group

    subroutine skip_group(cursor)
        type(mcp_cursor), intent(inout) :: cursor
        cursor%groups_skipped = cursor%groups_skipped + 1
        cursor%skipped_contributions = cursor%skipped_contributions + cursor%group%count
        cursor%bytes_skipped = cursor%bytes_skipped + cursor%group%bytes
    end subroutine skip_group

    subroutine ensure_buffer(buffer, length, status)
        type(mcp_buffer), intent(inout) :: buffer
        integer, intent(in) :: length
        integer, intent(out) :: status
        status = 0
        if (allocated(buffer%values)) then
            if (length <= size(buffer%values)) return
            deallocate(buffer%columns,buffer%positions,buffer%values)
        end if
        allocate(buffer%columns(length),buffer%positions(length),buffer%values(length),stat=status)
    end subroutine ensure_buffer

    subroutine read_segment(cursor, buffer, first, length, row_count)
        type(mcp_cursor), intent(inout) :: cursor
        type(mcp_buffer), intent(inout) :: buffer
        integer(int64), intent(in) :: first, row_count
        integer, intent(in) :: length
        integer :: status, n
        call ensure_buffer(buffer,length,status)
        if (status /= 0) call fail_input_io_mpi(status,cursor%unit,io_context//" allocate contribution buffer")
        call read_indexed_contributions(indexed_input,cursor%group,first,buffer%columns(:length), &
            buffer%positions(:length),buffer%values(:length))
        n = csf_counts_by_block(cursor%block)
        if (any(buffer%columns(:length) <= mpi_rank) .or. any(buffer%columns(:length) > n) .or. &
            any(buffer%positions(:length) < 1) .or. any(buffer%positions(:length) > row_count)) &
            call fail_input_io_mpi(0,cursor%unit,io_context//" contribution index outside rank block")
        if (any(mod(buffer%columns(:length)-1,mpi_process_count) /= mpi_rank)) &
            call fail_input_io_mpi(0,cursor%unit,io_context//" contribution column owned by another rank")
        if (first == 1) cursor%groups_read = cursor%groups_read + 1
        cursor%contributions = cursor%contributions + length
        cursor%bytes_read = cursor%bytes_read + contribution_bytes*length
    end subroutine read_segment
end module rmcdhf_mcp_io
