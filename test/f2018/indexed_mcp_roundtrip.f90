! Exercise bounded group reads, exact bits, empty ranks and 64-bit positions without a large matrix.
program indexed_mcp_roundtrip
    use iso_fortran_env, only: int64, real64
    use grasp_parameters, only: orbital_index_radix
    use mpi_state, only: mpi_init, mpi_finalize, mpi_rank, mpi_process_count, mpi_comm_rank, mpi_comm_size, mpi_comm_world
    use mpi_indexed_mcp, only: indexed_mcp, mcp_group, mcp_segment_size, open_indexed_mcp_writer, finish_indexed_mcp, &
        open_indexed_mcp_reader, write_indexed_sparse, start_indexed_channel, write_indexed_group, close_indexed_mcp, &
        read_indexed_sparse, read_indexed_group, read_indexed_contributions, contribution_bytes
    use write_mcp_sparsity_summary_mpi_interface, only: write_mcp_sparsity_summary_mpi
    implicit none (type, external)
    integer, parameter :: label = orbital_index_radix + 1
    integer, parameter :: entries = mcp_segment_size + 1
    integer :: status, i, columns, unused_columns(1), values_positions(1), segment_columns(mcp_segment_size), &
        segment_positions(mcp_segment_size)
    integer, allocatable :: indices(:), ends(:), rows(:)
    real(real64), allocatable :: weights(:)
    real(real64) :: last_value(1), segment_values(mcp_segment_size)
    integer(int64) :: position
    character(len=32) :: scenario, prefix
    type(indexed_mcp) :: file
    type(mcp_group) :: group

    call mpi_init(status)
    call mpi_comm_rank(mpi_comm_world,mpi_rank,status)
    call mpi_comm_size(mpi_comm_world,mpi_process_count,status)
    call get_command_argument(1,scenario)
    write(prefix,'(A,I0)') "probe-",mpi_rank
    select case (trim(scenario))
    case ("roundtrip")
        call open_indexed_mcp_writer(file,prefix,0,0,2,2,[1],["   0+   "],[0],0,repeat("a",64))
        columns = merge(1,0,mpi_rank==0)
        allocate(ends(columns),rows(columns))
        if (columns > 0) then
            ends = 1
            rows = 1
        end if
        call write_indexed_sparse(file,1,ends,rows)
        call start_indexed_channel(file,1)
        if (mpi_rank == 0) then
            allocate(indices(entries),weights(entries))
            indices = 1
            do i = 1, entries
                weights(i) = transfer(int(z'3ff0000000000000',int64)+i,0.0_real64)
            end do
            call write_indexed_group(file,label,indices,indices,weights)
            call write_indexed_group(file,label+1,indices(:1),indices(:1),weights(:1))
        end if
        call start_indexed_channel(file,2)
        call finish_indexed_mcp(file)

        call open_indexed_mcp_reader(file,prefix,1)
        call read_indexed_sparse(file,1,ends,rows)
        if (any(ends /= 1) .or. any(rows /= 1)) error stop "sparse round-trip"
        if (any(file%group_counts(2,:) /= 0)) error stop "empty channel"
        if (mpi_rank == 0) then
            if (file%group_counts(1,1) /= 2) error stop "channel group count"
            position = file%sparse_positions(1) + 8
            call read_indexed_group(file,1_int64,position,int(entries+1,int64),0_int64,group,.true.)
            if (group%count /= entries) error stop "group count"
            call read_indexed_contributions(file,group,1_int64,segment_columns,segment_positions,segment_values)
            if (any(segment_columns /= 1) .or. any(segment_positions /= 1)) error stop "group indices"
            if (any(transfer(segment_values,[0_int64],mcp_segment_size) /= &
                transfer(weights(:mcp_segment_size),[0_int64],mcp_segment_size))) error stop "group bits"
            call read_indexed_contributions(file,group,int(entries,int64),unused_columns,values_positions,last_value)
            if (transfer(last_value(1),0_int64) /= transfer(weights(entries),0_int64)) error stop "trailing segment bits"
        else
            if (file%group_count /= 0 .or. file%row_counts(1) /= 0) error stop "empty rank has data"
        end if
        call close_indexed_mcp(file)
    case ("large-offset")
        call open_indexed_mcp_reader(file,prefix,1)
        if (mpi_rank == 0) then
            position = file%sparse_positions(1) + 8 + contribution_bytes*(file%channel_counts(1,1)-1)
            if (position <= huge(0)) error stop "large-offset fixture below limit"
            call read_indexed_group(file,2_int64,position,1_int64,int(label,int64),group,.true.)
            call read_indexed_contributions(file,group,1_int64,unused_columns,values_positions,last_value)
            if (transfer(last_value(1),0_int64) /= int(z'3ff0000000000001',int64)) error stop "64-bit offset read"
        end if
        call close_indexed_mcp(file)
    case ("sparsity-count")
        call write_mcp_sparsity_summary_mpi(.false.,1500000000,100000)
    case default
        error stop "unknown round-trip scenario"
    end select
    call mpi_finalize(status)
end program indexed_mcp_roundtrip
