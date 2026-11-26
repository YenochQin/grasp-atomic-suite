! Exercise the production sparse writer on an occupied and an empty MPI rank.
program sparse_empty_rank
    use iso_fortran_env, only: int64
    use orbital_csf_basis, only: csf_count
    use mpi_state, only: mpi_rank, mpi_process_count, mpi_comm_world, mpi_init, mpi_comm_rank, mpi_comm_size, mpi_finalize
    use mpi_indexed_mcp, only: indexed_mcp, open_indexed_mcp_writer, start_indexed_channel, finish_indexed_mcp, &
        open_indexed_mcp_reader, read_indexed_sparse
    use write_mcp_sparsity_interface, only: write_mcp_sparsity
    use mcp_generation_options, only: mcp_temporary_units
    implicit none (type, external)
    type(indexed_mcp) :: shard
    integer :: status, local_size
    integer, allocatable :: ends(:), rows(:)
    character(16) :: prefix
    logical :: printing

    call mpi_init(status)
    call mpi_comm_rank(mpi_comm_world,mpi_rank,status)
    call mpi_comm_size(mpi_comm_world,mpi_process_count,status)
    write(prefix,'(A,I0)') "sparse-rank-",mpi_rank
    allocate(mcp_temporary_units(30:30))
    csf_count = 1
    local_size = merge(1, 0, mpi_rank == 0)
    printing = .false.
    open(newunit=mcp_temporary_units(30),status="scratch",form="unformatted")
    write(mcp_temporary_units(30)) "MCP",1
    if (local_size == 1) write(mcp_temporary_units(30)) 1,1,1
    rewind(mcp_temporary_units(30))
    call open_indexed_mcp_writer(shard,prefix,0,0,2,1,[1],["   0+   "],[0],0,repeat("a",64))
    call write_mcp_sparsity(check_callback,local_size,printing,1,mpi_rank,mpi_process_count,shard)
    close(mcp_temporary_units(30))
    call start_indexed_channel(shard,1)
    call start_indexed_channel(shard,2)
    call finish_indexed_mcp(shard)
    call open_indexed_mcp_reader(shard,prefix,1)
    if (shard%row_counts(1) /= int(local_size,int64)) error stop "Sparse count"
    allocate(ends(local_size),rows(local_size))
    call read_indexed_sparse(shard,1,ends,rows)
    if (any(ends /= 1) .or. any(rows /= 1)) error stop "Sparse owned column"
    close(shard%unit,status="delete")
    call mpi_finalize(status)
    print *, "PASS indexed sparse occupied and empty ranks"
contains
    subroutine check_callback(lprint,nnonz,n)
        logical, intent(in) :: lprint
        integer, intent(in) :: nnonz,n
        if (lprint .or. nnonz /= local_size .or. n /= csf_count) error stop "Sparse callback"
    end subroutine check_callback
end program sparse_empty_rank
