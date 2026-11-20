! Exercise the production memory loader, including its label-count sentinels.
program mcp_memory_io
    use, intrinsic :: iso_fortran_env, only: real64
    use rmcdhf_mcp_memory, only: lab_32, ncontr_32, coeff_32, indx_32, iclmn_32
    use cache_mcp_files_interface, only: cache_mcp_files
    use memory_io_injection, only: mcp_input_units, inject, injected_case, injected_file, injected_rank, injected_ranks, &
       & rewrite_input
#ifdef MCP_MPI
    use mpi_state, only: mpi_barrier, mpi_comm_rank, mpi_comm_size, mpi_comm_world, mpi_finalize, mpi_init, &
        mpi_process_count, mpi_rank
#else
    use serial_mpi_context, only: mpi_process_count, mpi_rank
#endif
    implicit none (type, external)
    character(40) :: scenario
    character(80) :: filename
    integer :: status
    logical :: corrupt
#ifdef MCP_MPI
    call MPI_Init(status)
    call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, status)
    call MPI_Comm_size(MPI_COMM_WORLD, mpi_process_count, status)
#else
    mpi_rank = 0
    mpi_process_count = 1
#endif
    call get_command_argument(1, scenario)
    corrupt = trim(scenario) /= "valid" .and. trim(scenario) /= "empty"
#ifdef MCP_MPI
    ! Only the last rank is damaged; other ranks must not hang in a collective.
    corrupt = corrupt .and. mpi_rank == mpi_process_count-1
#endif
    write(filename, "(A,I0,A)") "mcp-rank-", mpi_rank, ".32"
    inject = corrupt .and. index(scenario, "load-") == 1
    injected_case = scenario
    injected_file = filename
    injected_rank = mpi_rank
    injected_ranks = mpi_process_count
    allocate(mcp_input_units(32:32))
    open(newunit=mcp_input_units(32), file=trim(filename), form="unformatted", status="replace")
    write(mcp_input_units(32)) "MCP", 1, mpi_rank, mpi_process_count
    write(mcp_input_units(32)) 2, 1, 1
    write(mcp_input_units(32)) .false., 1, .false.
    if (.not. (corrupt .and. scenario == "block-header")) then
        if (scenario == "empty") then
            write(mcp_input_units(32)) "MCP", 1, 1, 0
        else
            write(mcp_input_units(32)) "MCP", 1, 1, 1
        end if
    end if
    if (.not. (corrupt .and. (scenario == "block-header" .or. scenario == "label"))) then
        if (scenario == "empty") then
            write(mcp_input_units(32)) 0, 0
        else
            if (corrupt .and. scenario == "short-label") then
                write(mcp_input_units(32)) 42
            else
                write(mcp_input_units(32)) 42, 1
            end if
            if (.not. (corrupt .and. scenario == "payload")) then
                if (corrupt .and. scenario == "short-payload") then
                    write(mcp_input_units(32)) 1
                else
                    write(mcp_input_units(32)) 1, 1, 0.25_real64
                end if
                if (.not. (corrupt .and. scenario == "sentinel")) write(mcp_input_units(32)) 0, 0
            end if
        end if
    end if
    rewind(mcp_input_units(32))
    read(mcp_input_units(32))
    read(mcp_input_units(32))
    read(mcp_input_units(32))
    call cache_mcp_files(1, 32)
    if (corrupt) error stop "Corrupt MCP was accepted"
    if (scenario == "empty") then
        if (size(lab_32) /= 1 .or. any(lab_32 /= 0)) error stop "Empty block labels"
        if (size(coeff_32) /= 0) error stop "Empty block coefficients"
    else
        if (size(lab_32) /= 2 .or. any(lab_32 /= [42,0])) error stop "MCP labels"
        if (any(ncontr_32 /= [1,0])) error stop "MCP contribution counts"
        if (size(coeff_32) /= 1 .or. coeff_32(1) /= 0.25_real64) error stop "MCP coefficient"
        if (indx_32(1) /= 1 .or. iclmn_32(1) /= 1) error stop "MCP indices"
    end if
    close(mcp_input_units(32), status="delete")
#ifdef MCP_MPI
    call MPI_Barrier(MPI_COMM_WORLD, status)
    call MPI_Finalize(status)
#endif
    print *, "PASS MCP memory input"
end program mcp_memory_io
