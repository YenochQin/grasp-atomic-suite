! Fault injection between the production loader's counting and loading passes.
module memory_io_injection
    use mcp_file_metadata, only: mcp_input_units
    implicit none (type, external)
    logical :: inject = .false.
    character(40) :: injected_case
    character(80) :: injected_file
    integer :: injected_rank, injected_ranks
contains
    subroutine rewrite_input()
        use, intrinsic :: iso_fortran_env, only: real64
        logical :: omit
        if (.not.inject) return
        inject = .false.
        close(mcp_input_units(32))
        open(newunit=mcp_input_units(32), file=trim(injected_file), status="replace", form="unformatted")
        write(mcp_input_units(32)) "MCP", 1, injected_rank, injected_ranks
        write(mcp_input_units(32)) 2, 1, 1
        write(mcp_input_units(32)) .false., 1, .false.
        omit = injected_case == "load-block-header"
        if (.not.omit) write(mcp_input_units(32)) "MCP", 1, 1, 1
        omit = omit .or. injected_case == "load-label"
        if (.not.omit) then
            if (injected_case == "load-short-label") then
                write(mcp_input_units(32)) 42
            else
                write(mcp_input_units(32)) 42, 1
            end if
        end if
        omit = omit .or. injected_case == "load-payload"
        if (.not.omit) then
            if (injected_case == "load-short-payload") then
                write(mcp_input_units(32)) 1
            else
                write(mcp_input_units(32)) 1, 1, 0.25_real64
            end if
        end if
        omit = omit .or. injected_case == "load-sentinel"
        if (.not.omit) write(mcp_input_units(32)) 0, 0
        flush(mcp_input_units(32))
    end subroutine rewrite_input
end module memory_io_injection

! Test-only allocation seam: preserve the pointer/shape contract and change
! the on-disk input after the final channel-32 allocation, before REWIND.
! The header-reader tests separately exercise the production allocator.
module array_memory_management
    use memory_io_injection, only: rewrite_input
    implicit none (type, external)
    interface alloc
        module procedure fixture_alloc_i, fixture_alloc_il, fixture_alloc_r, fixture_alloc_rl
    end interface
contains
    subroutine fixture_alloc_i(p, n, var, caller)
        integer, pointer :: p(:)
        integer, intent(in) :: n
        character(*), intent(in) :: var, caller
        allocate(p(n))
        p = 0
    end subroutine fixture_alloc_i
    subroutine fixture_alloc_il(p, n, var, caller)
        use, intrinsic :: iso_fortran_env, only: int64
        integer, pointer :: p(:)
        integer(int64), intent(in) :: n
        character(*), intent(in) :: var, caller
        allocate(p(n))
        p = 0
    end subroutine fixture_alloc_il
    subroutine fixture_alloc_r(p, n, var, caller)
        use, intrinsic :: iso_fortran_env, only: real64
        real(real64), pointer :: p(:)
        integer, intent(in) :: n
        character(*), intent(in) :: var, caller
        allocate(p(n))
        p = 0
        if (var == "COEFF_32") call rewrite_input()
    end subroutine fixture_alloc_r
    subroutine fixture_alloc_rl(p, n, var, caller)
        use, intrinsic :: iso_fortran_env, only: real64, int64
        real(real64), pointer :: p(:)
        integer(int64), intent(in) :: n
        character(*), intent(in) :: var, caller
        allocate(p(n))
        p = 0
        if (var == "COEFF_32") call rewrite_input()
    end subroutine fixture_alloc_rl
end module array_memory_management

! Exercise the production memory loader, including its label-count sentinels.
program mcp_memory_io
    use, intrinsic :: iso_fortran_env, only: real64
    use rmcdhf_mem_c, only: lab_32, ncontr_32, coeff_32, indx_32, iclmn_32
    use setmcp_mem_I, only: setmcp_mem
    use memory_io_injection, only: mcp_input_units, inject, injected_case, injected_file, injected_rank, injected_ranks, &
       & rewrite_input
#ifdef MCP_MPI
    use mpi_state, only: mpi_barrier, mpi_comm_rank, mpi_comm_size, mpi_comm_world, mpi_finalize, mpi_init, &
       & mpi_process_count, mpi_rank
#else
    use mpi_s, only: mpi_process_count, mpi_rank
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
    call setmcp_mem(1, 32)
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

! Keep deterministic memory-reporting seams; no host-memory assumption is needed.
function systemmem() result(value)
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)
    real(real64) :: value
    value = 1.0_real64
end function systemmem

function systemfreemem() result(value)
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)
    real(real64) :: value
    value = 1.0_real64
end function systemfreemem
