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
