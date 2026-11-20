! Real SETMCP and disk SETCOF entry points with controlled angular weights.
program mcp_reader_io
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters, only: orbital_index_radix
    use orbital_csf_basis, only: orbital_kappa, orbital_two_j
    use hamiltonian_symmetry_blocks, only: symmetry_block_count, csf_counts_by_block
    use symmetry_block_offsets, only: csf_offsets_by_block
    use scf_iteration_state, only: ucf, nddim, nxdim, nydim, ndcof, da
    use setmcp_I, only: setmcp
#ifdef MCP_MEM
    use rmcdhf_mem_c, only: lab_31, ncontr_31, coeff_31
#else
    use setcof_I, only: setcof
#endif
#ifdef MCP_MPI
    use mpi_state, only: mpi_barrier, mpi_comm_rank, mpi_comm_size, mpi_comm_world, mpi_finalize, mpi_init, &
       & mpi_process_count, mpi_rank
#else
    use mpi_s, only: mpi_process_count, mpi_rank
#endif
    implicit none (type, external)
    character(40) :: scenario, prefix
    character(80) :: filename
    character(8) :: labels(1)
    integer :: unit, status, ncore, rank_in_file, packed_label, block_id, coefficient_count
    logical :: corrupt, omit
#ifdef MCP_MPI
    call MPI_Init(status)
    call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, status)
    call MPI_Comm_size(MPI_COMM_WORLD, mpi_process_count, status)
#else
    mpi_rank = 0
    mpi_process_count = 1
#endif
    call get_command_argument(1, scenario)
    corrupt = scenario /= "valid" .and. scenario /= "empty"
#ifdef MCP_MPI
    corrupt = corrupt .and. mpi_rank == mpi_process_count-1
#endif
    write(prefix, "(A,I0)") "mcp-rank-", mpi_rank
    packed_label = orbital_index_radix+1
    do unit=30,32
        write(filename, "(A,A,I0)") trim(prefix), ".", unit
        if (corrupt .and. scenario == "missing" .and. unit == 30) cycle
        open(unit, file=trim(filename), form="unformatted", status="replace")
        omit = .false.
        if (unit == 30) then
            omit = corrupt .and. scenario == "global-header"
            if (.not.omit) then
                if (corrupt .and. scenario=="short-global-header") then
                    write(unit) 0
                else
                    write(unit) 0, 1, 0
                end if
            end if
            omit = omit .or. (corrupt .and. scenario == "block-sizes")
            if (.not.omit) write(unit) mpi_process_count
            omit = omit .or. (corrupt .and. scenario == "block-labels")
            if (.not.omit) write(unit) "  0+    "
        end if
        omit = omit .or. (corrupt .and. unit==31 .and. scenario=="channel-header")
        if (.not.omit) then
            rank_in_file = mpi_rank
            if (corrupt .and. unit==31 .and. scenario=="rank") rank_in_file=mpi_rank+1
            if (corrupt .and. unit==31 .and. scenario=="magic") then
                write(unit) "BAD", 1, rank_in_file, mpi_process_count
            else
                write(unit) "MCP", 1, rank_in_file, mpi_process_count
            end if
        end if
        omit = omit .or. (corrupt .and. unit==31 .and. scenario=="basis")
        if (.not.omit) write(unit) 2, mpi_process_count, 1
        omit = omit .or. (corrupt .and. unit==31 .and. scenario=="flags")
        if (.not.omit) write(unit) .false., 1, .false.
        if (.not.omit) then
            block_id = 1
            if (corrupt .and. ((unit==30 .and. scenario=="sparse-block-id") .or. &
                              (unit==31 .and. scenario=="block-id"))) block_id = 2
            if (unit==30) then
                write(unit) "MCP", block_id, mpi_process_count
                write(unit) 1
                write(unit) 1, mpi_rank+1
            else
                coefficient_count = 0
                if (unit==31 .and. scenario/="empty") coefficient_count = 1
                write(unit) "MCP", block_id, mpi_process_count, coefficient_count
                if (unit==31 .and. scenario/="empty") then
                    omit = corrupt .and. scenario=="label"
                    if (.not.omit) then
                        if (corrupt .and. scenario=="short-label") then
                            write(unit) packed_label
                        else
                            write(unit) packed_label, 1
                        end if
                    end if
                    omit = omit .or. (corrupt .and. scenario=="payload")
                    if (.not.omit) then
                        if (corrupt .and. scenario=="short-payload") then
                            write(unit) mpi_rank+1
                        else
                            write(unit) mpi_rank+1, 1, 0.25_real64
                        end if
                    end if
                    omit = omit .or. (corrupt .and. scenario=="sentinel")
                end if
                if (.not.omit) write(unit) 0, 0
            end if
        end if
        close(unit)
    end do
    call setmcp(ncore, 1, labels, trim(prefix))
#ifndef MCP_MEM
    orbital_kappa(1) = -1
    orbital_two_j(1) = 1
    ucf(1) = 1.0_real64
    allocate(csf_offsets_by_block(1))
    csf_offsets_by_block = 0
    nddim = 0
    nxdim = 0
    nydim = 0
    call setcof(.true., 1)
#endif
    if (corrupt) error stop "Corrupt MCP was accepted"
    if (ncore /= 0 .or. symmetry_block_count /= 1) error stop "Global MCP header"
    if (csf_counts_by_block(1) /= mpi_process_count) error stop "Block sizes"
    if (labels(1) /= "  0+    ") error stop "Block labels"
#ifdef MCP_MEM
    if (scenario=="empty") then
        if (size(coeff_31)/=0 .or. any(lab_31/=0)) error stop "Empty cached coefficients"
    else
        if (any(lab_31/=[packed_label,0]) .or. any(ncontr_31/=[1,0])) error stop "Cached labels"
        if (size(coeff_31)/=1 .or. coeff_31(1)/=0.25_real64) error stop "Cached coefficients"
    end if
#else
    if (scenario=="empty") then
        if (ndcof/=0) error stop "Empty contraction"
    else
        if (ndcof/=1 .or. da(1)/=0.125_real64) error stop "Contracted coefficient"
    end if
#endif
    do unit=30,32
        close(unit, status="delete")
    end do
#ifdef MCP_MPI
    call MPI_Barrier(MPI_COMM_WORLD, status)
    call MPI_Finalize(status)
#endif
    print *, "PASS MCP reader input"
end program mcp_reader_io

function dsubrs(eol, i, j, jblock) result(value)
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)
    logical, intent(in) :: eol
    integer, intent(in) :: i, j
    integer :: jblock
    real(real64) :: value
    value = 1.0_real64
end function dsubrs

function fco(k, ir, ia, ib) result(value)
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)
    integer, intent(in) :: k, ia, ib
    integer :: ir
    real(real64) :: value
    value = 0.0_real64
end function fco

function gco(k, ir, ia, ib) result(value)
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)
    integer :: k, ir
    integer, intent(in) :: ia, ib
    real(real64) :: value
    value = 0.0_real64
end function gco

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
