! Real SETMCP and disk SETCOF entry points with controlled angular weights.
program mcp_reader_io
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters, only: orbital_index_radix
    use orbital_csf_basis, only: orbital_kappa, orbital_two_j
    use hamiltonian_symmetry_blocks, only: symmetry_block_count, csf_counts_by_block
    use symmetry_block_offsets, only: csf_offsets_by_block
    use scf_iteration_state, only: orbital_occupation_weights, direct_coefficient_capacity, exchange_coefficient_capacity, &
        y_potential_coefficient_capacity, direct_coefficient_count, direct_coefficients
    use connect_mcp_files_interface, only: connect_mcp_files
#ifdef MCP_MEM
    use rmcdhf_mcp_memory, only: lab_31, ncontr_31, coeff_31
#else
    use assemble_orbital_potential_coefficients_interface, only: assemble_orbital_potential_coefficients
#endif
#ifdef MCP_MPI
    use mpi_state, only: mpi_barrier, mpi_comm_rank, mpi_comm_size, mpi_comm_world, mpi_finalize, mpi_init, &
        mpi_process_count, mpi_rank
#else
    use serial_mpi_context, only: mpi_process_count, mpi_rank
#endif
    use mcp_file_metadata, only: mcp_input_units
    implicit none (type, external)
    character(40) :: scenario, prefix
    character(80) :: filename
    character(8) :: labels(1)
    integer :: fixture_output_unit, channel_suffix, status, ncore, rank_in_file, packed_label, block_id, coefficient_count
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
    do channel_suffix=30,32
        write(filename, "(A,A,I0)") trim(prefix), ".", channel_suffix
        if (corrupt .and. scenario == "missing" .and. channel_suffix == 30) cycle
        open(newunit=fixture_output_unit, file=trim(filename), form="unformatted", status="replace")
        omit = .false.
        if (channel_suffix == 30) then
            omit = corrupt .and. scenario == "global-header"
            if (.not.omit) then
                if (corrupt .and. scenario=="short-global-header") then
                    write(fixture_output_unit) 0
                else
                    write(fixture_output_unit) 0, 1, 0
                end if
            end if
            omit = omit .or. (corrupt .and. scenario == "block-sizes")
            if (.not.omit) write(fixture_output_unit) mpi_process_count
            omit = omit .or. (corrupt .and. scenario == "block-labels")
            if (.not.omit) write(fixture_output_unit) "  0+    "
        end if
        omit = omit .or. (corrupt .and. channel_suffix==31 .and. scenario=="channel-header")
        if (.not.omit) then
            rank_in_file = mpi_rank
            if (corrupt .and. channel_suffix==31 .and. scenario=="rank") rank_in_file=mpi_rank+1
            if (corrupt .and. channel_suffix==31 .and. scenario=="magic") then
                write(fixture_output_unit) "BAD", 1, rank_in_file, mpi_process_count
            else
                write(fixture_output_unit) "MCP", 1, rank_in_file, mpi_process_count
            end if
        end if
        omit = omit .or. (corrupt .and. channel_suffix==31 .and. scenario=="basis")
        if (.not.omit) write(fixture_output_unit) 2, mpi_process_count, 1
        omit = omit .or. (corrupt .and. channel_suffix==31 .and. scenario=="flags")
        if (.not.omit) write(fixture_output_unit) .false., 1, .false.
        if (.not.omit) then
            block_id = 1
            if (corrupt .and. ((channel_suffix==30 .and. scenario=="sparse-block-id") .or. &
                              (channel_suffix==31 .and. scenario=="block-id"))) block_id = 2
            if (channel_suffix==30) then
                write(fixture_output_unit) "MCP", block_id, mpi_process_count
                write(fixture_output_unit) 1
                write(fixture_output_unit) 1, mpi_rank+1
            else
                coefficient_count = 0
                if (channel_suffix==31 .and. scenario/="empty") coefficient_count = 1
                write(fixture_output_unit) "MCP", block_id, mpi_process_count, coefficient_count
                if (channel_suffix==31 .and. scenario/="empty") then
                    omit = corrupt .and. scenario=="label"
                    if (.not.omit) then
                        if (corrupt .and. scenario=="short-label") then
                            write(fixture_output_unit) packed_label
                        else
                            write(fixture_output_unit) packed_label, 1
                        end if
                    end if
                    omit = omit .or. (corrupt .and. scenario=="payload")
                    if (.not.omit) then
                        if (corrupt .and. scenario=="short-payload") then
                            write(fixture_output_unit) mpi_rank+1
                        else
                            write(fixture_output_unit) mpi_rank+1, 1, 0.25_real64
                        end if
                    end if
                    omit = omit .or. (corrupt .and. scenario=="sentinel")
                end if
                if (.not.omit) write(fixture_output_unit) 0, 0
            end if
        end if
        close(fixture_output_unit)
    end do
    call connect_mcp_files(ncore, 1, labels, trim(prefix))
#ifndef MCP_MEM
    orbital_kappa(1) = -1
    orbital_two_j(1) = 1
    orbital_occupation_weights(1) = 1.0_real64
    allocate(csf_offsets_by_block(1))
    csf_offsets_by_block = 0
    direct_coefficient_capacity = 0
    exchange_coefficient_capacity = 0
    y_potential_coefficient_capacity = 0
    call assemble_orbital_potential_coefficients(.true., 1)
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
        if (direct_coefficient_count/=0) error stop "Empty contraction"
    else
        if (direct_coefficient_count/=1 .or. direct_coefficients(1)/=0.125_real64) error stop "Contracted coefficient"
    end if
#endif
    do channel_suffix=30,32
        close(mcp_input_units(channel_suffix), status="delete")
    end do
#ifdef MCP_MPI
    call MPI_Barrier(MPI_COMM_WORLD, status)
    call MPI_Finalize(status)
#endif
    print *, "PASS MCP reader input"
end program mcp_reader_io
