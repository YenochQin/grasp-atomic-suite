program ci_damping
    use ci_fixture, only: shifted_eigenvalues, csf_mixing_coefficients, solves
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use scf_damping_factors, only: mixing_damping_factors
    use atomic_calculation_state, only: requested_level_indices, selected_level_count, highest_requested_level_index
    use debug_output_settings, only: ldbpg
    use hamiltonian_symmetry_blocks, only: symmetry_block_count, csf_counts_by_block, eigenvector_counts_by_block, &
        highest_selected_level_by_block
    use sparse_hamiltonian_storage, only: hamiltonian_values, matrix_partner_indices, stored_segment_end_offsets
    use mcp_file_metadata, only: max_mcp_file_multipole_rank, mcp_input_units
    use calculation_file_units, only: mixing_output_unit
    use orbital_csf_basis, only: csf_count, orbital_count
    use symmetry_block_offsets, only: total_csf_count, mixing_coefficient_storage_size, csf_offsets_by_block, &
        eigenvalue_offsets_by_block, mixing_offsets_by_block
    use block_mean_energies, only: mean_energy_by_block
    use level_symmetry_metadata, only: packed_block_symmetries
#ifdef CI_MPI
    use mpi_state, only: processor_name, mpi_error_code, mpi_status_buffer, processor_name_length, mpi_rank, &
        mpi_process_count, mpi_finalize, initialize_mpi_runtime
    use solve_hamiltonian_blocks_mpi_interface, only: solve_hamiltonian_blocks_mpi
#else
    use serial_mpi_context, only: mpi_rank, mpi_process_count
    use solve_hamiltonian_blocks_interface, only: solve_hamiltonian_blocks
#endif
#ifdef CI_MEM
    use rmcdhf_mcp_memory, only: imin1_30, imax1_30, imin2_30, imax2_30
#endif
#ifdef CI_INDEXED
    use mpi_mcp_calculation_input, only: indexed_input
#endif
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)
    integer :: channel_suffix, i, records, local_size
    real(real64), parameter :: tol = 1e-12_real64
#ifdef CI_MPI
    call initialize_mpi_runtime(mpi_rank,mpi_process_count,processor_name,processor_name_length)
#else
    mpi_rank = 0
    mpi_process_count = 1
#endif
    symmetry_block_count=1; csf_count=2; total_csf_count=2; orbital_count=1; selected_level_count=2; &
        highest_requested_level_index=2; &
        mixing_coefficient_storage_size=8; max_mcp_file_multipole_rank=1
    ldbpg=.false.; packed_block_symmetries=1
    allocate(csf_counts_by_block(1), eigenvector_counts_by_block(1), highest_selected_level_by_block(1), &
        csf_offsets_by_block(1), eigenvalue_offsets_by_block(1), mixing_offsets_by_block(1), mean_energy_by_block(1))
    allocate(requested_level_indices(2), mixing_damping_factors(2), csf_mixing_coefficients(8), shifted_eigenvalues(2))
    csf_counts_by_block=2; eigenvector_counts_by_block=2; highest_selected_level_by_block=2; csf_offsets_by_block=0; &
        eigenvalue_offsets_by_block=0; mixing_offsets_by_block=2
    requested_level_indices=[1,2]; mixing_damping_factors=[0.0_real64,0.5_real64]; csf_mixing_coefficients=99.0_real64
    nullify(hamiltonian_values,matrix_partner_indices,stored_segment_end_offsets)
    local_size = count([(mod(i-1,mpi_process_count)==mpi_rank, i=1,2)])
#ifdef CI_MEM
    allocate(imin1_30(1),imax1_30(1),imin2_30(1),imax2_30(1))
    imin1_30=0; imax1_30=2; imin2_30=0; imax2_30=local_size
#elif defined(CI_INDEXED)
    allocate(indexed_input%sparse_positions(1),indexed_input%row_counts(1),indexed_input%column_counts(1))
    indexed_input%sparse_positions = 1
    indexed_input%row_counts = local_size
    indexed_input%column_counts = local_size
    open(newunit=indexed_input%unit,status="scratch",access="stream",form="unformatted")
    write(indexed_input%unit) [(i,i=1,local_size)],[(mpi_rank+1+(i-1)*mpi_process_count,i=1,local_size)]
#else
    allocate(mcp_input_units(30:33))
    do channel_suffix=30,33
        open(newunit=mcp_input_units(channel_suffix),status="scratch",form="unformatted")
        records=3
        if(channel_suffix==30) records=6
        do i=1,records
            write(mcp_input_units(channel_suffix)) 0
        end do
        if(channel_suffix==30) then
            write(mcp_input_units(channel_suffix)) "MCP",1,2
            write(mcp_input_units(channel_suffix)) local_size
            write(mcp_input_units(channel_suffix)) [(i,i=1,local_size)], [(mpi_rank+1+(i-1)*mpi_process_count,i=1,local_size)]
        end if
    end do
#endif
    if(mpi_rank==0) then
        open(newunit=mixing_output_unit,status="scratch",form="unformatted")
        write(mixing_output_unit) "G92MIX"
        write(mixing_output_unit) 2,2,1,2,8,1
    end if
    call matrix_step(.true.)
    ! With fcheck=all, the production allocation/zeroing must include IENDC(0).
    if (maxval(abs(csf_mixing_coefficients(3:6)-[1.0_real64,0.0_real64,0.0_real64,1.0_real64]))>tol) &
        error stop "FIRST must bypass damping"
    call matrix_step(.false.)
    if(.not.all(ieee_is_finite(csf_mixing_coefficients))) error stop "nonfinite damped roots"
    if (maxval(abs(csf_mixing_coefficients(3:6)-[0.6_real64,0.8_real64,0.8_real64,-0.6_real64]))>tol) &
        error stop "DDOT damping branch differs from analytic roots"
    if(abs(dot_product(csf_mixing_coefficients(3:4),csf_mixing_coefficients(5:6)))>tol) error stop "orthogonality"
    if(abs(sum(csf_mixing_coefficients(3:4)**2)-1.0_real64)>tol .or. &
        abs(sum(csf_mixing_coefficients(5:6)**2)-1.0_real64)>tol) &
        error stop "normalization"
    if (any(csf_mixing_coefficients(1:2) /= 99.0_real64) .or. &
        any(csf_mixing_coefficients(7:8) /= 99.0_real64)) error stop "block offset guards"
    mixing_damping_factors=0.0_real64
    call matrix_step(.false.)
    if (maxval(abs(csf_mixing_coefficients(3:6)-[0.6_real64,0.8_real64,-0.8_real64,0.6_real64]))>tol) &
        error stop "zero damping control"
    if(solves/=3) error stop "fixture did not execute three production calls"
    if(mpi_rank==0) then
        close(mixing_output_unit)
        print *, "PASS: production FIRST / nonzero multi-root DDOT / zero damping"
    end if
#if defined(CI_INDEXED)
    close(indexed_input%unit)
#elif !defined(CI_MEM)
    do channel_suffix=30,33
        close(mcp_input_units(channel_suffix))
    end do
#endif
#ifdef CI_MPI
    call MPI_Finalize(mpi_error_code)
#endif
contains
    subroutine matrix_step(first)
        logical,intent(in)::first
#ifdef CI_MPI
        call solve_hamiltonian_blocks_mpi(first)
#else
        call solve_hamiltonian_blocks(first)
#endif
end subroutine matrix_step
end program ci_damping
