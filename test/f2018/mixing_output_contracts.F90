program mixing_output_contracts
    use, intrinsic :: iso_fortran_env, only: output_unit, real64
    use atomic_calculation_state, only: electron_count
    use orbital_csf_basis, only: orbital_count
    use hamiltonian_symmetry_blocks, only: symmetry_block_count, csf_counts_by_block, eigenvector_counts_by_block
    use level_symmetry_metadata, only: packed_block_symmetries
    use csf_descriptor_reader, only: loaded_csf_basis_id
    use input_error_handling, only: abort_input
    use mixing_file_format, only: open_mixing_stream, write_mixing_block, close_mixing_stream
#ifdef CHECK_MPI
    use mpi_state, only: mpi_init, mpi_comm_rank, mpi_comm_world, mpi_barrier, mpi_finalize
    use mpi_io_errors, only: abort_input_mpi
#endif
    implicit none (type, external)
    character(24) :: scenario
#ifdef CHECK_MPI
    integer :: rank, mpi_status
#endif

    call get_command_argument(1,scenario)
#ifdef CHECK_MPI
    call mpi_init(mpi_status)
    call mpi_comm_rank(mpi_comm_world,rank,mpi_status)
    if (rank == 0) call exercise(abort_input_mpi)
    call mpi_barrier(mpi_comm_world,mpi_status)
    call mpi_finalize(mpi_status)
#else
    call exercise()
#endif
contains
    subroutine exercise(abort_handler)
        procedure(abort_input), optional :: abort_handler
        integer :: unit, status

        electron_count = 2
        orbital_count = 2
        loaded_csf_basis_id = repeat("a",64)
        symmetry_block_count = 2
        allocate(csf_counts_by_block(0:2), eigenvector_counts_by_block(2))
        csf_counts_by_block = [0,2,2]
        eigenvector_counts_by_block = [1,1]
        packed_block_symmetries(:2) = [1,-3]
        call open_mixing_stream("fixture.mix",unit,status)
        if (status /= 0) error stop "cannot open mixing fixture"
        if (scenario == "order") call write_mixing_block(unit,2,[1],[-1.0_real64],[0.0_real64,1.0_real64],abort_handler)
        call write_mixing_block(unit,1,[2],[-2.0_real64],[1.0_real64,0.0_real64],abort_handler)
        if (scenario == "incomplete") call close_mixing_stream(unit,abort_handler)
        call write_mixing_block(unit,2,[1],[-1.0_real64],[0.0_real64,1.0_real64],abort_handler)
        if (scenario == "limited-rewrite") then
            print *, "First pass complete"
            flush(output_unit)
            read(*,*)
        end if
        call write_mixing_block(unit,1,[2],[-3.0_real64],[1.0_real64,0.0_real64],abort_handler)
        call write_mixing_block(unit,2,[1],[-1.5_real64],[0.0_real64,1.0_real64],abort_handler)
        call close_mixing_stream(unit,abort_handler)
    end subroutine exercise
end program mixing_output_contracts
