! Exercise the RCI block solver with a shifted tridiagonal Hamiltonian of known spectrum.
program rci_storage
    use, intrinsic :: iso_fortran_env, only: int8, int64, real64, output_unit
    use calculation_file_units, only: summary_output_unit
    use ci_restart_units, only: restart_matrix_unit
    use csf_coupling_data, only: csf_angular_codes
    use diagonalize_hamiltonian_block_interface, only: diagonalize_hamiltonian_block
    use eigenstate_expansions, only: hamiltonian_energy_shift, shifted_eigenvalues, csf_mixing_coefficients
    use grasp_parameters, only: csf_angular_code_count, csf_total_j_slot, csf_parity_slot
    use orbital_csf_basis, only: csf_count
    use rci_matrix_file_position, only: nposition
    use selected_level_indices, only: active_level_count, highest_active_level_index, active_level_indices
    use sparse_hamiltonian_storage, only: stored_matrix_element_count
#ifdef CHECK_MPI
    use mpi_state, only: mpi_init, mpi_finalize, mpi_comm_rank, mpi_comm_size, mpi_comm_world, &
        mpi_rank, mpi_process_count, mpi_allreduce, mpi_integer8, mpi_max
#endif
    implicit none (type, external)

    integer, parameter :: order = 12, leading_order = 3, roots(*) = [1, 3]
    real(real64), parameter :: diagonal = -4.0_real64, coupling = -0.5_real64
    real(real64) :: matrix(order,order), angle, energy, matrix_norm, phase
    real(real64) :: expected(order), vectors(order,size(roots)), values(order)
    real(real64) :: energy_error, residual, coefficient_error, gram_error
    integer :: column, row, first, indices(order), entries, root, j_code, parity, memory_limit, read_status
#ifdef CHECK_MPI
    integer :: status
    integer(int64) :: max_rank_entries
#else
    integer, parameter :: mpi_rank = 0, mpi_process_count = 1
#endif
    character(16) :: storage
    character(80) :: summary_line
    character(:), allocatable :: storage_message

#ifdef CHECK_MPI
    call mpi_init(status)
    call mpi_comm_rank(mpi_comm_world, mpi_rank, status)
    call mpi_comm_size(mpi_comm_world, mpi_process_count, status)
#endif
    call get_command_argument(1, storage)
    memory_limit = huge(memory_limit)
    select case (storage)
    case ("disk")
        memory_limit = 0
        storage_message = "matrix stored on disc"
    case ("sparse")
#ifdef CHECK_MPI
        storage_message = "sparse matrix stored in memory"
#else
        storage_message = "matrix stored in sparse representation"
#endif
    case ("dense")
        storage_message = "full matrix stored"
    case default
        error stop "Unknown RCI storage test"
    end select

    csf_count = order
    allocate(csf_angular_codes(csf_angular_code_count,order), source=0_int8)
    csf_angular_codes(csf_total_j_slot,:) = 1
    csf_angular_codes(csf_parity_slot,:) = 1
    active_level_count = size(roots)
    highest_active_level_index = maxval(roots)
    allocate(active_level_indices(size(roots)), source=roots)
    hamiltonian_energy_shift = diagonal
    matrix = 0.0_real64
    matrix(1,1) = diagonal
    do column = 2, order
        matrix(column,column) = diagonal
        matrix(column-1,column) = coupling
        matrix(column,column-1) = coupling
    end do
    matrix_norm = norm2(reshape(matrix, [order**2]))

    open(newunit=summary_output_unit, status="scratch", form="formatted")
    open(newunit=restart_matrix_unit, status="scratch", form="unformatted")
    write(restart_matrix_unit) "preceding block"
    nposition = 1
    write(restart_matrix_unit) order, order, mpi_rank, mpi_process_count
    stored_matrix_element_count = 0_int64
    do column = mpi_rank + 1, order, mpi_process_count
        first = max(1, column-1)
        ! Store explicit zeros to exercise the packed dense path with the same operator.
        if (storage == "dense") first = 1
        entries = column - first + 1
        values(:entries) = matrix(first:column,column)
        indices(:entries) = [(row, row=first,column)]
        write(restart_matrix_unit) entries, 0.0_real64, values(:entries), indices(:entries)
        stored_matrix_element_count = stored_matrix_element_count + entries
    end do
#ifdef CHECK_MPI
    call mpi_allreduce(stored_matrix_element_count, max_rank_entries, 1, mpi_integer8, mpi_max, mpi_comm_world, status)
    call diagonalize_hamiltonian_block(j_code, parity, max_rank_entries, leading_order, memory_limit)
#else
    call diagonalize_hamiltonian_block(j_code, parity, leading_order, memory_limit)
#endif
    if (mpi_rank == 0) then
        rewind(summary_output_unit)
        do
            read(summary_output_unit, "(a)", iostat=read_status) summary_line
            if (read_status /= 0) error stop "Expected RCI storage branch absent from solver summary"
            if (index(summary_line, storage_message) > 0) exit
        end do
    end if
    if (j_code /= 1 .or. parity /= 1) error stop "Wrong RCI symmetry"
    vectors = reshape(csf_mixing_coefficients, shape(vectors))
    do root = 1, size(roots)
        angle = roots(root)*acos(-1.0_real64)/(order+1)
        expected = [(sin(row*angle), row=1,order)] * sqrt(2.0_real64/(order+1))
        energy = shifted_eigenvalues(root) + hamiltonian_energy_shift
        energy_error = abs(energy - (diagonal + 2*coupling*cos(angle)))
        phase = sign(1.0_real64, dot_product(vectors(:,root),expected))
        coefficient_error = maxval(abs(vectors(:,root) - phase*expected))
        residual = norm2(matmul(matrix,vectors(:,root))-energy*vectors(:,root)) / &
            ((matrix_norm+abs(energy))*norm2(vectors(:,root)))
        if (mpi_rank == 0) write(output_unit, '(a,a,i0,3(a,es10.3))') trim(storage), " root ", roots(root), &
            " energy error ", energy_error, " residual ", residual, " coefficient error ", coefficient_error
        if (.not. (energy_error <= 1.0e-10_real64 .and. coefficient_error <= 1.0e-8_real64 .and. &
            residual <= 1.0e-10_real64)) error stop "Wrong RCI eigenpair"
    end do
    gram_error = maxval(abs(matmul(transpose(vectors),vectors) - reshape([1,0,0,1],[2,2])))
    if (.not. (gram_error <= 1.0e-12_real64)) error stop "RCI eigenvectors are not orthonormal"
    close(restart_matrix_unit)
    close(summary_output_unit)
#ifdef CHECK_MPI
    call mpi_finalize(status)
#endif
    if (mpi_rank == 0) write(output_unit,'(a,es10.3)') "PASS RCI storage, orthogonality error ", gram_error
end program rci_storage
