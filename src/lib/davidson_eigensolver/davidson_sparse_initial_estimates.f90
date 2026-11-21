module davidson_sparse_initial_estimates
    implicit none (type, external)

    ! Module interface
    private
    public :: estimate_davidson_sparse_eigenpairs
contains
    subroutine estimate_davidson_sparse_eigenpairs(full_matrix_order, block_order, initial_vector_count, stored_matrix_values, &
        column_end_offsets, stored_row_indices, initial_estimates, csf_block_indices, target_block)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        ! Diagonalize a leading principal submatrix of the requested sparse symmetry block for initial guesses.
        use, intrinsic :: iso_fortran_env, only: real64
        use real_vector_initialization, only: initialize_real_vector
        use blas_lapack_interfaces, only: dcopy, dspevx
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: full_matrix_order
        integer, intent(in) :: block_order
        integer  :: initial_vector_count
        integer, intent(in) :: target_block
        integer, intent(in) :: column_end_offsets(0:*)
        integer, intent(in) :: stored_row_indices(*)
        integer, intent(in) :: csf_block_indices(*)
        real(real64), intent(in) :: stored_matrix_values(*)
        real(real64)  :: initial_estimates(*)

        ! Local declarations
        integer :: estimate_matrix_order, block_column, column_index, storage_index, row_index_shift, &
            computed_pair_count, lapack_status, status_code, block_row

        integer, allocatable :: lapack_integer_work(:), failed_vector_indices(:), block_to_full_index(:)
        real(real64), allocatable :: packed_estimate_matrix(:), initial_eigenvalues(:), &
            initial_block_vectors(:), lapack_real_work(:)

        estimate_matrix_order = min(1000,block_order)
        allocate (packed_estimate_matrix(estimate_matrix_order*(estimate_matrix_order + 1)/2),stat=status_code)
        if (status_code/=0) call mem_fail(6,estimate_matrix_order*(estimate_matrix_order + 1)/2,"iniest::ap", &
            status_code)

        call initialize_real_vector (estimate_matrix_order*(estimate_matrix_order + 1)/2, 0.0_real64, packed_estimate_matrix, 1)

        allocate (initial_eigenvalues(estimate_matrix_order),stat=status_code)
        if (status_code/=0) call mem_fail(6,estimate_matrix_order,"iniest::eigval",status_code)

        allocate (initial_block_vectors(estimate_matrix_order*initial_vector_count),stat=status_code)
        if (status_code/=0) call mem_fail(6,estimate_matrix_order*initial_vector_count,"iniest::vec",status_code)

        allocate (lapack_real_work(8*estimate_matrix_order),stat=status_code)
        if (status_code/=0) call mem_fail(6,8*estimate_matrix_order,"iniest::work",status_code)

        allocate (lapack_integer_work(5*estimate_matrix_order),stat=status_code)
        if (status_code/=0) call mem_fail(6,5*estimate_matrix_order,"iniest::iwork",status_code)

        allocate (failed_vector_indices(estimate_matrix_order),stat=status_code)
        if (status_code/=0) call mem_fail(6,estimate_matrix_order,"iniest::IFAIL",status_code)

        allocate (block_to_full_index(estimate_matrix_order),stat=status_code)
        if (status_code/=0) call mem_fail(6,estimate_matrix_order,"iniest::jsh",status_code)

        ! Extract the block columns in lower packed storage and retain the full-space index map.
        block_column = 0

        do column_index = 1, full_matrix_order
            if (block_column >= estimate_matrix_order) exit
            if (csf_block_indices(column_index) /= target_block) cycle
            block_column = block_column + 1
            block_to_full_index(block_column) = column_index
            storage_index = column_end_offsets(column_index-1) + 1
            row_index_shift = column_index - block_column
            do
                packed_estimate_matrix(stored_row_indices(storage_index)-row_index_shift+(block_column-1) &
                    *(2*estimate_matrix_order-block_column)/2) = stored_matrix_values(storage_index)
                storage_index = storage_index + 1
                if (storage_index > column_end_offsets(column_index)) exit
                ! Translate stored row indices by counting indices outside the selected block.
                row_index_shift = row_index_shift + count(csf_block_indices(stored_row_indices(storage_index &
                    -1):stored_row_indices(storage_index))/=target_block)
                if (stored_row_indices(storage_index) - row_index_shift > estimate_matrix_order) exit
            end do
        end do
        call dspevx ("Vectors also", "In a range", "Lower triangular", estimate_matrix_order, &
            packed_estimate_matrix, -1.0_real64, -1.0_real64, 1, initial_vector_count, 0.0_real64, &
            computed_pair_count, initial_eigenvalues, initial_block_vectors, estimate_matrix_order, &
            lapack_real_work, lapack_integer_work, failed_vector_indices, lapack_status)
        status_code = -abs(lapack_status)
        if (status_code /= 0) write (output_unit, *) "estimate_davidson_sparse_eigenpairs ierr =", status_code

        ! Scatter the small eigenvectors into zero-filled full-space columns.
        call initialize_real_vector (full_matrix_order*initial_vector_count, 0.0_real64, initial_estimates, 1)
        storage_index = 0
        do column_index = 1, initial_vector_count
            do block_row = 1, estimate_matrix_order
                initial_estimates(storage_index+block_to_full_index(block_row)) &
                    = initial_block_vectors((column_index - 1)*estimate_matrix_order + block_row)
            end do
            storage_index = storage_index + full_matrix_order
        end do
        ! Append eigenvalues after the stacked full-space eigenvectors.
        call dcopy (initial_vector_count, initial_eigenvalues, 1, &
            initial_estimates(initial_vector_count*full_matrix_order+1), 1)
        deallocate(packed_estimate_matrix)
        deallocate(initial_eigenvalues)
        deallocate(initial_block_vectors)
        deallocate(lapack_real_work)
        deallocate(lapack_integer_work)
        deallocate(failed_vector_indices)
        deallocate(block_to_full_index)

        return
    contains
        subroutine mem_fail(output_unit, element_count, variable_name, allocation_status)
            ! Arguments
            integer, intent(in) :: output_unit, element_count, allocation_status
            character(*), intent(in) :: variable_name
            write(output_unit, "(A,A,A,I0,A,I0)") "Unable to allocate ", variable_name, &
                ": elements=", element_count, ", status=", allocation_status
            error stop "INIEST allocation failed"
        end subroutine mem_fail
    end subroutine estimate_davidson_sparse_eigenpairs
end module davidson_sparse_initial_estimates
