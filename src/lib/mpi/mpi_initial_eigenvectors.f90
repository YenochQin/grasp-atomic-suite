module mpi_initial_eigenvectors
    implicit none (type, external)

    ! Module interface
    private
    public :: estimate_mpi_initial_eigenvectors
contains
    ! Build Davidson starting vectors from the leading principal submatrix of the Hamiltonian.
    subroutine estimate_mpi_initial_eigenvectors (submatrix_order_limit, csf_count, initial_vector_count, basis_vectors, &
        hamiltonian_values, stored_segment_end_offsets, matrix_partner_indices)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use array_memory_management, only: alloc, dalloc
        use mpi_state, only: mpi_error_code, mpi_process_count, mpi_rank

        use blas_lapack_interfaces, only: dcopy, dspevx
        use real_vector_initialization, only: initialize_real_vector
        use mpi_utilities, only: allreduce_mpi_real_sum
        use blas_lapack_interfaces, only: dlamch
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: submatrix_order_limit
        integer, intent(in) :: csf_count
        integer  :: initial_vector_count
        integer, intent(in) :: stored_segment_end_offsets(0:*)
        integer, intent(in) :: matrix_partner_indices(*)
        real(real64)  :: basis_vectors(*)
        real(real64), intent(in) :: hamiltonian_values(*)

        ! Local declarations
        integer :: submatrix_order, sparse_segment_offset, column_index, packed_column_offset, stored_element_index
        integer :: computed_pair_count, lapack_status

        integer, pointer :: lapack_integer_workspace(:),failed_vector_indices(:)
        real(real64), pointer :: packed_submatrix(:), initial_eigenvalues(:),submatrix_eigenvectors(:), lapack_real_workspace(:)
        real(real64) :: lapack_abs_tolerance

        submatrix_order = min(submatrix_order_limit,csf_count)
        call alloc (packed_submatrix, (submatrix_order*(submatrix_order + 1))/2, "packed_submatrix", &
            "estimate_mpi_initial_eigenvectors")
        call initialize_real_vector ((submatrix_order*(submatrix_order + 1))/2, 0.0_real64, packed_submatrix, 1)

        ! Expand the rank-local column segments into LAPACK upper-triangle packed storage.
        sparse_segment_offset = 0
        do column_index = mpi_rank + 1, submatrix_order, mpi_process_count
            packed_column_offset = (column_index*(column_index - 1))/2
            do stored_element_index = sparse_segment_offset + 1, stored_segment_end_offsets(column_index)
                packed_submatrix(matrix_partner_indices(stored_element_index)+packed_column_offset) = &
                    hamiltonian_values(stored_element_index)
            end do
            sparse_segment_offset = stored_segment_end_offsets(column_index)
        end do

        ! Sum the packed contributions so each rank diagonalizes the same submatrix.
        call allreduce_mpi_real_sum (packed_submatrix, (submatrix_order*(submatrix_order+1))/2)

        call alloc(initial_eigenvalues, submatrix_order,     "initial_eigenvalues", "estimate_mpi_initial_eigenvectors")
        call alloc (submatrix_eigenvectors,   submatrix_order*initial_vector_count, "submatrix_eigenvectors", &
            "estimate_mpi_initial_eigenvectors")
        call alloc (lapack_real_workspace,  8*submatrix_order,   "lapack_real_workspace",   "estimate_mpi_initial_eigenvectors")
        call alloc (lapack_integer_workspace, 5*submatrix_order,   "lapack_integer_workspace",  "estimate_mpi_initial_eigenvectors")
        call alloc (failed_vector_indices, submatrix_order,     "failed_vector_indices",  "estimate_mpi_initial_eigenvectors")

        lapack_abs_tolerance = 2*dlamch("S")

        ! Select the first initial_vector_count eigenpairs in ascending eigenvalue order.
        call dspevx ("Vectors also", "In a range", "Upper triangular", submatrix_order, packed_submatrix, &
            -1.0_real64, -1.0_real64, 1, initial_vector_count, lapack_abs_tolerance, computed_pair_count, &
            initial_eigenvalues, submatrix_eigenvectors, submatrix_order, lapack_real_workspace, &
            lapack_integer_workspace, failed_vector_indices, lapack_status)
        mpi_error_code = -abs(lapack_status)

        ! Zero-pad each starting vector to the full CSF dimension.
        ! Append the selected eigenvalues after the vector columns in basis_vectors.
        call initialize_real_vector (csf_count*initial_vector_count, 0.0_real64, basis_vectors, 1)

        do column_index = 1, initial_vector_count
            call dcopy (submatrix_order, &
                submatrix_eigenvectors(submatrix_order*(column_index-1)+1:submatrix_order*(column_index-1)+ &
                (submatrix_order)), &
                1, basis_vectors(csf_count*(column_index-1)+1), 1)
        end do

        call dcopy (initial_vector_count, initial_eigenvalues, 1, basis_vectors(initial_vector_count*csf_count+1), 1)

        call dalloc (packed_submatrix,     "packed_submatrix",     "estimate_mpi_initial_eigenvectors")
        call dalloc (initial_eigenvalues, "initial_eigenvalues", "estimate_mpi_initial_eigenvectors")
        call dalloc (submatrix_eigenvectors,    "submatrix_eigenvectors",    "estimate_mpi_initial_eigenvectors")
        call dalloc (lapack_real_workspace,   "lapack_real_workspace",   "estimate_mpi_initial_eigenvectors")
        call dalloc (lapack_integer_workspace,  "lapack_integer_workspace",  "estimate_mpi_initial_eigenvectors")
        call dalloc (failed_vector_indices,  "failed_vector_indices",  "estimate_mpi_initial_eigenvectors")

        return
    end subroutine estimate_mpi_initial_eigenvectors
end module mpi_initial_eigenvectors
