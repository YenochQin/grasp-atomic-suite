! Diagonalize the leading Hamiltonian submatrix to supply starting vectors and eigenvalues.
module sparse_eigensolver_initial_basis
    implicit none (type, external)

    ! Module interface
    private
    public :: estimate_sparse_initial_basis
contains
    subroutine estimate_sparse_initial_basis(maximum_submatrix_order, csf_count, requested_eigenvector_count, &
        & start_vectors_and_eigenvalues, &
        & sparse_matrix_values, column_end_offsets, row_indices)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use real_vector_initialization, only: initialize_real_vector
        use array_memory_management, only: alloc, dalloc
        use blas_lapack_interfaces, only: dcopy, dspevx
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: maximum_submatrix_order
        integer, intent(in) :: csf_count
        integer, intent(in) :: requested_eigenvector_count
        integer, intent(in) :: column_end_offsets(0:*)
        integer, intent(in) :: row_indices(*)
        real(real64) :: start_vectors_and_eigenvalues(*)
        real(real64), intent(in) :: sparse_matrix_values(*)

        ! Local declarations
        integer :: submatrix_order, previous_sparse_column_end, column_index, packed_column_offset, sparse_entry_index, &
            & found_eigenvector_count, lapack_status, negative_lapack_status

        integer, pointer :: integer_workspace(:),failed_eigenvector_indices(:)
        real(real64), pointer :: packed_submatrix(:), eigenvalues(:),eigenvectors(:), real_workspace(:)

        ! Convert the leading sparse block to LAPACK upper-triangle packed storage.
        submatrix_order = min(maximum_submatrix_order,csf_count)
        call alloc (packed_submatrix, (submatrix_order*(submatrix_order + 1))/2, "packed_submatrix", &
            "estimate_sparse_initial_basis")
        call initialize_real_vector ((submatrix_order*(submatrix_order + 1))/2, 0.0_real64, packed_submatrix, 1)

        previous_sparse_column_end = 0
        do column_index = 1, submatrix_order
            packed_column_offset = (column_index*(column_index - 1))/2
            do sparse_entry_index = previous_sparse_column_end + 1, column_end_offsets(column_index)
                packed_submatrix(row_indices(sparse_entry_index)+packed_column_offset) = sparse_matrix_values(sparse_entry_index)
            end do
            previous_sparse_column_end = column_end_offsets(column_index)
        end do

        call alloc(eigenvalues, submatrix_order, "eigenvalues", "estimate_sparse_initial_basis")
        call alloc (eigenvectors, submatrix_order*requested_eigenvector_count, "eigenvectors", "estimate_sparse_initial_basis")
        call alloc (real_workspace, 8*submatrix_order, "real_workspace", "estimate_sparse_initial_basis")
        call alloc (integer_workspace, 8*submatrix_order, "integer_workspace", "estimate_sparse_initial_basis")
        call alloc (failed_eigenvector_indices, submatrix_order, "failed_eigenvector_indices", "estimate_sparse_initial_basis")

        call dspevx ("Vectors also", "In a range", "Upper triangular", submatrix_order, packed_submatrix, &
            -1.0_real64, -1.0_real64, 1, requested_eigenvector_count, 0.0_real64, found_eigenvector_count, eigenvalues, &
            & eigenvectors, submatrix_order, real_workspace, integer_workspace, failed_eigenvector_indices, &
            lapack_status)
        negative_lapack_status = -abs(lapack_status)

        call initialize_real_vector (csf_count*requested_eigenvector_count, 0.0_real64, start_vectors_and_eigenvalues, 1)

        do column_index = 1, requested_eigenvector_count
            call dcopy (submatrix_order, eigenvectors(submatrix_order*(column_index-1)+1:submatrix_order*(column_index-1)+ &
                & (submatrix_order)), 1, start_vectors_and_eigenvalues(csf_count*(column_index-1)+1), 1)
        end do

        call dcopy (requested_eigenvector_count, eigenvalues, 1, &
            & start_vectors_and_eigenvalues(requested_eigenvector_count*csf_count+1), 1)

        call dalloc (packed_submatrix, "packed_submatrix", "estimate_sparse_initial_basis")
        call dalloc (eigenvalues, "eigenvalues", "estimate_sparse_initial_basis")
        call dalloc (eigenvectors, "eigenvectors", "estimate_sparse_initial_basis")
        call dalloc (real_workspace, "real_workspace", "estimate_sparse_initial_basis")
        call dalloc (integer_workspace, "integer_workspace", "estimate_sparse_initial_basis")
        call dalloc (failed_eigenvector_indices, "failed_eigenvector_indices", "estimate_sparse_initial_basis")

        return
    end subroutine estimate_sparse_initial_basis
end module sparse_eigensolver_initial_basis
