module davidson_initial_subspace
    implicit none (type, external)

    ! Module interface
    private
    public :: initialize_davidson_subspace
contains
    subroutine initialize_davidson_subspace(operation_code, request_slots, matrix_order, constraint_count, initial_vector_count, &
        minimum_initial_basis_size, subspace_limit, solve_high_end, orthogonalization_work, &
        basis_with_constraints, basis_vectors, operator_basis, projected_matrix)
        ! Module dependencies
        ! Extend the supplied basis with orthogonalized matrix-vector products and build its projection.
        use, intrinsic :: iso_fortran_env, only: real64
        use real_vector_initialization, only: initialize_real_vector
        use davidson_orthonormalization, only: orthonormalize_davidson_vectors
        use davidson_projected_matrix, only: update_davidson_projected_matrix
        use blas_lapack_interfaces, only: dcopy, dscal
        implicit none (type, external)

        ! Arguments
        integer, intent(inout) :: operation_code
        integer, intent(in) :: matrix_order
        integer, intent(in) :: constraint_count
        integer, intent(inout) :: initial_vector_count
        integer, intent(in) :: minimum_initial_basis_size
        integer, intent(in) :: subspace_limit
        logical, intent(in) :: solve_high_end
        integer, intent(out) :: request_slots(*)
        real(real64)  :: orthogonalization_work(*)
        ! The combined array starts with the constraints and overlaps basis_vectors after them.
        real(real64)  :: basis_with_constraints(matrix_order*(constraint_count + subspace_limit))
        real(real64)  :: basis_vectors(matrix_order*subspace_limit)
        real(real64)  :: operator_basis(matrix_order*subspace_limit)
        real(real64)  :: projected_matrix(*)

        ! Local declarations
        integer :: first_batch_vector, last_basis_vector, batch_vector_count, subspace_size

        integer :: target_basis_size

        save first_batch_vector, last_basis_vector, batch_vector_count, subspace_size
        target_basis_size = min(minimum_initial_basis_size, subspace_limit)
        if (operation_code /= 2) then
            if (initial_vector_count == 0) then
                call initialize_real_vector (matrix_order, 1.0_real64/sqrt(dble(matrix_order)), basis_vectors, 1)
                initial_vector_count = 1
            end if
            first_batch_vector = 1
            last_basis_vector = initial_vector_count
            batch_vector_count = initial_vector_count
        else if (last_basis_vector < target_basis_size) then
            ! Reuse the products from the previous batch to generate new trial vectors.
            batch_vector_count = min(target_basis_size - last_basis_vector, &
                last_basis_vector - first_batch_vector + 1)
            call dcopy (matrix_order*batch_vector_count, operator_basis((first_batch_vector-1)*matrix_order+1), 1, &
                basis_vectors(1+last_basis_vector*matrix_order), 1)
            call orthonormalize_davidson_vectors (matrix_order, constraint_count + last_basis_vector, batch_vector_count, &
                orthogonalization_work, basis_with_constraints)
            first_batch_vector = last_basis_vector + 1
            last_basis_vector = last_basis_vector + batch_vector_count
        else
            initial_vector_count = last_basis_vector
            ! Apply the spectral sign to the completed products before forming the first projection.
            if (solve_high_end) call dscal (matrix_order*initial_vector_count, -1.0_real64, operator_basis, 1)
            subspace_size = 0
            call update_davidson_projected_matrix (matrix_order, subspace_limit, subspace_size, initial_vector_count, &
                & basis_vectors, &
                operator_basis, projected_matrix)

            operation_code = 0
            return
        end if

        operation_code = 2
        request_slots(1) = batch_vector_count
        request_slots(2) = first_batch_vector
        request_slots(3) = first_batch_vector
    end subroutine initialize_davidson_subspace
end module davidson_initial_subspace
