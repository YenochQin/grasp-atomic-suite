module davidson_subspace_restart
    implicit none (type, external)

    ! Module interface
    private
    public :: restart_davidson_subspace
contains
    subroutine restart_davidson_subspace(matrix_order, retained_vector_count, subspace_size, row_workspace, basis_vectors, &
        operator_basis, projected_matrix, ritz_coefficients, ritz_values)
        ! Module dependencies
        ! Restart with the leading Ritz vectors and transform their operator products in the same way.
        use, intrinsic :: iso_fortran_env, only: real64
        use davidson_basis_transform, only: transform_davidson_basis
        use real_vector_initialization, only: initialize_real_vector
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: matrix_order
        integer, intent(in) :: retained_vector_count
        integer, intent(inout) :: subspace_size
        real(real64)  :: row_workspace(retained_vector_count)
        real(real64)  :: basis_vectors(matrix_order*subspace_size)
        real(real64)  :: operator_basis(matrix_order*subspace_size)
        real(real64), intent(out) :: projected_matrix((subspace_size*(subspace_size + 1))/2)
        real(real64)  :: ritz_coefficients(subspace_size*retained_vector_count)
        real(real64) , intent(in) :: ritz_values(subspace_size)

        ! Local declarations
        integer :: packed_column_offset, coefficient_column_offset, pair_index
        call transform_davidson_basis (matrix_order, subspace_size, retained_vector_count, ritz_coefficients, row_workspace, &
            basis_vectors)
        call transform_davidson_basis (matrix_order, subspace_size, retained_vector_count, ritz_coefficients, row_workspace, &
            operator_basis)
        ! In the restarted basis the projection is diagonal and its coefficient matrix is the identity.
        call initialize_real_vector ((retained_vector_count*(retained_vector_count + 1))/2, 0.0_real64, projected_matrix, 1)
        call initialize_real_vector (retained_vector_count*retained_vector_count, 0.0_real64, ritz_coefficients, 1)
        packed_column_offset = 0
        coefficient_column_offset = 0
        do pair_index = 1, retained_vector_count
            projected_matrix(packed_column_offset+pair_index) = ritz_values(pair_index)
            ritz_coefficients(coefficient_column_offset+pair_index) = 1
            coefficient_column_offset = coefficient_column_offset + retained_vector_count
            packed_column_offset = packed_column_offset + pair_index
        end do

        subspace_size = retained_vector_count

        return
    end subroutine restart_davidson_subspace
end module davidson_subspace_restart
