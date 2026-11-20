module davidson_projected_matrix
    implicit none (type, external)

    ! Module interface
    private
    public :: adds
contains
    subroutine adds(matrix_order, subspace_limit, subspace_size, new_vector_count, basis_vectors, &
        operator_basis, projected_matrix)
        ! Module dependencies
        ! Append projection columns from the stored operator products in packed upper-triangular form.
        use, intrinsic :: iso_fortran_env, only: real64
        use blas_lapack_interfaces, only: ddot
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: matrix_order
        integer , intent(in) :: subspace_limit
        integer , intent(in) :: subspace_size
        integer , intent(in) :: new_vector_count
        real(real64)  :: basis_vectors(matrix_order*subspace_limit)
        real(real64)  :: operator_basis(matrix_order*subspace_limit)
        real(real64), intent(inout) :: projected_matrix(subspace_limit*(subspace_limit + 1)/2)

        ! Local declarations
        integer :: operator_column_start, packed_column_offset, new_vector_index, basis_column_start, &
            basis_vector_index
        real(real64) :: projected_element
        operator_column_start = subspace_size*matrix_order + 1
        ! The one-based packed index for projection element (i,j) is i + j*(j-1)/2 for i <= j.
        ! The existing subspace occupies the first subspace_size*(subspace_size+1)/2 entries.
        packed_column_offset = subspace_size*(subspace_size + 1)/2
        do new_vector_index = 1, new_vector_count
            basis_column_start = 1
            do basis_vector_index = 1, subspace_size + new_vector_index
                projected_element = ddot(matrix_order,basis_vectors(basis_column_start),1, &
                    operator_basis(operator_column_start),1)
                projected_matrix(packed_column_offset+basis_vector_index) = projected_element
                basis_column_start = basis_column_start + matrix_order
            end do
            packed_column_offset = packed_column_offset + subspace_size + new_vector_index
            operator_column_start = operator_column_start + matrix_order
        end do

        return
    end subroutine adds
end module davidson_projected_matrix
