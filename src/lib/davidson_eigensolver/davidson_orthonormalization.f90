module davidson_orthonormalization
    implicit none (type, external)

    ! Module interface
    private
    public :: orthonormalize_davidson_vectors
contains
    subroutine orthonormalize_davidson_vectors(vector_length, existing_vector_count, new_vector_count, projection_coefficients, &
        vectors)
        ! Module dependencies
        ! Apply two passes of modified Gram-Schmidt to the appended columns.
        use, intrinsic :: iso_fortran_env, only: real64
        use blas_lapack_interfaces, only: ddot, dgemv, dscal
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: vector_length
        integer , intent(in) :: existing_vector_count
        integer, intent(in) :: new_vector_count
        real(real64), intent(inout) :: projection_coefficients(new_vector_count)
        real(real64), intent(inout) :: vectors((existing_vector_count + new_vector_count)*vector_length)

        ! Local declarations
        integer :: pass_index, new_vectors_start, reference_vector_start, reference_vector_index, &
            target_vector_start, target_vector_index, component_offset
        real(real64) :: vector_norm

        do pass_index = 1, 2

            new_vectors_start = existing_vector_count*vector_length + 1
            reference_vector_start = 1
            ! Remove projections onto the existing constraints and basis columns.
            do reference_vector_index = 1, existing_vector_count
                target_vector_start = new_vectors_start
                call dgemv ("T", vector_length, new_vector_count, 1.0_real64, vectors(target_vector_start), &
                    vector_length, vectors(reference_vector_start), 1, 0.0_real64, projection_coefficients, 1)
                do target_vector_index = 1, new_vector_count
                    do component_offset = 0, vector_length - 1
                        vectors(target_vector_start+component_offset) &
                            = vectors(target_vector_start+component_offset) &
                            - projection_coefficients(target_vector_index) &
                            *vectors(reference_vector_start+component_offset)
                    end do
                    target_vector_start = target_vector_start + vector_length
                end do
                reference_vector_start = reference_vector_start + vector_length
            end do
            ! Normalize each new column before projecting the subsequent new columns against it.
            do reference_vector_index = 1, new_vector_count
                target_vector_start = reference_vector_start + vector_length
                vector_norm = ddot(vector_length,vectors(reference_vector_start),1,vectors(reference_vector_start),1)
                vector_norm = sqrt(vector_norm)
                ! A small first-pass residual can still carry an independent direction.
                ! On pass two the input had unit norm, so an O(n*epsilon) remainder
                ! is at projection roundoff. Replace it before normalization.
                if (vector_norm <= tiny(vector_norm) .or. &
                    (pass_index == 2 .and. vector_norm <= vector_length*epsilon(vector_norm))) then
                    call complete_orthogonal_basis(vector_length, existing_vector_count + reference_vector_index - 1, &
                        vectors(1:reference_vector_start-1), vectors(reference_vector_start:target_vector_start-1))
                    vector_norm = norm2(vectors(reference_vector_start:target_vector_start-1))
                end if
                call dscal (vector_length, 1/vector_norm, vectors(reference_vector_start), 1)
                if (reference_vector_index < new_vector_count) then
                    call dgemv ("T", vector_length, new_vector_count - reference_vector_index, 1.0_real64, &
                        vectors(target_vector_start), vector_length, vectors(reference_vector_start), 1, 0.0_real64, &
                        projection_coefficients, 1)
                end if
                do target_vector_index = reference_vector_index + 1, new_vector_count
                    do component_offset = 0, vector_length - 1
                        vectors(target_vector_start+component_offset) &
                            = vectors(target_vector_start+component_offset) &
                            - projection_coefficients(target_vector_index-reference_vector_index) &
                            *vectors(reference_vector_start+component_offset)
                    end do
                    target_vector_start = target_vector_start + vector_length
                end do
                reference_vector_start = reference_vector_start + vector_length
            end do
        end do

        return
    end subroutine orthonormalize_davidson_vectors

    subroutine complete_orthogonal_basis(vector_length, basis_count, basis, vector)
        use, intrinsic :: iso_fortran_env, only: real64
        use blas_lapack_interfaces, only: daxpy, ddot
        implicit none (type, external)
        integer, intent(in) :: vector_length, basis_count
        real(real64), intent(in) :: basis(vector_length,basis_count)
        real(real64), intent(out) :: vector(vector_length)
        integer :: column, coordinate, pass_index
        real(real64) :: projection

        ! The coordinate with the smallest diagonal of Q Q^T has the largest
        ! orthogonal complement. Reuse the output column to accumulate this diagonal.
        vector = 0.0_real64
        do column = 1, basis_count
            vector = vector + basis(:,column)**2
        end do
        coordinate = minloc(vector,dim=1)
        vector = 0.0_real64
        vector(coordinate) = 1.0_real64
        do pass_index = 1, 2
            do column = 1, basis_count
                projection = ddot(vector_length, basis(:,column), 1, vector, 1)
                call daxpy(vector_length, -projection, basis(:,column), 1, vector, 1)
            end do
        end do
    end subroutine complete_orthogonal_basis
end module davidson_orthonormalization
