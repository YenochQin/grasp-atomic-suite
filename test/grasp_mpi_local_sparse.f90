! Check the real MPI sparse kernels against a dense matrix.  Only the
! collective is replaced: SPICMV returns each local contribution for summing
! here; INIEST's local packed matrix is checked before supplying the global
! packed matrix to LAPACK.  This tests storage semantics without launching MPI.
module sparse_reduction_check
    use, intrinsic :: iso_fortran_env, only: real64
    use mpi_state, only: mpi_rank, mpi_process_count
    implicit none (type, external)
    private :: real64
    integer, parameter :: matrix_size = 8
    real(real64) :: dense(matrix_size, matrix_size)
    integer :: packed_size = 0
    real(real64) :: packed_error = 0.0_real64
contains
    subroutine check_reduction(x, n)
        integer, intent(in) :: n
        real(real64), intent(inout) :: x(n)
        integer :: i, j, k
        real(real64) :: expected_local

        if (packed_size == 0) return
        if (n /= packed_size * (packed_size + 1) / 2) error stop "packed size"
        k = 0
        do j = 1, packed_size
            do i = 1, j
                k = k + 1
                expected_local = 0.0_real64
                if (mod(j - 1, mpi_process_count) == mpi_rank) expected_local = dense(i, j)
                packed_error = max(packed_error, abs(x(k) - expected_local))
                x(k) = dense(i, j)
            end do
        end do
    end subroutine check_reduction
end module sparse_reduction_check

module mpi_utilities
    implicit none (type, external)
    private
    public :: allreduce_mpi_real_sum
contains

subroutine allreduce_mpi_real_sum(x, n)
    use sparse_reduction_check, only: check_reduction
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)
    integer, intent(in) :: n
    real(real64), intent(inout) :: x(n)
    call check_reduction(x, n)
end subroutine allreduce_mpi_real_sum
end module mpi_utilities



program grasp_mpi_local_sparse
    use sparse_reduction_check, only: check_reduction, dense, matrix_size, mpi_rank, mpi_process_count, packed_error, &
       & packed_size
    use sparse_hamiltonian_storage, only: hamiltonian_values, matrix_partner_indices, stored_segment_end_offsets, &
        stored_matrix_element_count
    use mpi_state, only: mpi_error_code
    use mpi_sparse_matrix_product, only: multiply_mpi_sparse_matrix
    use mpi_initial_eigenvectors, only: estimate_mpi_initial_eigenvectors
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)
    integer, parameter :: rhs_count = 3, root_count = 2
    integer, parameter :: rank_counts(*) = [1, 2, 3, 4, 46]
    real(real64) :: rhs(matrix_size, rhs_count), local(matrix_size, rhs_count)
    real(real64) :: total(matrix_size, rhs_count), product_error, eigen_error
    real(real64) :: basis((matrix_size + 1) * root_count), vector(matrix_size)
    real(real64) :: eigenvalue
    integer :: i, j, k, count_index, mode, nmax, niv, root, failures

    dense = 0.0_real64
    do i = 1, matrix_size
        dense(i, i) = 2.0_real64 * i
    end do
    ! Unequal column lengths, diagonals only in some columns, and distant
    ! couplings expose reuse of entries from an earlier local column.
    dense(1, 2) = 0.3_real64
    dense(1, 4) = -0.2_real64
    dense(3, 4) = 0.4_real64
    dense(2, 5) = 0.1_real64
    dense(1, 7) = 0.6_real64
    dense(4, 7) = -0.3_real64
    dense(3, 8) = 0.2_real64
    do j = 1, matrix_size
        do i = 1, j - 1
            dense(j, i) = dense(i, j)
        end do
        do k = 1, rhs_count
            rhs(j, k) = real(j - 2 * k, real64) / 7.0_real64
        end do
    end do

    allocate(stored_segment_end_offsets(0:matrix_size), matrix_partner_indices(matrix_size**2), hamiltonian_values(matrix_size**2))
    failures = 0
    do count_index = 1, size(rank_counts)
        mpi_process_count = rank_counts(count_index)
        total = 0.0_real64
        packed_error = 0.0_real64
        eigen_error = 0.0_real64
        do mpi_rank = 0, mpi_process_count - 1
            ! MATRIXmpi reads endpoints only for the columns owned by this
            ! rank.  Its EMT/IROW contain just those columns, packed together.
            stored_segment_end_offsets = 0
            hamiltonian_values = 0.0_real64
            matrix_partner_indices = 1
            k = 0
            do j = mpi_rank + 1, matrix_size, mpi_process_count
                do i = 1, j
                    if (dense(i, j) == 0.0_real64) cycle
                    k = k + 1
                    hamiltonian_values(k) = dense(i, j)
                    matrix_partner_indices(k) = i
                end do
                stored_segment_end_offsets(j) = k
            end do
            stored_matrix_element_count = k

            packed_size = 0
            call multiply_mpi_sparse_matrix(matrix_size, rhs_count, rhs, local)
            total = total + local

            ! Full and truncated initial subspaces, including ranks with no
            ! columns, exercise the INIEST caller's two size regimes.
            do mode = 1, 2
                nmax = matrix_size
                if (mode == 2) nmax = 5
                packed_size = nmax
                niv = root_count
                call estimate_mpi_initial_eigenvectors(nmax, matrix_size, niv, basis, hamiltonian_values, &
                    & stored_segment_end_offsets, &
                    matrix_partner_indices)
                if (mpi_error_code /= 0) error stop "INIEST LAPACK failure"
                do root = 1, root_count
                    vector = basis((root - 1) * matrix_size + 1 : root * matrix_size)
                    eigenvalue = basis(root_count * matrix_size + root)
                    eigen_error = max(eigen_error, maxval(abs( &
                        matmul(dense(:nmax, :nmax), vector(:nmax)) - eigenvalue * vector(:nmax))))
                    eigen_error = max(eigen_error, abs(dot_product(vector, vector) - 1.0_real64))
                    if (nmax < matrix_size) eigen_error = max(eigen_error, maxval(abs(vector(nmax + 1:))))
                end do
            end do
        end do
        product_error = maxval(abs(total - matmul(dense, rhs)))
        write(*, "(A,I0,3(A,ES12.4))") "partitions=", mpi_process_count, &
            " matmul_error=", product_error, " packed_error=", packed_error, " eigen_residual=", eigen_error
        if (max(product_error, packed_error, eigen_error) > 1.0e-10_real64) failures = failures + 1
    end do
    deallocate(stored_segment_end_offsets, matrix_partner_indices, hamiltonian_values)
    if (failures > 0) error stop "rank-local sparse storage mismatch"
end program grasp_mpi_local_sparse
