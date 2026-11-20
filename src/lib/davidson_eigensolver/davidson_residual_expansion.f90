module davidson_residual_expansion
    implicit none (type, external)
    private
    public :: newvec
contains
      subroutine newvec(matrix_order, extreme_pair_count, subspace_limit, expansion_block_size, subspace_size, &
         residual_tolerance, candidate_count, pending_pair_indices, ritz_coefficients, ritz_values, &
         previous_ritz_values, operator_basis, basis_vectors, convergence_flags, residual_norms, &
         preconditioner_shifts, targets_converged)
      ! Form residuals in priority order and prepare up to one block of correction vectors.
      use, intrinsic :: iso_fortran_env, only: real64
      use blas_lapack_interfaces, only: daxpy, ddot, dgemv
      implicit none (type, external)
      integer, intent(in) :: matrix_order
      integer , intent(in) :: extreme_pair_count
      integer , intent(in) :: subspace_limit
      integer , intent(in) :: expansion_block_size
      integer, intent(in) :: subspace_size
      integer , intent(inout) :: candidate_count
      real(real64) , intent(in) :: residual_tolerance
      logical , intent(out) :: targets_converged
      integer , intent(inout) :: pending_pair_indices(extreme_pair_count)
      integer, intent(inout) :: convergence_flags(extreme_pair_count)
      real(real64)  :: ritz_coefficients(subspace_limit*extreme_pair_count)
      real(real64) , intent(in) :: ritz_values(subspace_limit)
      real(real64) , intent(in) :: previous_ritz_values(extreme_pair_count)
      real(real64)  :: operator_basis(matrix_order*subspace_limit)
      real(real64)  :: basis_vectors(matrix_order*subspace_limit)
      real(real64) , intent(inout) :: residual_norms(subspace_limit)
      real(real64) , intent(out) :: preconditioner_shifts(subspace_limit)
      integer :: expansion_start, added_vector_count, converged_candidate_count, expanded_size_limit, &
         expansion_vector_start, candidate_position, pair_index
      real(real64) :: squared_residual_norm, lower_ritz_gap, nearest_ritz_gap, upper_ritz_gap, &
         ritz_value_change, residual_norm, residual_shift
      targets_converged = .false.
      expansion_start = subspace_size*matrix_order + 1
      added_vector_count = 0
      converged_candidate_count = 0
      expanded_size_limit = min(subspace_limit,expansion_block_size + subspace_size)
      expansion_vector_start = expansion_start
      do candidate_position = 1, candidate_count
         pair_index = pending_pair_indices(candidate_position)

         call dgemv ("N", matrix_order, subspace_size, 1.0_real64, basis_vectors, matrix_order, &
            ritz_coefficients((pair_index-1)*subspace_size+1), 1, 0.0_real64, &
            basis_vectors(expansion_vector_start), 1)
         call dgemv ("N", matrix_order, subspace_size, 1.0_real64, operator_basis, matrix_order, &
            ritz_coefficients((pair_index-1)*subspace_size+1), 1, 0.0_real64, &
            operator_basis(expansion_vector_start), 1)
         call daxpy (matrix_order, (-ritz_values(pair_index)),basis_vectors(expansion_vector_start), 1, &
            operator_basis(expansion_vector_start), 1)
         squared_residual_norm = ddot(matrix_order,operator_basis(expansion_vector_start),1, &
            operator_basis(expansion_vector_start),1)
         residual_norms(pair_index) = sqrt(squared_residual_norm)

         ! Skip a converged candidate without consuming an expansion slot.
         if (residual_norms(pair_index) < residual_tolerance) then
            converged_candidate_count = converged_candidate_count + 1
            convergence_flags(pair_index) = 1
            if (converged_candidate_count < candidate_count) cycle
            targets_converged = .true.
            return
         else
            ! Estimate the local spectral gap from adjacent Ritz values.
            if (pair_index == 1) then
               lower_ritz_gap = 1.0e99_real64
               nearest_ritz_gap = abs(ritz_values(2)-ritz_values(1))
            else
               lower_ritz_gap = abs(ritz_values(pair_index)-ritz_values(pair_index-1))
               upper_ritz_gap = abs(ritz_values(pair_index+1)-ritz_values(pair_index))
               nearest_ritz_gap = min(lower_ritz_gap,upper_ritz_gap)
            end if
            ritz_value_change = abs(previous_ritz_values(pair_index)-ritz_values(pair_index))
            residual_norm = residual_norms(pair_index)

            if (nearest_ritz_gap > residual_norm) then
               residual_shift = sqrt(ritz_value_change*residual_norm)
            else
               residual_shift = min(ritz_value_change,sqrt(ritz_value_change))
            end if

            ! Shift the residual along its Ritz vector. The separate diagonal-preconditioner shift is zero.
            preconditioner_shifts(pair_index) = 0.0_real64
            call daxpy (matrix_order, (-residual_shift), basis_vectors(expansion_vector_start), 1, &
               operator_basis(expansion_vector_start), 1)

            added_vector_count = added_vector_count + 1
            pending_pair_indices(added_vector_count) = pair_index
            if (added_vector_count + subspace_size == expanded_size_limit) exit
            expansion_vector_start = expansion_vector_start + matrix_order
         end if
      end do

      candidate_count = added_vector_count

      return
      end subroutine newvec
end module davidson_residual_expansion
