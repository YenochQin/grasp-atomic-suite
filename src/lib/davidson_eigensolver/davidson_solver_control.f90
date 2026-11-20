! Copyright (c) by Andreas Stathopoulos and Charlotte F. Fischer
module davidson_solver_control
    implicit none (type, external)
    private
    public :: dvdson
contains
      subroutine dvdson(operation_code, request_slots, matrix_order, subspace_limit, constraint_count, &
         lowest_pair_index, highest_pair_index, selected_pair_indices, initial_vector_count, &
         expansion_block_size, eigenvalue_tolerance, coefficient_tolerance, residual_tolerance, &
         iteration_limit, real_workspace, real_workspace_size, integer_workspace, integer_workspace_size, &
         solve_high_end, iteration_count, solver_status)
      ! Set up the target pairs and work arrays, then dispatch reverse-communication requests.
      use, intrinsic :: iso_fortran_env, only: real64
      use davidson_initial_subspace, only: initdvd
      use davidson_iteration, only: dvdrvr
      use blas_lapack_interfaces, only: dcopy, dscal
      implicit none (type, external)
      integer, intent(inout) :: operation_code
      integer, intent(in) :: matrix_order
      integer, intent(in) :: subspace_limit
      integer, intent(in) :: constraint_count
      integer, intent(inout) :: lowest_pair_index
      integer, intent(inout) :: highest_pair_index
      integer, intent(inout) :: initial_vector_count
      integer, intent(in) :: expansion_block_size
      integer  :: iteration_limit
      integer, intent(in) :: real_workspace_size
      integer, intent(in) :: integer_workspace_size
      integer, intent(inout) :: iteration_count
      integer, intent(inout) :: solver_status
      real(real64)  :: eigenvalue_tolerance
      real(real64)  :: coefficient_tolerance
      real(real64)  :: residual_tolerance
      logical, intent(inout) :: solve_high_end
      integer, intent(inout) :: request_slots(*)
      integer, intent(inout) :: selected_pair_indices(subspace_limit)
      integer, intent(inout) :: integer_workspace(integer_workspace_size)
      real(real64), intent(inout) :: real_workspace(real_workspace_size)
      real(kind(0.0_real64)) :: unused_dot_product
      integer :: selected_pair_count, selection_position, extreme_pair_count, constrained_basis_start, &
         basis_start, ritz_values_start, operator_basis_start, projected_matrix_start, &
         lapack_packed_copy_start, ritz_coefficients_start, real_scratch_start, previous_ritz_values_start, &
         integer_scratch_start, pending_indices_start, convergence_flags_start
      save
      ! Codes 1 and 3 resume the iteration. Code 2 resumes initialization.
      if (operation_code==1 .or. operation_code==3) then
         goto 200
      else if (operation_code == 2) then
         goto 100
      end if
      solver_status = 0
      if (subspace_limit > matrix_order) solver_status = solver_status + 1
      if (subspace_limit <= 0) solver_status = solver_status + 2

      solve_high_end = .false.

      ! Use an explicit index list when the range start is outside the spectrum.
      if (lowest_pair_index <= 0 .or. lowest_pair_index > matrix_order) then
         if (selected_pair_indices(1) <= 0) then
            solver_status = solver_status + 4
         else
            selected_pair_count = 1
            lowest_pair_index = selected_pair_indices(1)
            highest_pair_index = selected_pair_indices(1)
            do selection_position = 2, subspace_limit
               if (selected_pair_indices(selection_position) <= 0) exit
               lowest_pair_index = min(lowest_pair_index,selected_pair_indices(selection_position))
               highest_pair_index = max(highest_pair_index,selected_pair_indices(selection_position))
               selected_pair_count = selected_pair_count + 1
            end do
            if (highest_pair_index > matrix_order) solver_status = solver_status + 8
         end if
      else
         if (highest_pair_index > matrix_order) solver_status = solver_status + 8
         selected_pair_count = highest_pair_index - lowest_pair_index + 1
         if (selected_pair_count <= 0) solver_status = solver_status + 16
         if (selected_pair_count > subspace_limit) then
            solver_status = solver_status + 32
         else
            do selection_position = 1, selected_pair_count
               selected_pair_indices(selection_position) = lowest_pair_index + selection_position - 1
            end do
         end if
      end if

      if (solver_status /= 0) return

      extreme_pair_count = highest_pair_index
      ! Map upper-end indices to the low end of the negated spectrum.
      if (lowest_pair_index + highest_pair_index - 1 > matrix_order) then
         solve_high_end = .true.
         extreme_pair_count = matrix_order - lowest_pair_index + 1
         selection_position = 1
         if (selected_pair_count > 0) then
            selected_pair_indices(:selected_pair_count) = matrix_order &
               - selected_pair_indices(:selected_pair_count) + 1
            selection_position = selected_pair_count + 1
         end if
      end if
      if (selected_pair_count > extreme_pair_count) solver_status = solver_status + 64
      if (extreme_pair_count > subspace_limit .or. &
         extreme_pair_count == subspace_limit .and. extreme_pair_count /= matrix_order) &
         solver_status = solver_status + 128
      if (expansion_block_size < 1 .or. expansion_block_size > selected_pair_count) &
         solver_status = solver_status + 256

      if (real_workspace_size < subspace_limit*(2*matrix_order + subspace_limit + (extreme_pair_count + 10)) &
         + extreme_pair_count + constraint_count*matrix_order .or. &
         integer_workspace_size < 6*subspace_limit + extreme_pair_count) solver_status = solver_status + 512

      if (initial_vector_count > subspace_limit) then
         write (6, *) "WARNING: Too many initial estimates.?"
         write (6, *) "The routine needs at most:", subspace_limit
         solver_status = solver_status + 4096
      else if (initial_vector_count < extreme_pair_count) then
         write (6, *) "WARNING: Not enough initial estimates"
         write (6, *) extreme_pair_count - initial_vector_count, " Lanczos vectors will be added"
      end if

      if (solver_status /= 0) return
      ! Partition the work arrays into the basis, operator products, projection and solver scratch.
      constrained_basis_start = 1
      basis_start = constrained_basis_start + matrix_order*constraint_count
      ritz_values_start = basis_start + matrix_order*subspace_limit
      operator_basis_start = ritz_values_start + subspace_limit
      projected_matrix_start = operator_basis_start + matrix_order*subspace_limit
      lapack_packed_copy_start  = projected_matrix_start      + subspace_limit*(subspace_limit + 1)/2
      ritz_coefficients_start   = lapack_packed_copy_start  + subspace_limit*(subspace_limit + 1)/2
      real_scratch_start  = ritz_coefficients_start   + subspace_limit*subspace_limit
      previous_ritz_values_start = real_scratch_start  + 8*subspace_limit
      integer_scratch_start = 1
      pending_indices_start = integer_scratch_start + 5*subspace_limit
      convergence_flags_start = pending_indices_start + subspace_limit

  100 continue
      call initdvd (operation_code, request_slots, matrix_order, constraint_count, initial_vector_count, &
         extreme_pair_count + 1, subspace_limit, solve_high_end, real_workspace(real_scratch_start), &
         real_workspace(constrained_basis_start), real_workspace(basis_start), &
         real_workspace(operator_basis_start), real_workspace(projected_matrix_start))
      ! Convert the basis and product column numbers to one-based workspace offsets.
      if (operation_code == 2) then
         request_slots(2) = basis_start + (request_slots(2)-1)*matrix_order
         request_slots(3) = operator_basis_start + (request_slots(3)-1)*matrix_order
         return
      end if
      iteration_count = 1

  200 continue
      call dvdrvr (operation_code, request_slots, matrix_order, solve_high_end, subspace_limit, &
         expansion_block_size, constraint_count, extreme_pair_count, initial_vector_count, selected_pair_count, &
         selected_pair_indices, eigenvalue_tolerance, coefficient_tolerance, residual_tolerance, &
         iteration_limit, real_workspace(ritz_values_start), real_workspace(basis_start), &
         real_workspace(constrained_basis_start), real_workspace(operator_basis_start), &
         real_workspace(projected_matrix_start), real_workspace(lapack_packed_copy_start), &
         real_workspace(ritz_coefficients_start), real_workspace(real_scratch_start), &
         integer_workspace(integer_scratch_start), integer_workspace(pending_indices_start), &
         integer_workspace(convergence_flags_start), real_workspace(previous_ritz_values_start), &
         iteration_count, solver_status)
      ! Preconditioning reads the residual from the product workspace and writes a new basis vector.
      if (operation_code == 1) then
         request_slots(2) = operator_basis_start + (request_slots(2)-1)*matrix_order
         request_slots(3) = basis_start + (request_slots(3)-1)*matrix_order
         request_slots(4) = ritz_values_start
         request_slots(5) = pending_indices_start
         request_slots(6) = real_scratch_start
         request_slots(7) = real_scratch_start + subspace_limit
         return
      else if (operation_code == 3) then
         request_slots(2) = basis_start + (request_slots(2)-1)*matrix_order
         request_slots(3) = operator_basis_start + (request_slots(3)-1)*matrix_order
         return
      end if

      ! Return Ritz vectors, eigenvalues, eigenvalue changes and residual norms in contiguous blocks.
      if (solve_high_end) call dscal (extreme_pair_count, -1.0_real64, real_workspace(ritz_values_start), 1)
      call dcopy (extreme_pair_count, real_workspace(ritz_values_start), 1, &
         real_workspace(basis_start+matrix_order*extreme_pair_count), 1)
      call dcopy (extreme_pair_count, real_workspace(previous_ritz_values_start), 1, &
         real_workspace(basis_start+(matrix_order+1)*extreme_pair_count), 1)
      call dcopy (extreme_pair_count, real_workspace(real_scratch_start), 1, &
         real_workspace(basis_start+(matrix_order+2)*extreme_pair_count), 1)
      operation_code = 0
      return
      end subroutine dvdson
end module davidson_solver_control
