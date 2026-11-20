module davidson_convergence
    implicit none (type, external)
    private
    public :: tstsel
contains
      logical function tstsel (subspace_size, extreme_pair_count, selected_pair_count, selected_pair_indices, &
         ritz_coefficients, ritz_values, convergence_flags, eigenvalue_tolerance, coefficient_tolerance, &
         recent_coefficient_maxima, unsorted_pair_indices, previous_ritz_values, recent_or_pending_count, &
         pending_pair_indices) result (targets_converged)
      ! Check the selected Ritz values and rank unconverged pairs by their newest basis coefficients.
      use, intrinsic :: iso_fortran_env, only: real64
      use blas_lapack_interfaces, only: idamax
      implicit none (type, external)
      integer , intent(in) :: subspace_size
      integer , intent(in) :: extreme_pair_count
      integer , intent(in) :: selected_pair_count
      integer , intent(inout) :: recent_or_pending_count
      real(real64) , intent(in) :: eigenvalue_tolerance
      real(real64) , intent(in) :: coefficient_tolerance
      integer , intent(in) :: selected_pair_indices(selected_pair_count)
      integer, intent(inout) :: convergence_flags(extreme_pair_count)
      integer , intent(inout) :: unsorted_pair_indices(selected_pair_count)
      integer , intent(out) :: pending_pair_indices(selected_pair_count)
      real(real64) , intent(in) :: ritz_coefficients(subspace_size*extreme_pair_count)
      real(real64) , intent(in) :: ritz_values(extreme_pair_count)
      real(real64), intent(inout) :: recent_coefficient_maxima(selected_pair_count)
      real(real64) , intent(in) :: previous_ritz_values(extreme_pair_count)
      integer :: unconverged_value_count, target_position, pair_index, pending_count, last_coefficient_index, &
         coefficient_offset, relative_max_index, saved_pair_index
      real(real64) :: max_recent_coefficient, saved_coefficient_maximum
      logical :: coefficients_converged
      coefficients_converged = .false.
      unconverged_value_count = 0
      do target_position = 1, selected_pair_count
         pair_index = selected_pair_indices(target_position)
         if (abs(previous_ritz_values(pair_index)-ritz_values(pair_index)) < eigenvalue_tolerance) cycle
         unconverged_value_count = unconverged_value_count + 1
      end do
      if (unconverged_value_count == 0) then
         targets_converged = .true.
         return
      end if
      ! Only the coefficients on the most recently appended columns enter this test.
      pending_count = 0
      do target_position = 1, selected_pair_count
         last_coefficient_index = subspace_size*selected_pair_indices(target_position)
         max_recent_coefficient = abs(ritz_coefficients(last_coefficient_index))
         do coefficient_offset = 1, recent_or_pending_count - 1
            max_recent_coefficient = max(max_recent_coefficient, &
               abs(ritz_coefficients(last_coefficient_index-coefficient_offset)))
         end do
         if (max_recent_coefficient < coefficient_tolerance) then
            convergence_flags(selected_pair_indices(target_position)) = 1
         else
            convergence_flags(selected_pair_indices(target_position)) = 0
            pending_count = pending_count + 1
            unsorted_pair_indices(pending_count) = selected_pair_indices(target_position)
            recent_coefficient_maxima(pending_count) = max_recent_coefficient
         end if
      end do

      recent_or_pending_count = pending_count
      if (recent_or_pending_count == 0) then
         coefficients_converged = .true.
      end if
      ! Order the pending pairs by descending coefficient magnitude for residual evaluation.
      do target_position = 1, recent_or_pending_count
         relative_max_index = idamax(recent_or_pending_count - target_position + 1, &
            recent_coefficient_maxima(target_position),1)
         pending_pair_indices(target_position) = unsorted_pair_indices(relative_max_index+target_position-1)

         saved_coefficient_maximum = recent_coefficient_maxima(relative_max_index+target_position-1)
         recent_coefficient_maxima(relative_max_index+target_position-1) &
            = recent_coefficient_maxima(target_position)
         recent_coefficient_maxima(target_position) = saved_coefficient_maximum
         saved_pair_index = unsorted_pair_indices(relative_max_index+target_position-1)
         unsorted_pair_indices(relative_max_index+target_position-1) = unsorted_pair_indices(target_position)
         unsorted_pair_indices(target_position) = saved_pair_index
      end do

      targets_converged = coefficients_converged
      return
      end function tstsel
end module davidson_convergence
