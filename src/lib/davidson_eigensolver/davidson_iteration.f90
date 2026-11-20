module davidson_iteration
    implicit none (type, external)

    ! Module interface
    private
    public :: iterate_davidson_solver
contains
    subroutine iterate_davidson_solver(operation_code, request_slots, matrix_order, solve_high_end, subspace_limit, &
        expansion_block_size, constraint_count, extreme_pair_count, initial_vector_count, selected_pair_count, &
        selected_pair_indices, eigenvalue_tolerance, coefficient_tolerance, residual_tolerance, &
        iteration_limit, ritz_values, basis_vectors, basis_with_constraints, operator_basis, projected_matrix, &
        lapack_packed_copy, ritz_coefficients, real_scratch, integer_scratch, pending_pair_indices, &
        convergence_flags, previous_ritz_values, iteration_count, solver_status)
        ! Module dependencies
        ! Solve the projected eigenproblem and expand or restart the Davidson subspace.
        use, intrinsic :: iso_fortran_env, only: real64
        use blas_lapack_interfaces, only: slamch
        use davidson_subspace_restart, only: restart_davidson_subspace
        use davidson_convergence, only: check_davidson_convergence
        use davidson_residual_expansion, only: expand_davidson_residuals
        use davidson_orthonormalization, only: orthonormalize_davidson_vectors
        use davidson_projected_matrix, only: update_davidson_projected_matrix
        use davidson_basis_transform, only: transform_davidson_basis
        use blas_lapack_interfaces, only: daxpy, dcopy, ddot, dscal, dspevx
        implicit none (type, external)

        ! Arguments
        integer, intent(inout) :: operation_code
        integer, intent(in) :: matrix_order
        integer, intent(in) :: subspace_limit
        integer  :: expansion_block_size
        integer, intent(in) :: constraint_count
        integer, intent(in) :: extreme_pair_count
        integer, intent(in) :: initial_vector_count
        integer, intent(in) :: selected_pair_count
        integer, intent(in) :: iteration_limit
        integer, intent(inout) :: iteration_count
        integer, intent(inout) :: solver_status
        real(real64)  :: eigenvalue_tolerance
        real(real64)  :: coefficient_tolerance
        real(real64)  :: residual_tolerance
        logical, intent(in) :: solve_high_end
        integer, intent(out) :: request_slots(*)
        integer  :: selected_pair_indices(selected_pair_count)
        integer  :: integer_scratch(5*subspace_limit)
        integer, intent(inout) :: pending_pair_indices(subspace_limit)
        integer, intent(inout) :: convergence_flags(extreme_pair_count + 1)
        real(real64), intent(inout) :: ritz_values(subspace_limit)
        real(real64)  :: basis_vectors(matrix_order*subspace_limit)
        real(real64)  :: basis_with_constraints(matrix_order*subspace_limit + constraint_count*matrix_order)
        real(real64)  :: operator_basis(matrix_order*subspace_limit)
        real(real64)  :: projected_matrix(subspace_limit*(subspace_limit + 1)/2)
        real(real64)  :: lapack_packed_copy(subspace_limit*(subspace_limit + 1)/2)
        real(real64)  :: ritz_coefficients(subspace_limit*subspace_limit)
        real(real64), intent(inout) :: real_scratch(8*subspace_limit)
        real(real64), intent(inout) :: previous_ritz_values(extreme_pair_count + 1)

        ! Local declarations
        integer :: restart_size, pair_position, subspace_size, expansion_count, requested_ritz_count, &
            computed_ritz_count, lapack_status, expansion_start
        real(real64) :: lapack_abs_tolerance
        logical :: first_iteration, targets_converged, iteration_limit_reached

        save first_iteration, targets_converged, restart_size, pair_position, subspace_size, expansion_count, &
            requested_ritz_count, lapack_abs_tolerance, computed_ritz_count, lapack_status, expansion_start
        iteration_limit_reached = .false.
        select case (operation_code)
        case (1)
            ! Restore the working Ritz values, then orthonormalize the correction vectors.
            if (solve_high_end) call dscal (extreme_pair_count, -1.0_real64, ritz_values, 1)
            do pair_position = 1, expansion_count
                ritz_values(pending_pair_indices(pair_position)) = ritz_values(pending_pair_indices(pair_position)) &
                    + real_scratch(subspace_limit+pending_pair_indices(pair_position))
            end do
            expansion_start = subspace_size*matrix_order + 1

            call orthonormalize_davidson_vectors (matrix_order, constraint_count + subspace_size, expansion_count, &
                real_scratch(subspace_limit+1), basis_with_constraints)
            request_slots(1) = expansion_count
            request_slots(2) = subspace_size + 1
            request_slots(3) = subspace_size + 1
            operation_code = 3
            return
        case (3)
            ! Apply the spectral sign to the new products and append the projection columns.
            if (solve_high_end) call dscal (matrix_order*expansion_count, -1.0_real64, &
                operator_basis(expansion_start), 1)
            call update_davidson_projected_matrix (matrix_order, subspace_limit, subspace_size, expansion_count, basis_vectors, &
                & operator_basis, &
                projected_matrix)

            subspace_size = subspace_size + expansion_count
            iteration_count = iteration_count + 1

            iteration_limit_reached = iteration_count > iteration_limit
            if (iteration_limit_reached) then
                solver_status = solver_status + 2048
                iteration_count = iteration_count - 1
                subspace_size = subspace_size - expansion_count
            end if
        case default
            pair_position = 1
            if (extreme_pair_count > 0) then
                ritz_values(:extreme_pair_count) = 1.0e30_real64
                convergence_flags(:extreme_pair_count) = 0
                pair_position = extreme_pair_count + 1
            end if
            first_iteration = .true.
            subspace_size = initial_vector_count
            expansion_count = subspace_size
            ! Leave at least three columns available after restarting.
            restart_size = max(subspace_limit/2,2*extreme_pair_count)
            restart_size = min(restart_size,subspace_limit - 3)
        end select

        if (.not. iteration_limit_reached) then
            ! The correction step uses neighboring Ritz values to estimate each target's spectral gap.
            requested_ritz_count = min(subspace_size,restart_size)
            lapack_abs_tolerance = slamch("S")
            call dcopy (requested_ritz_count, ritz_values, 1, previous_ritz_values, 1)
            ! DSPEVX overwrites its packed input, so solve on a copy of the projection.
            call dcopy ((subspace_size*(subspace_size + 1))/2, projected_matrix, 1, lapack_packed_copy, 1)
            call dspevx ("Vectors also", "In a range", "Upper triangular", subspace_size, lapack_packed_copy, &
                -1.0_real64, -1.0_real64, 1, requested_ritz_count, lapack_abs_tolerance, computed_ritz_count, &
                ritz_values, ritz_coefficients, subspace_size, real_scratch, integer_scratch, pending_pair_indices, &
                lapack_status)
            solver_status = -abs(lapack_status)
            if (solver_status /= 0) then
                operation_code = 0
                return
            end if
            targets_converged = check_davidson_convergence(subspace_size, extreme_pair_count, selected_pair_count, &
                & selected_pair_indices, &
                ritz_coefficients, ritz_values, convergence_flags, eigenvalue_tolerance, coefficient_tolerance, &
                real_scratch, integer_scratch, previous_ritz_values, expansion_count, pending_pair_indices)
            if (.not. targets_converged .and. subspace_size < matrix_order) then
                if (subspace_size == subspace_limit) then
                    call restart_davidson_subspace (matrix_order, restart_size, subspace_size, real_scratch, basis_vectors, &
                        operator_basis, projected_matrix, ritz_coefficients, ritz_values)
                end if
                call expand_davidson_residuals (matrix_order, extreme_pair_count, subspace_limit, expansion_block_size, &
                    & subspace_size, &
                    residual_tolerance, expansion_count, pending_pair_indices, ritz_coefficients, ritz_values, &
                    previous_ritz_values, operator_basis, basis_vectors, convergence_flags, real_scratch, &
                    real_scratch(subspace_limit+1), targets_converged)

                if (.not. targets_converged) then
                    ! Prepare the Ritz values used by the preconditioner for the original matrix.
                    do pair_position = 1, expansion_count
                        ritz_values(pending_pair_indices(pair_position)) = &
                            ritz_values(pending_pair_indices(pair_position)) &
                            - real_scratch(subspace_limit+pending_pair_indices(pair_position))
                    end do
                    if (solve_high_end) call dscal (extreme_pair_count, -1.0_real64, ritz_values, 1)

                    request_slots(1) = expansion_count
                    request_slots(2) = subspace_size + 1
                    request_slots(3) = subspace_size + 1
                    operation_code = 1
                    return
                end if
            end if
        end if

        do pair_position = 1, extreme_pair_count
            previous_ritz_values(pair_position) = abs(previous_ritz_values(pair_position) &
                - ritz_values(pair_position))
        end do

        ! Form full-space Ritz vectors and their images with the same coefficients.
        call transform_davidson_basis (matrix_order, subspace_size, extreme_pair_count, ritz_coefficients, real_scratch, &
            basis_vectors)
        call transform_davidson_basis (matrix_order, subspace_size, extreme_pair_count, ritz_coefficients, real_scratch, &
            operator_basis)
        ! Use the first column after the returned products as scratch for the final residuals.
        do pair_position = 1, extreme_pair_count
            call dcopy (matrix_order, operator_basis((pair_position-1)*matrix_order+1), 1, &
                operator_basis(extreme_pair_count*matrix_order+1), 1)
            call daxpy (matrix_order, (-ritz_values(pair_position)), &
                basis_vectors((pair_position-1)*matrix_order+1), 1, &
                operator_basis(extreme_pair_count*matrix_order+1), 1)
            real_scratch(pair_position) = ddot(matrix_order,operator_basis(extreme_pair_count*matrix_order+1),1, &
                operator_basis(extreme_pair_count*matrix_order+1),1)
            real_scratch(pair_position) = sqrt(real_scratch(pair_position))
        end do
        operation_code = 0
        return
    end subroutine iterate_davidson_solver
end module davidson_iteration
