module davidson_eigensolver
    implicit none (type, external)

    ! Module interface
    private
    public :: gdvd
    public :: dvd_operator

    ! Procedure interfaces
    ! Apply the original matrix to a block of column vectors.
    abstract interface
        subroutine dvd_operator(matrix_order, vector_count, input_vectors, output_vectors)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: matrix_order, vector_count
            real(real64), intent(in) :: input_vectors(matrix_order,vector_count)
            real(real64), intent(inout) :: output_vectors(matrix_order,vector_count)
        end subroutine dvd_operator
    end interface
contains

    subroutine gdvd(matrix_vector_callback, matrix_order, subspace_limit, matrix_diagonal, lowest_pair_index, &
        highest_pair_index, selected_pair_indices, initial_vector_count, expansion_block_size, &
        eigenvalue_tolerance, coefficient_tolerance, residual_tolerance, orthogonality_tolerance, &
        iteration_limit, real_workspace, real_workspace_size, integer_workspace, integer_workspace_size, &
        solve_high_end, iteration_count, matrix_vector_product_count, solver_status)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use, intrinsic :: iso_fortran_env, only: real64

        ! Service the solver requests with a matrix callback or diagonal preconditioning.
        use davidson_solver_control, only: dvdson
        implicit none (type, external)

        ! Arguments
        procedure(dvd_operator) :: matrix_vector_callback
        integer, intent(in) :: matrix_order
        integer, intent(in) :: subspace_limit
        integer  :: lowest_pair_index
        integer  :: highest_pair_index
        integer, intent(inout) :: initial_vector_count
        integer  :: expansion_block_size
        integer  :: iteration_limit
        integer, intent(in) :: real_workspace_size
        integer, intent(in) :: integer_workspace_size
        integer  :: iteration_count
        integer, intent(out) :: matrix_vector_product_count
        integer, intent(out) :: solver_status
        real(real64)  :: eigenvalue_tolerance
        real(real64)  :: coefficient_tolerance
        real(real64)  :: residual_tolerance
        ! Retain this interface parameter. Orthonormalization uses two fixed passes.
        real(real64)  :: orthogonality_tolerance
        logical  :: solve_high_end
        integer  :: selected_pair_indices(subspace_limit)
        integer, intent(inout) :: integer_workspace(integer_workspace_size)
        real(real64), intent(inout) :: matrix_diagonal(matrix_order)
        ! The workspace starts with the input guesses. On return its leading blocks contain
        ! Ritz vectors, eigenvalues, eigenvalue changes and residual norms.
        real(real64), intent(inout) :: real_workspace(real_workspace_size)

        ! Local declarations
        integer :: request_slots(7)
        integer :: constraint_count, operation_code, vector_count, input_start, output_start, ritz_values_start, &
            pending_indices_start, residual_norms_start, preconditioner_shifts_start, vector_offset, vector_index, &
            pair_offset, component_index
        real(real64) :: preconditioning_value, stored_shift, residual_norm, preconditioner_shift

        ! The public wrapper starts with an unconstrained, caller-supplied basis.
        constraint_count = 0
        operation_code = 0

        if (initial_vector_count == 0) then
            write (output_unit, *) "GDVD Error : No initial estimate!!!"
            solver_status = -1000
            return
        end if

        matrix_vector_product_count = 0
        do
            call dvdson (operation_code, request_slots, matrix_order, subspace_limit, constraint_count, &
                lowest_pair_index, highest_pair_index, selected_pair_indices, initial_vector_count, &
                expansion_block_size, eigenvalue_tolerance, coefficient_tolerance, residual_tolerance, &
                iteration_limit, real_workspace, real_workspace_size, integer_workspace, integer_workspace_size, &
                solve_high_end, iteration_count, solver_status)
            vector_count = request_slots(1)
            input_start = request_slots(2)
            output_start = request_slots(3)
            ritz_values_start = request_slots(4)
            pending_indices_start = request_slots(5)
            residual_norms_start = request_slots(6)
            preconditioner_shifts_start = request_slots(7)

            ! Precondition the corrected residuals with the shifted diagonal of the original matrix.
            if (operation_code == 1) then
                vector_offset = 0
                do vector_index = 1, vector_count
                    pair_offset = integer_workspace(pending_indices_start+vector_index-1) - 1
                    preconditioning_value = real_workspace(ritz_values_start+pair_offset)
                    stored_shift = real_workspace(preconditioner_shifts_start+pair_offset)
                    residual_norm = real_workspace(residual_norms_start+pair_offset)
                    preconditioner_shift = real_workspace(residual_norms_start+subspace_limit+pair_offset)
                    matrix_diagonal(:matrix_order) = matrix_diagonal(:matrix_order) - preconditioning_value
                    do component_index = 1, matrix_order
                        if (abs(matrix_diagonal(component_index)) > 1.0e-6_real64) then
                            real_workspace(output_start + vector_offset + component_index - 1) = &
                                real_workspace(input_start + vector_offset + component_index - 1) &
                                / matrix_diagonal(component_index)
                        else
                            real_workspace(output_start + vector_offset + component_index - 1) = &
                                real_workspace(input_start + vector_offset + component_index - 1)*1.0e6_real64
                        end if
                    end do
                    matrix_diagonal(:matrix_order) = matrix_diagonal(:matrix_order) + preconditioning_value

                    vector_offset = vector_offset + matrix_order
                end do

            ! Initialization and iteration share the same matrix-vector callback.
            else if (operation_code==2 .or. operation_code==3) then
                call matrix_vector_callback (matrix_order, vector_count, real_workspace(input_start), &
                    real_workspace(output_start))
                matrix_vector_product_count = matrix_vector_product_count + vector_count
            else
                exit
            end if
        end do
    end subroutine gdvd
end module davidson_eigensolver
