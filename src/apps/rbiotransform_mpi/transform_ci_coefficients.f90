! Counter-transform configuration mixing coefficients after the orbital transformation.
module transform_ci_coefficients_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: transform_ci_coefficients
contains

    ! Calculate the action of the  operator
    !
    !  T(L) = T(nl) T(n-1 L ) .... T(1 l)
    !
    ! where
    !
    !  T(i l ) = ( sum(k=0,2l+1) t(i l ) ** k ) (t_ilil ** E(il,il))
    !  t(i l ) = sum(j.lt.i) t(j,i) E(lj,li)
    !  E(lj,li) total symmetric shell excitation operator
    !
    ! Modification of the CITRA routine for MCHF
    !
    ! Input
    !
    ! CIIN      : Input CI vectors ( Destroyed in the process )
    ! NCSCF     : Length of CI expansion
    ! NCIV      : Number of CI vectors
    ! L         : L value of shells in excitations
    ! NSHL      : Number of shells with this L
    ! T         : List of coefficients
    ! NIN       : Number of inactive shells
    !
    ! Output
    !
    ! CIOUT : List of output CI vectors
    subroutine transform_ci_coefficients(input_mixing_coefficients,csf_count,level_count,symmetry_index,subshell_count, &
                & transformation_matrix,inactive_subshell_count,global_trace_level, &
        & transformed_mixing_coefficients,transfer_power_coefficients)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use biorthogonal_subshell_mapping, only: active_orbitals_by_group
        use orbital_csf_basis, only: orbital_kappa
        use mpi_f08, only: mpi_count_kind
        use mpi_element_counts, only: checked_mpi_count
        use mpi_state, only: mpi_allreduce, mpi_comm_world, mpi_double_precision, mpi_error_code, mpi_sum

        ! Procedure interfaces
        use write_matrix_interface, only: write_matrix
        use scale_vector_interface, only: scale_vector
        use apply_occupation_factor_interface, only: apply_occupation_factor
        use apply_subshell_transfer_interface, only: apply_subshell_transfer
        use combine_vectors_interface, only: combine_vectors
        use copy_vector_interface, only: copy_vector
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: csf_count
        integer, intent(in) :: level_count
        integer, intent(in) :: symmetry_index
        integer, intent(in) :: subshell_count
        integer, intent(in) :: inactive_subshell_count
        integer, intent(in) :: global_trace_level
        real(real64)  :: input_mixing_coefficients(csf_count,level_count)
        real(real64), intent(in) :: transformation_matrix(subshell_count,subshell_count)
        real(real64)  :: transformed_mixing_coefficients(csf_count,level_count)
        real(real64)  :: transfer_power_coefficients(csf_count,level_count)

        ! Local declarations
        integer :: coefficient_payload_count
        integer :: local_trace_level, trace_level, inactive_subshell, subshell_capacity, subshell_index, transfer_order
        real(real64) :: inactive_shell_factor, diagonal_element, inverse_factorial
        real(real64)  :: scrtmp(csf_count,level_count)

        coefficient_payload_count = checked_mpi_count(size(input_mixing_coefficients,kind=mpi_count_kind), &
            "transform_ci_coefficients")
        local_trace_level = 0
        trace_level = max(global_trace_level,local_trace_level)
        trace_level = 0

        if (trace_level >= 10) then
            write (output_unit, *)
            write (output_unit, *) " ***************"
            write (output_unit, *) " Entering CITRAG"
            write (output_unit, *) " ***************"
            write (output_unit, *)
        end if
        if (trace_level >= 100) then
            write (output_unit, *) " Input CI vectors "
            call write_matrix (input_mixing_coefficients, csf_count, level_count, csf_count, level_count)
            write (output_unit, *) " Transformation matrix T"
            call write_matrix (transformation_matrix, subshell_count, subshell_count, subshell_count, subshell_count)
        end if
        ! . Factor from inactive shells
        if (inactive_subshell_count /= 0) then
            inactive_shell_factor = 1.0_real64
            do inactive_subshell = 1, inactive_subshell_count
                inactive_shell_factor = inactive_shell_factor*transformation_matrix(inactive_subshell,inactive_subshell)
            end do
            ! IPOT = 2*(2*L+1)  (number of m_lm_s. This should be replaced
            !                    by (2j+1) corresponding to L)
            subshell_capacity = 2*iabs(orbital_kappa(active_orbitals_by_group(symmetry_index,inactive_subshell)))
            inactive_shell_factor = inactive_shell_factor**subshell_capacity
            call scale_vector (input_mixing_coefficients, inactive_shell_factor, coefficient_payload_count)
        end if
        if (inactive_subshell_count == subshell_count) call copy_vector (input_mixing_coefficients, &
                & transformed_mixing_coefficients, &
            & coefficient_payload_count)
        do subshell_index = inactive_subshell_count + 1, subshell_count
            if (trace_level >= 100) write (output_unit, *) " Loop I,L = ", subshell_index, symmetry_index
            ! . The diagonal contribution
            diagonal_element = transformation_matrix(subshell_index,subshell_index)
            call apply_occupation_factor (input_mixing_coefficients, csf_count, level_count, subshell_index, symmetry_index, &
                & diagonal_element, &
                & transformed_mixing_coefficients, global_trace_level)

            call copy_vector (transformed_mixing_coefficients, transfer_power_coefficients, coefficient_payload_count)
            ! .  Off diagonal contributions
            inverse_factorial = 1.0_real64
            do transfer_order = 1, 2*iabs(orbital_kappa(active_orbitals_by_group(symmetry_index,subshell_index)))
                if (trace_level >= 100) write (output_unit, *) " Loop N = ", transfer_order
                ! T ** (N-1) is supposed to be in SCR, copy to CIIN
                ! and apply S
                call copy_vector (transfer_power_coefficients, input_mixing_coefficients, coefficient_payload_count)
                call apply_subshell_transfer (input_mixing_coefficients,csf_count,level_count,subshell_index,symmetry_index, &
                & transformation_matrix(1,subshell_index),subshell_count, &
                    & scrtmp,global_trace_level)
                call mpi_allreduce(scrtmp,transfer_power_coefficients,coefficient_payload_count,         &
                    mpi_double_precision, mpi_sum, mpi_comm_world, mpi_error_code)
                inverse_factorial = inverse_factorial/float(transfer_order)
                call combine_vectors (transformed_mixing_coefficients, transformed_mixing_coefficients, &
                & transfer_power_coefficients, 1.0_real64, inverse_factorial, &
                    & coefficient_payload_count)
            end do
            call copy_vector (transformed_mixing_coefficients, input_mixing_coefficients, coefficient_payload_count)
        end do
        if (trace_level >= 100) then
            write (output_unit, *) " Output CI vectors L = ", symmetry_index
            call write_matrix (transformed_mixing_coefficients, csf_count, level_count, csf_count, level_count)
        end if
        return
    end subroutine transform_ci_coefficients
end module transform_ci_coefficients_interface
