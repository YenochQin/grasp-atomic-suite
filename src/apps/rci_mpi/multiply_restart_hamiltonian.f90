! Apply the Hamiltonian represented by the CI restart records to trial vectors.
module multiply_restart_hamiltonian_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: multiply_restart_hamiltonian
contains

    ! Matrix-matrix product: C = AB.  A  sparse  representation of the
    ! lower triangle of the  (NxN)  matrix  A  is read from the disk
    !
    ! This is an adaptation of  Andreas Stathopulos'  routine  SPSBMV,
    ! and is specific to GRASP2 derivatives.
    !
    ! F A Parpia and A Stathopoulos. Last revision: 13 Oct 1992
    ! Block Version by Xinghong He. Last revision: 18 Jun 1998
    subroutine multiply_restart_hamiltonian(n, m, b, c)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use eigenstate_expansions, only: hamiltonian_energy_shift
        use ci_restart_units, only: restart_matrix_unit
        use rci_matrix_file_position, only: nposition
        use mpi_state,           only:  mpi_rank, mpi_process_count
        use real_vector_initialization, only: initialize_real_vector
        use file_positioning, only: position_calculation_file
        use sparse_dot_and_scatter, only: sparse_dot_scatter
        use mpi_utilities, only: allreduce_mpi_real_sum
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: n
        integer, intent(in) :: m
        real(real64), intent(in) :: b(n, m)
        real(real64), intent(inout) :: c(n, m)

        ! Local declarations
        integer :: stored_row_indices(n)
        integer ::               csf_count, file_csf_count, file_first_order_cutoff, file_mpi_rank, file_mpi_process_count, &
            & stored_column_index&
            , stored_entry_count, stored_entry_index, rhs_vector_index
        real(real64) :: stored_column_values(n)
        real(real64) :: stored_column_energy, output_with_diagonal, off_diagonal_dot_product
        ! !...nposition+1 is the current position of the .res file
        ! !...It is set in matrix and used in maneig, spodmv
        write (output_unit, *) "Calling multiply_restart_hamiltonian..."
        csf_count = n

        ! Initialise the result matrix; note that this is specific to the
        ! data structure of DVDSON

        call initialize_real_vector (n*m, 0.0_real64, c, 1)

        ! ...moved from maneig before "CALL GDVD (SPODMV..."
        call position_calculation_file (0, restart_matrix_unit, nposition)

        read (restart_matrix_unit) file_csf_count, file_first_order_cutoff, file_mpi_rank, file_mpi_process_count
        if (csf_count/=file_csf_count .or. mpi_rank/=file_mpi_rank .or. file_mpi_process_count/=mpi_process_count) stop &
            "multiply_restart_hamiltonian: restart file CSF count mismatch"

        do stored_column_index = mpi_rank + 1, n, mpi_process_count
            read (restart_matrix_unit) stored_entry_count, stored_column_energy, (stored_column_values(stored_entry_index), &
                & stored_entry_index=1,stored_entry_count), (stored_row_indices(stored_entry_index),stored_entry_index=1, &
                & stored_entry_count)
            do rhs_vector_index = 1, m
                output_with_diagonal = c(stored_column_index,rhs_vector_index) + (stored_column_values(stored_entry_count)- &
                    & hamiltonian_energy_shift)*b(stored_column_index,rhs_vector_index)
                call sparse_dot_scatter (stored_entry_count - 1, b(1:n,rhs_vector_index), c(1:n,rhs_vector_index), &
                    & stored_row_indices(1:n), &
                    stored_column_values(1:n), b(stored_column_index,rhs_vector_index), off_diagonal_dot_product)
                c(stored_column_index,rhs_vector_index) = output_with_diagonal + off_diagonal_dot_product
            end do
        end do

        call allreduce_mpi_real_sum (c, n*m)

        return
    end subroutine multiply_restart_hamiltonian
end module multiply_restart_hamiltonian_interface
