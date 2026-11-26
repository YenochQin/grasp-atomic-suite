! Apply the stored sparse symmetric Hamiltonian to trial vectors.
module multiply_sparse_hamiltonian_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: multiply_sparse_hamiltonian
contains

    !   Matrix-matrix product: C = AB.  The lower triangle of the  (NxN)
    !   matrix is assumed available in packed form in the array EMT. The
    !   matrices B and C are (NxM).
    !
    !   This is an adaptation of  Andreas Stathopulos'  routine  TRSBMV,
    !   and is specific to GRASP2 derivatives.
    !
    ! F A Parpia and A Stathopoulos. Last revision: 09 Oct 1992
    ! Block version by Xinghong He. Last revision: 18 Jun 1998
    subroutine multiply_sparse_hamiltonian(n, m, b, c)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use sparse_hamiltonian_storage,          only: hamiltonian_values
        use mpi_state, only: mpi_process_count, mpi_rank
        use real_vector_initialization, only: initialize_real_vector

        ! Procedure interfaces
        use reduce_sparse_matrix_product_interface, only: reduce_sparse_matrix_product

        ! Module dependencies
        use mpi_f08, only: mpi_count_kind
        use mpi_element_counts, only: checked_mpi_count
        use mpi_utilities, only: allreduce_mpi_real_sum
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: n
        integer, intent(in) :: m
        real(real64), intent(in) :: b(n, m)
        real(real64), intent(inout) :: c(n, m)

        ! Local declarations
        integer :: value_count
        integer      :: ibeg, iend, icol, nelc, iv
        real(real64) :: diag, dl

        ! Initialise the result matrix; note that this is specific to the
        ! data structure of DVDSON --- there is no overdimensioning

        value_count = checked_mpi_count(size(c,kind=mpi_count_kind), "multiply_sparse_hamiltonian")
        call initialize_real_vector (value_count, 0.0_real64, c, 1)

        ibeg = 1
        iend = 0
        do icol = mpi_rank + 1, n, mpi_process_count
            iend = iend + icol
            nelc = iend - ibeg + 1
            do iv = 1, m
                diag = c(icol,iv) + hamiltonian_values(iend)*b(icol,iv)
                call reduce_sparse_matrix_product (nelc - 1, b(1:n,iv), c(1:n,iv), &
                    hamiltonian_values(ibeg:iend), b(icol,iv), dl)
                c(icol,iv) = diag + dl
            end do
            ibeg = iend + 1
        end do

        call allreduce_mpi_real_sum (c,value_count)

        return
    end subroutine multiply_sparse_hamiltonian
end module multiply_sparse_hamiltonian_interface
