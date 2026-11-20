module assemble_hamiltonian_block_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine assemble_hamiltonian_block (symmetry_block_index, mpi_rank, mpi_process_count, total_contribution_count, &
            & total_secondary_contribution_count)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: int64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: symmetry_block_index
            integer, intent(in) :: mpi_rank
            integer, intent(in) :: mpi_process_count
            integer(int64) :: total_contribution_count(20), total_secondary_contribution_count(20)
        end subroutine assemble_hamiltonian_block
    end interface
end module assemble_hamiltonian_block_interface
