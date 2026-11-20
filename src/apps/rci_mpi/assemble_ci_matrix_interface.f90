module assemble_ci_matrix_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine assemble_ci_matrix (atwinv, symmetry_block_index, mpi_rank, mpi_process_count, elsto, irestart, slf_en)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            real(real64) :: atwinv
            integer, intent(in) :: symmetry_block_index
            integer, intent(in) :: mpi_rank
            integer, intent(in) :: mpi_process_count
            real(real64), intent(out) :: elsto
            integer, intent(out) :: irestart
            real(real64) :: slf_en(*)
        end subroutine assemble_ci_matrix
    end interface
end module assemble_ci_matrix_interface
