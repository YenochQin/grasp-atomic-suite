module assemble_hamiltonian_block_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine assemble_hamiltonian_block (mpi_rank, mpi_process_count, symmetry_block_index, elsto,icstrt, nelmntt,  &
            atwinv,slf_en)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64, int64
            implicit none (type, external)

            ! Arguments
            integer(int64)       :: nelmntt
            integer             :: symmetry_block_index, icstrt
            integer, intent(in) :: mpi_rank, mpi_process_count
            real(real64) :: elsto, atwinv
            real(real64) :: slf_en(*)
        end subroutine assemble_hamiltonian_block
    end interface
end module assemble_hamiltonian_block_interface
