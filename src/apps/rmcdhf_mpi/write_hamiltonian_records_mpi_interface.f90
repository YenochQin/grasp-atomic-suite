module write_hamiltonian_records_mpi_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_hamiltonian_records_mpi (mpi_rank, mpi_process_count, block_csf_count)
            implicit none (type, external)

            ! Arguments
            integer :: mpi_rank
            !VAST...Dummy argument MYID is not referenced in this routine.
            integer :: mpi_process_count
            !VAST...Dummy argument NPROCS is not referenced in this routine.
            integer, intent(in) :: block_csf_count
        end subroutine write_hamiltonian_records_mpi
    end interface
end module write_hamiltonian_records_mpi_interface
