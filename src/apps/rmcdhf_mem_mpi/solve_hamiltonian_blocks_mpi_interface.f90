module solve_hamiltonian_blocks_mpi_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine solve_hamiltonian_blocks_mpi (dvdfirst)
            implicit none (type, external)

            ! Arguments
            logical, intent(in) ::  dvdfirst
        end subroutine solve_hamiltonian_blocks_mpi
    end interface
end module solve_hamiltonian_blocks_mpi_interface
