module solve_hamiltonian_blocks_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine solve_hamiltonian_blocks (dvdfirst)
            implicit none (type, external)

            ! Arguments
            logical, intent(in) ::  dvdfirst
        end subroutine solve_hamiltonian_blocks
    end interface
end module solve_hamiltonian_blocks_interface
