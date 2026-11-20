module solve_hamiltonian_blocks_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine solve_hamiltonian_blocks (unused_core_subshell_count, max_orbital_two_j)
            implicit none (type, external)

            ! Arguments
            integer,  intent(in):: unused_core_subshell_count, max_orbital_two_j
        end subroutine solve_hamiltonian_blocks
    end interface
end module solve_hamiltonian_blocks_interface
