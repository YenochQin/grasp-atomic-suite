module diagonalize_hamiltonian_block_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine diagonalize_hamiltonian_block (iatjpo, iaspar)
            implicit none (type, external)

            ! Arguments
            integer, intent(out) :: iatjpo
            integer, intent(out) :: iaspar
        end subroutine diagonalize_hamiltonian_block
    end interface
end module diagonalize_hamiltonian_block_interface
