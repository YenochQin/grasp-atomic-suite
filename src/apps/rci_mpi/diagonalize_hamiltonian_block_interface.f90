module diagonalize_hamiltonian_block_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine diagonalize_hamiltonian_block (iatjpo, iaspar, nelmnt_a)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: int64
            implicit none (type, external)

            ! Arguments
            integer, intent(out) :: iatjpo
            integer, intent(out) :: iaspar
            integer(int64)        :: nelmnt_a
        end subroutine diagonalize_hamiltonian_block
    end interface
end module diagonalize_hamiltonian_block_interface
