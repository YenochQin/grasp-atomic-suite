module diagonalize_hamiltonian_block_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine diagonalize_hamiltonian_block (dvdfirst,print_details, symmetry_block_index, ncfpat, ncminpat, nevecpat, &
            & total_csf_count)
            implicit none (type, external)

            ! Arguments
            logical, intent(in) ::  dvdfirst
            logical :: print_details
            !VAST...Dummy argument LPRINT is not referenced in this routine.
            integer, intent(in) :: symmetry_block_index
            integer, intent(in) :: ncfpat
            integer, intent(in) :: ncminpat
            integer, intent(in) :: nevecpat
            integer :: total_csf_count
        end subroutine diagonalize_hamiltonian_block
    end interface
end module diagonalize_hamiltonian_block_interface
