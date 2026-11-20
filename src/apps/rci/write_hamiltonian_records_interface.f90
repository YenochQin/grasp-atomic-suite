module write_hamiltonian_records_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_hamiltonian_records (restart_matrix_unit)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: restart_matrix_unit
        end subroutine write_hamiltonian_records
    end interface
end module write_hamiltonian_records_interface
