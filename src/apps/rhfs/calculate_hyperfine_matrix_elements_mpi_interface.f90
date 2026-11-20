module calculate_hyperfine_matrix_elements_mpi_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine calculate_hyperfine_matrix_elements_mpi
            implicit none (type, external)
        end subroutine calculate_hyperfine_matrix_elements_mpi
    end interface
end module calculate_hyperfine_matrix_elements_mpi_interface
