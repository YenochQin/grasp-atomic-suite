module apply_lagrange_constraints_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine apply_lagrange_constraints (orbital_index, mpi_process_count)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: orbital_index
            integer, intent(in) :: mpi_process_count
        end subroutine apply_lagrange_constraints
    end interface
end module apply_lagrange_constraints_interface
