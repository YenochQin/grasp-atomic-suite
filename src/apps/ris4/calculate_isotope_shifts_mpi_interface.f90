module calculate_isotope_shifts_mpi_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine calculate_isotope_shifts_mpi (name)
            implicit none (type, external)

            ! Arguments
            character(len=24), intent(in) :: name
        end subroutine calculate_isotope_shifts_mpi
    end interface
end module calculate_isotope_shifts_mpi_interface
