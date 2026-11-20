module build_orbital_potential_mpi_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine build_orbital_potential_mpi (eol, j, npts)
            implicit none (type, external)

            ! Arguments
            logical :: eol
            integer :: j
            integer :: npts
        end subroutine build_orbital_potential_mpi
    end interface
end module build_orbital_potential_mpi_interface
