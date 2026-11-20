module improve_orbital_mpi_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine improve_orbital_mpi (extended_optimal_level, orbital_index, lsort, dampmx)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            logical, intent(in) :: extended_optimal_level
            integer, intent(in) :: orbital_index
            logical :: lsort
            real(real64), intent(inout) :: dampmx
        end subroutine improve_orbital_mpi
    end interface
end module improve_orbital_mpi_interface
