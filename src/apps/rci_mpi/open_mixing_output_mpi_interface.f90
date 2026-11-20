module open_mixing_output_mpi_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine open_mixing_output_mpi (name, idblk)
            implicit none (type, external)

            ! Arguments
            character (len = *), intent(in) :: name
            character (len = 8) :: idblk(*)
        end subroutine open_mixing_output_mpi
    end interface
end module open_mixing_output_mpi_interface
