module connect_ci_restart_file_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine connect_ci_restart_file (isotope_file_path, radial_file_path, block_labels)
            implicit none (type, external)

            ! Arguments
            character (len = *) :: isotope_file_path
            character (len = *) :: radial_file_path
            character (len = 8) :: block_labels(*)
        end subroutine connect_ci_restart_file
    end interface
end module connect_ci_restart_file_interface
