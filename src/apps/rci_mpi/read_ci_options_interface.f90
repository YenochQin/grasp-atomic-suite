module read_ci_options_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine read_ci_options (isotope_file_path, radial_file_path, block_labels)
            implicit none (type, external)

            ! Arguments
            character(len=*):: isotope_file_path, radial_file_path
            character(len=8) ::  block_labels(*)
        end subroutine read_ci_options
    end interface
end module read_ci_options_interface
