module read_scf_options_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine read_scf_options (extended_optimal_level, block_labels, isotope_file_path, radial_file_path, &
            & nuclear_polarization_file_path)
            implicit none (type, external)

            ! Arguments
            logical, intent(out) :: extended_optimal_level
            character (len = 8) :: block_labels(*)
            character (len = *) :: isotope_file_path
            character (len = *) :: radial_file_path
            character (len = *), intent(in) :: nuclear_polarization_file_path
        end subroutine read_scf_options
    end interface
end module read_scf_options_interface
