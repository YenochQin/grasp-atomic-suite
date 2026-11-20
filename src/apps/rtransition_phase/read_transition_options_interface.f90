module read_transition_options_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine read_transition_options (isotope_file_path,name)
            implicit none (type, external)

            ! Arguments
            character(len=*):: isotope_file_path
            character (len = 24) :: name(2)
        end subroutine read_transition_options
    end interface
end module read_transition_options_interface
