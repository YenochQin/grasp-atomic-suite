module read_calculation_options_interface
    ! Module dependencies
    use csf_descriptor_reader, only: descriptor_source
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine read_calculation_options(source, status)
            ! Module dependencies
            import descriptor_source
            implicit none (type, external)

            ! Arguments
            type(descriptor_source), intent(in) :: source
            integer, intent(out) :: status
        end subroutine read_calculation_options
    end interface
end module read_calculation_options_interface
