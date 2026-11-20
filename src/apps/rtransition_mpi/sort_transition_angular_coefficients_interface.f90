module sort_transition_angular_coefficients_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine sort_transition_angular_coefficients (name, coefficient_scratch_unit, transition_angular_unit, print_details, &
            & transition_type_index, iblki, &
            & iblkf)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(in) :: name(2)
            integer, intent(in) :: coefficient_scratch_unit
            integer, intent(inout) :: transition_angular_unit
            logical, intent(in) :: print_details
            integer :: transition_type_index
            integer, intent(in) :: iblki
            integer, intent(in) :: iblkf
        end subroutine sort_transition_angular_coefficients
    end interface
end module sort_transition_angular_coefficients_interface
