module calculate_oscillator_strengths_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine calculate_oscillator_strengths (name)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(inout) :: name(2)
        end subroutine calculate_oscillator_strengths
    end interface
end module calculate_oscillator_strengths_interface
