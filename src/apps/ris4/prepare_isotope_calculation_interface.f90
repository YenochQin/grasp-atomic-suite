module prepare_isotope_calculation_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine prepare_isotope_calculation (name)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(in) :: name
        end subroutine prepare_isotope_calculation
    end interface
end module prepare_isotope_calculation_interface
