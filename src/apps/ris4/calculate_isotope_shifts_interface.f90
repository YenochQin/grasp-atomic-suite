module calculate_isotope_shifts_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Gediminas Gaigalas 11/18/19
        subroutine calculate_isotope_shifts(name)
            implicit none (type, external)

            ! Arguments
            character(len=24), intent(in) :: name
        end subroutine calculate_isotope_shifts
    end interface
end module calculate_isotope_shifts_interface
