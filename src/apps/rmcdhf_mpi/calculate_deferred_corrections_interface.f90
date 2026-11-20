module calculate_deferred_corrections_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine calculate_deferred_corrections (orbital_index)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: orbital_index
        end subroutine calculate_deferred_corrections
    end interface
end module calculate_deferred_corrections_interface
