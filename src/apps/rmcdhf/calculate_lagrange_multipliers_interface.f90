module calculate_lagrange_multipliers_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine calculate_lagrange_multipliers (extended_optimal_level)
            implicit none (type, external)

            ! Arguments
            logical :: extended_optimal_level
            !...This routine performs I/O.
        end subroutine calculate_lagrange_multipliers
    end interface
end module calculate_lagrange_multipliers_interface
