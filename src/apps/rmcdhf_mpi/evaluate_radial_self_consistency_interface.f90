module evaluate_radial_self_consistency_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine evaluate_radial_self_consistency (orbital_index)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: orbital_index
        end subroutine evaluate_radial_self_consistency
    end interface
end module evaluate_radial_self_consistency_interface
