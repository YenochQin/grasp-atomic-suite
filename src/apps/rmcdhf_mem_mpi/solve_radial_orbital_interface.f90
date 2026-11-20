module solve_radial_orbital_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine solve_radial_orbital (orbital_index, solve_failed, orbital_sign_reversal_count, matching_grid_index, nnp)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: orbital_index
            logical, intent(out) :: solve_failed
            integer, intent(inout) :: orbital_sign_reversal_count
            integer, intent(out) :: matching_grid_index
            integer, intent(out) :: nnp
        end subroutine solve_radial_orbital
    end interface
end module solve_radial_orbital_interface
