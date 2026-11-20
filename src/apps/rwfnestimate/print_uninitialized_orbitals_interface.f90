module print_uninitialized_orbitals_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine print_uninitialized_orbitals (all)
            ! Module dependencies
            use grasp_parameters, only: max_orbitals
            implicit none (type, external)

            ! Arguments
            logical, intent(out) :: all
        end subroutine print_uninitialized_orbitals
    end interface
end module print_uninitialized_orbitals_interface
