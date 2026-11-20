module set_hydrogenic_screening_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine set_hydrogenic_screening (core_subshell_count)
            ! Module dependencies
            use grasp_parameters, only: max_orbitals
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: core_subshell_count
        end subroutine set_hydrogenic_screening
    end interface
end module set_hydrogenic_screening_interface
