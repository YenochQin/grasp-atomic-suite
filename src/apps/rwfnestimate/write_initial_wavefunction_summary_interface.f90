module write_initial_wavefunction_summary_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_initial_wavefunction_summary (nunit)
            ! Module dependencies
            use grasp_parameters, only: radial_workspace_capacity, max_orbitals
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: nunit
        end subroutine write_initial_wavefunction_summary
    end interface
end module write_initial_wavefunction_summary_interface
