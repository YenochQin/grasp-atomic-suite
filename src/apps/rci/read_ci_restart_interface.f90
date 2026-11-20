module read_ci_restart_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine read_ci_restart
            implicit none (type, external)
        end subroutine read_ci_restart
    end interface
end module read_ci_restart_interface
