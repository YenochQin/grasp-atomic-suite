module build_vacuum_polarization_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine build_vacuum_polarization
            implicit none (type, external)
        end subroutine build_vacuum_polarization
    end interface
end module build_vacuum_polarization_interface
