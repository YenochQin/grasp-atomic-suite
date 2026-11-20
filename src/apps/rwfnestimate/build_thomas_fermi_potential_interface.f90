module build_thomas_fermi_potential_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine build_thomas_fermi_potential
            implicit none (type, external)
        end subroutine build_thomas_fermi_potential
    end interface
end module build_thomas_fermi_potential_interface
