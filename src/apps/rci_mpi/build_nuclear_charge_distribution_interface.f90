module build_nuclear_charge_distribution_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine build_nuclear_charge_distribution
            implicit none (type, external)
        end subroutine build_nuclear_charge_distribution
    end interface
end module build_nuclear_charge_distribution_interface
