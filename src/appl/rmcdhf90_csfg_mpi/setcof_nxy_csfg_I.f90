module setcof_nxy_csfg_I
    implicit none
    interface
        subroutine setcof_nxy_csfg(iroff, icoff, vcoeff)
            use vast_kind_param, only: double, byte
            integer, intent(in) :: iroff, icoff
            real(double), intent(in) :: vcoeff
        end subroutine
    end interface
end module
