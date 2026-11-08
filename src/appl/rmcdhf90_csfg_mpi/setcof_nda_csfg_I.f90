module setcof_nda_csfg_I
    implicit none
    interface
        subroutine setcof_nda_csfg(iroff, icoff, ia, ib, tcoeff)
            use vast_kind_param, only: double, byte
            integer, intent(in) :: iroff, icoff, ia, ib
            real(double), intent(in) :: tcoeff
        end subroutine
    end interface
end module
