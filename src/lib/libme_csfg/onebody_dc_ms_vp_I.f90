module onebody_dc_ms_vp_I
    implicit none
    interface
        subroutine onebody_dc_ms_vp(ia, ib, atwinv, tcoeff, emttmp)
            use vast_kind_param, only: double, byte
            integer :: ia, ib
            real(double) :: atwinv, tcoeff, emttmp
        end subroutine
    end interface
end module
