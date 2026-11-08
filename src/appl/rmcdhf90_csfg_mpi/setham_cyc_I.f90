module setham_cyc_I
    implicit none
    interface
        subroutine setham_cyc(icc, irr, ic, ir)
            use vast_kind_param, only: double, byte
            integer, intent(in) :: icc, irr, ic, ir
        end subroutine
    end interface
end module
