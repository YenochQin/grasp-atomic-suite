module matrixblock2_I
    implicit none
    interface
        subroutine matrixblock2(ic, ir, ncoec, incor, nctec, inc2, nmcbp, ncore, elsto)
            use vast_kind_param, only: double, byte
            integer :: ic, ir, ncoec, incor, nctec, inc2, nmcbp, ncore
            real(double) :: elsto
        end subroutine
    end interface
end module
