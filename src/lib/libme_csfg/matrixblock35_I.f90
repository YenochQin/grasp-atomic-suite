module matrixblock35_I
    implicit none
    interface
        subroutine matrixblock35(icsav, irsav, ncoec, incor, nctec, inc2, nmcbp, ncore, elsto)
            use vast_kind_param, only: double, byte
            integer :: icsav, irsav, ncoec, incor, nctec, inc2, nmcbp, ncore
            real(double) :: elsto
        end subroutine
    end interface
end module
