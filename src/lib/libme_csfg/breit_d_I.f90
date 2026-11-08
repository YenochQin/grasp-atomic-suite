module breit_d_I
    implicit none
    interface
        subroutine breit_d(ic, ir, iv, vcoeff, ncore, elsto, emttmp)
            use vast_kind_param, only: double, byte
            integer :: ic, ir, iv, ncore
            real(double) :: vcoeff, elsto, emttmp
        end subroutine
    end interface
end module
