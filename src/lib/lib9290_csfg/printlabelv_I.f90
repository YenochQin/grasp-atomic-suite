module printlabelv_I
    implicit none
    interface
        subroutine printlabelv(iflag, icw, irw, ntxx, mtype, ia, ib, tshell1)
            use vast_kind_param, only: double, byte
            integer :: iflag, icw, irw, ntxx, mtype, ia, ib
            real(double) :: tshell1
        end subroutine
    end interface
end module
