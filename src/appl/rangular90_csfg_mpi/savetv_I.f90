module savetv_I
    implicit none
    interface
        subroutine savetv(iflag, icw, irw, iunitf, mtype, iaa, ibb, tshell1)
            use vast_kind_param, only: double, byte
            integer :: iflag, icw, irw, iunitf, mtype, iaa, ibb
            real(double) :: tshell1
        end subroutine
    end interface
end module
