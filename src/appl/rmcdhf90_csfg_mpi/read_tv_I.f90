module read_tv_I
    implicit none
    interface
        subroutine read_tv(iflag, icw, irw, iunitf, mtype, iaa, ibb, tshell1)
            use vast_kind_param, only: double, byte
            integer, intent(in) :: iflag, icw, irw, iunitf, mtype
            integer :: iaa, ibb
            real(double) :: tshell1
        end subroutine
    end interface
end module
