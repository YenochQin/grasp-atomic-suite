module sortoffdv_I
    implicit none
    interface
        subroutine sortoffdv(nv, vcoeff, labcdk, labp)
            use vast_kind_param, only: double, byte
            integer, intent(in) :: nv
            real(double) :: vcoeff(nv)
            integer(byte) :: labcdk(5,nv)
            integer :: labp(nv)
        end subroutine
    end interface
end module
