module exchange_data_I
    implicit none
    interface
        subroutine exchange_data(nv, i, j, vcoeff, labcdk, labp, isort)
            use vast_kind_param, only: double, byte
            integer, intent(in) :: nv, i, j
            real(double) :: vcoeff(nv)
            integer(byte) :: labcdk(5,nv)
            integer :: labp(nv), isort(nv)
        end subroutine
    end interface
end module
