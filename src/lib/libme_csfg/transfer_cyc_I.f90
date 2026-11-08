module transfer_cyc_I
    implicit none
    interface
        subroutine transfer_cyc(ic, ir)
            use vast_kind_param, only: double, byte
            integer :: ic, ir
        end subroutine
    end interface
end module
