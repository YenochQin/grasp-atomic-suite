module transfer_csfg_I
    implicit none
    interface
        subroutine transfer_csfg(ic, ir)
            use vast_kind_param, only: double, byte
            integer :: ic, ir
        end subroutine
    end interface
end module
