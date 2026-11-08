module setarrs_csfg_I
    implicit none
    interface
        subroutine setarrs_csfg(jblock)
            use vast_kind_param, only: double, byte
            integer :: jblock
        end subroutine
    end interface
end module
