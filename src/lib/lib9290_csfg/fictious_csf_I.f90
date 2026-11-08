module fictious_csf_I
    implicit none
    interface
        subroutine fictious_csf(intype, icsfn, icsfo, iorbn1, iorbo1, iorbn2, iorbo2)
            use vast_kind_param, only: double, byte
            integer :: intype, icsfn, icsfo, iorbn1, iorbo1, iorbn2, iorbo2
        end subroutine
    end interface
end module
