module findtype_I
    implicit none
    interface
        subroutine findtype(ncfsmall, ncfgen, ncftot)
            use vast_kind_param, only: double, byte
            integer :: ncfsmall, ncfgen, ncftot
        end subroutine
    end interface
end module
