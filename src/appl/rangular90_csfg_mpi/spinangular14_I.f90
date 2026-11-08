module spinangular14_I
    implicit none
    interface
        subroutine spinangular14(jasav, jbsav, iunitf, npos)
            use vast_kind_param, only: double, byte
            integer :: jasav, jbsav, iunitf, npos
        end subroutine
    end interface
end module
