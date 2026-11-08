module spinangular2_I
    implicit none
    interface
        subroutine spinangular2(jasav,jbsav,iunitf,npos,jcsfasav,jcsfbsav)
            use vast_kind_param, only: double, byte
            integer :: jasav,jbsav,iunitf,npos,jcsfasav,jcsfbsav
        end subroutine
    end interface
end module
