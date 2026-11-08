module twobody_dc_I
    implicit none
    interface
        subroutine twobody_dc(tegral)
            use vast_kind_param, only: double, byte
            real(double) :: tegral
        end subroutine
    end interface
end module
