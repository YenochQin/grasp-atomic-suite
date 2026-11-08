module twobody_dc_sms_I
    implicit none
    interface
        subroutine twobody_dc_sms(atwinv, vcoeff, emttmp)
            use vast_kind_param, only: double, byte
            real(double) :: atwinv, vcoeff, emttmp
        end subroutine
    end interface
end module
