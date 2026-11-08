module csfg_expand_findtype_I
    implicit none
    interface
        subroutine csfg_expand_findtype(name, file_rl, file_out, iflag)
            use vast_kind_param, only: double, byte
            character(*), intent(in) :: name, file_rl, file_out
            integer, intent(in) :: iflag
        end subroutine
    end interface
end module
