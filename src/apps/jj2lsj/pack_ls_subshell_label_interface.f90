module pack_ls_subshell_label_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine pack_ls_subshell_label (m, el, q, couple, str)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: m
            character (len = 3), intent(in) :: el(*)
            integer, intent(in) :: q(*)
            character (len = 3), intent(in) :: couple(*)
            character (len = 164), intent(out) :: str
        end subroutine pack_ls_subshell_label
    end interface
end module pack_ls_subshell_label_interface
