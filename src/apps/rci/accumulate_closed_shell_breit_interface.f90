module accumulate_closed_shell_breit_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine accumulate_closed_shell_breit (bra_csf_index, ket_csf_index, ja1, ipca, jb1)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: bra_csf_index
            integer, intent(in) :: ket_csf_index
            integer, intent(in) :: ja1
            integer, intent(in) :: ipca
            integer, intent(in) :: jb1
        end subroutine accumulate_closed_shell_breit
    end interface
end module accumulate_closed_shell_breit_interface
