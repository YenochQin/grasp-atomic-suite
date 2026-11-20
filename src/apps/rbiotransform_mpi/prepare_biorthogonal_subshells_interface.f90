module prepare_biorthogonal_subshells_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine prepare_biorthogonal_subshells (ntestg, ncore1, ncore2)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: ntestg
            integer, intent(in) :: ncore1
            integer, intent(in) :: ncore2
        end subroutine prepare_biorthogonal_subshells
    end interface
end module prepare_biorthogonal_subshells_interface
