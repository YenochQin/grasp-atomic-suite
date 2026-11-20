module evaluate_breit_hartree_z_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine evaluate_breit_hartree_z (k, i, j)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: k
            integer :: i
            !VAST...Dummy argument I is not referenced in this routine.
            integer :: j
        end subroutine evaluate_breit_hartree_z
    end interface
end module evaluate_breit_hartree_z_interface
