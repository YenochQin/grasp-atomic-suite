module orthogonalize_same_kappa_orbitals_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine orthogonalize_same_kappa_orbitals (nw, jp, lsort)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: nw
            integer, intent(in) :: jp
            logical, intent(in) :: lsort
        end subroutine orthogonalize_same_kappa_orbitals
    end interface
end module orthogonalize_same_kappa_orbitals_interface
