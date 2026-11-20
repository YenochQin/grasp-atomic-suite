module orthogonalize_to_more_consistent_orbitals_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine orthogonalize_to_more_consistent_orbitals (orbital_index, orbital_sign_reversal_count)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: orbital_index
            integer, intent(inout) :: orbital_sign_reversal_count
        end subroutine orthogonalize_to_more_consistent_orbitals
    end interface
end module orthogonalize_to_more_consistent_orbitals_interface
