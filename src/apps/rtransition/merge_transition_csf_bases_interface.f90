module merge_transition_csf_bases_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine merge_transition_csf_bases (state_name)
            implicit none (type, external)

            ! Arguments
            character (len = 24) :: state_name(2)
        end subroutine merge_transition_csf_bases
    end interface
end module merge_transition_csf_bases_interface
