module right_csf_two_j_plus_one_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        integer function right_csf_two_j_plus_one (csf_index)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: csf_index
        end function right_csf_two_j_plus_one
    end interface
end module right_csf_two_j_plus_one_interface
