module right_csf_intermediate_coupling_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        integer function right_csf_intermediate_coupling (loc, csf_index)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: loc
            integer, intent(in) :: csf_index
        end function right_csf_intermediate_coupling
    end interface
end module right_csf_intermediate_coupling_interface
