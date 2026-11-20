module basis_csf_two_j_plus_one_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        integer function basis_csf_two_j_plus_one (csf_index)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: csf_index
            !VAST.../ORB2/ NCF(IN)
            !VAST.../IOUNIT/ ISTDE(IN)
            !VAST...Calls: JCUPA
            !...This routine performs I/O.
        end function basis_csf_two_j_plus_one
    end interface
end module basis_csf_two_j_plus_one_interface
