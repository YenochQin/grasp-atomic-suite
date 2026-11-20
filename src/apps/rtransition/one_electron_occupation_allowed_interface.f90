module one_electron_occupation_allowed_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        integer function one_electron_occupation_allowed (ja, jb)
            implicit none (type, external)

            ! Arguments
            integer :: ja
            integer :: jb
        end function one_electron_occupation_allowed
    end interface
end module one_electron_occupation_allowed_interface
