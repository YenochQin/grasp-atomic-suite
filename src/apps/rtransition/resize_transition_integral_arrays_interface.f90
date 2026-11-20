module resize_transition_integral_arrays_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine resize_transition_integral_arrays(isldr, isldr1, xsldr, ntdim, imode)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, pointer :: isldr(:), isldr1(:)
            real(real64), pointer :: xsldr(:)
            integer, intent(inout) :: ntdim
            integer, intent(in) :: imode
            !...This routine performs I/O.
        end subroutine resize_transition_integral_arrays
    end interface
end module resize_transition_integral_arrays_interface
