module resize_transition_coefficient_arrays_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine resize_transition_coefficient_arrays(iptr, isldr, isldr1, xsldr, nmdim, imode)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, pointer :: iptr(:), isldr(:), isldr1(:)
            real(real64), pointer :: xsldr(:)
            integer, intent(inout) :: nmdim
            integer, intent(in) :: imode
            !...This routine performs I/O.
        end subroutine resize_transition_coefficient_arrays
    end interface
end module resize_transition_coefficient_arrays_interface
