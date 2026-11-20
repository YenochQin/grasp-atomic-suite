module calculate_subshell_overlap_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine calculate_subshell_overlap (s, nshli, nshlf)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: nshli
            integer, intent(in) :: nshlf
            real(real64), intent(out) :: s(nshli, nshlf)
        end subroutine calculate_subshell_overlap
    end interface
end module calculate_subshell_overlap_interface
