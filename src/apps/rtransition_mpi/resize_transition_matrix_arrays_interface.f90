module resize_transition_matrix_arrays_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine resize_transition_matrix_arrays(jja, jjb, hb1, hb2, hc1, hc2, hm1, &
            hm2, lab, nptr, nsdim, imode)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer  :: nsdim
            integer, intent(in) :: imode
            integer, pointer :: jja(:), jjb(:), lab(:), nptr(:)
            real(real64), pointer :: hb1(:), hb2(:), hc1(:), hc2(:), hm1(:), hm2(:)
        end subroutine resize_transition_matrix_arrays
    end interface
end module resize_transition_matrix_arrays_interface
