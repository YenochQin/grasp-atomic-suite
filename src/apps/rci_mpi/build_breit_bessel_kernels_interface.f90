module build_breit_bessel_kernels_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine build_breit_bessel_kernels (orbital_a, orbital_b, ik, iw, k)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: orbital_a
            integer, intent(in) :: orbital_b
            integer, intent(in) :: ik
            integer, intent(in) :: iw
            integer, intent(in) :: k
        end subroutine build_breit_bessel_kernels
    end interface
end module build_breit_bessel_kernels_interface
