module lookup_breit_radial_integral_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: lookup_breit_radial_integral
contains

    !***********************************************************************
    !                                                                      *
    real(real64) function lookup_breit_radial_integral (breit_integral_type, orbital_a, orbital_b, orbital_c, orbital_d, &
        & multipole_rank)
        ! Module dependencies
        !                                                                      *
        !   Computes integrals for the transverse photon interaction.          *
        !                                                                      *
        !   Call(s) to: [RCI92]: BESSEL, BRRA.                                 *
        !                                                                      *
        !                                           Last update: 09 Oct 1992   *
        !                                                                      *
        !***********************************************************************
        use grasp_parameters,   only: orbital_index_radix
        use breit_bessel_cache_keys, only: breit_bessel_kernel_keys

        ! Procedure interfaces
        use build_breit_bessel_kernels_interface, only: build_breit_bessel_kernels
        use breit_radial_integral_interface, only: breit_radial_integral

        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: breit_integral_type
        integer, intent(in) :: orbital_a
        integer, intent(in) :: orbital_b
        integer, intent(in) :: orbital_c
        integer, intent(in) :: orbital_d
        integer, intent(in) :: multipole_rank

        ! Local constants
        integer, parameter :: orbital_label_radix_squared = orbital_index_radix*orbital_index_radix

        ! Local declarations
        integer :: packed_orbital_pair, second_pair_lower_rank_key, second_pair_upper_rank_key, first_pair_lower_rank_key, &
            & first_pair_upper_rank_key
        logical :: kernels_cached
        !-----------------------------------------------
        !
        select case (breit_integral_type)
        case (2, 6)
            !        Types 2 and 6 use Bessel functions at k-1 and k+1.
            kernels_cached = .false.
            if (breit_integral_type == 2) then
                if (orbital_a/=orbital_c .or. orbital_b/=orbital_d .or. orbital_a/=orbital_d .or. orbital_c/=orbital_b) then
                    !
                    packed_orbital_pair = max(orbital_c,orbital_d) + orbital_index_radix*min(orbital_c,orbital_d)
                    second_pair_lower_rank_key = packed_orbital_pair + orbital_label_radix_squared*(multipole_rank - 1)
                    second_pair_upper_rank_key = packed_orbital_pair + orbital_label_radix_squared*(multipole_rank + 1)
                    packed_orbital_pair = max(orbital_a,orbital_b) + orbital_index_radix*min(orbital_a,orbital_b)
                    first_pair_lower_rank_key = packed_orbital_pair + orbital_label_radix_squared*(multipole_rank - 1)
                    first_pair_upper_rank_key = packed_orbital_pair + orbital_label_radix_squared*(multipole_rank + 1)
                    kernels_cached = (second_pair_lower_rank_key==breit_bessel_kernel_keys(1,2) .and. &
                        second_pair_upper_rank_key==breit_bessel_kernel_keys(2,2) .and. &
                        & first_pair_lower_rank_key==breit_bessel_kernel_keys(1,1)&
                        .and. first_pair_upper_rank_key==breit_bessel_kernel_keys(2,1)) .or. &
                        (second_pair_lower_rank_key==breit_bessel_kernel_keys(1,1) .and. &
                            second_pair_upper_rank_key==breit_bessel_kernel_keys(2,1) .and. &
                            & first_pair_lower_rank_key==breit_bessel_kernel_keys(1,2)&
                        .and. first_pair_upper_rank_key==breit_bessel_kernel_keys(2,2))
                    if (.not. kernels_cached) then
                        call build_breit_bessel_kernels (orbital_c, orbital_d, 1, 2, multipole_rank - 1)
                        call build_breit_bessel_kernels (orbital_c, orbital_d, 2, 2, multipole_rank + 1)
                    end if
                end if
                !
            end if
            if (.not. kernels_cached) then
                call build_breit_bessel_kernels (orbital_a, orbital_b, 1, 1, multipole_rank - 1)
                call build_breit_bessel_kernels (orbital_a, orbital_b, 2, 1, multipole_rank + 1)
            end if
        case (3)
            !        Type 3 uses the Bessel functions for each nonzero frequency.
            if (orbital_a /= orbital_b) call build_breit_bessel_kernels (orbital_a, orbital_b, 1, 1, multipole_rank)
            if (orbital_c /= orbital_d) call build_breit_bessel_kernels (orbital_c, orbital_d, 1, 2, multipole_rank)
        case (4)
        !        Type 4 needs no Bessel kernels.
        case default
            !        Types 1 and 5 use the Bessel functions at k. Other types
            !        follow the type 1 path, as in the original dispatch.
            if (breit_integral_type /= 5) then
                if (orbital_a/=orbital_c .or. orbital_b/=orbital_d .or. orbital_a/=orbital_d .or. orbital_c==orbital_b) call &
                    & build_breit_bessel_kernels ( &
                    & orbital_c, orbital_d, 1, 2&
                    , multipole_rank)
            end if
            call build_breit_bessel_kernels (orbital_a, orbital_b, 1, 1, multipole_rank)
        end select
        !
        !   Compute the integral
        !
        lookup_breit_radial_integral = &
            & breit_radial_integral(breit_integral_type,orbital_a,orbital_b,orbital_c,orbital_d,multipole_rank)
        !
        return
    end function lookup_breit_radial_integral
end module lookup_breit_radial_integral_interface
