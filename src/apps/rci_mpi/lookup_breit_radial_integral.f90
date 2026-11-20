!***********************************************************************
!                                                                      *
real(real64) function lookup_breit_radial_integral (itype, orbital_a, orbital_b, orbital_c, orbital_d, k)
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
    use breit_bessel_cache_keys, only: keep

    ! Procedure interfaces
    use build_breit_bessel_kernels_interface, only: build_breit_bessel_kernels
    use breit_radial_integral_interface, only: breit_radial_integral

    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: itype
    integer, intent(in) :: orbital_a
    integer, intent(in) :: orbital_b
    integer, intent(in) :: orbital_c
    integer, intent(in) :: orbital_d
    integer, intent(in) :: k

    ! Local constants
    integer, parameter :: orbital_label_radix = orbital_index_radix
    integer, parameter :: key2 = orbital_label_radix*orbital_label_radix

    ! Local declarations
    integer :: icod, icod1, icod2, icod3, icod4
    logical :: kernels_cached
    !-----------------------------------------------
    !
    select case (itype)
    case (2, 6)
        !        Types 2 and 6 use Bessel functions at k-1 and k+1.
        kernels_cached = .false.
        if (itype == 2) then
            if (orbital_a/=orbital_c .or. orbital_b/=orbital_d .or. orbital_a/=orbital_d .or. orbital_c/=orbital_b) then
                !
                icod = max(orbital_c,orbital_d) + orbital_label_radix*min(orbital_c,orbital_d)
                icod1 = icod + key2*(k - 1)
                icod2 = icod + key2*(k + 1)
                icod = max(orbital_a,orbital_b) + orbital_label_radix*min(orbital_a,orbital_b)
                icod3 = icod + key2*(k - 1)
                icod4 = icod + key2*(k + 1)
                kernels_cached = (icod1==keep(1,2) .and. icod2==keep(2,2) .and. icod3==keep(1,1)&
                    .and. icod4==keep(2,1)) .or. &
                    (icod1==keep(1,1) .and. icod2==keep(2,1) .and. icod3==keep(1,2)&
                    .and. icod4==keep(2,2))
                if (.not. kernels_cached) then
                    call build_breit_bessel_kernels (orbital_c, orbital_d, 1, 2, k - 1)
                    call build_breit_bessel_kernels (orbital_c, orbital_d, 2, 2, k + 1)
                end if
            end if
            !
        end if
        if (.not. kernels_cached) then
            call build_breit_bessel_kernels (orbital_a, orbital_b, 1, 1, k - 1)
            call build_breit_bessel_kernels (orbital_a, orbital_b, 2, 1, k + 1)
        end if
    case (3)
        !        Type 3 uses the Bessel functions for each nonzero frequency.
        if (orbital_a /= orbital_b) call build_breit_bessel_kernels (orbital_a, orbital_b, 1, 1, k)
        if (orbital_c /= orbital_d) call build_breit_bessel_kernels (orbital_c, orbital_d, 1, 2, k)
    case (4)
    !        Type 4 needs no Bessel kernels.
    case default
        !        Types 1 and 5 use the Bessel functions at k. Other types
        !        follow the type 1 path, as in the original dispatch.
        if (itype /= 5) then
            if (orbital_a/=orbital_c .or. orbital_b/=orbital_d .or. orbital_a/=orbital_d .or. orbital_c==orbital_b) call &
                & build_breit_bessel_kernels ( &
                & orbital_c, orbital_d, 1, 2&
                , k)
        end if
        call build_breit_bessel_kernels (orbital_a, orbital_b, 1, 1, k)
    end select
    !
    !   Compute the integral
    !
    lookup_breit_radial_integral = breit_radial_integral(itype,orbital_a,orbital_b,orbital_c,orbital_d,k)
    !
    return
end function lookup_breit_radial_integral
