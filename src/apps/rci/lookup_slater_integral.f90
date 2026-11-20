!***********************************************************************
!                                                                      *
subroutine lookup_slater_integral(orbital_a, orbital_b, orbital_c, orbital_d, multipole_rank, tegral)
    ! Module dependencies
    !                                                                      *
    !                         k                                            *
    !   This routine returns R (abcd) integrals.                           *
    !                                                                      *
    !   Written by Per Jonsson                                             *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use distributed_slater_integrals, only: indteirk, valteirk, rk_local_count
    use orbital_csf_basis, only: orbital_count
    use coulomb_rank_offsets, only: kmax, kstart
    implicit none (type, external)

    ! Arguments
    integer, intent(inout) :: orbital_a, orbital_b, orbital_c, orbital_d
    integer, intent(in) :: multipole_rank
    real(real64), intent(out) :: tegral

    ! Local declarations
    integer :: key, iswap, index, jl, ju, jm, loc
    logical :: found, first
    !-----------------------------------------------
    !
    key = orbital_count + 1
    !
    !   Ensure that the indices are in `canonical' order
    !   Compute the composite (packed) index
    !


    if (orbital_a > orbital_c) then
        iswap = orbital_c
        orbital_c = orbital_a
        orbital_a = iswap
    end if
    if (orbital_b > orbital_d) then
        iswap = orbital_d
        orbital_d = orbital_b
        orbital_b = iswap
    end if
    if (orbital_a > orbital_b) then
        iswap = orbital_b
        orbital_b = orbital_a
        orbital_a = iswap
        iswap = orbital_d
        orbital_d = orbital_c
        orbital_c = iswap
    end if

    index = ((orbital_a*key + orbital_b)*key + orbital_c)*key + orbital_d
    !
    jl = kstart(multipole_rank)
    ju = kstart(multipole_rank+1) - 1

    if (index<indteirk(jl) .or. index>indteirk(ju)) then
        write (*, *) "Something wrong in rkintc"
        stop
    end if
    !
    !   The index is within the range of the indices stored; search
    !   for it in the list of indices
    !
    do while (ju - jl > 1)
        jm = (ju + jl)/2
        if (indteirk(jm) > index) then
            ju = jm
        else
            jl = jm
        end if
    end do
    !
    !   The range is bracketed to the extent possible
    !
    if (index == indteirk(ju)) then
        loc = ju
    else if (index == indteirk(jl)) then
        loc = jl
    else
        write (*, *) multipole_rank, orbital_a, orbital_b, orbital_c, orbital_d, index
        write (*, *) "Rkintc Integral not found"
        stop
    end if
    !
    !   Return the value of the integral
    !   from storage

    tegral = valteirk(loc)
    !
    return
end subroutine lookup_slater_integral
