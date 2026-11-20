!***********************************************************************
!                                                                      *
subroutine breit_radial_integral_1 (orbital_a,orbital_b,orbital_c,orbital_d,k,radial_integral_value)
    ! Module dependencies
    !                                                                      *
    !   Returns integrals for the transverse photon interaction.           *
    !                                                                      *
    !   Written by Per Jonsson                   Octaober 2014             *
    !                                                                      *
    !***********************************************************************
    !...Translated by Charlotte Froese Fischer
    !                       Gediminas Gaigalas  10/05/17
    use, intrinsic :: iso_fortran_env, only: real64
    use breit_integral_tables, only: ndtpa, ntpi, indtp1, indtp2, indtp3, indtp4, indtp5, indtp6, &
        valtp1, valtp2, valtp3, valtp4, valtp5, valtp6, first
    use orbital_csf_basis,           only: orbital_count
    use breit_rank_offsets, only: kmax, kstartbreit1, kstartbreit2
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: orbital_a, orbital_b, orbital_c, orbital_d, k
    real(real64), intent(out) :: radial_integral_value

    ! Local declarations
    integer ::  index, loc, ju, jl, jm, key
    !-----------------------------------------------
    !
    key = orbital_count + 1
    index = ((orbital_a*key+orbital_b)*key+orbital_c)*key+orbital_d
    !
    jl = kstartbreit1(k)
    ju = kstartbreit1(k+1) - 1
    !
    if (index < indtp1(jl).or.index > indtp1(ju)) then
        write(*,*) "Something wrong in brint1"
        stop
    end if
    !
    !   The index is within the range of the indices stored; search
    !   for it in the list of indices
    !
    do while (ju-jl > 1)
        jm = (ju+jl)/2
        if (indtp1(jm) > index) then
            ju = jm
        else
            jl = jm
        end if
    end do
    !
    !   The range is bracketed to the extent possible
    !
    if (index == indtp1(ju)) then
        loc = ju
    else if (index == indtp1(jl)) then
        loc = jl
    else
        write(*,*) k,orbital_a,orbital_b,orbital_c,orbital_d,index
        write(*,*) "Brint1 Integral not found"
        stop
    end if
    !
    !   Return the value of the integral
    !   from storage

    radial_integral_value = valtp1(loc)
    !
    return
end subroutine breit_radial_integral_1
