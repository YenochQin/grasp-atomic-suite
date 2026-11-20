!***********************************************************************
!                                                                      *
subroutine breit_radial_integral_4 (orbital_a,orbital_b,orbital_c,orbital_d,radial_integral_order,radial_integral_value)
    ! Module dependencies
    !                                                                      *
    !   Returns integrals for the transverse photon interaction.           *
    !      Integrals are stored in ordered lists. If the integral cannot   *
    !   be read from a list, it is computed by calling BRINTF.             *
    !                                                                      *
    !   Observ that it is not possible to use more than 100 orbitals when  *
    !   the Breit interaction is included. If 100 orbitals are used then   *
    !   NU <= 19.                                                          *
    !                                                                      *
    !   Call(s) to: [LIB92]: ALLOC, RALLOC.                                *
    !               [RCI92]: BRINTF                                        *
    !                                                                      *
    !   Written by Farid A Parpia               Last update: 09 Oct 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use array_memory_management, only: alloc, ralloc
    use breit_integral_tables, only: ndtpa, ntpi, indtp1, indtp2, indtp3, indtp4, indtp5, indtp6, &
        valtp1, valtp2, valtp3, valtp4, valtp5, valtp6, first
    use orbital_csf_basis,           only: csf_count, orbital_count, csf_subshell_occupations

    ! Procedure interfaces
    use lookup_breit_radial_integral_interface, only: lookup_breit_radial_integral
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: orbital_a, orbital_b, orbital_c, orbital_d, radial_integral_order
    real(real64), intent(out) :: radial_integral_value

    ! Local declarations
    logical :: found
    integer :: key, index, loc, ju, jl, jm, newsiz, i
    !-----------------------------------------------
    !
    key = orbital_count + 1
    !
    !   Compute the integral label
    !
    index = (((radial_integral_order*key+orbital_a)*key+orbital_b)*key+orbital_c)*key+orbital_d
    !
    if (.not. first(4)) then
        !
        !   This branch is executed on all entries except the first
        !
        if (index > indtp4(ntpi(4))) then
            !
            !   The index is greater than the largest stored
            !
            found = .false.
            loc = ntpi(4)
        !
        else if (index < indtp4(1)) then
            !
            !   The index is less than the smallest stored
            !
            found = .false.
            loc = 0
        !
        else
            !
            !   The index is within the range of the indices stored; search
            !   for it in the list of indices
            !
            ju = ntpi(4)
            jl = 1
            do while (ju-jl > 1)
                jm = (ju+jl)/2
                if (indtp4(jm) > index) then
                    ju = jm
                else
                    jl = jm
                end if
            end do
            !
            !   The range is bracketed to the extent possible
            !
            if (index == indtp4(ju)) then
                !
                found = .true.
                loc = ju
            !
            else if (index == indtp4(jl)) then
                !
                found = .true.
                loc = jl
            !
            else
                !
                found = .false.
                loc = jl
                !
            end if
            !
            !
        end if
        !
        if (found) then
            !
            !   Found the index in the list; return the value of the integral
            !   from storage
            !
            radial_integral_value = valtp4(loc)
        !
        else
            !
            !   Index not found; compute the integral
            !
            radial_integral_value = lookup_breit_radial_integral (4,orbital_a,orbital_b,orbital_c,orbital_d,radial_integral_order)
            !
            !   Increment the integral counter
            !
            ntpi(4) = ntpi(4)+1
            !
            !   Increase array length by half the present length if the latter
            !   is inadequate to store another pair
            !
            !
            if (ntpi(4) > ndtpa(4)) then
                newsiz = ndtpa(4)+ndtpa(4)/2
                call ralloc (indtp4,newsiz,"INDTP4", "BRINT4")
                call ralloc (valtp4,newsiz,"VALTP4", "BRINT4")
                ndtpa(4)= newsiz
            end if
            do i = ntpi(4),loc+2,-1
                indtp4(i) = indtp4(i-1)
                valtp4(i) = valtp4(i-1)
            end do
            !
            !   Put the new index and value into storage
            !
            indtp4(loc+1) = index
            valtp4(loc+1) = radial_integral_value
            !
        end if
    !
    else
        !
        !   This branch is executed only once per type of integral
        !
        first(4) = .false.
        !
        !   Designate the initial storage for arrays INDTPx and VALTPx;
        !   Array NDTPA stores the array dimensions
        !
        ndtpa(4) = 10000
        !
        !   Compute the integral's value
        !
        radial_integral_value = lookup_breit_radial_integral (4,orbital_a,orbital_b,orbital_c,orbital_d,radial_integral_order)
        !
        !   Initialise the integral counter
        !
        ntpi(4) = 1
        !
        !   Store the integral and its value
        !
        call alloc (indtp4,ndtpa(4),"INDTP4", "BRINT4")
        call alloc (valtp4,ndtpa(4),"VALTP4", "BRINT4")
        indtp4(1) = index
        valtp4(1) = radial_integral_value
        !
    end if
    !
    return
end subroutine breit_radial_integral_4
