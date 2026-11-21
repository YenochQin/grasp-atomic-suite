module lookup_slater_integral_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: lookup_slater_integral
contains

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
        use distributed_slater_integrals, only: local_slater_integral_indices, local_slater_integral_values
        use orbital_csf_basis, only: orbital_count
        use coulomb_rank_offsets, only: coulomb_integral_offsets_by_rank
        implicit none (type, external)

        ! Arguments
        integer, intent(inout) :: orbital_a, orbital_b, orbital_c, orbital_d
        integer, intent(inout) :: multipole_rank
        real(real64), intent(out) :: tegral

        ! Local declarations
        integer :: packing_radix, orbital_swap_temp, packed_integral_index, search_lower_index, search_upper_index, &
            & search_midpoint_index, table_position
        !-----------------------------------------------
        !
        packing_radix = orbital_count + 1
        !
        !   Ensure that the indices are in `canonical' order
        !   Compute the composite (packed) index
        !


        if (orbital_a > orbital_c) then
            orbital_swap_temp = orbital_c
            orbital_c = orbital_a
            orbital_a = orbital_swap_temp
        end if
        if (orbital_b > orbital_d) then
            orbital_swap_temp = orbital_d
            orbital_d = orbital_b
            orbital_b = orbital_swap_temp
        end if
        if (orbital_a > orbital_b) then
            orbital_swap_temp = orbital_b
            orbital_b = orbital_a
            orbital_a = orbital_swap_temp
            orbital_swap_temp = orbital_d
            orbital_d = orbital_c
            orbital_c = orbital_swap_temp
        end if

        packed_integral_index = ((orbital_a*packing_radix + orbital_b)*packing_radix + orbital_c)*packing_radix + orbital_d
        !
        search_lower_index = coulomb_integral_offsets_by_rank(multipole_rank)
        search_upper_index = coulomb_integral_offsets_by_rank(multipole_rank+1) - 1
        if (packed_integral_index<local_slater_integral_indices(search_lower_index) .or. &
            packed_integral_index>local_slater_integral_indices(search_upper_index)) then
            write (*, *) "Something wrong in lookup_slater_integral"
            stop
        end if
        !
        !   The index is within the range of the indices stored; search
        !   for it in the list of indices
        !
        do while (search_upper_index - search_lower_index > 1)
            search_midpoint_index = (search_upper_index + search_lower_index)/2
            if (local_slater_integral_indices(search_midpoint_index) > packed_integral_index) then
                search_upper_index = search_midpoint_index
            else
                search_lower_index = search_midpoint_index
            end if
        end do
        !
        !   The range is bracketed to the extent possible
        !
        if (packed_integral_index == local_slater_integral_indices(search_upper_index)) then
            table_position = search_upper_index
        else if (packed_integral_index == local_slater_integral_indices(search_lower_index)) then
            table_position = search_lower_index
        else
            write (*, *) multipole_rank, orbital_a, orbital_b, orbital_c, orbital_d, packed_integral_index
            write (*, *) "lookup_slater_integral Integral not found"
            stop
        end if
        !
        !   Return the value of the integral
        !   from storage

        tegral = local_slater_integral_values(table_position)
        !
        return
    end subroutine lookup_slater_integral
end module lookup_slater_integral_interface
