module one_electron_integral_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: one_electron_integral
contains

    !***********************************************************************
    !                                                                      *
    subroutine one_electron_integral (first_orbital,second_orbital,radial_integral_value)
        ! Module dependencies
        !                                                                      *
        !   This routine returns I(ab) integrals.                              *
        !      A list of previously computed integrals is maintained. If the   *
        !   required  integral is not already  available it is computed  and   *
        !   stored in the list.                                                *
        !      Note that indices IA and IB may be interchanged on output.      *
        !                                                                      *
        !   Call(s) to: [LIB92]: ALLOC, RALLOC, RINTI                          *
        !                                                                      *
        !   Written by Farid A Parpia             Last revision: 06 Oct 1992   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: orbital_index_radix
        use array_memory_management, only: alloc, ralloc
        use coulomb_integral_table, only: first_coulomb_integral_lookup, coulomb_integral_indices, coulomb_integral_count, &
            coulomb_integral_capacity, coulomb_integral_values
        use one_electron_radial_integral, only: integrate_one_electron_radial
        implicit none (type, external)

        ! Arguments
        integer, intent(inout) :: first_orbital, second_orbital
        real(real64), intent(out) :: radial_integral_value

        ! Local declarations
        logical ::  found
        integer :: i, packed_integral_index, table_position, search_upper_index, search_lower_index, search_midpoint_index, &
            new_integral_capacity, orbital_swap_temp
        !-----------------------------------------------
        !
        !   Ensure that the indices are in `canonical' order
        !
        if (first_orbital > second_orbital) then
            orbital_swap_temp = second_orbital
            second_orbital = first_orbital
            first_orbital = orbital_swap_temp
        end if
        !
        !   Compute the composite (packed) index
        !
        packed_integral_index = first_orbital*orbital_index_radix+second_orbital
        !
        if (.not. first_coulomb_integral_lookup) then
            !
            !   This branch is executed on all entries except the first
            !
            if (packed_integral_index > coulomb_integral_indices(coulomb_integral_count)) then
                !
                !   The index is greater than the largest stored
                !
                found = .false.
                table_position = coulomb_integral_count
            !
            else if (packed_integral_index < coulomb_integral_indices(1)) then
                !
                !   The index is less than the smallest stored
                !
                found = .false.
                table_position = 0
            !
            else
                !
                !   The index is within the range of the indices stored; search
                !   for it in the list of indices
                !
                search_upper_index = coulomb_integral_count
                search_lower_index = 1
                !
                do while (search_upper_index-search_lower_index > 1)
                    !
                    search_midpoint_index = (search_upper_index+search_lower_index)/2
                    if (coulomb_integral_indices(search_midpoint_index) > packed_integral_index) then
                        search_upper_index = search_midpoint_index
                    else
                        search_lower_index = search_midpoint_index
                    end if
                end do
                !
                !   The range is bracketed to the extent possible
                !
                if (packed_integral_index == coulomb_integral_indices(search_upper_index)) then
                    !
                    found = .true.
                    table_position = search_upper_index
                !
                else if (packed_integral_index == coulomb_integral_indices(search_lower_index)) then
                    !
                    found = .true.
                    table_position = search_lower_index
                !
                else
                    !
                    found = .false.
                    table_position = search_lower_index
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
                radial_integral_value = coulomb_integral_values(table_position)
            !
            else
                !
                !   Index not found; compute the integral
                !
                radial_integral_value = integrate_one_electron_radial (first_orbital,second_orbital,0)
                !
                !   Increment the integral counter
                !
                coulomb_integral_count = coulomb_integral_count+1
                !
                !   Increase array length by half the present length if the latter
                !   is inadequate to store another pair
                !
                if (coulomb_integral_count > coulomb_integral_capacity) then
                    new_integral_capacity = coulomb_integral_capacity+coulomb_integral_capacity/2
                    call ralloc (coulomb_integral_indices,new_integral_capacity,"INDOEI","IABINT" )
                    call ralloc (coulomb_integral_values,new_integral_capacity,"VALOEI","IABINT")
                    coulomb_integral_capacity = new_integral_capacity
                end if
                !
                do i = coulomb_integral_count,table_position+2,-1
                    coulomb_integral_indices(i) = coulomb_integral_indices(i-1)
                    coulomb_integral_values(i) = coulomb_integral_values(i-1)
                end do
                !
                !   Put the new index and value into storage
                !
                coulomb_integral_indices(table_position+1) = packed_integral_index
                coulomb_integral_values(table_position+1) = radial_integral_value
                !
            end if
        !
        else
            !
            !   This branch is executed only once
            !
            first_coulomb_integral_lookup = .false.
            !
            !   Allocate the initial storage for arrays INDOEI and VALOEI;
            !   NDCOEA stores the array dimension
            !
            coulomb_integral_capacity = 2
            call alloc (coulomb_integral_indices,coulomb_integral_capacity,"INDOEI","IABINT" )
            call alloc (coulomb_integral_values,coulomb_integral_capacity,"VALOEI","IABINT")
            !
            !   Compute the integral's value
            !
            radial_integral_value = integrate_one_electron_radial (first_orbital,second_orbital,0)
            !
            !   Initialise the integral counter
            !
            coulomb_integral_count = 1
            !
            !   Store the integral and its value
            !
            coulomb_integral_indices(1) = packed_integral_index
            coulomb_integral_values(1) = radial_integral_value
            !
        end if
        !
        return
    end subroutine one_electron_integral
end module one_electron_integral_interface
