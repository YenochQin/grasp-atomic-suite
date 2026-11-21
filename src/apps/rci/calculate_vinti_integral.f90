module calculate_vinti_integral_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_vinti_integral
contains

    !***********************************************************************
    !                                                                      *
    subroutine calculate_vinti_integral (first_orbital,second_orbital,radial_integral_value)
        ! Module dependencies
        !                                                                      *
        !   This routine returns Vinti integrals.                              *
        !      A list of previously computed integrals is maintained. If the   *
        !   required  integral is not already  available it is computed  and   *
        !   stored in the list.                                                *
        !                                                                      *
        !   Call(s) to: [LIB92]: ALLOC, RALLOC, VINTI.                         *
        !                                                                      *
        !   Written by Farid A Parpia             Last revision: 06 Jun 1992   *
        !                                                                      *
        !***********************************************************************
        !   Converted by CFF
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: orbital_index_radix
        use array_memory_management, only: alloc, ralloc
        use vinti_integral_table, only: first_vinti_lookup, vinti_integral_indices, vinti_integral_capacity, &
            vinti_integral_count, vinti_integral_values

        ! Procedure interfaces
        use vinti_radial_integral_interface, only: vinti_radial_integral
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: first_orbital, second_orbital
        real(real64), intent(out) :: radial_integral_value

        ! Local constants

        ! Local declarations
        logical :: found
        integer :: packed_integral_index, table_position, search_upper_index, search_lower_index, search_midpoint_index, &
            & new_integral_capacity, i
        !------------------------------------------------
        !
        !   Compute the composite (packed) index
        !
        packed_integral_index = first_orbital*orbital_index_radix+second_orbital
        !
        if (.not. first_vinti_lookup) then
            !
            !   This branch is executed on all entries except the first
            !
            if (packed_integral_index > vinti_integral_indices(vinti_integral_count)) then
                !
                !   The index is greater than the largest stored
                !
                found = .false.
                table_position = vinti_integral_count
            !
            else if (packed_integral_index < vinti_integral_indices(1)) then
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
                search_upper_index = vinti_integral_count
                search_lower_index = 1
                !
                do while (search_upper_index-search_lower_index > 1)
                    !
                    search_midpoint_index = (search_upper_index+search_lower_index)/2
                    if (vinti_integral_indices(search_midpoint_index) > packed_integral_index) then
                        search_upper_index = search_midpoint_index
                    else
                        search_lower_index = search_midpoint_index
                    end if
                end do
                !
                !   The range is bracketed to the extent possible
                !
                if (packed_integral_index == vinti_integral_indices(search_upper_index)) then
                    !
                    found = .true.
                    table_position = search_upper_index
                !
                else if (packed_integral_index == vinti_integral_indices(search_lower_index)) then
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
                radial_integral_value = vinti_integral_values(table_position)
            !
            else
                !
                !   Index not found; compute the integral
                !
                radial_integral_value = vinti_radial_integral (first_orbital,second_orbital)
                !
                !   Increment the integral counter
                !
                vinti_integral_count = vinti_integral_count+1
                !
                !   Increase array length by half the present length if the latter
                !   is inadequate to store another pair
                !
                if (vinti_integral_count > vinti_integral_capacity) then
                    new_integral_capacity = vinti_integral_capacity+vinti_integral_capacity/2
                    call ralloc (vinti_integral_indices,new_integral_capacity, "vinti_integral_indices", "calculate_vinti_integral")
                    call ralloc (vinti_integral_values,new_integral_capacity, "vinti_integral_values", "calculate_vinti_integral")
                    vinti_integral_capacity = new_integral_capacity
                end if
                !
                do i = vinti_integral_count,table_position+2,-1
                    vinti_integral_indices(i) = vinti_integral_indices(i-1)
                    vinti_integral_values(i) = vinti_integral_values(i-1)
                end do
                !
                !   Put the new index and value into storage
                !
                vinti_integral_indices(table_position+1) = packed_integral_index
                vinti_integral_values(table_position+1) = radial_integral_value
                !
            end if
        !
        else
            !
            !   This branch is executed only once
            !
            first_vinti_lookup = .false.
            !
            !   Allocate the initial storage for arrays INDTEI and VALTEI;
            !   NDVIN stores the array dimension
            !
            vinti_integral_capacity = 2
            call alloc (vinti_integral_indices,vinti_integral_capacity, "vinti_integral_indices", "calculate_vinti_integral")
            call alloc (vinti_integral_values,vinti_integral_capacity, "vinti_integral_values", "calculate_vinti_integral")
            !
            !   Compute the integral's value
            !
            radial_integral_value = vinti_radial_integral (first_orbital,second_orbital)
            !
            !   Initialise the integral counter
            !
            vinti_integral_count = 1
            !
            !   Store the integral and its value
            !
            vinti_integral_indices(1) = packed_integral_index
            vinti_integral_values(1) = radial_integral_value
            !
        end if
        !
        return
    end subroutine calculate_vinti_integral
end module calculate_vinti_integral_interface
