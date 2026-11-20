module calculate_vacuum_polarization_integral_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_vacuum_polarization_integral
contains

    !***********************************************************************
    !                                                                      *
    subroutine calculate_vacuum_polarization_integral (first_orbital,second_orbital,radial_integral_value)
        ! Module dependencies
        !                                                                      *
        !   This  routine computes nuclear vacuum polarization integrals.      *
        !                                                                      *
        !   Call(s) to: [LIB92]: ALLOC, RALLOC.                                *
        !               [RCI92]: VPINTF.                                       *
        !                                                                      *
        !                                           Last update: 09 Oct 1992   *
        !                                                                      *
        !***********************************************************************
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: orbital_index_radix
        use array_memory_management, only: alloc, ralloc
        use vacuum_polarization_integral_table, only: vacuum_polarization_integral_indices, &
            vacuum_polarization_integral_values, first_vacuum_polarization_lookup, vacuum_polarization_integral_capacity, &
            vacuum_polarization_integral_count

        ! Procedure interfaces
        use vacuum_polarization_integral_interface, only: vacuum_polarization_integral
        implicit none (type, external)

        ! Arguments
        integer, intent(inout) :: first_orbital, second_orbital
        real(real64), intent(out) :: radial_integral_value

        ! Local constants

        ! Local declarations
        logical :: found
        integer :: orbital_swap_temp, packed_integral_index, table_position, search_upper_index, search_lower_index, &
            & search_midpoint_index, new_integral_capacity, i
        !------------------------------------------------
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
        if (.not. first_vacuum_polarization_lookup) then
            !
            !   This branch is executed on all entries except the first
            !
            if (packed_integral_index > vacuum_polarization_integral_indices(vacuum_polarization_integral_count)) then
                !
                !   The index is greater than the largest stored
                !
                found = .false.
                table_position = vacuum_polarization_integral_count
            !
            else if (packed_integral_index < vacuum_polarization_integral_indices(1)) then
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
                search_upper_index = vacuum_polarization_integral_count
                search_lower_index = 1
                !
                do while (search_upper_index-search_lower_index > 1)
                    !
                    search_midpoint_index = (search_upper_index+search_lower_index)/2
                    if (vacuum_polarization_integral_indices(search_midpoint_index) > packed_integral_index) then
                        search_upper_index = search_midpoint_index
                    else
                        search_lower_index = search_midpoint_index
                    end if
                end do
                !
                !   The range is bracketed to the extent possible
                !
                if (packed_integral_index == vacuum_polarization_integral_indices(search_upper_index)) then
                    !
                    found = .true.
                    table_position = search_upper_index
                !
                else if (packed_integral_index == vacuum_polarization_integral_indices(search_lower_index)) then
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
                radial_integral_value = vacuum_polarization_integral_values(table_position)
            !
            else
                !
                !   Index not found; compute the integral
                !
                radial_integral_value = vacuum_polarization_integral (first_orbital,second_orbital)
                !
                !   Increment the integral counter
                !
                vacuum_polarization_integral_count = vacuum_polarization_integral_count+1
                !
                !   Increase array length by half the present length if the latter
                !   is inadequate to store another pair
                !
                if (vacuum_polarization_integral_count > vacuum_polarization_integral_capacity) then
                    new_integral_capacity = vacuum_polarization_integral_capacity+vacuum_polarization_integral_capacity/2
                    call ralloc (vacuum_polarization_integral_indices,new_integral_capacity,"INDOEI", "VPINT")
                    call ralloc (vacuum_polarization_integral_values,new_integral_capacity,"VALOEI", "VPINT")
                    vacuum_polarization_integral_capacity = new_integral_capacity
                end if
                !
                do i = vacuum_polarization_integral_count,table_position+2,-1
                    vacuum_polarization_integral_indices(i) = vacuum_polarization_integral_indices(i-1)
                    vacuum_polarization_integral_values(i) = vacuum_polarization_integral_values(i-1)
                end do
                !
                !   Put the new index and value into storage
                !
                vacuum_polarization_integral_indices(table_position+1) = packed_integral_index
                vacuum_polarization_integral_values(table_position+1) = radial_integral_value
                !
            end if
        !
        else
            !
            !   This branch is executed only once
            !
            first_vacuum_polarization_lookup = .false.
            !
            !   Allocate the initial storage for arrays INDOEI and VALOEI;
            !   NDVPA stores the array dimension
            !
            vacuum_polarization_integral_capacity = 2
            call alloc (vacuum_polarization_integral_indices,vacuum_polarization_integral_capacity,"INDOEI", "VPINT")
            call alloc (vacuum_polarization_integral_values,vacuum_polarization_integral_capacity,"VALOEI", "VPINT")
            !
            !   Compute the integral's value
            !
            radial_integral_value = vacuum_polarization_integral (first_orbital,second_orbital)
            !
            !   Initialise the integral counter
            !
            vacuum_polarization_integral_count = 1
            !
            !   Store the integral and its value
            !
            vacuum_polarization_integral_indices(1) = packed_integral_index
            vacuum_polarization_integral_values(1) = radial_integral_value
            !
        end if
        !
        return
    end subroutine calculate_vacuum_polarization_integral
end module calculate_vacuum_polarization_integral_interface
