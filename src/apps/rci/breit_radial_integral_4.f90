module breit_radial_integral_4_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: breit_radial_integral_4
contains

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
        use breit_integral_tables, only: breit_integral_capacities, breit_integral_counts, breit_integral_indices_4, &
            breit_integral_values_4, first_breit_integral_lookup
        use orbital_csf_basis, only: orbital_count

        ! Procedure interfaces
        use lookup_breit_radial_integral_interface, only: lookup_breit_radial_integral
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_a, orbital_b, orbital_c, orbital_d, radial_integral_order
        real(real64), intent(out) :: radial_integral_value

        ! Local declarations
        logical :: found
        integer :: packing_radix, packed_integral_index, table_position, search_upper_index, search_lower_index, &
            & search_midpoint_index, new_integral_capacity, i
        !-----------------------------------------------
        !
        packing_radix = orbital_count + 1
        !
        !   Compute the integral label
        !
        packed_integral_index = &
            & (((radial_integral_order*packing_radix+orbital_a)*packing_radix+orbital_b) &
            & *packing_radix+orbital_c)*packing_radix+orbital_d
        !
        if (.not. first_breit_integral_lookup(4)) then
            !
            !   This branch is executed on all entries except the first
            !
            if (packed_integral_index > breit_integral_indices_4(breit_integral_counts(4))) then
                !
                !   The index is greater than the largest stored
                !
                found = .false.
                table_position = breit_integral_counts(4)
            !
            else if (packed_integral_index < breit_integral_indices_4(1)) then
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
                search_upper_index = breit_integral_counts(4)
                search_lower_index = 1
                do while (search_upper_index-search_lower_index > 1)
                    search_midpoint_index = (search_upper_index+search_lower_index)/2
                    if (breit_integral_indices_4(search_midpoint_index) > packed_integral_index) then
                        search_upper_index = search_midpoint_index
                    else
                        search_lower_index = search_midpoint_index
                    end if
                end do
                !
                !   The range is bracketed to the extent possible
                !
                if (packed_integral_index == breit_integral_indices_4(search_upper_index)) then
                    !
                    found = .true.
                    table_position = search_upper_index
                !
                else if (packed_integral_index == breit_integral_indices_4(search_lower_index)) then
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
                radial_integral_value = breit_integral_values_4(table_position)
            !
            else
                !
                !   Index not found; compute the integral
                !
                radial_integral_value = lookup_breit_radial_integral (4,orbital_a,orbital_b,orbital_c,orbital_d, &
                    & radial_integral_order)
                !
                !   Increment the integral counter
                !
                breit_integral_counts(4) = breit_integral_counts(4)+1
                !
                !   Increase array length by half the present length if the latter
                !   is inadequate to store another pair
                !
                !
                if (breit_integral_counts(4) > breit_integral_capacities(4)) then
                    new_integral_capacity = breit_integral_capacities(4)+breit_integral_capacities(4)/2
                    call ralloc (breit_integral_indices_4,new_integral_capacity,"INDTP4", "BRINT4")
                    call ralloc (breit_integral_values_4,new_integral_capacity,"VALTP4", "BRINT4")
                    breit_integral_capacities(4)= new_integral_capacity
                end if
                do i = breit_integral_counts(4),table_position+2,-1
                    breit_integral_indices_4(i) = breit_integral_indices_4(i-1)
                    breit_integral_values_4(i) = breit_integral_values_4(i-1)
                end do
                !
                !   Put the new index and value into storage
                !
                breit_integral_indices_4(table_position+1) = packed_integral_index
                breit_integral_values_4(table_position+1) = radial_integral_value
                !
            end if
        !
        else
            !
            !   This branch is executed only once per type of integral
            !
            first_breit_integral_lookup(4) = .false.
            !
            !   Designate the initial storage for arrays INDTPx and VALTPx;
            !   Array NDTPA stores the array dimensions
            !
            breit_integral_capacities(4) = 10000
            !
            !   Compute the integral's value
            !
            radial_integral_value = lookup_breit_radial_integral (4,orbital_a,orbital_b,orbital_c,orbital_d,radial_integral_order)
            !
            !   Initialise the integral counter
            !
            breit_integral_counts(4) = 1
            !
            !   Store the integral and its value
            !
            call alloc (breit_integral_indices_4,breit_integral_capacities(4),"INDTP4", "BRINT4")
            call alloc (breit_integral_values_4,breit_integral_capacities(4),"VALTP4", "BRINT4")
            breit_integral_indices_4(1) = packed_integral_index
            breit_integral_values_4(1) = radial_integral_value
            !
        end if
        !
        return
    end subroutine breit_radial_integral_4
end module breit_radial_integral_4_interface
