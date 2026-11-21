module breit_radial_integral_2_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: breit_radial_integral_2
contains

    !***********************************************************************
    !                                                                      *
    subroutine breit_radial_integral_2 (orbital_a,orbital_b,orbital_c,orbital_d,k,radial_integral_value)
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
        use breit_integral_tables, only: breit_integral_indices_2, breit_integral_values_2
        use orbital_csf_basis,           only: orbital_count
        use breit_rank_offsets, only: second_breit_integral_group_offsets
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_a, orbital_b, orbital_c, orbital_d, k
        real(real64), intent(out) :: radial_integral_value

        ! Local declarations
        integer :: packing_radix, packed_integral_index, table_position, search_upper_index, search_lower_index, &
            & search_midpoint_index
        !-----------------------------------------------

        packing_radix = orbital_count + 1

        packed_integral_index = ((orbital_a*packing_radix+orbital_b)*packing_radix+orbital_c)*packing_radix+orbital_d
        !
        search_lower_index = second_breit_integral_group_offsets(k)
        search_upper_index = second_breit_integral_group_offsets(k+1) - 1

        if (packed_integral_index < breit_integral_indices_2(search_lower_index) .or. &
            packed_integral_index > breit_integral_indices_2(search_upper_index)) then
            write(*,*) "Something wrong in breit_radial_integral_2"
            stop
        end if
        !
        !   The index is within the range of the indices stored; search
        !   for it in the list of indices
        !
        do while (search_upper_index-search_lower_index > 1)
            search_midpoint_index = (search_upper_index+search_lower_index)/2
            if (breit_integral_indices_2(search_midpoint_index) > packed_integral_index) then
                search_upper_index = search_midpoint_index
            else
                search_lower_index = search_midpoint_index
            end if
        end do
        !
        !   The range is bracketed to the extent possible
        !
        if (packed_integral_index == breit_integral_indices_2(search_upper_index)) then
            table_position = search_upper_index
        else if (packed_integral_index == breit_integral_indices_2(search_lower_index)) then
            table_position = search_lower_index
        else
            write(*,*) k,orbital_a,orbital_b,orbital_c,orbital_d,packed_integral_index
            write(*,*) "breit_radial_integral_2 Integral not found"
            stop
        end if
        !
        !   Return the value of the integral
        !   from storage

        radial_integral_value = breit_integral_values_2(table_position)
        !
        return
    end subroutine breit_radial_integral_2
end module breit_radial_integral_2_interface
