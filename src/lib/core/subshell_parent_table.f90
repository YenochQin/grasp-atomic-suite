! Map active subshell occupations to the table rows for their possible parent states.
module subshell_parent_table
    implicit none (type, external)

    ! Module interface
    private
    public :: lookup_subshell_parents
contains
    subroutine lookup_subshell_parents(orbital_indices, occupations, subshell_capacities, parent_table_rows)
        ! Module dependencies
        use subshell_term_table, only: nrows
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_indices(4)
        integer, intent(inout) :: occupations(4)
        integer, intent(in) :: subshell_capacities(4)
        integer, intent(inout) :: parent_table_rows(4)

        ! Local declarations
        integer :: table_population_index(4)
        integer :: parent_electron_count, parent_hole_count
        integer  :: active_shell_position
        if (orbital_indices(1) == orbital_indices(2)) occupations(1) = occupations(2) - 1
        if (orbital_indices(3) == orbital_indices(4)) occupations(3) = occupations(4) - 1
        do active_shell_position = 1, 4
            if (occupations(active_shell_position)<=0 .or. &
                & occupations(active_shell_position)>subshell_capacities(active_shell_position)) then
                write (*, '("lookup_subshell_parents: ",1i3," Electrons in shell ",1i3," with 2j+1 = ",1i3)') &
                    occupations(active_shell_position), orbital_indices(active_shell_position), &
                    & subshell_capacities(active_shell_position)
                stop
            end if
            parent_electron_count = occupations(active_shell_position) - 1
            parent_hole_count = subshell_capacities(active_shell_position) - parent_electron_count
            table_population_index(active_shell_position) = min(parent_electron_count,parent_hole_count) + 1
            if (table_population_index(active_shell_position) /= 1) then
                parent_table_rows(active_shell_position) = &
                    & (subshell_capacities(active_shell_position)*(subshell_capacities(active_shell_position)-2))/8 + &
                    & table_population_index(active_shell_position)
            else
                parent_table_rows(active_shell_position) = 1
            end if
            if (parent_table_rows(active_shell_position) <= nrows) cycle
            write (*, '("lookup_subshell_parents: Subshell term table capacity exceeded")')
            stop
        end do
        return
        return
    end subroutine lookup_subshell_parents
end module subshell_parent_table
