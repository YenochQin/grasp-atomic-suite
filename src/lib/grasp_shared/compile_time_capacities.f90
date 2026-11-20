! Defines the orbital packing radix and capacities used to dimension shared arrays.
module compile_time_capacities
    implicit none (type, external)

    ! Module interface
    private
    public :: orbital_index_radix, max_radial_points, radial_workspace_capacity, max_orbitals, max_orbitals_minus_one, &
        & max_orbitals_minus_two

    ! Module constants
    integer, parameter :: orbital_index_radix = 215
    integer, parameter :: max_radial_points = 590
    integer, parameter :: radial_workspace_capacity = max_radial_points + 10
    integer, parameter :: max_orbitals = 127
    integer, parameter :: max_orbitals_minus_one = max_orbitals - 1
    integer, parameter :: max_orbitals_minus_two = max_orbitals - 2
end module compile_time_capacities
