! Shared compile-time limits. Keep upstream values as the default configuration.
! Edit here (or use the grid patcher), then rebuild all libraries and programs.
module suite_parameters
    implicit none
    private
    public :: orbital_index_radix, max_radial_points, radial_workspace_capacity, max_orbitals, max_orbitals_minus_one, &
       & max_orbitals_minus_two

    integer, parameter :: orbital_index_radix = 215
    integer, parameter :: max_radial_points = 590
    integer, parameter :: radial_workspace_capacity = max_radial_points + 10
    integer, parameter :: max_orbitals = 127
    integer, parameter :: max_orbitals_minus_one = max_orbitals - 1
    integer, parameter :: max_orbitals_minus_two = max_orbitals - 2
end module suite_parameters
