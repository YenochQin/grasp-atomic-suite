module field_shift_radial_integral_vector_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: field_shift_radial_integral_vector
contains

    !***********************************************************************
    !                                                                      *
    subroutine field_shift_radial_integral_vector (first_orbital,second_orbital,orbital_density_profiles,nuclear_grid_point_count)
        ! Module dependencies
        !                                                                      *
        !   The value of DINT1VEC is an approximation to:                      *
        !                                                                      *
        !                                                                      *
        !      (4pi^)-1  r^-2 |P (r)*P (r) + Q (r)*Q (r) |                     *
        !                       I     J       I     J                          *
        !                                                                      *
        !                                                                      *
        !                                                                      *
        !   Written by Per Jonsson  and Jörgen Ekman                           *
        !   Last revision: 24 Dec 2013                                         *
        !                                                                      *
        !***********************************************************************
        !...Translated by Gediminas Gaigalas 11/18/19
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,    only: max_orbitals
        use atomic_calculation_state,            only: pi
        use radial_grid_state, only: radial_point_count, radial_grid
        use radial_wavefunctions, only: max_orbitals, radial_large_components, radial_small_components
        implicit none (type, external)

        ! Arguments and local declarations
        real(real64), intent(out) :: orbital_density_profiles(max_orbitals, max_orbitals, radial_point_count)
        integer, intent(in) :: first_orbital, second_orbital
        integer :: radial_index
        integer, intent(in) :: nuclear_grid_point_count
        !-----------------------------------------------
        !
        do radial_index = 2, nuclear_grid_point_count
            orbital_density_profiles(first_orbital,second_orbital,radial_index) = (radial_large_components(radial_index, &
                & first_orbital)* &
                & radial_large_components(radial_index, &
                & second_orbital)+radial_small_components(radial_index,first_orbital)* &
                & radial_small_components(radial_index,second_orbital))/          &
                (4.0_real64*pi*radial_grid(radial_index)*radial_grid(radial_index))
        end do
        !
        return
    end subroutine field_shift_radial_integral_vector
end module field_shift_radial_integral_vector_interface
