! Evaluate the near-origin radial density contribution used in the field shift.
module field_shift_radial_integral_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: field_shift_radial_integral
contains

    !   The value of RINTDENS is an approximation to:
    !
    !      (4pi^)-1  r^-2 |P (r)*P (r) + Q (r)*Q (r) | r -> 0
    !
    ! Written by Per Jonsson. Last revision: 24 Dec 1992
    real(real64) function field_shift_radial_integral (first_orbital, second_orbital)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state, only: pi
        use radial_grid_state, only: radial_grid
        use radial_wavefunctions, only: radial_large_components, radial_small_components

        ! Procedure interfaces
        use interpolate_electron_density_interface, only: interpolate_electron_density
        implicit none (type, external)

        ! Arguments
        integer  :: first_orbital, second_orbital

        ! Local declarations
        integer :: radial_index, sample_index, first_grid_index
        real(real64) :: sample_radii(3), sample_densities(3)
        real(real64) :: density_at_origin

        first_grid_index = 10

        sample_index = 1
        do radial_index = first_grid_index,first_grid_index+4,2
            sample_radii(sample_index) = radial_grid(radial_index)
            sample_densities(sample_index) = (radial_large_components(radial_index,first_orbital)* &
                & radial_large_components(radial_index, &
                & second_orbital) + radial_small_components(radial_index,first_orbital)* &
                & radial_small_components(radial_index,second_orbital))/(4.0_real64*pi*radial_grid(radial_index)* &
                    & radial_grid(radial_index))
            sample_index = sample_index + 1
        end do
        call interpolate_electron_density (sample_radii,sample_densities,density_at_origin)
        field_shift_radial_integral = density_at_origin
        return
    end function field_shift_radial_integral
end module field_shift_radial_integral_interface
