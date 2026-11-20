! Evaluate the double-orbital Slater radial integral using the Hartree Y function.
module slater_radial_integral
    implicit none (type, external)
    private
    public :: slater
contains
      real(real64) function slater (orbital_a, orbital_b, orbital_c, orbital_d, tensor_rank)
      use, intrinsic :: iso_fortran_env, only: real64
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use orbital_csf_basis, only: orbital_principal_n, orbital_symmetry_labels
      use radial_quadrature_workspace,          only: weighted_integrand => ta, hartree_y_values => tb, last_radial_index => mtp
      use radial_wavefunctions,          only: orbital_last_grid_index, radial_small_components, radial_large_components
      use hartree_yz_functions, only: yzk
      use radial_quadrature, only: quad
      implicit none (type, external)
      integer, intent(in) :: orbital_a
      integer, intent(in) :: orbital_b
      integer, intent(in) :: orbital_c
      integer, intent(in) :: orbital_d
      integer, intent(in) :: tensor_rank
      integer :: radial_index
      real(real64) :: slater_integral
      call yzk (tensor_rank, orbital_b, orbital_d)
      if (tensor_rank==0 .and. orbital_b==orbital_d) then
         last_radial_index = min(orbital_last_grid_index(orbital_a),orbital_last_grid_index(orbital_c))
      else
         last_radial_index = min(min(last_radial_index,orbital_last_grid_index(orbital_a)),orbital_last_grid_index(orbital_c))
      end if
      weighted_integrand(1) = 0.0_real64
      do radial_index = 2, last_radial_index
         weighted_integrand(radial_index) = (radial_large_components(radial_index, &
            & orbital_a)*radial_large_components(radial_index,orbital_c) + radial_small_components(radial_index,orbital_a)* &
            & radial_small_components(radial_index, &
               & orbital_c))*radial_derivative_over_radius(radial_index)*hartree_y_values(radial_index)
      end do
      call quad (slater_integral)
      slater = slater_integral
      if (ldbpr(10)) write (99, '(/,"  (",1i1,")",/," R   (",1i2,1a2,",",1i2,1a2,";",1i2,1a2,",",1i2,1a2,") ","= ",&
         &1p,d19.12)') tensor_rank, orbital_principal_n(orbital_a), orbital_symmetry_labels(orbital_a), &
         & orbital_principal_n(orbital_b), &
         & orbital_symmetry_labels(orbital_b), orbital_principal_n(orbital_c)&
         , orbital_symmetry_labels(orbital_c), orbital_principal_n(orbital_d), orbital_symmetry_labels(orbital_d), slater
      return
      return
      end function slater
end module slater_radial_integral
