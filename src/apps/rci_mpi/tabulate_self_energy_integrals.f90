module tabulate_self_energy_integrals_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: tabulate_self_energy_integrals
contains

    !   This  routine estimates the F(Z\alpha) function of self energy for *
    !   each orbital.
    !
    ! Modified from subroutine QED by Yu Zou, Last update: 13 Mar 2000
    subroutine tabulate_self_energy_integrals(slfint)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: max_orbitals, max_radial_points
        use atomic_calculation_state, only: max_radial_points, nuclear_charge, pi, speed_of_light_au
        use orbital_csf_basis, only: orbital_count, orbital_principal_n, orbital_kappa
        use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, radial_large_components, radial_small_components
        use qed_orbital_cutoffs, only: self_energy_cutoff_flag, self_energy_max_principal_n
        use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_integer, mpi_rank

        ! Procedure interfaces
        use radial_density_ratio_interface, only: radial_density_ratio
        use hydrogenic_self_energy_factor_interface, only: hydrogenic_self_energy_factor
        implicit none (type, external)

        ! Arguments
        real(real64), intent(out) :: slfint(max_orbitals)

        ! Local declarations
        integer :: maxiter, orbital_index, npj, kappa, last_radial_index, i, npjmax
        real(real64) :: ptemp(max_radial_points), qtemp(max_radial_points)
        real(real64) :: zeff, ratio, valu
        ! Pre-set tolerable number for iteration in finding effective
        ! nuclear charge.
        maxiter = 20
        if (mpi_rank == 0) then
            if (self_energy_cutoff_flag==1) then
                npjmax = self_energy_max_principal_n
            else
                npjmax = 8
            end if
        end if
        call mpi_bcast (npjmax,1,mpi_integer,0,mpi_comm_world,mpi_error_code)

        do orbital_index = 1, orbital_count
            npj = orbital_principal_n(orbital_index)
            if (npj <= npjmax) then
                ! Only orbitals with principal quantum number 8 or less can
                ! be treated by this section of code
                kappa = orbital_kappa(orbital_index)
                ! Begin by transferring the function to a temporary array
                last_radial_index = orbital_last_grid_index(orbital_index)
                ptemp(1) = 0.0_real64
                qtemp(1) = 0.0_real64
                do i = 2, last_radial_index
                    ptemp(i) = radial_large_components(i,orbital_index)
                    qtemp(i) = radial_small_components(i,orbital_index)
                end do
                zeff = nuclear_charge
                ratio = radial_density_ratio(ptemp,qtemp,last_radial_index,npj,kappa,zeff)
                valu = ratio*hydrogenic_self_energy_factor(npj,kappa,zeff)/dble(npj**3)
                slfint(orbital_index) = valu*zeff**4/(pi*speed_of_light_au**3)
            else
                ! The self-energy for orbitals with principal quantum number
                ! greater than 8 is set to zero
                slfint(orbital_index) = 0.0_real64
            end if
        end do
        ! Deallocate storage for the `generalised occupation numbers'
        return
    end subroutine tabulate_self_energy_integrals
end module tabulate_self_energy_integrals_interface
