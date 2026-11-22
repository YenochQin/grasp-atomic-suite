! Prepare the hydrogenic reference data and factors used in the QED self-energy estimate.
module build_qed_self_energy_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: build_qed_self_energy
contains

    !   This  routine estimates corrections to the  energy levels due to
    !   self-energy.
    !
    ! Last update: 30 Oct 1992
    ! Modified by Xinghong He. Last update: 24 Jun 1997
    subroutine build_qed_self_energy(jstate, slfint, ucf)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: max_orbitals, max_radial_points
        use atomic_calculation_state, only: max_radial_points, nuclear_charge, pi, speed_of_light_au
        use eigenstate_expansions, only: csf_mixing_coefficients
        use orbital_csf_basis, only: orbital_count, csf_count, orbital_principal_n, orbital_kappa
        use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, radial_large_components, radial_small_components
        use csf_quantum_queries, only: csf_occupation

        ! Procedure interfaces
        use radial_density_ratio_interface, only: radial_density_ratio
        use hydrogenic_self_energy_factor_interface, only: hydrogenic_self_energy_factor
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: jstate
        real(real64), intent(out) :: slfint(max_orbitals)
        real(real64), intent(out) :: ucf(:)

        ! Local declarations
        integer :: maxiter, orbital_index, i, csf_index, npj, kappa, last_radial_index
        real(real64) :: ptemp(max_radial_points), qtemp(max_radial_points)
        real(real64) :: ucfj, zeff, ratio, valu
        ! Pre-set tolerable number for iteration in finding effective
        ! nuclear charge.
        maxiter = 20
        ! Modified so that UCFJ describes the current eigenstate
        do orbital_index = 1, orbital_count
            ucfj = 0.0_real64
            i = jstate
            do csf_index = 1, csf_count
                ucfj = ucfj + dble(csf_occupation(orbital_index, &
                    & csf_index))*csf_mixing_coefficients(csf_index + (i - 1)*csf_count)**2
            end do

            ucf(orbital_index) = ucfj
            ! UCF(J) = UCFJ/DBLE (NCF)
        end do
        do orbital_index = 1, orbital_count
            npj = orbital_principal_n(orbital_index)
            if (npj <= 8) then
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
    end subroutine build_qed_self_energy
end module build_qed_self_energy_interface
