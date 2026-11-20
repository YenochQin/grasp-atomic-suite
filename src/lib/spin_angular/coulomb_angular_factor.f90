! Form the Coulomb angular factor from two reduced spherical harmonics.
module coulomb_angular_factor
    implicit none (type, external)

    ! Module interface
    private
    public :: coulom
contains
    subroutine coulom(abs_kappa_a,abs_kappa_b,abs_kappa_c,abs_kappa_d,orbital_l_a,orbital_l_b,orbital_l_c, &
        orbital_l_d,multipole_rank,coulomb_angular_factor)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only:  angular_tolerance => eps, zero
        use angular_triangle_and_parity, only: ittk
        use reduced_angular_matrix_element, only: cre
        implicit none (type, external)

        ! Arguments
        integer, intent(in)        :: abs_kappa_a,abs_kappa_b,abs_kappa_c,abs_kappa_d,orbital_l_a,orbital_l_b, &
            orbital_l_c,orbital_l_d,multipole_rank
        real(real64), intent(out)  :: coulomb_angular_factor

        ! Local declarations
        integer :: harmonic_rank, twice_phase_exponent
        coulomb_angular_factor=zero
        if(ittk(orbital_l_a,orbital_l_c,multipole_rank) == 0)return
        if(ittk(orbital_l_b,orbital_l_d,multipole_rank) == 0)return
        harmonic_rank=(2*multipole_rank+1)/2
        coulomb_angular_factor=cre (abs_kappa_a,harmonic_rank,abs_kappa_c)
        if(dabs(coulomb_angular_factor) < angular_tolerance)return
        coulomb_angular_factor=coulomb_angular_factor*cre (abs_kappa_b,harmonic_rank,abs_kappa_d)
        if(dabs(coulomb_angular_factor) < angular_tolerance)return
        twice_phase_exponent=orbital_l_c-2*multipole_rank-orbital_l_a+orbital_l_d-orbital_l_b
        if((twice_phase_exponent/4)*4 /= twice_phase_exponent)coulomb_angular_factor=-coulomb_angular_factor
        return
    end subroutine coulom
end module coulomb_angular_factor
