module second_spin_angular_factor_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: second_spin_angular_factor
contains

    !***********************************************************************
    !                                                                      *
    subroutine second_spin_angular_factor (spin_matrix_element_convention,i1,i2,spin_angular_factor)
        ! Module dependencies
        !                                                                      *
        !   This  routine computes                                             *
        !   [-kappa_a || sigma^(2) ||kappa_b] for convention 1                      *
        !   <-kappa_a || sigma^(2) ||kappa_b> for convention 2                      *
        !                                                                      *
        !                                                                      *
        !   Written by  G. Gaigalas                                            *
        !                                      Last revision:  11 Aprel 2019   *
        !                                                                      *
        !***********************************************************************
        !
        use, intrinsic :: iso_fortran_env, only: real64
        use orbital_csf_basis, only: orbital_kappa
        use racah_coefficients, only: evaluate_racah_coefficient
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: spin_matrix_element_convention, i1, i2
        real(real64), intent(out) :: spin_angular_factor

        ! Local declarations
        integer :: l1, l2, kap1, kap2, bra_two_j, ket_two_j, l2_tilda
        real(real64) :: rac
        !-----------------------------------------------
        !
        !   Set KAP1 and KAP2
        !
        kap1 = orbital_kappa(i1)
        kap2 = orbital_kappa(i2)
        !
        !   Determine the l quantum numbers
        !
        if (kap1 > 0) then
            l1 =  kap1
        else
            l1 = -kap1-1
        end if
        !
        if (kap2 > 0) then
            l2 =  kap2
        else
            l2 = -kap2-1
        end if
        !
        !   Determine the j quantum numbers and l_2 tilda
        !
        bra_two_j = iabs(orbital_kappa(i1))*2-1
        ket_two_j = iabs(orbital_kappa(i2))*2-1
        l2_tilda = ket_two_j - l2
        if (l2_tilda == l1) then
            !
            !   Determine the Racah W coefficients.
            !
            call evaluate_racah_coefficient (1,bra_two_j,1,ket_two_j,2*l1,2,rac)
            if (mod(1+bra_two_j+1+ket_two_j,4) == 0) rac=-rac
            spin_angular_factor =  rac*dsqrt(6.0_real64*(ket_two_j+1))
        else
            spin_angular_factor = 0.0_real64
        end if
        !
        !   Determine the phase factori
        !
        if(mod(2*l1+1+bra_two_j+2,4) == 0) spin_angular_factor=-spin_angular_factor
        if(spin_matrix_element_convention == 2) then
            spin_angular_factor = dsqrt(dble(bra_two_j+1))*spin_angular_factor
        end if
        return
    end subroutine second_spin_angular_factor
end module second_spin_angular_factor_interface
