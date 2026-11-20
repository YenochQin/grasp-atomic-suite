!***********************************************************************
!                                                                      *
subroutine second_spin_angular_factor (ipar,i1,i2,apart)
    ! Module dependencies
    !                                                                      *
    !   This  routine computes                                             *
    !   [-kappa_a || sigma^(2) ||kappa_b] if IPAR = 1                      *
    !   <-kappa_a || sigma^(2) ||kappa_b> if IPAR = 2                      *
    !                                                                      *
    !                                                                      *
    !   Written by  G. Gaigalas                                            *
    !                                      Last revision:  11 Aprel 2019   *
    !                                                                      *
    !***********************************************************************
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use orbital_csf_basis,           only: orbital_principal_n, orbital_kappa
    use racah_coefficients, only: dracah
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: ipar, i1, i2
    real(real64), intent(out) :: apart

    ! Local declarations
    integer :: l1, l2, kap1, kap2, j1, j2, l2_tilda
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
    j1 = iabs(orbital_kappa(i1))*2-1
    j2 = iabs(orbital_kappa(i2))*2-1
    l2_tilda = j2 - l2
    if (l2_tilda == l1) then
        !
        !   Determine the Racah W coefficients.
        !
        call dracah (1,j1,1,j2,2*l1,2,rac)
        if (mod(1+j1+1+j2,4) == 0) rac=-rac
        apart =  rac*dsqrt(6.0_real64*(j2+1))
    else
        apart = 0.0_real64
    end if
    !
    !   Determine the phase factori
    !
    if(mod(2*l1+1+j1+2,4) == 0) apart=-apart
    if(ipar == 2) then
        apart = dsqrt(dble(j1+1))*apart
    end if
    return
end subroutine second_spin_angular_factor
