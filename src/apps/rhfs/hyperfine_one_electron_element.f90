module hyperfine_one_electron_element_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: hyperfine_one_electron_element
contains

    !***********************************************************************
    !                                                                      *
    subroutine hyperfine_one_electron_element(bra_orbital_index, tensor_rank, ket_orbital_index, hyperfine_angular_factor, &
                & lande_g_angular_factor, lande_g_correction_angular_factor)
        ! Module dependencies
        !                                                                      *
        !   This  routine computes  the angular part  of the  reduced matrix   *
        !   elements of the magnetic and electric multipole operators as       *
        !   well as for th operators associated with the g_j factor            *
        !                                                                      *
        !   Calls to: [LIB92]: CLRX                                            *
        !                                                                      *
        !   Written by Per Jonsson                Last revision: 22 Oct 1999   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use orbital_csf_basis, only: orbital_kappa
        use relativistic_three_j, only: relativistic_three_j_value
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: bra_orbital_index
        integer, intent(in) :: tensor_rank
        integer, intent(in) :: ket_orbital_index
        real(real64), intent(out) :: hyperfine_angular_factor
        real(real64), intent(out) :: lande_g_angular_factor
        real(real64), intent(out) :: lande_g_correction_angular_factor

        ! Local declarations
        integer      :: bra_kappa, tensor_ket_kappa, bra_orbital_l, tensor_ket_orbital_l
        real(real64) :: angular_phase, phase_times_sqrt_ket_multiplicity
        !-----------------------------------------------
        !
        !
        !   Set KAP1 and KAP2
        !
        bra_kappa = orbital_kappa(bra_orbital_index)
        !
        if (mod(tensor_rank,2) == 1) then
            tensor_ket_kappa = -orbital_kappa(ket_orbital_index)
        else
            tensor_ket_kappa = orbital_kappa(ket_orbital_index)
        end if
        !
        !   Determine the l quantum numbers
        !
        if (bra_kappa > 0) then
            bra_orbital_l = bra_kappa
        else
            bra_orbital_l = (-bra_kappa) - 1
        end if
        !
        if (tensor_ket_kappa > 0) then
            tensor_ket_orbital_l = tensor_ket_kappa
        else
            tensor_ket_orbital_l = (-tensor_ket_kappa) - 1
        end if
        !
        if (mod(bra_orbital_l + tensor_rank + tensor_ket_orbital_l,2) == 0) then
            !
            !   Parity selection rule satisfied
            !
            !   Determine the phase factor
            !
            if (mod(bra_kappa + 1,2) == 0) then
                angular_phase = 1.0_real64
            else
                angular_phase = -1.0_real64
            end if
            !
            !   The other factor is \sqrt (2 j_2 + 1); since j = | \kappa | - 1/2,
            !   we have 2 j_2 + 1 = 2 | \kappa |; the factor \sqrt (2 j_2 + 1)
            !   has been accounted for in MCT
            !
            phase_times_sqrt_ket_multiplicity = angular_phase*sqrt(dble(2*abs(tensor_ket_kappa)))
            !
            if (mod(tensor_rank,2) == 1) then
                !
                !   These are for the magnetic multipole moments and the two operators of the g_j factor
                !
                hyperfine_angular_factor = (bra_kappa + orbital_kappa(ket_orbital_index))*relativistic_three_j_value(bra_kappa, &
                & tensor_rank,tensor_ket_kappa)*phase_times_sqrt_ket_multiplicity
                !
                lande_g_angular_factor = hyperfine_angular_factor
                lande_g_correction_angular_factor = -(bra_kappa + orbital_kappa(ket_orbital_index)-1.0_real64)* &
                & relativistic_three_j_value(bra_kappa,tensor_rank,tensor_ket_kappa)*phase_times_sqrt_ket_multiplicity
            else
                !
                !   These are for the electric multipole moments
                !
                hyperfine_angular_factor = relativistic_three_j_value(bra_kappa,tensor_rank,tensor_ket_kappa)* &
                & phase_times_sqrt_ket_multiplicity
                !
            end if
        !
        else
            !
            hyperfine_angular_factor = 0.0_real64
            lande_g_angular_factor = 0.0_real64
            lande_g_correction_angular_factor = 0.0_real64
            !
        end if
        !
        return
    end subroutine hyperfine_one_electron_element
end module hyperfine_one_electron_element_interface
