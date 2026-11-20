! Apply the phase and subshell degeneracy factors to the relativistic 3-j coefficient.
module reduced_angular_matrix_element
    implicit none (type, external)

    ! Module interface
    private
    public :: cre
contains
    real(real64) function cre (bra_kappa, tensor_rank, ket_kappa)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use relativistic_three_j, only: clrx
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: bra_kappa
        integer, intent(in) :: tensor_rank
        integer, intent(in) :: ket_kappa

        ! Local declarations
        integer :: abs_bra_kappa
        real(real64) :: subshell_degeneracy_product
        abs_bra_kappa = abs(bra_kappa)
        subshell_degeneracy_product = dble(4*abs_bra_kappa*iabs(ket_kappa))
        cre = sqrt(subshell_degeneracy_product)*clrx(bra_kappa,tensor_rank,ket_kappa)
        if (mod(abs_bra_kappa,2) == 1) cre = -cre
        return
    end function cre
end module reduced_angular_matrix_element
