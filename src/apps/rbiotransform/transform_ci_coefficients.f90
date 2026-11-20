!***********************************************************************
!                                                                      *
!                                                                      *
subroutine transform_ci_coefficients(input_mixing_coefficients, csf_count, level_count, l, subshell_count, t, nin, ntestg, &
    & transformed_mixing_coefficients, scr)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    !                                                                      *
    !   Calculate the action of the  operator                              *
    !                                                                      *
    !    T(L) = T(nl) T(n-1 L ) .... T(1 l)                                *
    !                                                                      *
    !   where                                                              *
    !                                                                      *
    !    T(i l ) = ( sum(k=0,2l+1) t(i l ) ** k ) (t_ilil ** E(il,il))     *
    !    t(i l ) = sum(j.lt.i) t(j,i) E(lj,li)                             *
    !    E(lj,li) total symmetric shell excitation operator                *
    !                                                                      *
    !   Modification of the CITRA routine for MCHF                         *
    !                                                                      *
    !   =====                                                              *
    !   Input                                                              *
    !   =====                                                              *
    !   CIIN      : Input CI vectors ( Destroyed in the process )          *
    !   NCSCF     : Length of CI expansion                                 *
    !   NCIV      : Number of CI vectors                                   *
    !   L         : L value of shells in excitations                       *
    !   NSHL      : Number of shells with this L                           *
    !   T         : List of coefficients                                   *
    !   NIN       : Number of inactive shells                              *
    !                                                                      *
    !   ======                                                             *
    !   Output                                                             *
    !   ======                                                             *
    !                                                                      *
    !   CIOUT : List of output CI vectors                                  *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use biorthogonal_subshell_mapping, only: max_orbitals, nlmax, nshlp, nshlpp
    use orbital_csf_basis, only: max_orbitals, orbital_kappa

    ! Procedure interfaces
    !
    use write_matrix_interface, only: write_matrix
    use scale_vector_interface, only: scale_vector
    use apply_occupation_factor_interface, only: apply_occupation_factor
    use apply_subshell_transfer_interface, only: apply_subshell_transfer
    use combine_vectors_interface, only: combine_vectors
    use copy_vector_interface, only: copy_vector
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: csf_count
    integer, intent(in) :: level_count
    integer, intent(in) :: l
    integer, intent(in) :: subshell_count
    integer, intent(in) :: nin
    integer, intent(in) :: ntestg
    real(real64)  :: input_mixing_coefficients(csf_count,level_count)
    real(real64), intent(in) :: t(subshell_count,subshell_count)
    real(real64)  :: transformed_mixing_coefficients(csf_count,level_count)
    real(real64)  :: scr(csf_count,level_count)

    ! Local declarations
    integer :: ntestl, ntest, iin, ipot, i, n
    real(real64) :: factor, tii, xnfaci
    !-----------------------------------------------
    !

    ntestl = 0
    ntest = max(ntestg,ntestl)
    ntest = 0

    !
    if (ntest >= 10) then
        write (output_unit, *)
        write (output_unit, *) " ***************"
        write (output_unit, *) " Entering CITRAG"
        write (output_unit, *) " ***************"
        write (output_unit, *)
    end if
    if (ntest >= 100) then
        write (output_unit, *) " Input CI vectors"
        call write_matrix (input_mixing_coefficients, csf_count, level_count, csf_count, level_count)
        write (output_unit, *) " Transformation matrix T"
        call write_matrix (t, subshell_count, subshell_count, subshell_count, subshell_count)
    end if
    !
    !. Factor from inactive shells
    !
    if (nin /= 0) then
        factor = 1.0_real64
        do iin = 1, nin
            factor = factor*t(iin,iin)
        end do
        !
        !       IPOT = 2*(2*L+1)  (number of m_lm_s. This should be replaced
        !                          by (2j+1) corresponding to L)
        !
        ipot = 2*iabs(orbital_kappa(nshlp(l,iin)))
        factor = factor**ipot
        call scale_vector (input_mixing_coefficients, factor, level_count*csf_count)
    end if
    if (nin == subshell_count) call copy_vector (input_mixing_coefficients, transformed_mixing_coefficients, &
        & level_count*csf_count)
    !
    do i = nin + 1, subshell_count
        if (ntest >= 100) write (output_unit, *) " Loop I,L = ", i, l
        !
        !. The diagonal contribution
        !
        tii = t(i,i)
        call apply_occupation_factor (input_mixing_coefficients, csf_count, level_count, i, l, tii, &
            & transformed_mixing_coefficients, ntestg)
        !        IF (LWORK2.LT.NCIV*NCSF) THEN
        !          WRITE(*,*) 'In CITRAG: Dimension of LWORK2 must be',
        !     &               'increased to at least',NCIV*NCSF
        !        ENDIF
        call copy_vector (transformed_mixing_coefficients, scr, level_count*csf_count)
        !
        !.  Off diagonal contributions
        !
        xnfaci = 1.0_real64
        do n = 1, 2*iabs(orbital_kappa(nshlp(l,i)))
            if (ntest >= 100) write (output_unit, *) " Loop N = ", n
            !
            !   T ** (N-1) is supposed to be in SCR, copy to CIIN
            !   and apply S
            !
            call copy_vector (scr, input_mixing_coefficients, level_count*csf_count)
            call apply_subshell_transfer (input_mixing_coefficients, csf_count, level_count, i, l, t(1,i), subshell_count, scr, &
                & ntestg)
            xnfaci = xnfaci/float(n)
            call combine_vectors (transformed_mixing_coefficients, transformed_mixing_coefficients, scr, 1.0_real64, xnfaci, &
                & level_count*csf_count)
        end do
        call copy_vector (transformed_mixing_coefficients, input_mixing_coefficients, level_count*csf_count)
        !
    end do
    !
    if (ntest >= 100) then
        write (output_unit, *) " Output CI vectors L = ", l
        call write_matrix (transformed_mixing_coefficients, csf_count, level_count, csf_count, level_count)
    end if
    !
    return
end subroutine transform_ci_coefficients
