!***********************************************************************
!                                                                      *
subroutine read_transition_angular_coefficients(operator_parity, transition_type_index, name, transition_angular_unit)
    ! Module dependencies
    !                                                                      *
    !   This routine loads  coefficients with parity and  rank specified   *
    !   by KP(JKP) into the arrays ISLDR and XSLDR.  IOPAR is the parity   *
    !   (+/- 1) and is determined from the sign of  KP(JKP).               *
    !                                                                      *
    !                                         Last revision: 28 Dec 1992   *
    !   Updated by Jacek Bieron               Last revision: 10 Mar 1994   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: orbital_index_radix, max_orbitals
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
    use first_order_csf_cutoffs, only: iccut
    use transition_off_diagonal_flags, only: noffd1, noffd2
    use orbital_csf_basis, only: max_orbitals, orbital_count
    use transition_matrix_workspace, only: nsdim, hb1, hb2, hc1, hc2, hm1, hm2, jja, jjb, lk, kk, ntdim, &
        xsldr, totc, totb, isldr, isldr1, nint, ninteg, nptr, lab, nkp, kp, ltc

    ! Procedure interfaces
    use resize_transition_matrix_arrays_interface, only: resize_transition_matrix_arrays
    use resize_transition_integral_arrays_interface, only: resize_transition_integral_arrays
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: transition_angular_unit
    integer  :: operator_parity
    integer, intent(in) :: transition_type_index
    character  :: name(2)*24

    ! Local declarations
    integer :: m, k, iblki, iblkf, i, labl, ncsf, j
    real(real64) :: shell_coefficients(max_orbitals)
    logical :: cache_exists
    !-----------------------------------------------
    !
    !   Read the data back as required by OSCL conventions
    !

    m = 0
    k = 0
    !
    read (transition_angular_unit) iblki, iblkf, orbital_count, nkp
    read (transition_angular_unit) nint
    !
    do i = 1, nint
        !
        read (transition_angular_unit) labl, ncsf
        !
        m = m + 1
        !bieron
        if (m >= nsdim) call resize_transition_matrix_arrays (jja, jjb, hb1, hb2, hc1, &
            hc2, hm1, hm2, lab, nptr, nsdim, 2)
        lab(m) = labl
        !
        !   Read configuration pairs and coefficients for this integral
        !
        do while (ncsf + k > ntdim)
            call resize_transition_integral_arrays (isldr, isldr1, xsldr, ntdim, 2)
        end do
        nptr(m) = k
        read (transition_angular_unit) (isldr(j + k),isldr1(j + k),xsldr(j + k),j=1,ncsf)
        !          write(*,*) (ISLDR(J+K),XSLDR(J+K),J = 1,NCSF)
        k = k + ncsf
        !
    end do
    !
    !   Close (and hence release) the scratch file
    !
    !     IF (AVAIL) THEN
    !
    !        CLOSE (unit=coefficient_scratch_unit,status="DELETE")
    !      ENDIF
    !
    nptr(m+1) = k
    ninteg = m
    !
    return
    !
    return
    !
end subroutine read_transition_angular_coefficients
