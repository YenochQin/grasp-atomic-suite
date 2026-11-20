!***********************************************************************
!                                                                      *
real(real64) function direct_angular_coefficient (multipole_rank, csf_index, orbital_a, orbital_b)
    ! Module dependencies
    !                                                                      *
    !   This routine evaluates a coefficient                               *
    !                                                                      *
    !                                K                                     *
    !                               f   (IA,IB)                            *
    !                                IR                                    *
    !                                                                      *
    !   Here  K  is the multipolarity, IR  is the sequence number of the   *
    !   configuration, and  IA  and  IB  are orbital  sequence  numbers.   *
    !   ( I P Grant,  B J McKenzie,  P H Norrington, D F Mayers, and N C   *
    !   Pyper,  Computer Phys Commun 21 (1980) 207-231, Eqs (6). )         *
    !                                                                      *
    !   Call(s) to: [LIB92]: CLRX, IQ.                                     *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford  Last revision: 21 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    ! . Two calls to IQ are physically inlined here to reduce the overhead.
    ! . In the inlined text, the checks for the integer range (ref.
    !    lib92/iq.f for details) are removed since FCO is called only by
    !    SETCOF and SETHAM where the input arguments always fall in the
    !    suitable range.
    !XHH 1997.03.05
    !
    !************************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use orbital_csf_basis, only: csf_subshell_occupations, orbital_kappa, orbital_two_j
    use relativistic_three_j, only: clrx
    use csf_quantum_queries, only: iq
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: multipole_rank
    integer  :: csf_index
    integer , intent(in) :: orbital_a
    integer , intent(in) :: orbital_b

    ! Local declarations
    integer :: iqa, iqf, kappa
    real(real64) :: fac
    !-----------------------------------------------
    !
    if (orbital_a == orbital_b) then
        !
        iqa = iq (orbital_a,csf_index)
        !         IQA = IBITS(IIQA((IA - 1)/4 + 1,IR),8*MOD(IA - 1,4),8)
        !
        if (multipole_rank == 0) then
            direct_angular_coefficient = dble((iqa*(iqa - 1))/2)
        else
            iqf = orbital_two_j(orbital_a) + 1
            if (iqa == iqf) then
                kappa = orbital_kappa(orbital_a)
                fac = clrx(kappa,multipole_rank,kappa)*dble(iqa)
                direct_angular_coefficient = -0.5_real64*fac*fac
            else
                direct_angular_coefficient = 0.0_real64
            end if
        end if
    !
    else
        !
        if (multipole_rank == 0) then
            direct_angular_coefficient = dble (iq (orbital_a,csf_index)*iq (orbital_b,csf_index))
        !            FCO = DBLE(IBITS(IIQA((IA - 1)/4 + 1,IR),8*MOD(IA - 1,4),8)*IBITS(&
        !               IIQA((IB - 1)/4 + 1,IR),8*MOD(IB - 1,4),8))
        else
            direct_angular_coefficient = 0.0_real64
        end if
        !
    end if
    return
    !*
    !  300 FORMAT (/'  ',1I2
    !     :        /' f    (',1I2,1A2,',',1I2,1A2,') = ',1PD21.14,
    !     :        /'  ',1I3/)
    !*
end function direct_angular_coefficient
