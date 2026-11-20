!***********************************************************************
!                                                                      *
real(real64) function exchange_angular_coefficient (multipole_rank, csf_index, orbital_a, orbital_b)
    ! Module dependencies
    !                                                                      *
    !   This routine evaluates a coefficient                               *
    !                                                                      *
    !                                K                                     *
    !                               g   (IA,IB)                            *
    !                                IR                                    *
    !                                                                      *
    !                                                                      *
    !   Here  K  is the multipolarity, IR  is the sequence number of the   *
    !   configuration, and  IA and IB are orbital sequence  numbers. See   *
    !   I P Grant,  B J McKenzie,  P H Norrington,  D F  Mayers, and N C   *
    !   Pyper, Computer Phys Commun 21 (1980) 207-231, Eq (7).             *
    !                                                                      *
    !   Call(s) to: [LIB92]: CLRX.                                         *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford  Last revision: 18 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use orbital_csf_basis, only: csf_subshell_occupations, orbital_kappa, orbital_two_j
    use relativistic_three_j, only: clrx
    use csf_quantum_queries, only: iq
    implicit none (type, external)

    ! Arguments
    integer  :: multipole_rank
    integer  :: csf_index
    integer , intent(in) :: orbital_a
    integer , intent(in) :: orbital_b

    ! Local declarations
    integer :: iqa, iqb
    real(real64) :: fac
    !      LOGICAL :: FULLA, FULLB
    !-----------------------------------------------
    !
    iqa = iq (orbital_a,csf_index)
    iqb = iq (orbital_b,csf_index)
    !GG      IQA = IBITS(IIQA((IA - 1)/4 + 1,IR),8*MOD(IA - 1,4),8)
    !GG      IQB = IBITS(IIQA((IB - 1)/4 + 1,IR),8*MOD(IB - 1,4),8)
    if (iqa==orbital_two_j(orbital_a) + 1 .or. iqb==orbital_two_j(orbital_b)+1) then
        fac = clrx(orbital_kappa(orbital_a),multipole_rank,orbital_kappa(orbital_b))
        exchange_angular_coefficient = -dble(iqa*iqb)*fac*fac
    else
        exchange_angular_coefficient = 0.0_real64
    end if

    return
end function exchange_angular_coefficient
