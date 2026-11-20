!***********************************************************************
!                                                                      *
integer function one_electron_occupation_allowed (ja, jb)
    ! Module dependencies
    !                                                                      *
    !   This routine is to check the occupation condition for one electron *
    !   operator.                                                          *
    !                                                                      *
    !   Call(s) to: [LIB92]: IQ.                                           *
    !                                                                      *
    !   Yu Zou                                Last revision: 8/16/00       *
    !                                                                      *
    !***********************************************************************
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use orbital_csf_basis, only: csf_count,orbital_count,csf_subshell_occupations
    use csf_quantum_queries, only: iq
    implicit none (type, external)

    ! Arguments
    integer  :: ja
    integer  :: jb

    ! Local declarations
    integer :: k, orbital_index, iqa, iqb
    !-----------------------------------------------
    !
    !
    one_electron_occupation_allowed = 0
    k = 0
    do orbital_index = 1, orbital_count
        iqa = iq(orbital_index,ja)
        iqb = iq(orbital_index,jb)
        if (iqa == iqb) cycle
        k = k + 1
        if (k > 2) return
        if (iabs(iqa - iqb) <= 1) cycle
        return
    end do
    if (k==2 .or. k==0) one_electron_occupation_allowed = 1
    return
end function one_electron_occupation_allowed
