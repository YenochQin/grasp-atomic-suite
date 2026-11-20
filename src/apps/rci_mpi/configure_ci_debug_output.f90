!***********************************************************************
!                                                                      *
subroutine configure_ci_debug_output
    ! Module dependencies
    !                                                                      *
    !   This subroutine sets the arrays that control debug printout from   *
    !   the radial and angular modules of the GRASP92 suite.               *
    !                                                                      *
    !                                                                      *
    !   Written by Farid A Parpia               Last update: 21 Oct 1992   *
    !                                                                      *
    !***********************************************************************
    !-----------------------------------------------
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    implicit none (type, external)
    !
    ldbpa = .false.
    ldbpg = .false.

    ldbpr = .false.

    return
end subroutine configure_ci_debug_output
