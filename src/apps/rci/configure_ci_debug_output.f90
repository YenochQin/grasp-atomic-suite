module configure_ci_debug_output_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: configure_ci_debug_output
contains

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
        use debug_output_settings, only: ldbpa, ldbpg, rmcdhf_debug_flags
        implicit none (type, external)
        !
        ldbpa = .false.
        ldbpg = .false.

        rmcdhf_debug_flags = .false.

        return
    end subroutine configure_ci_debug_output
end module configure_ci_debug_output_interface
