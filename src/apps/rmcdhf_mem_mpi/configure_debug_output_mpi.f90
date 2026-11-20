!***********************************************************************
!                                                                      *
subroutine configure_debug_output_mpi(debug_file_path)
    ! Module dependencies
    !                                                                      *
    !   This subroutine sets the arrays that control debug printout from   *
    !   the radial and angular modules of the GRASP92 suite.               *
    !                                                                      *
    !   Call(s) to: [LIB92]: OPENFL.                               *
    !                                                                      *
    !   Written by Farid A Parpia               Last update: 10 Dec 1992   *
    !   Modified bu Xinghong He                 Last update: 06 Jul 1998   *
    !                                                                      *
    !***********************************************************************
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_logical, mpi_rank

    ! Procedure interfaces
    use configure_debug_output_interface, only: configure_debug_output
    implicit none (type, external)

    ! Arguments
    character(len=*) , intent(in) :: debug_file_path
    !-----------------------------------------------
    !
    if (mpi_rank == 0) then
        call configure_debug_output (debug_file_path)
    end if

    call mpi_bcast (ldbpa, 5, mpi_logical, 0, mpi_comm_world, mpi_error_code)
    call mpi_bcast (ldbpg, 5, mpi_logical, 0, mpi_comm_world, mpi_error_code)
    call mpi_bcast (ldbpr,30, mpi_logical, 0, mpi_comm_world, mpi_error_code)

    return
end subroutine configure_debug_output_mpi
