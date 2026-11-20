!***********************************************************************
!                                                                      *
subroutine configure_debug_output_mpi (debug, fullname)
    ! Module dependencies
    !                                                                      *
    !   This routine calls setdbg to open the  .dbg  file and sets the     *
    !   arrays that control debug printout from the GENMCP program.        *
    !                                                                      *
    !   Written by Xinghong He                  Last update: 29 Jun 1998   *
    !                                                                      *
    !***********************************************************************
    use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_logical, mpi_rank

    ! Procedure interfaces
    use configure_debug_output_interface, only: configure_debug_output

    ! Module dependencies
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    implicit none (type, external)

    ! Arguments
    logical, intent(out) :: debug
    character(len=*) :: fullname
    !-----------------------------------------------
    !
    if (mpi_rank == 0) then
        call configure_debug_output (debug, fullname)
    end if
    call mpi_bcast (debug,1,mpi_logical,0,mpi_comm_world,mpi_error_code)
    call mpi_bcast (ldbpa,5,mpi_logical,0,mpi_comm_world,mpi_error_code)
    call mpi_bcast (ldbpg,5,mpi_logical,0,mpi_comm_world,mpi_error_code)
    return
end subroutine configure_debug_output_mpi
