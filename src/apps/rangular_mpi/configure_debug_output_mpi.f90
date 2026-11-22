! Configure diagnostic output consistently across MPI ranks.
module configure_debug_output_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: configure_debug_output_mpi
contains

    !   This routine calls setdbg to open the  .dbg  file and sets the
    !   arrays that control debug printout from the GENMCP program.
    !
    ! Written by Xinghong He. Last update: 29 Jun 1998
    subroutine configure_debug_output_mpi (debug, file_path)
        ! Module dependencies
        use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_logical, mpi_rank

        ! Procedure interfaces
        use configure_debug_output_interface, only: configure_debug_output

        ! Module dependencies
        use debug_output_settings, only: ldbpa, ldbpg
        implicit none (type, external)

        ! Arguments
        logical :: debug
        character(len=*) :: file_path

        if (mpi_rank == 0) then
            call configure_debug_output (debug, file_path)
        end if
        call mpi_bcast (debug,1,mpi_logical,0,mpi_comm_world,mpi_error_code)
        call mpi_bcast (ldbpa,5,mpi_logical,0,mpi_comm_world,mpi_error_code)
        call mpi_bcast (ldbpg,5,mpi_logical,0,mpi_comm_world,mpi_error_code)
        return
    end subroutine configure_debug_output_mpi
end module configure_debug_output_mpi_interface
