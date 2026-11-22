! Configure diagnostic output consistently across MPI ranks.
module configure_debug_output_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: configure_debug_output_mpi
contains

    !   This subroutine sets the arrays that control debug printout from
    !   the radial and angular modules of the GRASP92 suite.
    !
    ! Written by Farid A Parpia. Last update: 10 Dec 1992
    ! Modified bu Xinghong He. Last update: 06 Jul 1998
    subroutine configure_debug_output_mpi(debug_file_path)
        ! Module dependencies
        use debug_output_settings, only: ldbpa, ldbpg, rmcdhf_debug_flags
        use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_logical, mpi_rank

        ! Procedure interfaces
        use configure_debug_output_interface, only: configure_debug_output
        implicit none (type, external)

        ! Arguments
        character(len=*) , intent(in) :: debug_file_path

        if (mpi_rank == 0) then
            call configure_debug_output (debug_file_path)
        end if

        call mpi_bcast (ldbpa, 5, mpi_logical, 0, mpi_comm_world, mpi_error_code)
        call mpi_bcast (ldbpg, 5, mpi_logical, 0, mpi_comm_world, mpi_error_code)
        call mpi_bcast (rmcdhf_debug_flags,30, mpi_logical, 0, mpi_comm_world, mpi_error_code)

        return
    end subroutine configure_debug_output_mpi
end module configure_debug_output_mpi_interface
