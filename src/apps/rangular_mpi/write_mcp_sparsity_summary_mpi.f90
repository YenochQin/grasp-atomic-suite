module write_mcp_sparsity_summary_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_mcp_sparsity_summary_mpi
contains

    !***********************************************************************
    subroutine write_mcp_sparsity_summary_mpi(print_details,nonzero_element_count, csf_count)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        !***********************************************************************
        use mpi_state, only: mpi_comm_world, mpi_error_code, mpi_integer, mpi_rank, mpi_reduce, mpi_sum
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: nonzero_element_count
        integer, intent(in) :: csf_count
        logical, intent(in) :: print_details

        ! Local constants
        character(len=*), parameter :: routine_name ="write_mcp_sparsity_summary_mpi"

        ! Local declarations
        integer :: nnonz_a
        !-----------------------------------------------
        !
        call mpi_reduce (nonzero_element_count, nnonz_a, 1, mpi_integer, mpi_sum, 0,  &
            mpi_comm_world, mpi_error_code)
        if (mpi_rank == 0) &
            write(output_unit, *)" ... complete; density of non-zero elements of H(DC):", &
            nnonz_a, "/", (csf_count*(csf_count + 1))/2
        !
        !   Debug printout
        !
        if (print_details) then
            if (mpi_rank == 0) then
                call ensure_debug_output()
                write (angular_debug_unit, *)
                call ensure_debug_output()
                write (angular_debug_unit, *) "From ", routine_name, " :"
                call ensure_debug_output()
                write (angular_debug_unit, '(" Number of nonzero elements in H(DC): ",1i4)') nnonz_a
                write (output_unit, *) "This part not finished. See ", routine_name
                call ensure_debug_output()
                write (angular_debug_unit, *) "This part not finished. See ", routine_name
            end if
        end if
        return
    end subroutine write_mcp_sparsity_summary_mpi
end module write_mcp_sparsity_summary_mpi_interface
