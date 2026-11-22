module write_mcp_sparsity_summary_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_mcp_sparsity_summary
contains

    subroutine write_mcp_sparsity_summary(print_details,nonzero_element_count, csf_count)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: nonzero_element_count
        integer, intent(in) :: csf_count
        logical, intent(in) :: print_details

        ! Local constants
        character(len=*), parameter :: routine_name ="write_mcp_sparsity_summary"

        ! Local declarations
        integer :: nnonz_a
        nnonz_a = nonzero_element_count

        write(output_unit, *)" ... complete; density of non-zero elements of H(DC):", &
            nnonz_a, "/", (csf_count*(csf_count + 1))/2

        if (print_details) then
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

        return
    end subroutine write_mcp_sparsity_summary
end module write_mcp_sparsity_summary_interface
