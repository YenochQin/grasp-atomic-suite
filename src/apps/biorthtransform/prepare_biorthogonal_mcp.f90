! Prepare the angular coefficients used to transform mixing coefficients between orbital bases.
module prepare_biorthogonal_mcp_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: prepare_biorthogonal_mcp
contains

    !   one-particle coupling coefficients.
    !
    ! Written by Per Jonsson                                 June 1996
    subroutine prepare_biorthogonal_mcp(name, ic, global_trace_level, inpci, output_path)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit

        ! Procedure interfaces
        use load_merged_csf_basis_interface, only: load_merged_csf_basis

        ! Module dependencies
        use log_factorial_table, only: initialize_log_factorial_table

        ! Procedure interfaces
        use write_biorthogonal_mcp_interface, only: write_biorthogonal_mcp
        use read_biorthogonal_mcp_interface, only: read_biorthogonal_mcp
        implicit none (type, external)

        ! Arguments
        integer  :: ic
        integer  :: global_trace_level
        integer  :: inpci
        character, intent(in) :: name*24
        character(*), intent(in) :: output_path

        ! Local declarations
        integer :: unused_core_subshell_count

        write (output_unit, *)
        write (output_unit, *) "prepare_biorthogonal_mcp: Execution begins for ", name
        ! Open, check, load data from, and close the  csl  file
        call load_merged_csf_basis (name, unused_core_subshell_count,3)
        ! Set up the table of logarithms of factorials for use by
        ! angular modules
        call initialize_log_factorial_table
        ! Proceed with the generation of MCP coefficients
        call write_biorthogonal_mcp (name, ic, global_trace_level, inpci)
        call read_biorthogonal_mcp (name, ic, global_trace_level, inpci, output_path)

        write (output_unit, *)
        return
    end subroutine prepare_biorthogonal_mcp
end module prepare_biorthogonal_mcp_interface
