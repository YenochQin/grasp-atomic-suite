!***********************************************************************
!                                                                      *
subroutine prepare_biorthogonal_mcp(name, ic, ntestg, inpci)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit

    ! Procedure interfaces
    !                                                                      *
    !   Entry routine for GENMCP. Controls the computation of the          *
    !   one-particle coupling coefficients.                                *
    !                                                                      *
    !   Written by Per Jonsson                                 June 1996   *
    !                                                                      *
    !***********************************************************************
    use load_merged_csf_basis_interface, only: load_merged_csf_basis

    ! Module dependencies
    use log_factorial_table, only: factt

    ! Procedure interfaces
    use write_biorthogonal_mcp_interface, only: write_biorthogonal_mcp
    use read_biorthogonal_mcp_interface, only: read_biorthogonal_mcp
    implicit none (type, external)

    ! Arguments
    integer  :: ic
    integer  :: ntestg
    integer  :: inpci
    character  :: name*24

    ! Local declarations
    integer :: unused_core_subshell_count
    !-----------------------------------------------
    !
    write (output_unit, *)
    write (output_unit, *) "GENMCP: Execution begins for ", name
    !
    !   Open, check, load data from, and close the  csl  file
    !
    call load_merged_csf_basis (name, unused_core_subshell_count,3)
    !
    !   Set up the table of logarithms of factorials for use by
    !   angular modules
    !
    call factt
    !
    !   Proceed with the generation of MCP coefficients
    !
    call write_biorthogonal_mcp (name, ic, ntestg, inpci)
    call read_biorthogonal_mcp (name, ic, ntestg, inpci)
    !
    !   Print completion message
    !
    write (output_unit, *)
    !
    return
end subroutine prepare_biorthogonal_mcp
