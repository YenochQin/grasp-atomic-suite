module merge_transition_csf_bases_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: merge_transition_csf_bases
contains

    !***********************************************************************
    !                                                                      *
    subroutine merge_transition_csf_bases(state_name)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        !                                                                      *
        !   Entry routine for merging two csl lists                            *
        !                                                                      *
        !***********************************************************************

        ! Procedure interfaces
        use load_initial_transition_basis_interface, only: load_initial_transition_basis
        use load_final_transition_basis_interface, only: load_final_transition_basis
        use merge_transition_basis_records_interface, only: merge_transition_basis_records
        implicit none (type, external)

        ! Arguments
        character  :: state_name(2)*24

        ! Local declarations
        integer :: initial_core_subshell_count, final_core_subshell_count
        !-----------------------------------------------
        !

        write (output_unit, *)
        write (output_unit, *) "MRGCSL: Execution begins ..."
        !
        !   Load the first  .csl  file
        !
        call load_initial_transition_basis (initial_core_subshell_count, state_name(1))
        !
        !   Load the second  .csl  file
        !
        call load_final_transition_basis (final_core_subshell_count, state_name(2))
        !
        !   Merge the two  .csl  lists, observe that there may be doublets
        !   among the CSF's
        !
        call merge_transition_basis_records (state_name, initial_core_subshell_count, final_core_subshell_count)

        return
    end subroutine merge_transition_csf_bases
end module merge_transition_csf_bases_interface
