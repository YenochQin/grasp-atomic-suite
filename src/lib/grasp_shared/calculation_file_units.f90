! Units for files shared by the setup, computational, and reporting routines.
module calculation_file_units
    implicit none (type, external)

    ! Module state
    ! For the units initialized to -1, that value marks a file that has yet to be assigned.
    integer :: summary_output_unit = -1
    integer :: calculation_log_unit = -1
    integer :: angular_log_unit
    integer :: mixing_input_unit
    integer :: mixing_output_unit
    integer :: hyperfine_output_unit
    integer :: lsj_transition_output_unit = -1
    integer :: one_body_cache_unit = -1
    integer :: two_body_cache_unit = -1

contains

    subroutine ensure_summary_output()
        ! Some MPI ranks historically create their own summary on the first write.
        if (summary_output_unit == -1) then
            open(newunit=summary_output_unit, file="fort.24", status="unknown", form="formatted")
        end if
    end subroutine ensure_summary_output

    subroutine ensure_calculation_log()
        ! Standalone callers of the shared trace and input routines have no named log.
        if (calculation_log_unit == -1) then
            open(newunit=calculation_log_unit, file="fort.734", status="unknown", form="formatted")
        end if
    end subroutine ensure_calculation_log
end module calculation_file_units
