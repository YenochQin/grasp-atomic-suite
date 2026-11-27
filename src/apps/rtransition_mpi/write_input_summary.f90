module write_input_summary_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_input_summary
contains

    subroutine write_input_summary(name, inpci, ilbl)
        ! Module dependencies
        use calculation_file_units, only: lsj_transition_output_unit, summary_output_unit
        use calculation_state_files, only: transition_report_path
        use transition_lsj_labels, only: iopen_status1, iopen_status2
        implicit none (type, external)

        ! Arguments
        character(*), intent(in) :: name(2)
        integer, intent(in) :: inpci, ilbl

        ! Local declarations
        logical, parameter :: with_phase = .false.
        character(:), allocatable :: path

        if (ilbl == 0) then
            path = transition_report_path(name(1), name(2), inpci == 0, .false., with_phase)
            open(newunit=summary_output_unit, file=path, form="formatted", status="replace", action="write")
        else if (ilbl == 1 .and. iopen_status1 == 0 .and. iopen_status2 == 0) then
            if (lsj_transition_output_unit == -1) then
                path = transition_report_path(name(1), name(2), inpci == 0, .true., with_phase)
                open(newunit=lsj_transition_output_unit, file=path, form="formatted", status="replace", action="write")
            end if
        end if
    end subroutine write_input_summary
end module write_input_summary_interface
