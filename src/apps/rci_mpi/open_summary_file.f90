module open_summary_file_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: open_summary_file
contains

    ! Open the  .csum  file on stream 24.
    ! Xinghong He                                          10 Jun 1998
    subroutine open_summary_file(name)
        ! Module dependencies
        use calculation_file_units, only: summary_output_unit
        use file_open, only: open_calculation_file
        implicit none (type, external)

        ! Arguments
        character(len=*) , intent(in) :: name

        ! Local declarations
        integer :: open_status

        call open_calculation_file (summary_output_unit, trim(name)//".csum", "FORMATTED", "UNKNOWN", open_status)

        return
    end subroutine open_summary_file
end module open_summary_file_interface
