!***********************************************************************

subroutine open_summary_file(name)
    ! Module dependencies
    use calculation_file_units, only: summary_output_unit
    !
    !   Open the  .csum  file on stream 24.
    !   Xinghong He                                          10 Jun 1998
    !
    !***********************************************************************
    use file_open, only: openfl
    implicit none (type, external)

    ! Arguments
    character(len=*) , intent(in) :: name

    ! Local declarations
    integer :: open_status
    !-----------------------------------------------------------------------


    call openfl (summary_output_unit, trim(name)//".csum", "FORMATTED", "UNKNOWN", open_status)

    return
end subroutine open_summary_file
