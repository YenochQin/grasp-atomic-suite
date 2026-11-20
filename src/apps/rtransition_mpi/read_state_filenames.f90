module read_state_filenames_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_state_filenames
contains

    !***********************************************************************
    !                                                                      *
    subroutine read_state_filenames(name)
        ! Module dependencies
        use calculation_state_files, only: read_state_pair
        !                                                                      *
        !   Determines the name of the initial and final states                *
        !                                                                      *
        !   Written by Per Jonsson                                             *
        !                                                                      *
        !***********************************************************************
        implicit none (type, external)

        ! Arguments
        character(len=24), intent(out) :: name(2)
        call read_state_pair(name)

        return
    end subroutine read_state_filenames
end module read_state_filenames_interface
