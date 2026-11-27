! Resolve and check the files for the initial and final states.
module read_state_filenames_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_state_filenames
contains

    ! Determines the name of the initial and final states
    ! In addition this subroutine determines which J symmetries
    ! that are to be transformed
    !
    ! Written by Per Jonsson
    subroutine read_state_filenames(name)
        ! Module dependencies
        use input_error_handling, only: fail_input_io
        use calculation_state_files, only: read_state_pair
        use, intrinsic :: iso_fortran_env, only: output_unit
        use biorthogonal_transition_selection, only: itjqj1, itjqj2, jqj1, jqj2, ntrans
        use yes_no_input, only: read_yes_no
        implicit none (type, external)

        ! Arguments
        character(len=24), intent(out) :: name(2)

        ! Local declarations
        integer :: i
        logical :: yes
        ! Obtain the names of the initial and final state files
        call read_state_pair(name)

        ! Per april 2007
        ! Check if the initial and final states are identical.

        if (trim(name(1)) == trim(name(2))) then
            print *
            print *, " Initial and final states are identical and there is"
            print *, " no need for the biorthogonal transformation. Just  "
            print *, " copy name.wfn to name.bw and name.mix/name.cimix to name.(c)bm "
            print *, " and run rtransition.                               "
            print *

            print *, " Do you want to continue anyway ? "
            yes = read_yes_no ()
            if (yes.eqv..false.) stop
        end if

        ! end Per 2007
        write (output_unit, *) " Transformation of all J symmetries?"
        yes = read_yes_no()
        if (yes) then
            ntrans = 0
        else
            ntrans = 1
            write (output_unit, *) " Number of initial state J symmetries to be transformed"
            read (*, *) jqj1
            if (jqj1 < 0 .or. jqj1 > size(itjqj1)) &
                call fail_input_io(0, -1, "Invalid number of selected J symmetries")
            write (output_unit, *) " Give the J symmetries in the form 2*J"
            read (*, *) (itjqj1(i),i=1,jqj1)
            itjqj1(:jqj1) = itjqj1(:jqj1) + 1
            write (output_unit, *) " Number of final state J symmetries to be transformed"
            read (*, *) jqj2
            if (jqj2 < 0 .or. jqj2 > size(itjqj2)) &
                call fail_input_io(0, -1, "Invalid number of selected J symmetries")
            write (output_unit, *) " Give the J symmetries in the form 2*J"
            read (*, *) (itjqj2(i),i=1,jqj2)
            itjqj2(:jqj2) = itjqj2(:jqj2) + 1
        end if

        return
    end subroutine read_state_filenames
end module read_state_filenames_interface
