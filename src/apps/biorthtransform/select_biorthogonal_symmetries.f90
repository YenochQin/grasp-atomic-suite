! Select which J symmetries to transform for the validated state pair.
module select_biorthogonal_symmetries_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: select_biorthogonal_symmetries
contains

    ! Written by Per Jonsson
    subroutine select_biorthogonal_symmetries()
        ! Module dependencies
        use input_error_handling, only: fail_input_io
        use, intrinsic :: iso_fortran_env, only: output_unit
        use biorthogonal_transition_selection, only: itjqj1, itjqj2, jqj1, jqj2, ntrans
        use yes_no_input, only: read_yes_no
        implicit none (type, external)

        ! Local declarations
        integer :: i
        logical :: yes
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
    end subroutine select_biorthogonal_symmetries
end module select_biorthogonal_symmetries_interface
