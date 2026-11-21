!***********************************************************************
!                                                                      *
subroutine write_ci_eigenvectors
    ! Module dependencies
    use calculation_file_units, only: summary_output_unit
    !                                                                      *
    !   Routine for printing eigenvectors.                                 *
    !                                                                      *
    !   Call(s) to: [LIB92]: ALLOC, DALLOC.                                *
    !                                                                      *
    !   Written by Farid A Parpia             Last revision: 06 Oct 1992   *
    !                                                                      *
    !***********************************************************************
    use array_memory_management, only: alloc, dalloc
    use eigenstate_expansions, only: csf_mixing_coefficients
    use orbital_csf_basis, only: csf_count
    use selected_level_indices, only: active_level_indices, active_level_count
    implicit none (type, external)

    ! Local declarations
    integer :: i, j1, j2, nvex, jp, l, csf_index
    integer, pointer :: iloc(:)
    !-----------------------------------------------
    !
    !   Header
    !
    write (summary_output_unit,'(//" Eigenvectors:"/)')
    !
    !   Determine which eigenvectors are to be printed out
    !
    nvex = active_level_count
    call alloc (iloc,active_level_count,"ILOC", "write_ci_eigenvectors")
    do i = 1,active_level_count
        iloc(i) = active_level_indices(i)
    end do
    !
    !   There are eight columns across the width of a page
    !
    jp = 8
    !
    !   Set up for the first set of columns
    !
    j1 = 1
    j2 = jp
    !
    !   Loop over sets of columns
    !
    do
        if (j2 > nvex) j2 = nvex
        !
        !   Write out the column numbers; skip a line
        !
        write (summary_output_unit,'(1x,8(i8,7x))') (iloc(l),l = j1,j2)
        write (summary_output_unit,'(1x)')
        !
        !   Write out the rows
        !
        do csf_index = 1,csf_count
            !
            write (summary_output_unit,'(1x,1p,8d15.7)') (csf_mixing_coefficients(csf_index+(l-1)*csf_count),l = j1,j2)
            !
            !   Skip a line after every tenth row
            !
            if (mod(csf_index,10) == 0) write (summary_output_unit,'(1x)')
            !
        end do
        !
        !   Set up for the next set of columns
        !
        write (summary_output_unit,'(1x,1p,8d15.7)')
        if (j2 >= nvex) exit
        j1 = j1+jp
        j2 = j2+jp
    end do
    !
    !   Deallocate storage for ILOC
    !
    call dalloc (iloc, "ILOC", "write_ci_eigenvectors")
    !
    return
    !
    !
end subroutine write_ci_eigenvectors
