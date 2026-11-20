!***********************************************************************
!                                                                      *
subroutine open_biorthogonal_angular_cache(name, avail, kamax)
    ! Module dependencies
    use biorthogonal_subshell_layout, only: biorthogonal_angular_unit
    use, intrinsic :: iso_fortran_env, only: output_unit
    !                                                                      *
    !   Checks if the angular file name.T is available and appropriate     *
    !                                                                      *
    !   Written by Per Jonsson                      6 March 1997           *
    !                                                                      *
    !***********************************************************************
    use orbital_csf_basis,           only: csf_count, orbital_count
    implicit none (type, external)

    ! Arguments
    integer, intent(in)      :: kamax
    logical, intent(out)     :: avail
    character, intent(inout) :: name*24

    ! Local declarations
    integer :: j, ncfd, nwd, kamaxd
    logical :: found
    !-----------------------------------------------
    !
    j = index(name," ")
    inquire(file=name(1:j-1)//".TB", exist=found)
    if (.not.found) then
        write (output_unit, *) " Angular file not available"
        avail = .false.
        return
    else
        !
        !  Open the file and check if it is appropriate for the present case
        !
        open(newunit=biorthogonal_angular_unit, file=name(1:j-1)//".TB", status="OLD", form=&
            "UNFORMATTED", position="asis")
        rewind (biorthogonal_angular_unit)
        read (biorthogonal_angular_unit) ncfd, nwd, kamaxd
        if (.not.(ncfd==csf_count .and. nwd==orbital_count .and. kamaxd==kamax)) then
            write (output_unit, *) " Angular file not appropriate"
            avail = .false.
            return
        else
            write (output_unit, *) " Angular data read from file"
            avail = .true.
        end if
    end if
    return
end subroutine open_biorthogonal_angular_cache
