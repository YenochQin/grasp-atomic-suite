! Open the cached angular coefficients for transition matrix elements.
module open_transition_angular_cache_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: open_transition_angular_cache
contains

    ! Checks if the angular file name(1).name(2).T is available
    ! and appropriate
    !
    ! Written by Per Jonsson. 6 March 1997
    subroutine open_transition_angular_cache(name, avail, jkp, transition_angular_unit)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use orbital_csf_basis, only: orbital_count
        use transition_matrix_workspace, only: requested_multipole_count, requested_multipole_codes
        use mpi_state, only: mpi_process_count, mpi_rank
        implicit none (type, external)

        ! Arguments
        integer  :: jkp
        integer, intent(out) :: transition_angular_unit
        logical, intent(out) :: avail
        character, intent(inout) :: name(2)*24

        ! Local declarations
        integer :: j1,j2,initial_block_index,final_block_index,nwd,nkpd,nprocsd,myidd
        logical :: found
        character :: s(-9:9)*2

        s((-9)) = "-9"
        s((-8)) = "-8"
        s((-7)) = "-7"
        s((-6)) = "-6"
        s((-5)) = "-5"
        s((-4)) = "-4"
        s((-3)) = "-3"
        s((-2)) = "-2"
        s((-1)) = "-1"
        s(0) = "+0"
        s(1) = "+1"
        s(2) = "+2"
        s(3) = "+3"
        s(4) = "+4"
        s(5) = "+5"
        s(6) = "+6"
        s(7) = "+7"
        s(8) = "+8"
        s(9) = "+9"

        j1 = index(name(1)," ")
        j2 = index(name(2)," ")
        inquire(file=name(1)(1:j1-1)//"."//name(2)(1:j2-1)//"."//s(requested_multipole_codes(jkp))//"T"&
            , exist=found)
        if (.not.found) then
            write (output_unit, *)
            write (output_unit, *) " Angular file not available"
            avail = .false.
            return
        else
            ! Open the file and check if it is appropriate for the present case
            open(newunit=transition_angular_unit, file=name(1)(1:j1-1)//"."//name(2)(1:j2-1)//"."//s(&
                requested_multipole_codes(jkp))//"T", status="OLD", form="UNFORMATTED", position="asis")
            rewind (transition_angular_unit)
            read (transition_angular_unit) initial_block_index,final_block_index,nwd,nkpd,nprocsd,myidd
            if (.not.(nwd==orbital_count .and. nkpd==requested_multipole_count                          &
                .and.nprocsd==mpi_process_count.and.myidd==mpi_rank)) then
                write (output_unit, *) " Angular file not appropriate"
                avail = .false.
                close(transition_angular_unit, status="DELETE")
                return
            else
                rewind (transition_angular_unit)
                write (output_unit, *) " Angular data read from file"
                avail = .true.
            end if
        end if
        return
    end subroutine open_transition_angular_cache
end module open_transition_angular_cache_interface
