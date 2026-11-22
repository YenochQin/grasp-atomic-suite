! Open the saved angular coefficients for the biorthogonal transformation when they can be reused.
module open_biorthogonal_angular_cache_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: open_biorthogonal_angular_cache
contains

    ! Checks if the angular file name.T is available and appropriate
    !
    ! Written by Per Jonsson. 6 March 1997
    subroutine open_biorthogonal_angular_cache(name, matching_angular_cache_available, kappa_group_count)
        ! Module dependencies
        use biorthogonal_subshell_layout, only: biorthogonal_angular_unit
        use, intrinsic :: iso_fortran_env, only: output_unit
        use orbital_csf_basis,           only: csf_count, orbital_count
        implicit none (type, external)

        ! Arguments
        integer, intent(in)      :: kappa_group_count
        logical, intent(out)     :: matching_angular_cache_available
        character, intent(inout) :: name*24

        ! Local declarations
        integer :: j, cached_csf_count, cached_orbital_count, cached_kappa_group_count
        logical :: found

        j = index(name," ")
        inquire(file=name(1:j-1)//".TB", exist=found)
        if (.not.found) then
            write (output_unit, *) " Angular file not available"
            matching_angular_cache_available = .false.
            return
        else
            ! Open the file and check if it is appropriate for the present case
            open(newunit=biorthogonal_angular_unit, file=name(1:j-1)//".TB", status="OLD", form=&
                "UNFORMATTED", position="asis")
            rewind (biorthogonal_angular_unit)
            read (biorthogonal_angular_unit) cached_csf_count, cached_orbital_count, cached_kappa_group_count
            if (.not.(cached_csf_count==csf_count .and. cached_orbital_count==orbital_count .and. cached_kappa_group_count== &
                & kappa_group_count)) then
                write (output_unit, *) " Angular file not appropriate"
                matching_angular_cache_available = .false.
                return
            else
                write (output_unit, *) " Angular data read from file"
                matching_angular_cache_available = .true.
            end if
        end if
        return
    end subroutine open_biorthogonal_angular_cache
end module open_biorthogonal_angular_cache_interface
