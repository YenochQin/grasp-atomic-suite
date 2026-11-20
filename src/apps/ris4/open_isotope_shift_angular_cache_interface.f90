module open_isotope_shift_angular_cache_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Gediminas Gaigalas 11/18/19
        subroutine open_isotope_shift_angular_cache (name,avail,whichone,status,cache_unit)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(in) :: name
            logical, intent(out) :: avail
            integer, intent(in)  :: whichone
            integer, intent(out) :: status
            integer, intent(out) :: cache_unit
        end subroutine open_isotope_shift_angular_cache
    end interface
end module open_isotope_shift_angular_cache_interface
