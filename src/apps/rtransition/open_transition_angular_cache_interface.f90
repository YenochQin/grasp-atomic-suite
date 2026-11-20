module open_transition_angular_cache_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine open_transition_angular_cache (name, avail, jkp, transition_angular_unit)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(inout) :: name(2)
            logical, intent(out) :: avail
            integer :: jkp
            integer, intent(out) :: transition_angular_unit
        end subroutine open_transition_angular_cache
    end interface
end module open_transition_angular_cache_interface
