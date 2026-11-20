module open_biorthogonal_angular_cache_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine open_biorthogonal_angular_cache (name, avail, kamax)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(inout) :: name
            logical, intent(out) :: avail
            integer, intent(in) :: kamax
        end subroutine open_biorthogonal_angular_cache
    end interface
end module open_biorthogonal_angular_cache_interface
