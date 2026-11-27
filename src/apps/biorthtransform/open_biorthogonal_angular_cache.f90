! Open the saved angular coefficients for the biorthogonal transformation when they can be reused.
module open_biorthogonal_angular_cache_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: open_biorthogonal_angular_cache
contains

    subroutine open_biorthogonal_angular_cache(name, side, key, available)
        use biorthogonal_angular_cache, only: angular_cache_key, current_cache_key, open_cache_file
        use biorthogonal_subshell_layout, only: biorthogonal_angular_unit
        use calculation_input_defaults, only: angular_output_dump_flag
        character(*), intent(in) :: name
        integer, intent(in) :: side
        type(angular_cache_key), intent(out) :: key
        logical, intent(out) :: available

        key = current_cache_key(side, 1, 0)
        call open_cache_file(trim(name)//".biorth.angular.cache", key, angular_output_dump_flag == 1, &
            biorthogonal_angular_unit, available)
    end subroutine open_biorthogonal_angular_cache
end module open_biorthogonal_angular_cache_interface
