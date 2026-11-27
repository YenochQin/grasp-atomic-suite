! Open the angular record stream for one multipole in the local working directory.
module open_transition_angular_cache_interface
    use transition_angular_cache, only: transition_cache_key, open_cache
    implicit none (type, external)
    private
    public :: open_transition_angular_cache
contains
    subroutine open_transition_angular_cache(name, available, index, unit, key)
        use transition_matrix_workspace, only: requested_multipole_codes
        character(*), intent(in) :: name(2)
        integer, intent(in) :: index
        integer, intent(out) :: unit
        logical, intent(out) :: available
        type(transition_cache_key), intent(out) :: key

        call open_cache(name, requested_multipole_codes(index), 1, -1, unit, key, available)
    end subroutine open_transition_angular_cache
end module open_transition_angular_cache_interface
