! Identify completed biorthogonal angular-coefficient files in the rank's working directory.
module biorthogonal_angular_cache
    use csf_descriptor_reader, only: loaded_csf_basis_id
    use biorthogonal_subshell_layout, only: kappa_by_group, kappa_group_count, max_biorthogonal_kappa_groups
    use biorthogonal_transition_selection, only: ntrans, itjqj1, itjqj2, jqj1, jqj2
    use orbital_csf_basis, only: csf_count, orbital_count
    implicit none (type, external)
    private
    public :: angular_cache_key, current_cache_key, open_cache_file, finish_cache_file

    type :: angular_cache_key
        character(8) :: signature = "BIOTB001"
        character(64) :: basis
        integer :: dimensions(4), groups(max_biorthogonal_kappa_groups), selection(size(itjqj1)), selection_count
    end type

contains
    function current_cache_key(side, ranks, rank) result(key)
        integer, intent(in) :: side, ranks, rank
        type(angular_cache_key) :: key

        key%basis = loaded_csf_basis_id
        key%dimensions = [csf_count, orbital_count, ranks, rank]
        key%groups = kappa_by_group
        key%selection = 0
        key%selection_count = -1
        if (ntrans == 0) return
        if (side == 1) then
            key%selection_count = jqj1
            key%selection(:jqj1) = itjqj1(:jqj1)
        else
            key%selection_count = jqj2
            key%selection(:jqj2) = itjqj2(:jqj2)
        end if
    end function current_cache_key

    subroutine open_cache_file(path, key, save_cache, unit, available)
        character(*), intent(in) :: path
        type(angular_cache_key), intent(in) :: key
        logical, intent(in) :: save_cache
        integer, intent(out) :: unit
        logical, intent(out) :: available
        type(angular_cache_key) :: stored
        logical :: exists
        integer :: status, header(3)

        available = .false.
        inquire(file=path, exist=exists)
        if (exists) then
            open(newunit=unit,file=path,status="old",action="read",form="unformatted",position="append")
            backspace(unit,iostat=status)
            if (status == 0) read(unit,iostat=status) stored
            if (status == 0) then
                available = stored%signature == key%signature .and. stored%basis == key%basis .and. &
                    all(stored%dimensions == key%dimensions) .and. all(stored%groups == key%groups) .and. &
                    stored%selection_count == key%selection_count .and. all(stored%selection == key%selection)
            end if
            if (available) then
                rewind(unit)
                read(unit,iostat=status) header
                available = status == 0
                if (available) available = all(header == [key%dimensions(:2), kappa_group_count])
            end if
            if (available) then
                rewind(unit)
                write(*,*) "Angular data read from file: ", path
                return
            end if
            close(unit)
            write(*,*) "Rebuilding incompatible or incomplete angular cache: ", path
        end if
        if (save_cache) then
            open(newunit=unit,file=path,status="replace",form="unformatted")
        else
            open(newunit=unit,status="scratch",form="unformatted")
        end if
    end subroutine open_cache_file

    subroutine finish_cache_file(unit, key)
        integer, intent(in) :: unit
        type(angular_cache_key), intent(in) :: key
        write(unit) key
    end subroutine finish_cache_file
end module biorthogonal_angular_cache
