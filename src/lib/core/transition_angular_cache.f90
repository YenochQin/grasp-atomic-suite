! Cache one radiative multipole for an ordered pair of CSF bases.
module transition_angular_cache
    use calculation_input_defaults, only: angular_output_dump_flag
    use orbital_csf_basis, only: orbital_count
    use transition_basis_blocks, only: transition_csf_basis_ids
    use transition_matrix_workspace, only: requested_multipole_count
    use transition_off_diagonal_flags, only: different_j_only
    implicit none (type, external)
    private
    public :: transition_cache_key, open_cache

    type :: transition_cache_key
        character(8) :: version = "RADKT001"
        character(64) :: bases(2)
        integer :: multipole, different_j, ranks, rank
    end type

contains
    subroutine open_cache(names, multipole, ranks, rank, unit, key, available)
        character(*), intent(in) :: names(2)
        integer, intent(in) :: multipole, ranks, rank
        integer, intent(out) :: unit
        type(transition_cache_key), intent(out) :: key
        logical, intent(out) :: available
        type(transition_cache_key) :: stored
        character(16) :: initial_name_length, multipole_rank
        character :: multipole_kind
        character(:), allocatable :: path
        integer :: status, header(6), fields
        logical :: exists

        key%bases = transition_csf_basis_ids
        key%multipole = multipole
        key%different_j = merge(1, 0, different_j_only(multipole))
        key%ranks = ranks
        key%rank = rank
        ! Keep the length-prefixed cache path for existing caches.
        write(initial_name_length, "(I0)") len_trim(names(1))
        write(multipole_rank, "(I0)") abs(multipole)
        multipole_kind = merge("e", "m", sign(1, multipole) == (-1)**abs(multipole))
        path = trim(initial_name_length)//"."//trim(names(1))//"."//trim(names(2))// &
            ".rtransition."//multipole_kind//trim(multipole_rank)//".angular.cache"
        inquire(file=path, exist=exists)
        available = .false.
        if (exists) then
            open(newunit=unit, file=path, status="old", action="read", form="unformatted", position="append")
            backspace(unit, iostat=status)
            if (status == 0) read(unit, iostat=status) stored
            if (status == 0) then
                available = stored%version == key%version .and. all(stored%bases == key%bases) .and. &
                    stored%multipole == multipole .and. stored%different_j == key%different_j .and. &
                    stored%ranks == ranks .and. stored%rank == rank
            end if
            if (available) then
                ! MPI block headers contain two additional fields, including in single-rank runs.
                fields = merge(6, 4, rank >= 0)
                rewind(unit)
                read(unit, iostat=status) header(:fields)
                available = status == 0
                if (available) available = all(header(:4) == [1, 1, orbital_count, requested_multipole_count])
                if (available .and. rank >= 0) available = all(header(5:6) == [ranks, rank])
            end if
            if (available) then
                rewind(unit)
                write(*,*) "Angular data read from file: ", path
                return
            end if
            close(unit)
            write(*,*) "Rebuilding incompatible or incomplete transition cache: ", path
        end if
        if (angular_output_dump_flag == 1) then
            open(newunit=unit, file=path, status="replace", form="unformatted")
        else
            open(newunit=unit, status="scratch", form="unformatted")
        end if
    end subroutine open_cache
end module transition_angular_cache
