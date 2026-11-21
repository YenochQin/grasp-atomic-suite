! Parse individual level indices and inclusive ranges into the shared selection list.
module level_selection_parser
    implicit none (type, external)

    ! Module interface
    private
    public :: parse_level_selection
contains
    subroutine parse_level_selection(selected_level_count, maximum_level_index, input_line, parse_status)
        ! Module dependencies
        use array_memory_management, only: alloc, ralloc
        use atomic_calculation_state, only: requested_level_indices
        use standard_io_units, only: stderr_unit
        use integer_format, only: format_integer
        implicit none (type, external)

        ! Arguments
        integer, intent(inout) :: selected_level_count
        integer, intent(in) :: maximum_level_index
        integer, intent(out) :: parse_status
        character (len = *), intent(in) :: input_line

        ! Local declarations
        integer :: length_of_record, initial_level_count, allocated_level_capacity, range_start_position, token_start, &
            & character_index, token_end, &
            token_width, width_text_length, read_status, parsed_level_index, range_start_level, additional_level_count, &
            & range_offset
        character :: current_character, integer_format*7, width_text*3

        length_of_record = len_trim(input_line)

        initial_level_count = selected_level_count
        allocated_level_capacity = selected_level_count

        if (allocated_level_capacity == 0) then
            allocated_level_capacity = 1
            call alloc (requested_level_indices, allocated_level_capacity, "requested_level_indices", "parse_level_selection")
        end if
        range_start_position = 0
        character_index = 1
        do while (character_index <= length_of_record)
            current_character = input_line(character_index:character_index)
            if (current_character == " " .or. current_character == ",") then
                character_index = character_index + 1
                cycle
            end if

            token_start = character_index
            do while (character_index <= length_of_record)
                current_character = input_line(character_index:character_index)
                if (current_character == " " .or. current_character == "," .or. current_character == "-") exit
                character_index = character_index + 1
            end do
            token_end = character_index - 1
            token_width = token_end - token_start + 1
            call format_integer (token_width, width_text, width_text_length)
            integer_format = "(1i"//width_text(1:width_text_length)//")"
            read (input_line(token_start:token_end), integer_format, iostat=read_status) parsed_level_index
            if (read_status /= 0) then
                write (stderr_unit, *) "parse_level_selection: unable to decode "//input_line(token_start:token_end)//";"
                parse_status = -1
                return
            end if

            if (range_start_position == 0) then
                selected_level_count = selected_level_count + 1
                if (selected_level_count > allocated_level_capacity) then
                    call ralloc (requested_level_indices, selected_level_count, "requested_level_indices", "parse_level_selection")
                    allocated_level_capacity = selected_level_count
                end if

                if (parsed_level_index<1 .or. parsed_level_index>maximum_level_index) then
                    write (stderr_unit, *) "parse_level_selection: serial numbers must be", &
                        " in the range [1,", maximum_level_index, "];"
                    parse_status = -2
                    return
                end if

                requested_level_indices(selected_level_count) = parsed_level_index
                character_index = character_index + 1
                if (current_character == "-") range_start_position = selected_level_count
            else
                range_start_level = requested_level_indices(selected_level_count)
                additional_level_count = parsed_level_index - range_start_level

                if (additional_level_count < 0) then
                    write (stderr_unit, *) range_start_level, "-", parsed_level_index, " not allowed"
                    parse_status = -3
                    return
                end if

                selected_level_count = selected_level_count + additional_level_count
                if (selected_level_count > allocated_level_capacity) then
                    call ralloc (requested_level_indices, selected_level_count, "requested_level_indices", "parse_level_selection")
                    allocated_level_capacity = selected_level_count
                end if
                do range_offset = 1, additional_level_count
                    requested_level_indices(range_start_position+range_offset) = range_start_level + range_offset
                end do
                character_index = character_index + 1
                range_start_position = 0
            end if

        end do

        if (selected_level_count == initial_level_count) then
            parse_status = 1
            return
        end if
        if (selected_level_count /= allocated_level_capacity) then
            write (stderr_unit, *) "parse_level_selection: selected level count differs from allocated capacity"
            call ralloc (requested_level_indices, selected_level_count, "requested_level_indices", "parse_level_selection")
        end if

        parse_status = 0
        return
    end subroutine parse_level_selection
end module level_selection_parser
