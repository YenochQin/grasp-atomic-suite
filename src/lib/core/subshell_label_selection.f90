! Resolve a comma- or space-separated list of subshell labels and wildcards against the orbital basis.
module subshell_label_selection
    implicit none (type, external)

    ! Module interface
    private
    public :: select_subshell_labels
contains
    subroutine select_subshell_labels(matched_orbital_indices, matched_subshell_count)
        ! Module dependencies
        use calculation_file_units, only: calculation_log_unit, ensure_calculation_log
        use standard_io_units, only: stderr_unit
        implicit none (type, external)

        ! Arguments
        integer, intent(out) :: matched_subshell_count
        integer, intent(inout) :: matched_orbital_indices(*)

        ! Local declarations
        logical :: retry_input
        character(len=500) :: input_line

        retry_input = .false.
        do
            if (retry_input) write (stderr_unit, *) " redo ..."
            retry_input = .true.
            matched_subshell_count = 0
            read (*, "(A)") input_line

            call ensure_calculation_log()
            write(calculation_log_unit,"(a)") trim(input_line)
            if (parse_selection(input_line, matched_orbital_indices, matched_subshell_count)) exit
        end do
    end subroutine select_subshell_labels

    logical function parse_selection(input_line, matched_orbital_indices, matched_subshell_count) result(parsed)
        ! Module dependencies
        use orbital_csf_basis, only: orbital_count, orbital_principal_n, orbital_symmetry_labels
        use standard_io_units, only: stderr_unit
        use decimal_digit_check, only: is_decimal_digit
        use integer_format, only: format_integer
        implicit none (type, external)

        ! Arguments
        character(len=*), intent(in) :: input_line
        integer, intent(inout) :: matched_orbital_indices(*), matched_subshell_count

        ! Local declarations
        integer :: token_start, character_index, token_end, principal_digit_count, orbital_index, width_text_length, &
            & read_status, requested_principal_n, match_index
        logical :: subshell_matched, any_orbital_l, orbital_l_matches, any_principal_n, principal_n_matches, &
            & any_spin_component, spin_component_matches
        character :: requested_orbital_letter, requested_spin_suffix, current_character, width_text*2, integer_format*5

        parsed = .false.
        token_start = 0
        do character_index = 1, len(input_line)
            current_character = input_line(character_index:character_index)
            if (current_character/=" " .and. current_character/=",") then
                if (token_start == 0) token_start = character_index
            else
                if (token_start /= 0) then
                    token_end = character_index - 1
                    if (input_line(token_start:token_start) == "*") then
                        any_principal_n = .true.
                        token_start = min(token_start + 1,token_end)
                    else
                        any_principal_n = .false.
                        principal_digit_count = 0
                        do orbital_index = token_start, token_end
                            if (.not.is_decimal_digit(input_line(orbital_index:orbital_index))) cycle
                            principal_digit_count = principal_digit_count + 1
                        end do
                        if (principal_digit_count == 0) then
                            write (stderr_unit, *) "GETRSL: Unable to interpret ", &
                                "the principal quantum number;"
                            return
                        end if
                        call format_integer (principal_digit_count, width_text, width_text_length)
                        integer_format = "(1I"//width_text(1:width_text_length)//")"
                        read (input_line(token_start:token_start+principal_digit_count-1), integer_format, &
                            & iostat=read_status) requested_principal_n
                        if (read_status /= 0) then
                            write (stderr_unit, *) "GETRSL: Unable to interpret ", "string ", &
                                input_line(token_start:token_end-2), " as a principal quantum number"
                            return
                        end if
                        token_start = token_start + principal_digit_count
                    end if
                    requested_orbital_letter = input_line(token_start:token_start)
                    if (requested_orbital_letter == "*") then
                        any_orbital_l = .true.
                    else
                        any_orbital_l = .false.
                    end if
                    if (token_end > token_start) then
                        requested_spin_suffix = input_line(token_end:token_end)
                        if (requested_spin_suffix == "*") then
                            any_spin_component = .true.
                        else if (requested_spin_suffix == "-") then
                            any_spin_component = .false.
                        else
                            write (stderr_unit, *) "GETRSL: Unable to interpret ", "string ", &
                                requested_spin_suffix, " as a spin-orbit component indicator"
                            return
                        end if
                    else
                        if (any_orbital_l) then
                            any_spin_component = .true.
                        else
                            any_spin_component = .false.
                            requested_spin_suffix = " "
                        end if
                    end if
                    subshell_matched = .false.
                    do orbital_index = 1, orbital_count
                        principal_n_matches = any_principal_n .or. orbital_principal_n(orbital_index)==requested_principal_n
                        orbital_l_matches = any_orbital_l .or. &
                            requested_orbital_letter==orbital_symmetry_labels(orbital_index)(1:1)
                        spin_component_matches = any_spin_component .or. &
                            & requested_spin_suffix==orbital_symmetry_labels(orbital_index)(2:2)
                        if (.not.(principal_n_matches .and. orbital_l_matches .and. spin_component_matches)) cycle
                        do match_index = 1, matched_subshell_count
                            if (matched_orbital_indices(match_index) /= orbital_index) cycle
                            write (stderr_unit, *) "GETRSL: ", "Repeated subshell in list;"
                            return
                        end do
                        subshell_matched = .true.
                        matched_subshell_count = matched_subshell_count + 1
                        matched_orbital_indices(matched_subshell_count) = orbital_index
                    end do
                    if (.not.subshell_matched) then
                        write (stderr_unit, *) "GETRSL: Subshell not occupied as ", &
                            " according to CSL  File;"
                        return
                    end if
                    token_start = 0
                end if
            end if
        end do
        parsed = .true.
    end function parse_selection
end module subshell_label_selection
