!***********************************************************************
!                                                                      *
subroutine write_input_summary
    ! Module dependencies
    use calculation_file_units, only: ensure_summary_output, summary_output_unit
    !                                                                      *
    !***********************************************************************
    use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
    use atomic_calculation_state, only: electron_count
    use csf_block_cutoffs, only: iccut, iccutblk
    use mcp_generation_options, only: diagonal_only, max_coulomb_multipole_rank, use_first_order_restriction
    use orbital_csf_basis, only: csf_count, orbital_count
    use calendar_clock, only: calen
    use integer_format, only: convrt
    implicit none (type, external)

    ! Local declarations
    integer :: text_length, symmetry_block_index
    character :: time_text*10, date_text*8
    character :: record*256, cdata*26
    !-----------------------------------------------
    !
    !
    !   Get the date and time of day; make this information the
    !   header of the summary file
    !
    call calen (time_text, date_text)
    call ensure_summary_output()
    write (summary_output_unit, *) "GENMCP run at ", time_text, " on ", date_text, "."
    !
    !   Write out the basic dimensions of the electron cloud
    !
    write (summary_output_unit, *)
    call convrt (electron_count, record, text_length)
    write (summary_output_unit, *) "There are "//record(1:text_length)//" electrons in the cloud"
    call convrt (csf_count, record, text_length)
    write (summary_output_unit, *) " in "//record(1:text_length)//" relativistic CSFs"
    call convrt (orbital_count, record, text_length)
    write (summary_output_unit, *) " based on "//record(1:text_length)//" relativistic subshells."
    !
    !   If the CSFs are not treated uniformly, write out an
    !   informative message
    !
    if (diagonal_only) then
        write (summary_output_unit, *)
        write (summary_output_unit, *) "Only diagonal matrix elements are computed."
    else
        if (use_first_order_restriction) then
            do symmetry_block_index = 1, symmetry_block_count
                write (summary_output_unit, *)
                call convrt (iccut(symmetry_block_index), record, text_length)
                write (summary_output_unit, *) "CSFs 1--"//record(1:text_length)//" constitute"//&
                    " the zero-order space."
            end do
        end if
    end if
    !
    return
end subroutine write_input_summary
