! Report the leading CSF contributions to the selected atomic states.
module write_csf_weights_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_csf_weights
contains

    ! Print  the  weights of the largest five CSF contributors to each
    ! ASF.
    !
    !                                        Last updated: 21 Dec 1992
    !                                        Last updated: 24 Feb 1997
    ! Modified by G. Gaigalas. 05 Feb 2017
    subroutine write_csf_weights(write_standard_output)
        ! Module dependencies
        use calculation_file_units, only: summary_output_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state, only: requested_level_indices, selected_level_count
        use eigenstate_expansions, only: csf_mixing_coefficients
        use angular_symmetry_labels, only: format_j_label, parity_labels
        use level_symmetry_metadata,          only: packed_block_symmetries
        use standard_io_units, only: stdout_unit
        use hamiltonian_symmetry_blocks, only: csf_counts_by_block
        use symmetry_block_offsets, only: csf_offsets_by_block, eigenvalue_offsets_by_block, mixing_offsets_by_block
        use transition_block_indices, only: level_block_indices
        implicit none (type, external)

        ! Arguments
        logical, intent(in) :: write_standard_output

        ! Local declarations
        integer :: dominant_csf_indices(5)
        integer :: selected_level_index, symmetry_block_index, block_csf_count, csf_offset, eigenvalue_offset, &
            block_mixing_offset, level_mixing_offset, reported_coefficient_count, level_number, two_j_plus_one, &
            parity_label_index, i, j, csf_index_swap_temp, parity_index
        real(real64) :: coefficient_swap_temp, candidate_coefficient
        real(real64) :: dominant_csf_coefficients(5)

        if (write_standard_output) then
            write (stdout_unit, '(/,"Weights of major contributors to ASF:"/,/,&
                &"Block Level J Parity      CSF contributions"/)')
        else
            write (summary_output_unit, &
                '(/,"Weights of major contributors to ASF:"/,/,"Block Level J Parity      CSF contributions"/)')
        end if

        do selected_level_index = 1, selected_level_count
            ! loop over eigenvectors

            symmetry_block_index = level_block_indices(selected_level_index)
            block_csf_count = csf_counts_by_block(symmetry_block_index)
            csf_offset = csf_offsets_by_block(symmetry_block_index)
            eigenvalue_offset = eigenvalue_offsets_by_block(symmetry_block_index)
            block_mixing_offset = mixing_offsets_by_block(symmetry_block_index)
            level_mixing_offset = block_mixing_offset + (selected_level_index - eigenvalue_offset - 1)*block_csf_count
            reported_coefficient_count = min(5,block_csf_count)                       !Find maximum 5 ... within block
            level_number = requested_level_indices(selected_level_index)
            two_j_plus_one = iabs(packed_block_symmetries(symmetry_block_index))
            if(packed_block_symmetries(symmetry_block_index) >= 0) then
                parity_label_index = 2
            else
                parity_label_index = 1
            end if
            ! IVTJPO = IATJPO(IV)                     ! j-value related
            ! IVSPAR = IASPAR(IV)                     ! parity related

            do i = 1, reported_coefficient_count
                dominant_csf_coefficients(i) = csf_mixing_coefficients(level_mixing_offset + i)
                dominant_csf_indices(i) = i
            end do

            ! sort the first nelt in decreasing order
            do i = 1, reported_coefficient_count
                do j = i + 1, reported_coefficient_count
                    if (abs(dominant_csf_coefficients(j)) <= abs(dominant_csf_coefficients(i))) cycle
                    coefficient_swap_temp = dominant_csf_coefficients(i)
                    dominant_csf_coefficients(i) = dominant_csf_coefficients(j)
                    dominant_csf_coefficients(j) = coefficient_swap_temp
                    csf_index_swap_temp = dominant_csf_indices(i)
                    dominant_csf_indices(i) = dominant_csf_indices(j)
                    dominant_csf_indices(j) = csf_index_swap_temp
                end do
            end do

            do i = reported_coefficient_count + 1, block_csf_count
                candidate_coefficient = csf_mixing_coefficients(level_mixing_offset + i)
                if (candidate_coefficient==0.0_real64 .or. &
                    & abs(candidate_coefficient)<=abs(dominant_csf_coefficients(reported_coefficient_count))) cycle
                ! we have a non-zero value larger than the largest so far
                do j = 1, reported_coefficient_count
                    if (abs(candidate_coefficient) <= abs(dominant_csf_coefficients(j))) cycle
                    dominant_csf_coefficients(reported_coefficient_count:1+j:(-1)) = &
                        & dominant_csf_coefficients(reported_coefficient_count-1:j:(-1))
                    dominant_csf_indices(reported_coefficient_count:1+j:(-1)) = &
                        & dominant_csf_indices(reported_coefficient_count-1:j:(-1))
                    dominant_csf_coefficients(j) = candidate_coefficient
                    dominant_csf_indices(j) = i
                    exit
                end do
            end do

            ! IP = (IASPAR(IV) + 3)/2
            parity_index =  parity_label_index

            if (write_standard_output) then
                write (stdout_unit, '(i3,1x,i5,2x,2a,5(3x,f8.4))') symmetry_block_index, level_number, &
                    & format_j_label(two_j_plus_one), &
                    & parity_labels(parity_index), (dominant_csf_coefficients(i),i&
                    =1,reported_coefficient_count)
                write (stdout_unit, '(19x,5(3x,i8))') (dominant_csf_indices(i),i=1,reported_coefficient_count)
            else

                write (summary_output_unit, '(i3,1x,i5,2x,2a,5(3x,f8.4))') symmetry_block_index, level_number, &
                    & format_j_label(two_j_plus_one), parity_labels(parity_index), (dominant_csf_coefficients(i),i=1,&
                    reported_coefficient_count)
                write (summary_output_unit, '(19x,5(3x,i8))') (dominant_csf_indices(i),i=1,reported_coefficient_count)
            end if
        end do

        return
    end subroutine write_csf_weights
end module write_csf_weights_interface
