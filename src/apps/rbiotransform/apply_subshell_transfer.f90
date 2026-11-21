module apply_subshell_transfer_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: apply_subshell_transfer
contains

    !***********************************************************************
    !                                                                      *
    subroutine apply_subshell_transfer(input_mixing_coefficients,csf_count,level_count,subshell_index,symmetry_index, &
                & transformation_column,subshell_count, &
        & transformed_mixing_coefficients,global_trace_level)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        !                                                                      *
        !   Calculate the action of the  operator                              *
        !   Sum( j) T(jl) E(jl,il) on a set of vectors                         *
        !                                                                      *
        !   only j<i terms ! included
        !                                                                      *
        !   Modification of the TI1TV program for MCHF                         *
        !                                                                      *
        !   Written by Per Jonsson                                             *
        !                                                                      *
        !   =====                                                              *
        !   Input                                                              *
        !   =====                                                              *
        !   CIIN  : Input CI vectors                                           *
        !   NCSCF : Length of CI expansion                                     *
        !   NCIV  : Number of CI vectors                                       *
        !   I     : Shell number                                               *
        !   L     : L value of shells of this excitation                       *
        !   NSHL  : Number of shells with this L                               *
        !   T     : The integrals T(*,il)                                      *
        !                                                                      *
        !   ======                                                             *
        !   Output                                                             *
        !   ======                                                             *
        !                                                                      *
        !   CIOUT : List of output CI vectors                                  *
        !                                                                      *
        !***********************************************************************
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: orbital_index_radix
        use mcp_coefficient_storage, only: mcp_coefficient_values, mcp_packed_orbital_pairs, mcp_coefficient_end_offsets, &
            mcp_bra_csf_indices, mcp_ket_csf_indices, mcp_integral_count
        use biorthogonal_subshell_mapping, only: active_orbitals_by_group, active_orbital_positions_by_group

        ! Procedure interfaces
        use fill_vector_interface, only: fill_vector
        use write_matrix_interface, only: write_matrix
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: csf_count, level_count, subshell_count, subshell_index, symmetry_index, global_trace_level
        real(real64), intent(in) :: input_mixing_coefficients(csf_count, level_count)
        real(real64), intent(in) :: transformation_column(subshell_count)
        real(real64), intent(out) :: transformed_mixing_coefficients(csf_count, level_count)

        ! Local constants

        ! Local declarations
        integer :: local_trace_level, trace_level, coefficient_count, packed_pair_index, first_orbital_index, &
                & second_orbital_index, target_subshell_position, first_subshell_position, second_subshell_position, &
                & first_coefficient, &
            source_subshell, bra_csf_index, ket_csf_index, level_index, coefficient_index
        real(real64) :: angular_coefficient, transfer_weight
        !
        local_trace_level = 00
        trace_level = max(local_trace_level,global_trace_level)
        trace_level = 000
        !
        if(trace_level>=10) write(output_unit,*) " Entering TI1TV"

        call fill_vector(transformed_mixing_coefficients,0.0_real64,level_count*csf_count)
        !
        !.  Obtain address of first coupling coefficient for h(jl,il) :  IFIRST
        !.  Obtain number of        coupling coefficient for h(jl,il) :  NFOUND
        !
        !        NFOUND : Number of coefficients obtained
        !
        coefficient_count = 0
        do packed_pair_index = 1,mcp_integral_count
            first_orbital_index = mcp_packed_orbital_pairs(packed_pair_index)/orbital_index_radix
            second_orbital_index = mod(mcp_packed_orbital_pairs(packed_pair_index),orbital_index_radix)
            !
            !   Determine the ordernumber of the orbitals I, IA and IB within the
            !   shell with symmetry L
            !
            target_subshell_position = active_orbital_positions_by_group(symmetry_index, &
                active_orbitals_by_group(symmetry_index,subshell_index))
            first_subshell_position = active_orbital_positions_by_group(symmetry_index,first_orbital_index)
            second_subshell_position = active_orbital_positions_by_group(symmetry_index,second_orbital_index)

            if (target_subshell_position==second_subshell_position.and.first_subshell_position<second_subshell_position) then
                if (packed_pair_index==1) then
                    first_coefficient = 1
                else
                    first_coefficient = mcp_coefficient_end_offsets(packed_pair_index-1) + 1
                end if
                coefficient_count = mcp_coefficient_end_offsets(packed_pair_index) - first_coefficient + 1

                do coefficient_index = 1, coefficient_count
                    angular_coefficient = mcp_coefficient_values(first_coefficient-1+coefficient_index)
                    source_subshell = first_subshell_position
                    transfer_weight = transformation_column(source_subshell)*angular_coefficient
                    !ww            IF (NORDII.EQ.0.AND.NORDFF.EQ.0) THEN
                    !
                    !  Expression for normal orbital ordering
                    !
                    bra_csf_index = mcp_bra_csf_indices(first_coefficient-1+coefficient_index)
                    ket_csf_index = mcp_ket_csf_indices(first_coefficient-1+coefficient_index)
                    !ww            ELSEIF (NORDII.EQ.1.AND.NORDFF.EQ.1) THEN
                    !
                    !  Expression for reverse orbital ordering
                    !
                    !ww              IRIGHT = JANN(IFIRST-1+IELMNT)
                    !ww              ILEFT = JBNN(IFIRST-1+IELMNT)
                    !ww            ELSE
                    !ww              WRITE(*,*) 'SOMETHING WRONG'
                    !ww              STOP
                    !ww            ENDIF
                    do level_index =1, level_count
                        transformed_mixing_coefficients(bra_csf_index,level_index) =                                       &
                            &        transformed_mixing_coefficients(bra_csf_index,level_index) + transfer_weight* &
                & input_mixing_coefficients(ket_csf_index,level_index)
                    end do
                end do
            end if
        end do
        !
        if(trace_level>=100) then
            write(output_unit,*)
            write(output_unit,*) " Input and output vectors from apply_subshell_transfer I,L",subshell_index,symmetry_index
            call write_matrix(input_mixing_coefficients,csf_count,level_count,csf_count,level_count)
            write(output_unit,*)
            call write_matrix(transformed_mixing_coefficients,csf_count,level_count,csf_count,level_count)
            write(output_unit,*)
        end if
        !
        if(trace_level>=10) write(output_unit,*) " LEAVING TI1TV"

        return
    end subroutine apply_subshell_transfer
end module apply_subshell_transfer_interface
