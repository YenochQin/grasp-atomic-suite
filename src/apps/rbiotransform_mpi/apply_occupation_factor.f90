module apply_occupation_factor_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: apply_occupation_factor
contains

    !***********************************************************************
    !                                                                      *
    subroutine apply_occupation_factor(input_mixing_coefficients, csf_count, level_count, subshell_index, symmetry_index, &
                & diagonal_factor, &
        & transformed_mixing_coefficients, global_trace_level)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        !                                                                      *
        !   Calculates the action of the operator                              *
        !   Const ** E(li,li) on a set of vectors                              *
        !                                                                      *
        !   Adapted for GRASP, daughter of TIINI, born February 1996           *
        !                                                                      *
        !   =====                                                              *
        !   Input                                                              *
        !   =====                                                              *
        !   CIIN      : Input CI vectors                                       *
        !   NCSCF     : Length of CI expansion                                 *
        !   NCIV      : Number of CI vectors                                   *
        !   I         : Shell number                                           *
        !   L         : symmetry                                               *
        !   NSHLP(L,K): Gives the shell number as defined in getcsl for        *
        !               the K:th shell with symmetry L                         *
        !   CONST     : The constant                                           *
        !                                                                      *
        !   ======                                                             *
        !   Output                                                             *
        !   ======                                                             *
        !                                                                      *
        !   CIOUT : List of output CI vectors                                  *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: orbital_index_radix
        use mcp_coefficient_storage, only: mcp_coefficient_values, mcp_packed_orbital_pairs, mcp_coefficient_end_offsets, &
            mcp_bra_csf_indices, mcp_coefficient_count, mcp_integral_count
        use biorthogonal_subshell_mapping, only: active_orbitals_by_group
        use mpi_state, only: mpi_allreduce, mpi_comm_world, mpi_double_precision, mpi_error_code, mpi_sum

        ! Procedure interfaces
        use fill_vector_interface, only: fill_vector
        use write_matrix_interface, only: write_matrix
        use copy_vector_interface, only: copy_vector
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: csf_count
        integer, intent(in) :: level_count
        integer, intent(in) :: subshell_index
        integer, intent(in) :: symmetry_index
        integer, intent(in) :: global_trace_level
        real(real64), intent(in) :: diagonal_factor
        real(real64), intent(in)    :: input_mixing_coefficients(csf_count, level_count)
        real(real64), intent(inout) :: transformed_mixing_coefficients(csf_count, level_count)

        ! Local constants

        ! Local declarations
        integer :: local_trace_level, trace_level, coefficient_count, packed_pair_index, first_orbital_index, &
            second_orbital_index, first_coefficient, coefficient_index, occupation_number, bra_csf_index
        real(real64) :: ciouttmp(csf_count, level_count)
        real(real64) :: occupation_factor
        !-----------------------------------------------
        !
        !
        local_trace_level = 0
        trace_level = max(local_trace_level,global_trace_level)

        if (trace_level >= 10) write (output_unit, *) " Entering TIINI"

        call fill_vector (transformed_mixing_coefficients, 0.0_real64, csf_count*level_count)
        if (mcp_coefficient_count /= 0) then
            !
            !.  Obtain address of first coupling coefficient for h(il,il) :  IFIRST
            !.  Obtain number of        coupling coefficient for h(il,il) :  NFOUND
            !.
            !     NFOUND : Number of coefficients obtained
            !     IVAL   : actual RACAH coefficient <CSF(L)!E(il,il)!CSF(R)>
            !     ILEFT  = CSF(L) ?
            !
            coefficient_count = 0
            do packed_pair_index = 1, mcp_integral_count
                first_orbital_index = mcp_packed_orbital_pairs(packed_pair_index)/orbital_index_radix
                second_orbital_index = mod(mcp_packed_orbital_pairs(packed_pair_index),orbital_index_radix)
                if (active_orbitals_by_group(symmetry_index,subshell_index)/=first_orbital_index .or. &
                    first_orbital_index/=second_orbital_index) cycle
                if (packed_pair_index == 1) then
                    first_coefficient = 1
                else
                    first_coefficient = mcp_coefficient_end_offsets(packed_pair_index - 1) + 1
                end if
                coefficient_count = mcp_coefficient_end_offsets(packed_pair_index) - first_coefficient + 1
                exit
            end do

            do coefficient_index = 1, coefficient_count
                ! Bug 2011-08-18 Per Jonsson        IVAL = CNN(IFIRST-1+IELMNT)
                !GG         IVAL = IDNINT(CNN(IFIRST-1+IELMNT))
                occupation_number = mcp_coefficient_values(first_coefficient-1+coefficient_index)
                occupation_factor = diagonal_factor**occupation_number
                bra_csf_index = mcp_bra_csf_indices(first_coefficient - 1 + coefficient_index)
                transformed_mixing_coefficients(bra_csf_index,:level_count) = occupation_factor* &
                & input_mixing_coefficients(bra_csf_index,:level_count)
            end do
        end if
        call copy_vector(transformed_mixing_coefficients, ciouttmp, level_count*csf_count)
        !GG      call MPI_ALLREDUCE(CIOUTtmp(1,1), CIOUT(1,1), NCIV*NCSF,      &
        call mpi_allreduce(ciouttmp, transformed_mixing_coefficients, level_count*csf_count,      &
            mpi_double_precision, mpi_sum, mpi_comm_world, mpi_error_code)
        !
        !.  The previous provided us with all
        !   terms with nonvanishing occupation.
        !   For terms with vanishing occupation of il,
        !   just copy coefficients, since (x) ** 0 = 1
        !
        where (transformed_mixing_coefficients(:csf_count,:level_count) == 0.0_real64)
            transformed_mixing_coefficients(:csf_count,:level_count) = input_mixing_coefficients(:csf_count,:level_count)
        end where
        !
        if (trace_level >= 100) then
            write (output_unit, *)
            write (output_unit, *) " Input and output vectors from TIINI I,L", subshell_index, symmetry_index
            call write_matrix (input_mixing_coefficients, csf_count, level_count, csf_count, level_count)
            write (output_unit, *)
            call write_matrix (transformed_mixing_coefficients, csf_count, level_count, csf_count, level_count)
            write (output_unit, *)
        end if
        !
        if (trace_level >= 10) write (output_unit, *) " Leaving  TIINI"

        return
    end subroutine apply_occupation_factor
end module apply_occupation_factor_interface
