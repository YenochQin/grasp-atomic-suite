module compute_mcp_angular_coefficients_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: compute_mcp_angular_coefficients
contains

    !***********************************************************************
    !                                                                      *
    subroutine compute_mcp_angular_coefficients (symmetry_block_index, is_restart, mpi_rank, mpi_process_count, &
        & output_file_prefix, temporary_file_prefix)
        ! Module dependencies
        use mcp_generation_options, only: mcp_temporary_units
        use calculation_file_units, only: summary_output_unit
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        !                                                                      *
        !   This routine controls the computation  and storage of the values   *
        !   and all indices of the angular coefficients                        *
        !                                                                      *
        !                                       k                              *
        !                   T  (ab)            V  (abcd)                       *
        !                    rs                 rs                             *
        !                                                                      *
        !   k is the multipolarity of a two-particle Coulomb integral. a, b,   *
        !   c and d are orbital sequence numbers.  r and s are configuration   *
        !   state function indices.                                            *
        !                                                                      *
        !   Call(s) to: [LIB92]: ALCBUF, ALLOC, CONVRT, DALLOC, RKCO_GG,       *
        !                        TNSRJJ.                                       *
        !               [GENMCP]: FNDBEG, SETSDA, SORT.                        *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 28 Sep 1993   *
        !   Modified by C. Froese Fischer for block computation.               *
        !   Modified by Gediminas Gaigalas:                         Feb 2017   *
        !                               1) for new spin-angular integration,   *
        !                               2) for sorting in the memory.          *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use iso_fortran_env, only: int64
        use rangular_profile, only: stage_clock, record_stage
        use rangular_reuse, only: reuse_requested, begin_reuse_block, end_reuse_block, load_reused_pair
        use grasp_parameters,   only:  max_orbitals, orbital_index_radix
        use array_memory_management, only: alloc, dalloc, ralloc

        ! Procedure interfaces
        use write_mcp_sparsity_summary_interface, only: write_mcp_sparsity_summary
        use locate_mcp_restart_position_interface, only: locate_mcp_restart_position

        ! Module dependencies
        use angular_coefficient_buffer, only: manage_angular_coefficient_buffer
        use two_electron_angular_coefficients, only: calculate_two_electron_angular_coefficients

        ! Procedure interfaces
        use write_mcp_sparsity_interface, only: write_mcp_sparsity
        use sort_mcp_coefficients_interface, only: sort_mcp_coefficients
        use sort_mcp_coefficients_in_memory_interface, only: sort_mcp_coefficients_in_memory
        use try_allocate_sort_buffers_interface, only: try_allocate_sort_buffers

        ! Module dependencies
        use angular_coefficient_buffer_state,  only: angular_coefficient_count, angular_coefficient_labels, angular_coefficients
        use debug_output_settings,   only: ldbpa
        use calculation_input_defaults, only: nondefault_options_flag
        use csf_block_cutoffs,    only: active_csf_cutoffs
        use mcp_generation_options,     only: max_coulomb_multipole_rank, diagonal_only, use_first_order_restriction
        use orbital_csf_basis,     only: csf_count, csf_subshell_occupations
        use csf_coupling_data,    only: csf_subshell_quantum_codes, csf_coupling_quantum_codes
        use angular_coefficient_channels, only: &
            channel_coefficients, channel_csf_columns, channel_matrix_positions, channel_packed_orbital_labels
        use one_electron_scalar_coefficients, only: calculate_one_electron_scalar_coefficients

        ! Procedure interfaces
        use open_temporary_mcp_channels_interface, only: open_temporary_mcp_channels

        ! Module dependencies
        use closed_shell_mcp, only: calculate_closed_shell_mcp
        implicit none (type, external)

        ! Arguments
        integer  :: symmetry_block_index
        integer  :: mpi_rank
        integer  :: mpi_process_count
        logical, intent(in) :: is_restart
        character(len=*), intent(in) :: output_file_prefix, temporary_file_prefix

        ! Local constants
        integer, parameter :: orbital_index_radix_squared = orbital_index_radix*orbital_index_radix
        real(real64), parameter :: angular_coefficient_cutoff = 1.0e-10_real64
        ! GG begin
        !     2^(23)
        integer, parameter :: initial_sort_buffer_capacity = 8388608

        ! Local declarations
        !     2^(24)
        !        PARAMETER (NGGMAX = 16777216)
        !     2^(27)
        !        PARAMETER (NGGMAX   = 134217728)
        !        PARAMETER (NGGMAX = 536870912)
        ! GG end
        integer, pointer :: two_electron_counts_by_rank(:), mark(:)
        integer :: first_csf_column,first_csf_row,matrix_element_position,one_electron_coefficient_count,multipole_rank, &
            & csf_index,csf_row_index,orbital_a,orbital_b,packed_orbital_label,i,orbital_c, &
            orbital_d, orbital_swap_count, orbital_swap_temporary, first_orbital_pair_label, second_orbital_pair_label, &
                & distinct_radial_integral_count
        real(real64) :: one_electron_shell_coefficients(max_orbitals)
        real(real64) :: two_electron_angular_coefficient
        logical :: is_direct_monopole_integral, increment_matrix_position
        integer :: irez, sort_buffer_capacity, sort_capacity_per_channel, markt, requested_sort_capacity, jj
        integer(int64) :: stage_start_count
        !-----------------------------------------------
        !
        !   Allocate storage to the array that stores the lengths
        !   of the lists of V coefficients of each multipolarity
        !
        ! GG begin
        stage_start_count = stage_clock()
        irez = 1
        sort_buffer_capacity = initial_sort_buffer_capacity
        sort_capacity_per_channel = sort_buffer_capacity/(2+max_coulomb_multipole_rank)
        sort_buffer_capacity = sort_capacity_per_channel*(2+max_coulomb_multipole_rank)
        !GG      PRINT*, "The program alloc Memory for sorting:",IGGMAX
        call alloc(two_electron_counts_by_rank, 0, max_coulomb_multipole_rank, "LLISTV", "MCP")
        call alloc (mark, max_coulomb_multipole_rank + 1,"MARK", "MCP")
        call alloc (channel_coefficients,sort_capacity_per_channel,max_coulomb_multipole_rank+2,"TCOEFF","MCP")
        call alloc (channel_csf_columns,sort_capacity_per_channel,max_coulomb_multipole_rank+2,"TICLMN","MCP")
        call alloc (channel_matrix_positions,sort_capacity_per_channel,max_coulomb_multipole_rank+2,"TINDEX","MCP")
        call alloc (channel_packed_orbital_labels,sort_capacity_per_channel,max_coulomb_multipole_rank+2,"TLABEL","MCP")
        markt = 0
        mark = 0
        ! GG end
        !
        !   If this is a restart, determine the pair of CSFs with which
        !   the computation should begin.
        !   If not restart, initialize the counters for the total number of
        !   T coefficients and V coefficients of each multipolarity.
        !
        if (is_restart) then
            print *, "Not ready for RESTART"
            stop
            call locate_mcp_restart_position (first_csf_column, first_csf_row, matrix_element_position, &
                & one_electron_coefficient_count, &
                & two_electron_counts_by_rank)
        else
            first_csf_column = 1 + mpi_rank
            matrix_element_position = 0
            one_electron_coefficient_count = 0
            do multipole_rank = 0, max_coulomb_multipole_rank
                two_electron_counts_by_rank(multipole_rank) = 0
            end do
        end if
        !
        !   Allocate storage for the arrays in BUFFER
        !
        call manage_angular_coefficient_buffer (1)
        if (reuse_requested) call begin_reuse_block()
        call record_stage("coefficient_buffer_allocation", symmetry_block_index, stage_start_count)
        stage_start_count = stage_clock()
        !
        !   Write header to  .dbg  file if appropriate
        !
        if (ldbpa(2) .or. ldbpa(3)) then
            call ensure_debug_output()
            write (angular_debug_unit,'(/"From MCP:")')
        end if
        !
        !   JA and JB respectively refer to the initial and final states
        !   in the list of NCF configurations
        !
        do csf_index = first_csf_column, csf_count, mpi_process_count
            first_csf_row = 1
            do csf_row_index = first_csf_row, csf_index
                if (diagonal_only .or. (use_first_order_restriction .and. (csf_row_index > &
                    active_csf_cutoffs(symmetry_block_index))) ) then
                    if (csf_row_index/=csf_index) cycle
                end if
                !   LINCR is .TRUE. if npos is to be incremented by 1; there
                !   is always a diagonal element in each column
                if (csf_row_index /= csf_index) then
                    increment_matrix_position = .true.
                else
                    matrix_element_position = matrix_element_position + 1
                    increment_matrix_position = .false.
                end if
                if (reuse_requested) call load_reused_pair(csf_index,csf_row_index,orbital_a,orbital_b, &
                    & one_electron_shell_coefficients(1))
                if (csf_row_index /= csf_index) then
                    if (.not. reuse_requested) call calculate_one_electron_scalar_coefficients(csf_index, &
                        & csf_row_index,orbital_a,orbital_b,one_electron_shell_coefficients)
                    if (orbital_a /= 0 .and. orbital_a /= orbital_b .and.       &
                        abs (one_electron_shell_coefficients(1)) > angular_coefficient_cutoff) then
                        matrix_element_position = matrix_element_position + 1
                        increment_matrix_position = .false.
                        one_electron_coefficient_count = one_electron_coefficient_count + 1
                        packed_orbital_label = min (orbital_a,orbital_b) * orbital_index_radix + max (orbital_a,orbital_b)
                        ! GG begin
                        if(markt == 0) then
                            if(one_electron_coefficient_count > sort_capacity_per_channel) then
                                !   Increase array length by half the present length if the latter
                                !   is inadequate to store another pair
                                !
                                requested_sort_capacity = (sort_capacity_per_channel+sort_capacity_per_channel/2)*(2+ &
                                    & max_coulomb_multipole_rank)
                                call try_allocate_sort_buffers (requested_sort_capacity,irez)
                                if(irez == 1) then
                                    sort_buffer_capacity = requested_sort_capacity
                                    sort_capacity_per_channel = sort_buffer_capacity/(2+max_coulomb_multipole_rank)
                                    call ralloc                                    &
                                        (channel_coefficients,sort_capacity_per_channel,max_coulomb_multipole_rank+2,"TCOEFF","MCP")
                                    call ralloc                                    &
                                        (channel_csf_columns,sort_capacity_per_channel,max_coulomb_multipole_rank+2,"TICLMN","MCP")
                                    call ralloc                                    &
                                        (channel_matrix_positions,sort_capacity_per_channel,max_coulomb_multipole_rank+2, &
                                        & "TINDEX","MCP")
                                    call ralloc                                    &
                                        (channel_packed_orbital_labels,sort_capacity_per_channel,max_coulomb_multipole_rank+2, &
                                        & "TLABEL","MCP")
                                    channel_csf_columns(one_electron_coefficient_count,1) = csf_index
                                    channel_matrix_positions(one_electron_coefficient_count,1) = matrix_element_position
                                    channel_packed_orbital_labels(one_electron_coefficient_count,1) = packed_orbital_label
                                    channel_coefficients(one_electron_coefficient_count,1) = one_electron_shell_coefficients(1)
                                else
                                    !GG                        print*,"LLISTT =",LLISTT
                                    !GG                        print*,
                                    !GG     :             "The program switches to the disk version of
                                    !sorting"
                                    call open_temporary_mcp_channels (symmetry_block_index, 31, temporary_file_prefix)
                                    markt = 1
                                    do jj = 1, one_electron_coefficient_count-1
                                        write (mcp_temporary_units(31)) channel_csf_columns(jj,1), channel_matrix_positions(jj,1), &
                                            channel_packed_orbital_labels(jj,1), channel_coefficients(jj,1)
                                    end do
                                    write (mcp_temporary_units(31)) csf_index, matrix_element_position, packed_orbital_label, &
                                        one_electron_shell_coefficients(1)
                                end if
                            else
                                channel_csf_columns(one_electron_coefficient_count,1) = csf_index
                                channel_matrix_positions(one_electron_coefficient_count,1) = matrix_element_position
                                channel_packed_orbital_labels(one_electron_coefficient_count,1) = packed_orbital_label
                                channel_coefficients(one_electron_coefficient_count,1) = one_electron_shell_coefficients(1)
                            end if
                        else
                            write (mcp_temporary_units(31)) csf_index, matrix_element_position, packed_orbital_label, &
                                one_electron_shell_coefficients(1)
                        end if
                        ! GG end
                    end if
                end if
                !
                !   Call the MCP package to generate V coefficients; ac and bd
                !   are the density pairs. NVCOEF is initialized here but changed
                !   in /rkco/cor[d] via COMMON.
                !
                if (.not. reuse_requested) then
                    angular_coefficient_count = 0
                    call calculate_two_electron_angular_coefficients (csf_index, csf_row_index, calculate_closed_shell_mcp, 0, 1)
                end if
                do i = 1, angular_coefficient_count
                    two_electron_angular_coefficient = angular_coefficients(i)
                    if (abs (two_electron_angular_coefficient) > angular_coefficient_cutoff) then
                        orbital_a = angular_coefficient_labels(1,i)
                        orbital_b = angular_coefficient_labels(2,i)
                        orbital_c = angular_coefficient_labels(3,i)
                        orbital_d = angular_coefficient_labels(4,i)
                        multipole_rank  = angular_coefficient_labels(5,i)
                        is_direct_monopole_integral=(multipole_rank==0).and.(orbital_a==orbital_c).and.(orbital_b==orbital_d)
                        if (.not. is_direct_monopole_integral) then
                            ! Swap index to make sure IA <= IC, IB <= ID and record the number
                            ! of swaps
                            orbital_swap_count = 0
                            if (orbital_a > orbital_c) then
                                orbital_swap_temporary = orbital_c
                                orbital_c = orbital_a
                                orbital_a = orbital_swap_temporary
                                orbital_swap_count = orbital_swap_count + 1
                            end if
                            if (orbital_b > orbital_d) then
                                orbital_swap_temporary = orbital_d
                                orbital_d = orbital_b
                                orbital_b = orbital_swap_temporary
                                orbital_swap_count = orbital_swap_count + 1
                            end if
                            if (increment_matrix_position) then
                                matrix_element_position = matrix_element_position + 1
                                increment_matrix_position = .false.
                            end if
                            two_electron_counts_by_rank(multipole_rank) = two_electron_counts_by_rank(multipole_rank) + 1
                            first_orbital_pair_label = orbital_a * orbital_index_radix + orbital_c
                            second_orbital_pair_label = orbital_b * orbital_index_radix + orbital_d
                            if (first_orbital_pair_label < second_orbital_pair_label) then
                                packed_orbital_label = first_orbital_pair_label * orbital_index_radix_squared + &
                                & second_orbital_pair_label
                            else
                                packed_orbital_label = second_orbital_pair_label * orbital_index_radix_squared + &
                                & first_orbital_pair_label
                            end if
                            ! GG begin
                            if(mark(multipole_rank+1) == 0) then
                                if(two_electron_counts_by_rank(multipole_rank) > sort_capacity_per_channel) then
                                    !   Increase array length by half the present length if the latter
                                    !   is inadequate to store another pair
                                    !
                                    requested_sort_capacity = (sort_capacity_per_channel+sort_capacity_per_channel/2)*(2+ &
                                        & max_coulomb_multipole_rank)
                                    call try_allocate_sort_buffers (requested_sort_capacity,irez)
                                    if(irez == 1) then
                                        sort_buffer_capacity = requested_sort_capacity
                                        sort_capacity_per_channel = sort_buffer_capacity/(2+max_coulomb_multipole_rank)
                                        call ralloc                                    &
                                            (channel_coefficients,sort_capacity_per_channel,max_coulomb_multipole_rank+2, &
                                            & "TCOEFF","MCP")
                                        call ralloc                                    &
                                            (channel_csf_columns,sort_capacity_per_channel,max_coulomb_multipole_rank+2, &
                                            & "TICLMN","MCP")
                                        call ralloc                                    &
                                            (channel_matrix_positions,sort_capacity_per_channel,max_coulomb_multipole_rank+2, &
                                            & "TINDEX","MCP")
                                        call ralloc                                    &
                                            (channel_packed_orbital_labels,sort_capacity_per_channel,max_coulomb_multipole_rank+2, &
                                            & "TLABEL","MCP")
                                        channel_csf_columns(two_electron_counts_by_rank(multipole_rank),multipole_rank+2) = &
                                            & csf_index
                                        channel_matrix_positions(two_electron_counts_by_rank(multipole_rank),multipole_rank+2) = &
                                            & matrix_element_position
                                        channel_packed_orbital_labels(two_electron_counts_by_rank(multipole_rank), &
                                            & multipole_rank+2) = packed_orbital_label
                                        channel_coefficients(two_electron_counts_by_rank(multipole_rank),multipole_rank+2) = &
                                            & two_electron_angular_coefficient
                                    else
                                        !GG                        print*,"K=",K," LLISTV =",LLISTV(K)
                                        !GG                        print*,
                                        !GG     :             "The program switches to the disk version of
                                        !sorting"
                                        call open_temporary_mcp_channels (symmetry_block_index, 32+multipole_rank, &
                                            & temporary_file_prefix)
                                        mark(multipole_rank+1) = 1
                                        do jj = 1,sort_capacity_per_channel
                                            write (mcp_temporary_units(32+multipole_rank)) channel_csf_columns(jj, &
                                                & multipole_rank+2), &
                                                & channel_matrix_positions(jj,multipole_rank+2),  &
                                                channel_packed_orbital_labels(jj,multipole_rank+2),channel_coefficients(jj, &
                                                & multipole_rank+2)
                                        end do
                                        write (mcp_temporary_units(32+multipole_rank)) csf_index, matrix_element_position, &
                                & packed_orbital_label, &
                                            two_electron_angular_coefficient
                                    end if
                                else
                                    channel_csf_columns(two_electron_counts_by_rank(multipole_rank),multipole_rank+2) = csf_index
                                    channel_matrix_positions(two_electron_counts_by_rank(multipole_rank),multipole_rank+2) = &
                                        & matrix_element_position
                                    channel_packed_orbital_labels(two_electron_counts_by_rank(multipole_rank), &
                                        & multipole_rank+2) = packed_orbital_label
                                    channel_coefficients(two_electron_counts_by_rank(multipole_rank),multipole_rank+2) = &
                                        & two_electron_angular_coefficient
                                end if
                            else
                                write (mcp_temporary_units(32+multipole_rank)) csf_index, matrix_element_position, &
                                & packed_orbital_label, &
                                    two_electron_angular_coefficient
                            end if
                            ! GG end


                        end if
                    end if
                end do
                !
                !   All angular coefficients for this pair of CSFs have been
                !   generated; update file 30
                !
                if ((csf_row_index == csf_index) .or. (.not. increment_matrix_position))   &
                    write (mcp_temporary_units(30)) csf_index, csf_row_index, matrix_element_position
            end do
            if (mod(csf_index,100) == 0) then
                if (csf_index == 1 .or. csf_index == csf_count .or. mod (csf_index,100) == 0) then
                    print *, "Row ", csf_index, "; nnonz = ", matrix_element_position,";  block = ", symmetry_block_index
                end if
            end if
        end do
        call record_stage("pair_generation_and_buffering", symmetry_block_index, stage_start_count, &
            int(one_electron_coefficient_count,int64) + sum(int(two_electron_counts_by_rank,int64)), &
            & storage="memory_or_spill")
        if (reuse_requested) call end_reuse_block(symmetry_block_index)
        stage_start_count = stage_clock()
        !
        !   Deallocate storage that is no longer required
        !
        call dalloc (csf_subshell_occupations, "IQA", "MCP")
        call dalloc (csf_subshell_quantum_codes, "JQSA", "MCP")
        call dalloc (csf_coupling_quantum_codes, "JCUPA", "MCP")
        call manage_angular_coefficient_buffer (3)
        call record_stage("csf_and_pair_buffer_release", symmetry_block_index, stage_start_count)
        !
        !   Write out a report for this run
        !
        if (nondefault_options_flag /= 0 .and. mpi_rank == 0) then
            write (summary_output_unit,*)
            write (summary_output_unit,*) one_electron_coefficient_count, " T coefficients generated;"
            do multipole_rank = 0, max_coulomb_multipole_rank
                write (summary_output_unit,*) two_electron_counts_by_rank(multipole_rank)," V(k=",multipole_rank, &
                    ") coefficients generated;"
            end do
        end if
        !
        !   Set up sparse structure definition arrays in file 30
        !
        stage_start_count = stage_clock()
        rewind (mcp_temporary_units(30))
        call write_mcp_sparsity (write_mcp_sparsity_summary, matrix_element_position, ldbpa(4), symmetry_block_index, mpi_rank, &
            & mpi_process_count, output_file_prefix)
        close (mcp_temporary_units(30), status = "DELETE")
        call record_stage("sparse_structure_and_output", symmetry_block_index, stage_start_count, &
            & int(matrix_element_position,int64), storage="disk")
        !
        !   Sort MCP coefficients into integral-based lists
        !
        ! 1) the T coefficient
        stage_start_count = stage_clock()
        if(markt == 0) then
            call sort_mcp_coefficients_in_memory(31,sort_capacity_per_channel,(one_electron_coefficient_count), &
                & distinct_radial_integral_count, ldbpa(2), &
                & symmetry_block_index, output_file_prefix)
            call record_stage("sort_and_output", symmetry_block_index, stage_start_count, &
                & int(one_electron_coefficient_count,int64), 31, "memory")
        else
            rewind (mcp_temporary_units(31))
            call sort_mcp_coefficients (31, (one_electron_coefficient_count), distinct_radial_integral_count, ldbpa(2), &
                & symmetry_block_index, &
                & output_file_prefix)
            close (mcp_temporary_units(31), status = "DELETE")
            call record_stage("sort_and_output", symmetry_block_index, stage_start_count, &
                & int(one_electron_coefficient_count,int64), 31, "disk")
        end if
        if (nondefault_options_flag /= 0 .and. mpi_rank == 0)                   &
            write (summary_output_unit,*) distinct_radial_integral_count, " I(ab) integrals;"
        ! 2) the V coefficient
        do multipole_rank = 0, max_coulomb_multipole_rank
            stage_start_count = stage_clock()
            if(mark(multipole_rank+1) == 0) then
                call sort_mcp_coefficients_in_memory(multipole_rank+32,sort_capacity_per_channel, &
                    & (two_electron_counts_by_rank(multipole_rank)),distinct_radial_integral_count, &
                    & ldbpa(3),symmetry_block_index,output_file_prefix)
                call record_stage("sort_and_output", symmetry_block_index, stage_start_count, &
                    & int(two_electron_counts_by_rank(multipole_rank),int64), multipole_rank+32, "memory")
            else
                rewind (mcp_temporary_units(multipole_rank+32))
                call sort_mcp_coefficients (multipole_rank+32, (two_electron_counts_by_rank(multipole_rank)), &
                    & distinct_radial_integral_count, ldbpa(3), &
                    & symmetry_block_index, output_file_prefix)
                close (mcp_temporary_units(multipole_rank+32), status = "DELETE")
                call record_stage("sort_and_output", symmetry_block_index, stage_start_count, &
                    & int(two_electron_counts_by_rank(multipole_rank),int64), multipole_rank+32, "disk")
            end if
            if (nondefault_options_flag /= 0 .and. mpi_rank == 0)               &
                write (summary_output_unit,*) " k = ", multipole_rank, ": ", distinct_radial_integral_count, " Slater integrals;"
        end do
        call dalloc (mark, "MARK","MCP")
        call dalloc (channel_coefficients,"TCOEFF","MCP")
        call dalloc (channel_csf_columns,"TICLMN","MCP")
        call dalloc (channel_matrix_positions,"TINDEX","MCP")
        call dalloc (channel_packed_orbital_labels,"TLABLE","MCP")
        call dalloc(two_electron_counts_by_rank, "LLISTV", "MCP")
        return
    end subroutine compute_mcp_angular_coefficients
end module compute_mcp_angular_coefficients_interface
