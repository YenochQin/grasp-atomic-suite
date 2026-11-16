!***********************************************************************
!                                                                      *
      subroutine mcp (symmetry_block_index, is_restart, mpi_rank, mpi_process_count, output_file_prefix)
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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use iso_fortran_env, only: int64
      use rangular_profile, only: stage_clock, record_stage
      use rangular_reuse, only: reuse_requested, begin_reuse_block, end_reuse_block, load_reused_pair
      use parameter_def,   only:  max_orbitals, orbital_index_radix
      use memory_man
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use outsda_i, only: outsda
      use fndbeg_i, only: fndbeg
      use alcbuf_i, only: alcbuf
      use rkco_gg_i, only: rkco_gg
      use setsda_i, only: setsda
      use sort_i, only: sort
      use sortmem_i, only: sortmem
      use alloccheck_i, only: alloccheck
!-----------------------------------------------
!   C O M M O N  B L O C K S
!-----------------------------------------------
      use buffer_c,  only: angular_coefficient_count, angular_coefficient_labels, angular_coefficients
      use debug_c,   only: ldbpa
      use default_c, only: ndef
      use iccu_c,    only: iccut
      use mcp_c,     only: max_coulomb_multipole_rank, diagonal_only, use_first_order_restriction
      use orb_c,     only: csf_count, csf_subshell_occupations
      use stat_c,    only: csf_subshell_quantum_codes, csf_coupling_quantum_codes
      use sacoef_c, only: channel_coefficients, channel_csf_columns, channel_matrix_positions, channel_packed_orbital_labels, double
      use onescalar_i, only: onescalar
      use settmpgg_i, only: settmpgg
      use cord_i, only: cord
      implicit none (type, external)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: symmetry_block_index
      integer  :: mpi_rank
      integer  :: mpi_process_count
      logical, intent(in) :: is_restart
      character(len=*), intent(in) :: output_file_prefix
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: coefficient_label_radix = orbital_index_radix
      integer, parameter :: orbital_index_radix_squared = orbital_index_radix*orbital_index_radix
      real(real64), parameter :: angular_coefficient_cutoff = 1.0e-10_real64
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! GG begin
!     2^(23)
      integer, parameter :: initial_sort_buffer_capacity = 8388608
!     2^(24)
!        PARAMETER (NGGMAX = 16777216)
!     2^(27)
!        PARAMETER (NGGMAX   = 134217728)
!        PARAMETER (NGGMAX = 536870912)
! GG end
      integer, pointer :: two_electron_counts_by_rank(:), mark(:)
      integer :: first_csf_column,first_csf_row,matrix_element_position,one_electron_coefficient_count,multipole_rank, &
         & csf_column_index,csf_row_index,orbital_a,orbital_b,lab,i,orbital_c, &
         orbital_d, orbital_swap_count, orbital_swap_temporary, lac, lbd, ntgi
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
         call fndbeg (first_csf_column, first_csf_row, matrix_element_position, one_electron_coefficient_count, &
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
      call alcbuf (1)
      if (reuse_requested) call begin_reuse_block()
      call record_stage("coefficient_buffer_allocation", symmetry_block_index, stage_start_count)
      stage_start_count = stage_clock()
!
!   Write header to  .dbg  file if appropriate
!
      if (ldbpa(2) .or. ldbpa(3)) write (99,300)
!
!   JA and JB respectively refer to the initial and final states
!   in the list of NCF configurations
!
      do 5 csf_column_index = first_csf_column, csf_count, mpi_process_count
         first_csf_row = 1
         do 4 csf_row_index = first_csf_row, csf_column_index
            if (diagonal_only .or. (use_first_order_restriction .and. (csf_row_index > iccut(symmetry_block_index))) ) then
               if (csf_row_index/=csf_column_index) cycle
            end if
!   LINCR is .TRUE. if npos is to be incremented by 1; there
!   is always a diagonal element in each column
            if (csf_row_index /= csf_column_index) then
               increment_matrix_position = .true.
            else
               matrix_element_position = matrix_element_position + 1
               increment_matrix_position = .false.
            end if
            if (reuse_requested) call load_reused_pair(csf_column_index,csf_row_index,orbital_a,orbital_b, &
               & one_electron_shell_coefficients(1))
            if (csf_row_index /= csf_column_index) then
               if (.not. reuse_requested) call onescalar(csf_column_index, &
                  & csf_row_index,orbital_a,orbital_b,one_electron_shell_coefficients)
               if (orbital_a /= 0 .and. orbital_a /= orbital_b .and.       &
                      abs (one_electron_shell_coefficients(1)) > angular_coefficient_cutoff) then
                  matrix_element_position = matrix_element_position + 1
                  increment_matrix_position = .false.
                  one_electron_coefficient_count = one_electron_coefficient_count + 1
                  lab = min (orbital_a,orbital_b) * coefficient_label_radix + max (orbital_a,orbital_b)
! GG begin
                  if(markt == 0) then
                    if(one_electron_coefficient_count > sort_capacity_per_channel) then
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
                      requested_sort_capacity = (sort_capacity_per_channel+sort_capacity_per_channel/2)*(2+ &
                         & max_coulomb_multipole_rank)
                      call alloccheck (requested_sort_capacity,irez)
                      if(irez == 1) then
                        sort_buffer_capacity = requested_sort_capacity
                        sort_capacity_per_channel = sort_buffer_capacity/(2+max_coulomb_multipole_rank)
                        call ralloc                                    &
                             (channel_coefficients,sort_capacity_per_channel,max_coulomb_multipole_rank+2,"TCOEFF","MCP")
                        call ralloc                                    &
                             (channel_csf_columns,sort_capacity_per_channel,max_coulomb_multipole_rank+2,"TICLMN","MCP")
                        call ralloc                                    &
                             (channel_matrix_positions,sort_capacity_per_channel,max_coulomb_multipole_rank+2,"TINDEX","MCP")
                        call ralloc                                    &
                             (channel_packed_orbital_labels,sort_capacity_per_channel,max_coulomb_multipole_rank+2,"TLABEL","MCP")
                        channel_csf_columns(one_electron_coefficient_count,1) = csf_column_index
                        channel_matrix_positions(one_electron_coefficient_count,1) = matrix_element_position
                        channel_packed_orbital_labels(one_electron_coefficient_count,1) = lab
                        channel_coefficients(one_electron_coefficient_count,1) = one_electron_shell_coefficients(1)
                      else
!GG                        print*,"LLISTT =",LLISTT
!GG                        print*,
!GG     :             "The program switches to the disk version of
!sorting"
                        call settmpgg (symmetry_block_index, 31, "tmp")
                        markt = 1
                        do jj = 1, one_electron_coefficient_count-1
                          write (31) channel_csf_columns(jj,1), channel_matrix_positions(jj,1),       &
                                     channel_packed_orbital_labels(jj,1), channel_coefficients(jj,1)
                        end do
                        write (31) csf_column_index, matrix_element_position, lab, one_electron_shell_coefficients(1)
                      end if
                    else
                      channel_csf_columns(one_electron_coefficient_count,1) = csf_column_index
                      channel_matrix_positions(one_electron_coefficient_count,1) = matrix_element_position
                      channel_packed_orbital_labels(one_electron_coefficient_count,1) = lab
                      channel_coefficients(one_electron_coefficient_count,1) = one_electron_shell_coefficients(1)
                    end if
                  else
                    write (31) csf_column_index, matrix_element_position, lab, one_electron_shell_coefficients(1)
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
               call rkco_gg (csf_column_index, csf_row_index, cord, 0, 1)
            end if
            do 2 i = 1, angular_coefficient_count
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
                     lac = orbital_a * coefficient_label_radix + orbital_c
                     lbd = orbital_b * coefficient_label_radix + orbital_d
                     if (lac < lbd) then
                        lab = lac * orbital_index_radix_squared + lbd
                     else
                        lab = lbd * orbital_index_radix_squared + lac
                     end if
! GG begin
                  if(mark(multipole_rank+1) == 0) then
                    if(two_electron_counts_by_rank(multipole_rank) > sort_capacity_per_channel) then
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
                      requested_sort_capacity = (sort_capacity_per_channel+sort_capacity_per_channel/2)*(2+ &
                         & max_coulomb_multipole_rank)
                      call alloccheck (requested_sort_capacity,irez)
                      if(irez == 1) then
                        sort_buffer_capacity = requested_sort_capacity
                        sort_capacity_per_channel = sort_buffer_capacity/(2+max_coulomb_multipole_rank)
                        call ralloc                                    &
                             (channel_coefficients,sort_capacity_per_channel,max_coulomb_multipole_rank+2,"TCOEFF","MCP")
                        call ralloc                                    &
                             (channel_csf_columns,sort_capacity_per_channel,max_coulomb_multipole_rank+2,"TICLMN","MCP")
                        call ralloc                                    &
                             (channel_matrix_positions,sort_capacity_per_channel,max_coulomb_multipole_rank+2,"TINDEX","MCP")
                        call ralloc                                    &
                             (channel_packed_orbital_labels,sort_capacity_per_channel,max_coulomb_multipole_rank+2,"TLABEL","MCP")
                        channel_csf_columns(two_electron_counts_by_rank(multipole_rank),multipole_rank+2) = csf_column_index
                        channel_matrix_positions(two_electron_counts_by_rank(multipole_rank),multipole_rank+2) = &
                           & matrix_element_position
                        channel_packed_orbital_labels(two_electron_counts_by_rank(multipole_rank),multipole_rank+2) = lab
                        channel_coefficients(two_electron_counts_by_rank(multipole_rank),multipole_rank+2) = &
                           & two_electron_angular_coefficient
                      else
!GG                        print*,"K=",K," LLISTV =",LLISTV(K)
!GG                        print*,
!GG     :             "The program switches to the disk version of
!sorting"
                        call settmpgg (symmetry_block_index, 32+multipole_rank, "tmp")
                        mark(multipole_rank+1) = 1
                        do jj = 1,sort_capacity_per_channel
                          write (32+multipole_rank) channel_csf_columns(jj,multipole_rank+2), &
                             & channel_matrix_positions(jj,multipole_rank+2),  &
                                       channel_packed_orbital_labels(jj,multipole_rank+2),channel_coefficients(jj,multipole_rank+2)
                        end do
                        write (32+multipole_rank) csf_column_index, matrix_element_position, lab, two_electron_angular_coefficient
                      end if
                    else
                      channel_csf_columns(two_electron_counts_by_rank(multipole_rank),multipole_rank+2) = csf_column_index
                      channel_matrix_positions(two_electron_counts_by_rank(multipole_rank),multipole_rank+2) = &
                         & matrix_element_position
                      channel_packed_orbital_labels(two_electron_counts_by_rank(multipole_rank),multipole_rank+2) = lab
                      channel_coefficients(two_electron_counts_by_rank(multipole_rank),multipole_rank+2) = &
                         & two_electron_angular_coefficient
                    end if
                  else
                    write (32+multipole_rank) csf_column_index, matrix_element_position, lab, two_electron_angular_coefficient
                  end if
! GG end


                  end if
               end if
    2       continue
!
!   All angular coefficients for this pair of CSFs have been
!   generated; update file 30
!
            if ((csf_row_index == csf_column_index) .or. (.not. increment_matrix_position))   &
               write (30) csf_column_index, csf_row_index, matrix_element_position
    4    continue
      if (mod(csf_column_index,100) == 0) then
         if (csf_column_index == 1 .or. csf_column_index == csf_count .or. mod (csf_column_index,100) == 0) then
            print *, "Row ", csf_column_index, "; nnonz = ", matrix_element_position,";  block = ", symmetry_block_index
         end if
      end if
    5 continue
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
      call alcbuf (3)
      call record_stage("csf_and_pair_buffer_release", symmetry_block_index, stage_start_count)
!
!   Write out a report for this run
!
      if (ndef /= 0 .and. mpi_rank == 0) then
         write (24,*)
         write (24,*) one_electron_coefficient_count, " T coefficients generated;"
         do multipole_rank = 0, max_coulomb_multipole_rank
            write (24,*) two_electron_counts_by_rank(multipole_rank)," V(k=",multipole_rank,") coefficients generated;"
         end do
      end if
!
!   Set up sparse structure definition arrays in file 30
!
      stage_start_count = stage_clock()
      rewind (30)
      call setsda (outsda, matrix_element_position, ldbpa(4), symmetry_block_index, mpi_rank, mpi_process_count, output_file_prefix)
      close (30, status = "DELETE")
      call record_stage("sparse_structure_and_output", symmetry_block_index, stage_start_count, &
         & int(matrix_element_position,int64), storage="disk")
!
!   Sort MCP coefficients into integral-based lists
!
! 1) the T coefficient
      stage_start_count = stage_clock()
      if(markt == 0) then
         call sortmem(31,sort_capacity_per_channel,(one_electron_coefficient_count), ntgi, ldbpa(2), &
            & symmetry_block_index, output_file_prefix)
         call record_stage("sort_and_output", symmetry_block_index, stage_start_count, &
            & int(one_electron_coefficient_count,int64), 31, "memory")
      else
         rewind (31)
         call sort (31, (one_electron_coefficient_count), ntgi, ldbpa(2), symmetry_block_index, output_file_prefix)
         close (31, status = "DELETE")
         call record_stage("sort_and_output", symmetry_block_index, stage_start_count, &
            & int(one_electron_coefficient_count,int64), 31, "disk")
      end if
      if (ndef /= 0 .and. mpi_rank == 0)                   &
                       write (24,*) ntgi, " I(ab) integrals;"
! 2) the V coefficient
      do multipole_rank = 0, max_coulomb_multipole_rank
        stage_start_count = stage_clock()
        if(mark(multipole_rank+1) == 0) then
         call sortmem(multipole_rank+32,sort_capacity_per_channel,(two_electron_counts_by_rank(multipole_rank)),ntgi, &
            & ldbpa(3),symmetry_block_index,output_file_prefix)
         call record_stage("sort_and_output", symmetry_block_index, stage_start_count, &
            & int(two_electron_counts_by_rank(multipole_rank),int64), multipole_rank+32, "memory")
        else
          rewind (multipole_rank+32)
          call sort (multipole_rank+32, (two_electron_counts_by_rank(multipole_rank)), ntgi, ldbpa(3), &
             & symmetry_block_index, output_file_prefix)
          close (multipole_rank+32, status = "DELETE")
          call record_stage("sort_and_output", symmetry_block_index, stage_start_count, &
             & int(two_electron_counts_by_rank(multipole_rank),int64), multipole_rank+32, "disk")
        end if
          if (ndef /= 0 .and. mpi_rank == 0)               &
            write (24,*) " k = ", multipole_rank, ": ", ntgi, " Slater integrals;"
      end do
      call dalloc (mark, "MARK","MCP")
      call dalloc (channel_coefficients,"TCOEFF","MCP")
      call dalloc (channel_csf_columns,"TICLMN","MCP")
      call dalloc (channel_matrix_positions,"TINDEX","MCP")
      call dalloc (channel_packed_orbital_labels,"TLABLE","MCP")
      call dalloc(two_electron_counts_by_rank, "LLISTV", "MCP")
  300 format (/"From MCP:")
      return
      end subroutine mcp
