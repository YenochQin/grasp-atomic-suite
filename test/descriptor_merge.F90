program descriptor_merge
   use iso_fortran_env, only: int8
   use merge_transition_csf_bases_interface, only: merge_transition_csf_bases
   use orbital_csf_basis, only: orbital_count, csf_count, csf_subshell_occupations, orbital_principal_n, orbital_kappa
   use csf_coupling_data, only: csf_subshell_quantum_codes, csf_angular_codes
   use transition_basis_blocks, only: initial_transition_block_count, final_transition_block_count, &
       initial_cumulative_csf_counts, final_cumulative_csf_counts
   use csf_quantum_queries, only: csf_total_j_code, csf_parity_code
#ifdef DESCRIPTOR_MPI
   use mpi_state, only: mpi_init, mpi_comm_rank, mpi_comm_size, mpi_comm_world, mpi_rank, mpi_process_count, mpi_finalize
#endif
   implicit none (type, external)
#ifdef DESCRIPTOR_MPI
   character(len=128) :: names(2)
   integer :: status
#else
   character(len=24) :: names(2)
#endif
   integer(int8), allocatable :: expected(:,:), expected_angular(:,:)
   integer :: csf_index
   character(len=16) :: scenario
#ifdef DESCRIPTOR_MPI
   call mpi_init(status)
   call mpi_comm_rank(mpi_comm_world,mpi_rank,status)
   call mpi_comm_size(mpi_comm_world,mpi_process_count,status)
#endif
   call get_command_argument(1,scenario)
   names(1)="left"
   names(2)="right"
   call merge_transition_csf_bases(names)
   if (trim(scenario) == "multi") then
      if (orbital_count /= 6 .or. csf_count /= 5) error stop "Merged dimensions"
      if (any(orbital_principal_n(:6) /= [1,2,3,4,5,2])) error stop "Merged principal numbers"
      if (any(orbital_kappa(:6) /= [-1,-1,-1,-1,-1,-2])) error stop "Merged kappas"
      allocate(expected(6,5))
      expected=reshape(int([2,2,0,0,0,0, 2,0,1,0,0,1, 2,0,0,2,0,0, &
                            2,0,0,0,1,1, 2,1,0,0,0,1],int8),[6,5])
      allocate(expected_angular(6,5))
      expected_angular=reshape(int([1,1,1,1,1,1, 1,1,2,2,2,3, 1,1,1,1,1,1, &
                                    1,1,1,1,2,3, 1,2,2,2,2,3],int8),[6,5])
      if (any(final_cumulative_csf_counts(:2) /= [1,3])) error stop "Right block prefixes"
   else
      if (orbital_count /= 5 .or. csf_count /= 4) error stop "Merged dimensions"
      if (any(orbital_principal_n(:5) /= [1,2,3,4,2])) error stop "Merged principal numbers"
      if (any(orbital_kappa(:5) /= [-1,-1,-1,-1,-2])) error stop "Merged kappas"
      allocate(expected(5,4))
      expected=reshape(int([2,2,0,0,0, 2,0,1,0,1, 2,0,0,2,0, 2,1,0,0,1],int8),[5,4])
      allocate(expected_angular(5,4))
      expected_angular=reshape(int([1,1,1,1,1, 1,1,2,2,3, 1,1,1,1,1, 1,2,2,2,3],int8),[5,4])
      if (any(final_cumulative_csf_counts(:2) /= [1,2])) error stop "Right block prefixes"
   end if
   if (any(csf_angular_codes(:orbital_count,:) /= expected_angular)) error stop "Merged cumulative coupling mapping"
   do csf_index=1,csf_count
      if (csf_total_j_code(csf_index) /= expected_angular(orbital_count,csf_index)) error stop "Merged total J"
      if (csf_parity_code(csf_index) /= (-1)**int(expected(orbital_count,csf_index))) error stop "Merged parity"
   end do
   if (any(csf_subshell_occupations(:orbital_count,:) /= expected)) error stop "Merged occupation mapping"
   if (any(csf_subshell_quantum_codes(:orbital_count,1,:) /= merge(1_int8,0_int8,expected == 1))) &
      error stop "Merged seniority mapping"
   if (initial_transition_block_count /= 2 .or. final_transition_block_count /= 2) error stop "Lost input blocks"
   if (any(initial_cumulative_csf_counts(:2) /= [1,2])) error stop "Lost input block prefixes"
#ifdef DESCRIPTOR_MPI
   call mpi_finalize(status)
#endif
   print *, "PASS unequal orbital tables merged without text"
end program descriptor_merge
