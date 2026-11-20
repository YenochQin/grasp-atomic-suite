program descriptor_merge
   use iso_fortran_env, only: int8
   use merge_transition_csf_bases_interface, only: merge_transition_csf_bases
   use orbital_csf_basis, only: orbital_count, csf_count, csf_subshell_occupations, orbital_principal_n, orbital_kappa
   use csf_coupling_data, only: csf_subshell_quantum_codes, csf_coupling_quantum_codes
   use transition_basis_blocks, only: nblocki, nblockf, ncfi, ncff
   use grasp_parameters, only: max_orbitals
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
   integer(int8), allocatable :: expected(:,:)
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
      if (any(csf_coupling_quantum_codes(max_orbitals,:) /= [1,-3,1,-3,-3])) error stop "Merged symmetries"
      if (any(ncff(:2) /= [1,3])) error stop "Right block prefixes"
   else
      if (orbital_count /= 5 .or. csf_count /= 4) error stop "Merged dimensions"
      if (any(orbital_principal_n(:5) /= [1,2,3,4,2])) error stop "Merged principal numbers"
      if (any(orbital_kappa(:5) /= [-1,-1,-1,-1,-2])) error stop "Merged kappas"
      allocate(expected(5,4))
      expected=reshape(int([2,2,0,0,0, 2,0,1,0,1, 2,0,0,2,0, 2,1,0,0,1],int8),[5,4])
      if (any(csf_coupling_quantum_codes(max_orbitals,:) /= [1,-3,1,-3])) error stop "Merged symmetries"
      if (any(ncff(:2) /= [1,2])) error stop "Right block prefixes"
   end if
   if (any(csf_subshell_occupations(:orbital_count,:) /= expected)) error stop "Merged occupation mapping"
   if (any(csf_subshell_quantum_codes(:orbital_count,1,:) /= merge(1_int8,0_int8,expected == 1))) &
      error stop "Merged seniority mapping"
   if (nblocki /= 2 .or. nblockf /= 2) error stop "Lost input blocks"
   if (any(ncfi(:2) /= [1,2])) error stop "Lost input block prefixes"
#ifdef DESCRIPTOR_MPI
   call mpi_finalize(status)
#endif
   print *, "PASS unequal orbital tables merged without text"
end program descriptor_merge
