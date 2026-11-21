program descriptor_lifecycle
   use iso_fortran_env, only: int8
   use iso_c_binding, only: c_associated
   use csf_coupling_data, only: csf_angular_code_count, csf_total_j_slot, csf_parity_slot
   use grasp_parameters, only: max_orbitals
   use orbital_csf_basis, only: csf_subshell_occupations, csf_count, orbital_count, orbital_two_j
   use csf_coupling_data, only: csf_subshell_quantum_codes, csf_angular_codes
   use csf_pair_quantum_numbers, only: prepare_csf_pair_quantum_numbers
   use csf_pair_angular_state, only: bra_intermediate_couplings, ket_intermediate_couplings, peel_orbital_count
   use csf_quantum_queries, only: csf_intermediate_coupling, csf_total_j_code, csf_parity_code
   use csf_descriptor_reader, only: descriptor_source, open_descriptor_source, close_descriptor_source, read_descriptor_block
   implicit none (type, external)
   type(descriptor_source) :: first, second
   character(len=4096) :: filename
   character(len=16) :: scenario
   integer :: status, b, n
   integer(int8), allocatable :: occupations(:,:), other(:,:), quantum(:,:,:), couplings(:,:)
   call get_command_argument(1,filename)
   call get_command_argument(2,scenario)
   call open_descriptor_source(first,trim(filename)//".missing",50,status=status)
   if (status == 0 .or. first%is_open .or. c_associated(first%handle)) error stop "failed open retained state"
   call open_descriptor_source(first,trim(filename),50,status=status)
   if (status /= 0 .or. first%format_version /= 3) error stop "valid V3 open failed"
   call open_descriptor_source(second,trim(filename),50,status=status)
   if (status /= 0) error stop "second handle open failed"
   do b=1,first%block_count
      n=first%block_lengths(b)
      allocate(occupations(max_orbitals,n),other(max_orbitals,n),quantum(max_orbitals,3,n),couplings(csf_angular_code_count,n))
      call read_descriptor_block(first,b,occupations,quantum,couplings,status)
      if (status /= 0) error stop "full block failed"
      if (scenario == "pairs") call check_pairs()
      if (first%peel_count == 0) then
         if (sum(int(occupations(:,1))) /= first%electron_count) error stop "Closed-core electron count"
         if (any(quantum(:first%core_count,1:2,:) /= 0) .or. &
             any(quantum(:first%core_count,3,:) /= 1) .or. any(couplings(:first%core_count,:) /= 1) .or. &
             any(couplings(csf_total_j_slot,:) /= 1) .or. any(couplings(csf_parity_slot,:) /= 1)) &
             error stop "Closed-core quantum numbers"
      end if
      call read_descriptor_block(second,b,other,status=status)
      if (status /= 0 .or. any(other /= occupations)) error stop "occupation-only block differs"
      deallocate(occupations,other,quantum,couplings)
   end do
   call close_descriptor_source(first)
   call close_descriptor_source(first)
   call close_descriptor_source(second)
   call open_descriptor_source(first,trim(filename),50,status=status)
   if (status /= 0) error stop "reopen failed"
   allocate(occupations(max_orbitals,1))
   call read_descriptor_block(first,0,occupations,status=status)
   if (status == 0 .or. first%is_open .or. c_associated(first%handle)) error stop "bad block did not close source"
   if (any(occupations /= 0)) error stop "bad block left occupation data"
   print *, "PASS descriptor lifecycle and occupation-only equivalence"
contains
   subroutine check_pairs()
      orbital_count = first%orbital_count
      orbital_two_j(:orbital_count) = first%single_electron_two_j
      csf_count = 3
      if (orbital_count /= 4 .or. n /= 3) error stop "Unexpected analytic fixture dimensions"
      allocate(csf_subshell_occupations, source=occupations)
      allocate(csf_subshell_quantum_codes, source=quantum)
      allocate(csf_angular_codes, source=couplings)
      if (any(couplings(:4,:) /= reshape([1,2,1,2, 1,2,3,2, 1,1,2,2],[4,3]))) &
         error stop "Unexpected cumulative codes"
      if (csf_total_j_code(1) /= 2 .or. csf_parity_code(1) /= 1) error stop "Unexpected symmetry"
      if (csf_intermediate_coupling(1,2) /= 3 .or. csf_intermediate_coupling(2,2) /= 2) &
         error stop "Unexpected open-shell step"
      call prepare_csf_pair_quantum_numbers(1,3)
      if (peel_orbital_count /= 3) error stop "Unexpected pair peel count"
      if (any(bra_intermediate_couplings(:2) /= [1,2]) .or. any(ket_intermediate_couplings(:2) /= [2,2])) &
         error stop "Unexpected A0/B pair"
      call prepare_csf_pair_quantum_numbers(3,1)
      if (any(bra_intermediate_couplings(:2) /= [2,2]) .or. any(ket_intermediate_couplings(:2) /= [1,2])) &
         error stop "Unexpected B/A0 pair"
      call prepare_csf_pair_quantum_numbers(2,3)
      if (any(bra_intermediate_couplings(:2) /= [3,2]) .or. any(ket_intermediate_couplings(:2) /= [2,2])) &
         error stop "Unexpected A1/B pair"
      call prepare_csf_pair_quantum_numbers(3,3)
      if (peel_orbital_count /= 1) error stop "Unexpected B/B peel count"
      deallocate(csf_subshell_occupations,csf_subshell_quantum_codes,csf_angular_codes)
   end subroutine check_pairs
end program descriptor_lifecycle
