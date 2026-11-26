program descriptor_lifecycle
   use iso_fortran_env, only: int8
   use csf_coupling_data, only: csf_angular_code_count, csf_total_j_slot, csf_parity_slot
   use grasp_parameters, only: max_orbitals
   use orbital_csf_basis, only: csf_subshell_occupations, csf_count, orbital_count, orbital_two_j
   use csf_coupling_data, only: csf_subshell_quantum_codes, csf_angular_codes
   use csf_pair_quantum_numbers, only: prepare_csf_pair_quantum_numbers
   use csf_pair_angular_state, only: bra_intermediate_couplings, ket_intermediate_couplings, peel_orbital_count
   use csf_quantum_queries, only: csf_intermediate_coupling, csf_total_j_code, csf_parity_code
   use csf_descriptor_reader, only: descriptor_source, open_descriptor_source, close_descriptor_source, read_descriptor_block, &
      loaded_csf_basis_id
   use csf_full_basis_reader, only: load_csf_basis
   implicit none (type, external)
   type(descriptor_source) :: first, second
   character(len=4096) :: filename
   character(len=16) :: scenario
   integer :: status, b, n
   integer(int8), allocatable :: occupations(:,:), other(:,:), quantum(:,:,:), couplings(:,:)
   call get_command_argument(1,filename)
   call get_command_argument(2,scenario)
   if (scenario == "batches") then
      call check_batch_boundary()
      print *, "PASS descriptor batch boundary"
      stop
   end if
   call open_descriptor_source(first,trim(filename)//".missing",50,status=status)
   if (status == 0 .or. first%is_open()) error stop "failed open retained state"
   call open_descriptor_source(first,trim(filename),50,status=status)
   if (status /= 0 .or. first%format_version /= 3) error stop "valid V3 open failed"
   call open_descriptor_source(second,trim(filename),50,status=status)
   if (status /= 0) error stop "second handle open failed"
   do b=first%block_count,1,-1
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
      call read_descriptor_block(first,b,other,status=status)
      if (status /= 0 .or. any(other /= occupations)) error stop "rereading block differs"
      deallocate(occupations,other,quantum,couplings)
   end do
   call close_descriptor_source(first)
   call close_descriptor_source(first)
   call close_descriptor_source(second)
   call open_descriptor_source(first,trim(filename),50,status=status)
   if (status /= 0) error stop "reopen failed"
   allocate(occupations(max_orbitals,1))
   call read_descriptor_block(first,0,occupations,status=status)
   if (status == 0 .or. first%is_open()) error stop "bad block did not close source"
   if (first%basis_id /= "" .or. allocated(first%block_lengths)) error stop "bad block retained basis metadata"
   if (any(occupations /= 0)) error stop "bad block left occupation data"
   if (scenario == "dual") call check_distinct_bases()
   print *, "PASS descriptor lifecycle and occupation-only equivalence"
contains
   subroutine check_batch_boundary()
      integer(int8) :: expected(max_orbitals,3)

      call open_descriptor_source(first,trim(filename),50,status=status)
      if (status /= 0) error stop "open descriptor batch fixture"
      if (first%block_count /= 1 .or. first%block_lengths(1) /= 8193) error stop "unexpected descriptor batch fixture size"
      n = first%block_lengths(1)
      allocate(occupations(max_orbitals,n),quantum(max_orbitals,3,n),couplings(csf_angular_code_count,n))
      call read_descriptor_block(first,1,occupations,quantum,couplings,status)
      if (status /= 0) error stop "read descriptor batches"
      expected = 0_int8
      expected([1,2,3,4],1) = [2_int8,2_int8,1_int8,1_int8]
      expected([8,15,1,3],2) = [2_int8,2_int8,1_int8,1_int8]
      expected([8,15,1,4],3) = [2_int8,2_int8,1_int8,1_int8]
      if (any(occupations(:,[1,8192,8193]) /= expected)) error stop "descriptor batch boundary lost row order"
      if (any(quantum([1,4],3,8193) /= 2_int8) .or. &
         any(couplings([1,4],8193) /= [2_int8,1_int8])) error stop "descriptor final batch lost quantum codes"
      call close_descriptor_source(first)
   end subroutine check_batch_boundary

   subroutine check_distinct_bases()
      character(len=4096) :: alternate
      character(len=64) :: first_basis, second_basis
      integer :: first_csf_count, core_count
      integer(int8), allocatable :: first_row_block(:,:)

      call get_command_argument(3,alternate)
      call open_descriptor_source(first,trim(filename),50,status=status)
      if (status /= 0) error stop "reopen first basis"
      first_basis = first%basis_id
      first_csf_count = first%total_csf_count
      call open_descriptor_source(second,trim(alternate),50,status=status)
      if (status /= 0 .or. second%basis_id == first_basis) error stop "expected distinct basis handles"
      second_basis = second%basis_id

      allocate(first_row_block(max_orbitals,first%block_lengths(1)))
      call read_descriptor_block(first,1,first_row_block,status=status)
      if (status /= 0) error stop "read first basis"
      call read_descriptor_block(second,1,occupations,status=status)
      if (status /= 0 .or. occupations(1,1) /= 2 .or. any(occupations(2:,1) /= 0)) error stop "read second basis"
      call read_descriptor_block(second,0,occupations,status=status)
      if (status == 0 .or. second%is_open() .or. second%basis_id /= "") error stop "bad block left the second source live"
      allocate(other(max_orbitals,first%block_lengths(1)))
      call read_descriptor_block(first,1,other,status=status)
      if (status /= 0 .or. first%basis_id /= first_basis .or. any(other /= first_row_block)) &
         error stop "failure affected the other basis"
      call close_descriptor_source(first)
      call open_descriptor_source(second,trim(alternate),50,status=status)
      if (status /= 0 .or. second%basis_id /= second_basis) error stop "second basis did not recover"
      call close_descriptor_source(second)

      call load_csf_basis(trim(filename),core_count,status)
      if (status /= 0 .or. loaded_csf_basis_id /= first_basis) error stop "install first basis"
      if (size(csf_subshell_occupations,2) /= first_csf_count) error stop "wrong full basis extent"
      call load_csf_basis(trim(alternate),core_count,status)
      if (status /= 0 .or. loaded_csf_basis_id /= second_basis) error stop "replace with second basis"
      if (size(csf_subshell_occupations,2) /= 1 .or. size(csf_subshell_quantum_codes,3) /= 1 .or. &
         size(csf_angular_codes,2) /= 1) error stop "old allocation survived basis replacement"
      call load_csf_basis(trim(alternate)//".missing.parquet",core_count,status)
      if (status == 0 .or. loaded_csf_basis_id /= second_basis) error stop "failed replacement changed the basis identity"
      if (csf_subshell_occupations(1,1) /= 2 .or. any(csf_subshell_occupations(2:,1) /= 0)) &
         error stop "failed replacement changed the installed basis"
      deallocate(csf_subshell_occupations,csf_subshell_quantum_codes,csf_angular_codes)
   end subroutine check_distinct_bases

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
