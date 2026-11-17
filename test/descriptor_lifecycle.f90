program descriptor_lifecycle
   use iso_fortran_env, only: int8
   use iso_c_binding, only: c_associated
   use parameter_def, only: max_orbitals
   use csf_descriptor_input, only: descriptor_source, open_descriptor_source, close_descriptor_source, read_descriptor_block
   implicit none (type, external)
   type(descriptor_source) :: first, second
   character(len=4096) :: filename
   integer :: status, b, n
   integer(int8), allocatable :: occupations(:,:), other(:,:), quantum(:,:,:), couplings(:,:)
   call get_command_argument(1,filename)
   call open_descriptor_source(first,trim(filename)//".missing",50,status=status)
   if (status == 0 .or. first%is_open .or. c_associated(first%handle)) error stop "failed open retained state"
   call open_descriptor_source(first,trim(filename),50,status=status)
   if (status /= 0) error stop "valid open failed"
   call open_descriptor_source(second,trim(filename),50,status=status)
   if (status /= 0) error stop "second handle open failed"
   do b=1,first%block_count
      n=first%block_lengths(b)
      allocate(occupations(max_orbitals,n),other(max_orbitals,n),quantum(max_orbitals,3,n),couplings(max_orbitals,n))
      call read_descriptor_block(first,b,occupations,quantum,couplings,status)
      if (status /= 0) error stop "full block failed"
      if (first%peel_count == 0) then
         if (sum(int(occupations(:,1))) /= first%electron_count) error stop "Closed-core electron count"
         if (any(quantum(:first%core_count,1:2,:) /= 0) .or. &
             any(quantum(:first%core_count,3,:) /= 1) .or. any(couplings(max_orbitals,:) /= 1)) &
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
   print *, "PASS descriptor lifecycle and occupation-only equivalence"
end program descriptor_lifecycle
