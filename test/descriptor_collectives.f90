program descriptor_collectives
   use iso_fortran_env, only: int8
   use mpi_c, only: mpi_init, mpi_comm_rank, mpi_comm_world, mpi_finalize
   use csf_descriptor_collectives, only: allocate_descriptor_arrays, broadcast_descriptor_arrays, check_descriptor_status
   implicit none (type, external)
   integer(int8), pointer :: occupations(:,:),quantum(:,:,:),couplings(:,:)
   integer :: rank, status, i, j
   character(len=16) :: scenario
   call mpi_init(status)
   call mpi_comm_rank(mpi_comm_world,rank,status)
   call get_command_argument(1,scenario)
   status=0
   if (scenario == "error" .and. rank == 1) status=-1
   call check_descriptor_status(status)
   call allocate_descriptor_arrays(3,occupations,quantum,couplings)
   occupations=0_int8
   quantum=0_int8
   couplings=0_int8
   if (rank == 0) then
      do j=1,3
         do i=1,size(occupations,1)
            occupations(i,j)=int(mod(i+j,127),int8)
         end do
      end do
      couplings=-occupations
      quantum(:,1,:)=occupations
      quantum(:,2,:)=couplings
      quantum(:,3,:)=1_int8
   end if
   if (scenario == "shape" .and. rank == 1) then
      call broadcast_descriptor_arrays(occupations(:,:2),quantum,couplings,column_limit=1)
   else
      call broadcast_descriptor_arrays(occupations,quantum,couplings,column_limit=1)
   end if
   do j=1,3
      do i=1,size(occupations,1)
         if (occupations(i,j) /= int(mod(i+j,127),int8)) error stop "occupation broadcast"
      end do
   end do
   if (any(couplings /= -occupations) .or. any(quantum(:,1,:) /= occupations) .or. &
       any(quantum(:,2,:) /= couplings) .or. any(quantum(:,3,:) /= 1)) error stop "quantum broadcast"
   deallocate(occupations,quantum,couplings)
   call mpi_finalize(status)
   print *, "PASS one-column descriptor collectives"
end program descriptor_collectives
