program descriptor_collectives
   use iso_fortran_env, only: int8
   use csf_coupling_data, only: csf_total_j_slot, csf_parity_slot
   use grasp_parameters, only: max_orbitals
   use mpi_state, only: mpi_init, mpi_comm_rank, mpi_comm_world, mpi_finalize
   use mpi_csf_collectives, only: allocate_descriptor_arrays, broadcast_descriptor_arrays, require_collective_success
   implicit none (type, external)
   integer(int8), allocatable :: occupations(:,:),quantum(:,:,:),couplings(:,:)
   integer :: rank, status, i, j
   integer(int8) :: strided(2*max_orbitals,3), empty(max_orbitals,0)
   character(len=16) :: scenario
   call mpi_init(status)
   call mpi_comm_rank(mpi_comm_world,rank,status)
   call get_command_argument(1,scenario)
   status=0
   if (scenario == "error" .and. rank == 1) status=-1
   call require_collective_success(status)
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
      couplings(:max_orbitals,:)=-occupations
      couplings(csf_total_j_slot,:)=[1,2,127]
      couplings(csf_parity_slot,:)=[1,-1,-1]
      quantum(:,1,:)=occupations
      quantum(:,2,:)=couplings(:max_orbitals,:)
      quantum(:,3,:)=1_int8
   end if
   if (scenario == "shape" .and. rank == 1) then
      call broadcast_descriptor_arrays(occupations(:,:2),quantum,couplings,csf_batch_limit=1)
   else
      call broadcast_descriptor_arrays(occupations,quantum,couplings,csf_batch_limit=1)
   end if
   do j=1,3
      do i=1,size(occupations,1)
         if (occupations(i,j) /= int(mod(i+j,127),int8)) error stop "occupation broadcast"
      end do
   end do
   if (any(couplings(:max_orbitals,:) /= -occupations) .or. any(quantum(:,1,:) /= occupations) .or. &
       any(quantum(:,2,:) /= couplings(:max_orbitals,:)) .or. any(quantum(:,3,:) /= 1)) error stop "quantum broadcast"
   if (any(couplings(csf_total_j_slot,:) /= [1,2,127]) .or. any(couplings(csf_parity_slot,:) /= [1,-1,-1])) &
       error stop "symmetry broadcast"
   strided = -1_int8
   if (rank == 0) strided(1::2,:) = occupations
   call broadcast_descriptor_arrays(strided(1::2,:),csf_batch_limit=2)
   if (any(strided(1::2,:) /= occupations) .or. any(strided(2::2,:) /= -1_int8)) error stop "strided broadcast"
   call broadcast_descriptor_arrays(empty,csf_batch_limit=1)
   deallocate(occupations,quantum,couplings)
   call mpi_finalize(status)
   print *, "PASS one-column descriptor collectives"
end program descriptor_collectives
