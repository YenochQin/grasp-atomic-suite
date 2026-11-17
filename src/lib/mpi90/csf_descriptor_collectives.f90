! Collective failure and bounded column transfers for native descriptor data.
module csf_descriptor_collectives
   use iso_fortran_env, only: int8, int64
   use mpi
   use parameter_def, only: max_orbitals
   implicit none (type, external)
   private
   public :: check_descriptor_status, allocate_descriptor_arrays, broadcast_descriptor_arrays
contains
   subroutine check_descriptor_status(status)
      integer, intent(in) :: status
      integer :: failed, combined, ierr
      failed=merge(1,0,status /= 0)
      call MPI_Allreduce(failed,combined,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
      if (combined /= 0 .or. ierr /= MPI_SUCCESS) then
         call MPI_Abort(MPI_COMM_WORLD,1,ierr)
         error stop "Descriptor collective failed"
      end if
   end subroutine check_descriptor_status

   subroutine allocate_descriptor_arrays(rows,occupations,quantum,couplings)
      integer, intent(in) :: rows
      integer(int8), pointer, intent(out) :: occupations(:,:),quantum(:,:,:),couplings(:,:)
      integer :: status
      call check_descriptor_status(merge(0,1,rows > 0))
      allocate(occupations(max_orbitals,rows),stat=status)
      if (status == 0) allocate(quantum(max_orbitals,3,rows),stat=status)
      if (status == 0) allocate(couplings(max_orbitals,rows),stat=status)
      call check_descriptor_status(status)
   end subroutine allocate_descriptor_arrays

   subroutine broadcast_descriptor_arrays(occupations,quantum,couplings,column_limit)
      integer(int8), intent(inout) :: occupations(:,:)
      integer(int8), intent(inout), optional :: quantum(:,:,:),couplings(:,:)
      integer, intent(in), optional :: column_limit
      integer :: dimensions(4), root_dimensions(4), first, last, columns, ierr, status
      integer(int64) :: width
      dimensions=[size(occupations,1),size(occupations,2),merge(1,0,present(quantum)),merge(1,0,present(couplings))]
      root_dimensions=dimensions
      call MPI_Bcast(root_dimensions,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      status=merge(0,1,all(root_dimensions == dimensions))
      if (present(quantum)) then
         if (any(shape(quantum) /= [dimensions(1),3,dimensions(2)])) status=1
      end if
      if (present(couplings)) then
         if (any(shape(couplings) /= dimensions(:2))) status=1
      end if
      if (dimensions(1) <= 0) status=1
      call check_descriptor_status(status)
      width=int(dimensions(1),int64)
      if (present(quantum)) width=3_int64*width
      columns=int(min(8192_int64,int(huge(0),int64)/width))
      if (present(column_limit)) columns=min(columns,column_limit)
      call check_descriptor_status(merge(0,1,columns > 0))
      ! All ranks must use the same segmentation, including optional test limits.
      root_dimensions(1)=columns
      call MPI_Bcast(root_dimensions,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call check_descriptor_status(merge(0,1,columns == root_dimensions(1)))
      first=1
      do while (first <= dimensions(2))
         last=int(min(int(dimensions(2),int64),int(first,int64)+columns-1))
         call MPI_Bcast(occupations(:,first:last),dimensions(1)*(last-first+1),MPI_INTEGER1,0,MPI_COMM_WORLD,ierr)
         if (ierr /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD,1,ierr)
         if (present(quantum)) then
            call MPI_Bcast(quantum(:,:,first:last),3*dimensions(1)*(last-first+1), &
                           MPI_INTEGER1,0,MPI_COMM_WORLD,ierr)
            if (ierr /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD,1,ierr)
         end if
         if (present(couplings)) then
            call MPI_Bcast(couplings(:,first:last),dimensions(1)*(last-first+1),MPI_INTEGER1,0,MPI_COMM_WORLD,ierr)
            if (ierr /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD,1,ierr)
         end if
         if (last == dimensions(2)) exit
         first=last+1
      end do
   end subroutine broadcast_descriptor_arrays
end module csf_descriptor_collectives
