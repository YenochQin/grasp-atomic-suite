module rmcdhf_descriptor_mpi
   use csf_descriptor_collectives, only: broadcast_descriptor_arrays, check_all => check_descriptor_status
   use mpi
   use mpi_c, only: mpi_rank
   use rmcdhf_descriptor_input, only: load_rmcdhf_descriptors
   use parameter_def, only: max_orbitals
   use orb_c, only: csf_count, orbital_count, csf_subshell_occupations, orbital_principal_n, orbital_kappa, &
                    orbital_l, orbital_two_j, orbital_symmetry_labels
   use def_c, only: electron_count
   use core_c, only: closed_core_subshell_count
   use hblock_c, only: symmetry_block_count, csf_counts_by_block
   use syma_c, only: jpgg, nblk0
   implicit none (type, external)
   private
   public :: load_rmcdhf_descriptors_mpi
contains


   subroutine load_rmcdhf_descriptors_mpi(name,ncore,labels)
      character(len=*), intent(in) :: name
      integer, intent(out) :: ncore
      character(len=8), intent(in) :: labels(:)
      integer :: status, ierr, info(5), lengths(nblk0)
      character(len=8) :: root_labels(nblk0)
      status = 0
      if (mpi_rank == 0) call load_rmcdhf_descriptors(name,ncore,labels,status)
      call check_all(status)
      if (mpi_rank == 0) then
         info = [csf_count,orbital_count,electron_count,symmetry_block_count,ncore]
         lengths = 0
         lengths(:symmetry_block_count) = csf_counts_by_block(1:symmetry_block_count)
         root_labels = ""
         root_labels(:symmetry_block_count) = labels
      end if
      call MPI_Bcast(info,5,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if (any(info /= [csf_count,orbital_count,electron_count,symmetry_block_count,closed_core_subshell_count])) status=1
      call check_all(status)
      call MPI_Bcast(lengths,nblk0,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(root_labels,8*nblk0,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      if (any(lengths(:symmetry_block_count) /= csf_counts_by_block(1:symmetry_block_count))) status=1
      if (any(root_labels(:symmetry_block_count) /= labels)) status=1
      call check_all(status)
      ncore=info(5)
      if (mpi_rank /= 0) allocate(csf_subshell_occupations(max_orbitals,csf_count),stat=status)
      call check_all(status)
      call MPI_Bcast(orbital_principal_n,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(orbital_kappa,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(orbital_l,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(orbital_two_j,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(orbital_symmetry_labels,2*orbital_count,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(jpgg,nblk0,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call broadcast_descriptor_arrays(csf_subshell_occupations)
   end subroutine load_rmcdhf_descriptors_mpi
end module rmcdhf_descriptor_mpi
