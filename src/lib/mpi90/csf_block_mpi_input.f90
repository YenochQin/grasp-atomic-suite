! Root-owned RCI stream and replicated block metadata/arrays.
module csf_block_mpi_input
   use csf_block_input, only: open_csf_blocks, load_csf_block
   use csf_descriptor_collectives, only: check_descriptor_status, broadcast_descriptor_arrays
   use mpi_c, only: mpi_rank, mpi_bcast, mpi_integer, mpi_character, mpi_comm_world
   use orb_c, only: csf_count, orbital_count, orbital_principal_n, orbital_kappa, orbital_l, orbital_two_j, &
                    orbital_symmetry_labels, csf_subshell_occupations
   use stat_c, only: csf_subshell_quantum_codes, csf_coupling_quantum_codes
   use def_c, only: electron_count
   use hblock_c, only: symmetry_block_count, csf_counts_by_block
   implicit none (type, external)
   private
   public :: open_csf_blocks_mpi, load_csf_block_mpi
contains
   subroutine open_csf_blocks_mpi(name,ncore,labels)
      character(len=*), intent(in) :: name
      integer, intent(out) :: ncore
      character(len=8), intent(out) :: labels(:)
      integer :: status, ierr, info(5)
      status=0
      if (mpi_rank == 0) call open_csf_blocks(name,ncore,labels,status)
      call check_descriptor_status(status)
      if (mpi_rank == 0) info=[csf_count,orbital_count,ncore,electron_count,symmetry_block_count]
      call mpi_bcast(info,5,mpi_integer,0,mpi_comm_world,ierr)
      call check_descriptor_status(ierr)
      csf_count=info(1)
      orbital_count=info(2)
      ncore=info(3)
      electron_count=info(4)
      symmetry_block_count=info(5)
      call check_descriptor_status(merge(0,1,size(labels) >= symmetry_block_count))
      if (mpi_rank /= 0) then
         if (associated(csf_counts_by_block)) deallocate(csf_counts_by_block)
         allocate(csf_counts_by_block(0:symmetry_block_count),stat=status)
         labels=""
      end if
      call check_descriptor_status(status)
      call mpi_bcast(csf_counts_by_block,symmetry_block_count+1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(labels,8*symmetry_block_count,mpi_character,0,mpi_comm_world,ierr)
      call mpi_bcast(orbital_principal_n,orbital_count,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(orbital_kappa,orbital_count,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(orbital_l,orbital_count,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(orbital_two_j,orbital_count,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(orbital_symmetry_labels,2*orbital_count,mpi_character,0,mpi_comm_world,ierr)
      call check_descriptor_status(ierr)
   end subroutine open_csf_blocks_mpi

   subroutine load_csf_block_mpi(block_index)
      integer, intent(in) :: block_index
      integer :: status
      status=0
      if (mpi_rank == 0) call load_csf_block(block_index,status)
      call check_descriptor_status(status)
      call broadcast_descriptor_arrays(csf_subshell_occupations,csf_subshell_quantum_codes,csf_coupling_quantum_codes)
   end subroutine load_csf_block_mpi
end module csf_block_mpi_input
