module csf_basis_mpi
   use csf_descriptor_collectives, only: broadcast_descriptor_arrays, check_all => check_descriptor_status
   use mpi
   use mpi_c, only: mpi_rank
   use csf_basis_input, only: load_csf_basis
   use parameter_def, only: max_orbitals
   use orb_c, only: csf_count, orbital_count, csf_subshell_occupations, orbital_principal_n, orbital_kappa, &
                    orbital_l, orbital_two_j, orbital_symmetry_labels
   use def_c, only: electron_count
   use stat_c, only: csf_subshell_quantum_codes, csf_coupling_quantum_codes
   use blk_c, only: nblock, ncfblk
   implicit none (type, external)
   private
   public :: load_csf_basis_mpi
contains


   subroutine load_csf_basis_mpi(name,ncore)
      character(len=*), intent(in) :: name
      integer, intent(out) :: ncore
      integer :: status, ierr, info(5)
      status = 0
      if (mpi_rank == 0) call load_csf_basis(name,ncore,status)
      call check_all(status)
      if (mpi_rank == 0) info=[csf_count,orbital_count,electron_count,nblock,ncore]
      call MPI_Bcast(info,5,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      csf_count=info(1)
      orbital_count=info(2)
      electron_count=info(3)
      nblock=info(4)
      ncore=info(5)
      call MPI_Bcast(ncfblk,size(ncfblk),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if (mpi_rank /= 0) then
         allocate(csf_subshell_occupations(max_orbitals,csf_count),stat=status)
         if (status == 0) allocate(csf_subshell_quantum_codes(max_orbitals,3,csf_count),stat=status)
         if (status == 0) allocate(csf_coupling_quantum_codes(max_orbitals,csf_count),stat=status)
      end if
      call check_all(status)
      call MPI_Bcast(orbital_principal_n,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(orbital_kappa,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(orbital_l,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(orbital_two_j,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(orbital_symmetry_labels,2*orbital_count,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call broadcast_descriptor_arrays(csf_subshell_occupations,csf_subshell_quantum_codes,csf_coupling_quantum_codes)
   end subroutine load_csf_basis_mpi
end module csf_basis_mpi
