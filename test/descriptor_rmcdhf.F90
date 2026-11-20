program descriptor_rmcdhf
   use closed_core_subshells, only: closed_core_subshell_count
   use, intrinsic :: iso_fortran_env, only: standard_output_unit => output_unit
   use orbital_csf_basis, only: csf_count, orbital_count, csf_subshell_occupations
   use csf_coupling_data, only: csf_subshell_quantum_codes, csf_coupling_quantum_codes
   use hamiltonian_symmetry_blocks, only: csf_counts_by_block, symmetry_block_count
   use level_symmetry_metadata, only: jpgg, nblk0
   use standard_io_units, only: stderr_unit
   use setmcp_i, only: setmcp
#ifdef DESCRIPTOR_MPI
   use mpi_state, only: mpi_init, mpi_comm_rank, mpi_comm_size, mpi_comm_world, mpi_rank, mpi_process_count, mpi_finalize
   use setcslmpi_i, only: setcslmpi
#else
   use mpi_s, only: mpi_rank, mpi_process_count
   use setcsl_i, only: setcsl
#endif
   implicit none (type, external)
   integer :: ncore, status, output_unit
   character(len=8) :: labels(nblk0)
   character(len=3) :: rank_text
   character(len=64) :: mcp_prefix, output_path
   nullify(csf_subshell_quantum_codes,csf_coupling_quantum_codes)
   stderr_unit=standard_output_unit
#ifdef DESCRIPTOR_MPI
   call mpi_init(status)
   call mpi_comm_rank(mpi_comm_world,mpi_rank,status)
   call mpi_comm_size(mpi_comm_world,mpi_process_count,status)
   write(rank_text,'(I3.3)') mpi_rank
   mcp_prefix=rank_text//"/mcp"
#else
   mpi_rank=0
   mpi_process_count=1
   rank_text="000"
   mcp_prefix="mcp"
#endif
   call setmcp(closed_core_subshell_count,nblk0,labels,trim(mcp_prefix))
#ifdef DESCRIPTOR_MPI
   call setcslmpi("rcsf_descriptors.parquet",ncore,labels)
#else
   call setcsl("rcsf_descriptors.parquet",ncore,labels)
#endif
   if (ncore /= 1 .or. csf_count /= 10 .or. orbital_count /= 7) error stop "basis dimensions"
   if (symmetry_block_count /= 3 .or. lbound(csf_counts_by_block,1) /= 0) error stop "block bounds"
   if (any(csf_counts_by_block /= [0,5,2,3])) error stop "MCP block lengths overwritten"
   if (any(jpgg(:3) /= [1,5,-3])) error stop "block symmetry"
   if (associated(csf_subshell_quantum_codes) .or. associated(csf_coupling_quantum_codes)) &
      error stop "full coupling allocation in RMCDHF"
   output_path="occupations."//rank_text//".bin"
   open(newunit=output_unit,file=trim(output_path),access="stream",form="unformatted",status="replace")
   write(output_unit) csf_subshell_occupations(:orbital_count,:)
   close(output_unit)
#ifdef DESCRIPTOR_MPI
   call mpi_finalize(status)
#endif
   print *, "PASS RMCDHF descriptor contract"
end program descriptor_rmcdhf
