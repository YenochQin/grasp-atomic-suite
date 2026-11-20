program descriptor_weights
   use, intrinsic :: iso_fortran_env, only: real64, output_unit
   use hamiltonian_symmetry_blocks, only: symmetry_block_count, csf_counts_by_block
   use level_symmetry_metadata, only: packed_block_symmetries
   use csf_coupling_data, only: csf_coupling_quantum_codes
   use standard_io_units, only: stdin_unit, stderr_unit
   use read_al_level_weights_interface, only: read_al_level_weights
   implicit none (type, external)
   real(real64) :: weights(3)
   integer :: input_unit
   nullify(csf_coupling_quantum_codes)
   symmetry_block_count=2
   allocate(csf_counts_by_block(0:2))
   csf_counts_by_block=[0,2,1]
   packed_block_symmetries=0
   packed_block_symmetries(:2)=[1,-3]
   open(newunit=input_unit,status="scratch",form="formatted")
   write(input_unit,*) 5
   rewind(input_unit)
   stdin_unit=input_unit
   stderr_unit=output_unit
   call read_al_level_weights(3,weights)
   if (maxval(abs(weights-[0.2_real64,0.2_real64,0.6_real64])) > 1.0e-15_real64) &
      error stop "Standard weights disagree with 2J+1"
   if (associated(csf_coupling_quantum_codes)) error stop "Unexpected coupling allocation"
   close(input_unit)
   print *, "PASS standard weights without coupling arrays"
end program descriptor_weights
