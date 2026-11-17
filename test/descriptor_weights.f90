program descriptor_weights
   use iso_fortran_env, only: real64
   use hblock_c, only: symmetry_block_count, csf_counts_by_block
   use syma_c, only: jpgg
   use stat_c, only: csf_coupling_quantum_codes
   use iounit_c, only: stdin_unit, stderr_unit
   use getaldwt_i, only: getaldwt
   implicit none (type, external)
   real(real64) :: weights(3)
   integer :: input_unit
   nullify(csf_coupling_quantum_codes)
   symmetry_block_count=2
   allocate(csf_counts_by_block(0:2))
   csf_counts_by_block=[0,2,1]
   jpgg=0
   jpgg(:2)=[1,-3]
   open(newunit=input_unit,status="scratch",form="formatted")
   write(input_unit,*) 5
   rewind(input_unit)
   stdin_unit=input_unit
   stderr_unit=6
   call getaldwt(3,weights)
   if (maxval(abs(weights-[0.2_real64,0.2_real64,0.6_real64])) > 1.0e-15_real64) &
      error stop "Standard weights disagree with 2J+1"
   if (associated(csf_coupling_quantum_codes)) error stop "Unexpected coupling allocation"
   close(input_unit)
   print *, "PASS standard weights without coupling arrays"
end program descriptor_weights
