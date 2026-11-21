! Compare both loaders directly, before any angular coefficient calculation.
program csf_descriptor_verify
   use iso_fortran_env, only: int8, int64, error_unit
   use grasp_parameters, only: max_orbitals
   use csf_coupling_data, only: csf_angular_code_count, csf_total_j_slot, csf_parity_slot
   use orbital_csf_basis, only: csf_subshell_occupations, csf_count, orbital_count, orbital_principal_n, orbital_kappa, &
       orbital_l, orbital_two_j, orbital_symmetry_labels
   use csf_text_oracle_coupling_data, only: csf_subshell_quantum_codes, csf_coupling_quantum_codes
   use hamiltonian_symmetry_blocks, only: symmetry_block_count, csf_counts_by_block
   use atomic_calculation_state, only: electron_count
   use standard_io_units, only: stderr_unit
   use csf_descriptor_reader, only: descriptor_source, open_descriptor_source, close_descriptor_source, &
                                   read_descriptor_block, descriptor_block_label
   use cslh_i, only: cslh
   use lodcsh2_i, only: lodcsh2
   implicit none (type, external)
   integer, parameter :: maximum_blocks = 50
   type(descriptor_source) :: source
   integer(int8), allocatable :: descriptor_occupations(:,:), descriptor_quantum_numbers(:,:,:), descriptor_couplings(:,:)
   integer :: core_count, block_index, shell_index, csf_index, quantum_channel
   integer :: argument_length, csf_input_unit, cumulative, open_shells, occupation, signed_symmetry
   integer(int64) :: verified_values
   character(len=:), allocatable :: binary_filename, text_filename
   character(len=8) :: block_labels(maximum_blocks)

   if (command_argument_count() /= 2) error stop "Usage: csf_descriptor_verify DESCRIPTORS.parquet ORIGINAL.c"
   call get_command_argument(1, length=argument_length)
   allocate(character(len=argument_length) :: binary_filename)
   call get_command_argument(1, value=binary_filename)
   call get_command_argument(2, length=argument_length)
   allocate(character(len=argument_length) :: text_filename)
   call get_command_argument(2, value=text_filename)
   stderr_unit = error_unit
   call open_descriptor_source(source, binary_filename, maximum_blocks)
   call cslh(text_filename, core_count, maximum_blocks, block_labels, csf_input_unit)
   call require_equal(core_count == source%core_count, "core count")
   call require_equal(orbital_count == source%orbital_count, "orbital count")
   call require_equal(csf_count == source%total_csf_count, "total CSF count")
   call require_equal(electron_count == source%electron_count, "electron count")
   call require_equal(symmetry_block_count == source%block_count, "block count")
   call require_equal(all(csf_counts_by_block(:symmetry_block_count) == source%block_lengths), "block lengths")
   call require_equal(all(orbital_principal_n(:orbital_count) == source%principal_numbers), "principal quantum numbers")
   call require_equal(all(orbital_kappa(:orbital_count) == source%relativistic_kappas), "relativistic kappas")
   call require_equal(all(orbital_l(:orbital_count) == source%orbital_momenta), "orbital angular momenta")
   call require_equal(all(orbital_two_j(:orbital_count) == source%single_electron_two_j), "single-electron 2j")
   call require_equal(all(orbital_symmetry_labels(:orbital_count) == source%orbital_labels), "orbital labels/order")
   verified_values = 0_int64
   do block_index = 1, symmetry_block_count
      call require_equal(block_labels(block_index) == descriptor_block_label(source%block_two_j(block_index), &
                         source%block_parities(block_index)), "block J/parity label/order")
      csf_count = csf_counts_by_block(block_index)
      allocate(csf_subshell_occupations(max_orbitals,csf_count), &
               csf_subshell_quantum_codes(max_orbitals,3,csf_count), &
               csf_coupling_quantum_codes(max_orbitals,csf_count))
      allocate(descriptor_occupations(max_orbitals,csf_count), &
               descriptor_quantum_numbers(max_orbitals,3,csf_count), &
               descriptor_couplings(csf_angular_code_count,csf_count))
      call lodcsh2(csf_input_unit, core_count, block_index)
      call read_descriptor_block(source, block_index, descriptor_occupations, descriptor_quantum_numbers, &
                                 descriptor_couplings)
      do csf_index = 1, csf_count
         cumulative = 1
         open_shells = 0
         do shell_index = 1, max_orbitals
            call compare_value(int(csf_subshell_occupations(shell_index,csf_index)), &
                               int(descriptor_occupations(shell_index,csf_index)), "occupation", 0)
            do quantum_channel = 1, 3
               call compare_value(int(csf_subshell_quantum_codes(shell_index,quantum_channel,csf_index)), &
                                  int(descriptor_quantum_numbers(shell_index,quantum_channel,csf_index)), &
                                  "subshell quantum number (v,w,2J+1)", quantum_channel)
            end do
            if (shell_index <= orbital_count) then
               occupation = csf_subshell_occupations(shell_index,csf_index)
               if (occupation > 0 .and. occupation < orbital_two_j(shell_index)+1) then
                  open_shells = open_shells + 1
                  if (open_shells == 1) then
                     cumulative = csf_subshell_quantum_codes(shell_index,3,csf_index)
                  else
                     cumulative = csf_coupling_quantum_codes(open_shells-1,csf_index)
                  end if
               end if
            else
               cumulative = 0
            end if
            call compare_value(cumulative, int(descriptor_couplings(shell_index,csf_index)), "cumulative coupling", 0)
         end do
         signed_symmetry = csf_coupling_quantum_codes(max_orbitals,csf_index)
         call compare_value(abs(signed_symmetry),int(descriptor_couplings(csf_total_j_slot,csf_index)), "total J", 0)
         call compare_value(sign(1,signed_symmetry),int(descriptor_couplings(csf_parity_slot,csf_index)), "parity", 0)
      end do
      verified_values = verified_values + (5_int64*max_orbitals+2)*csf_count
      write(*, "(A,I0,A,I0)") "Verified block ", block_index, ", CSFs=", csf_count
      deallocate(csf_subshell_occupations, csf_subshell_quantum_codes, csf_coupling_quantum_codes)
      nullify(csf_subshell_occupations, csf_subshell_quantum_codes, csf_coupling_quantum_codes)
      deallocate(descriptor_occupations, descriptor_quantum_numbers, descriptor_couplings)
   end do
   close(csf_input_unit)
   call close_descriptor_source(source)
   deallocate(csf_counts_by_block)
   nullify(csf_counts_by_block)
   write(*, "(A,I0)") "Verified array values: ", verified_values
   ! The wrapper requires this marker: legacy text readers can STOP with exit 0.
   write(*, "(A)") "CSF_DESCRIPTOR_VERIFICATION_SUCCESS"
contains
   subroutine require_equal(condition, field_name)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: field_name
      if (.not. condition) then
         write(error_unit, "(A)") "Descriptor/text header mismatch: "//field_name
         error stop 1
      end if
   end subroutine require_equal

   subroutine compare_value(text_value, descriptor_value, field_name, channel)
      integer, intent(in) :: text_value, descriptor_value, channel
      character(len=*), intent(in) :: field_name
      if (text_value == descriptor_value) return
      write(error_unit, "(A)") "Descriptor/text array mismatch: "//field_name
      write(error_unit, "(A,6(I0,1X))") "block, CSF, shell/slot, channel, text, descriptor: ", &
                                       block_index, csf_index, shell_index, channel, text_value, descriptor_value
      error stop 1
   end subroutine compare_value
end program csf_descriptor_verify
