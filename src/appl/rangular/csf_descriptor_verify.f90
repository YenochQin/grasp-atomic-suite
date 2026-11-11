! Compare both loaders directly, before any angular coefficient calculation.
program csf_descriptor_verify
   use iso_fortran_env, only: int8, int64, error_unit
   use parameter_def, only: orbital_capacity => NNNW
   use orb_C, only: occupations => IQA, current_csf_count => NCF, orbital_count => NW, &
                    principal_numbers => NP, relativistic_kappas => NAK, orbital_momenta => NKL, &
                    single_electron_two_j => NKJ, orbital_labels => NH
   use stat_C, only: subshell_quantum_numbers => JQSA, couplings => JCUPA
   use hblock_C, only: block_count => NBLOCK, block_lengths => NCFBLK
   use def_C, only: electron_count => NELEC
   use iounit_C, only: standard_error_unit => ISTDE
   use csf_descriptor_input, only: descriptor_source, open_descriptor_source, close_descriptor_source, &
                                   read_descriptor_block, descriptor_block_label
   use cslh_I, only: cslh
   use lodcsh2_I, only: lodcsh2
   implicit none
   integer, parameter :: maximum_blocks = 50
   type(descriptor_source) :: source
   integer(int8), allocatable :: descriptor_occupations(:,:), descriptor_quantum_numbers(:,:,:), descriptor_couplings(:,:)
   integer :: core_count, block_index, shell_index, csf_index, quantum_channel
   integer :: argument_length
   integer(int64) :: verified_values
   character(len=:), allocatable :: binary_filename, text_filename
   character(len=8) :: block_labels(maximum_blocks)

   if (command_argument_count() /= 2) error stop 'Usage: csf_descriptor_verify DESCRIPTORS.bin ORIGINAL.c'
   call get_command_argument(1, length=argument_length)
   allocate(character(len=argument_length) :: binary_filename)
   call get_command_argument(1, value=binary_filename)
   call get_command_argument(2, length=argument_length)
   allocate(character(len=argument_length) :: text_filename)
   call get_command_argument(2, value=text_filename)
   standard_error_unit = error_unit
   call open_descriptor_source(source, binary_filename, maximum_blocks)
   call cslh(text_filename, core_count, maximum_blocks, block_labels)
   call require_equal(core_count == source%core_count, 'core count')
   call require_equal(orbital_count == source%orbital_count, 'orbital count')
   call require_equal(current_csf_count == source%total_csf_count, 'total CSF count')
   call require_equal(electron_count == source%electron_count, 'electron count')
   call require_equal(block_count == source%block_count, 'block count')
   call require_equal(all(block_lengths(:block_count) == source%block_lengths), 'block lengths')
   call require_equal(all(principal_numbers(:orbital_count) == source%principal_numbers), 'principal quantum numbers')
   call require_equal(all(relativistic_kappas(:orbital_count) == source%relativistic_kappas), 'relativistic kappas')
   call require_equal(all(orbital_momenta(:orbital_count) == source%orbital_momenta), 'orbital angular momenta')
   call require_equal(all(single_electron_two_j(:orbital_count) == source%single_electron_two_j), 'single-electron 2j')
   call require_equal(all(orbital_labels(:orbital_count) == source%orbital_labels), 'orbital labels/order')
   verified_values = 0_int64
   do block_index = 1, block_count
      call require_equal(block_labels(block_index) == descriptor_block_label(source%block_two_j(block_index), &
                         source%block_parities(block_index)), 'block J/parity label/order')
      current_csf_count = block_lengths(block_index)
      allocate(occupations(orbital_capacity,current_csf_count), &
               subshell_quantum_numbers(orbital_capacity,3,current_csf_count), &
               couplings(orbital_capacity,current_csf_count))
      allocate(descriptor_occupations(orbital_capacity,current_csf_count), &
               descriptor_quantum_numbers(orbital_capacity,3,current_csf_count), &
               descriptor_couplings(orbital_capacity,current_csf_count))
      call lodcsh2(21, core_count, block_index)
      call read_descriptor_block(source, block_index, descriptor_occupations, descriptor_quantum_numbers, &
                                 descriptor_couplings)
      do csf_index = 1, current_csf_count
         do shell_index = 1, orbital_capacity
            call compare_value(int(occupations(shell_index,csf_index)), &
                               int(descriptor_occupations(shell_index,csf_index)), 'occupation', 0)
            do quantum_channel = 1, 3
               call compare_value(int(subshell_quantum_numbers(shell_index,quantum_channel,csf_index)), &
                                  int(descriptor_quantum_numbers(shell_index,quantum_channel,csf_index)), &
                                  'subshell quantum number (v,w,2J+1)', quantum_channel)
            end do
            call compare_value(int(couplings(shell_index,csf_index)), &
                               int(descriptor_couplings(shell_index,csf_index)), 'coupling / total J-parity', 0)
         end do
      end do
      verified_values = verified_values + 5_int64*orbital_capacity*current_csf_count
      write(*, '(A,I0,A,I0)') 'Verified block ', block_index, ', CSFs=', current_csf_count
      deallocate(occupations, subshell_quantum_numbers, couplings)
      nullify(occupations, subshell_quantum_numbers, couplings)
      deallocate(descriptor_occupations, descriptor_quantum_numbers, descriptor_couplings)
   end do
   close(21)
   call close_descriptor_source(source)
   deallocate(block_lengths)
   nullify(block_lengths)
   write(*, '(A,I0)') 'Verified array values: ', verified_values
   ! The wrapper requires this marker: legacy text readers can STOP with exit 0.
   write(*, '(A)') 'CSF_DESCRIPTOR_VERIFICATION_SUCCESS'
contains
   subroutine require_equal(condition, field_name)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: field_name
      if (.not. condition) then
         write(error_unit, '(A)') 'Descriptor/text header mismatch: '//field_name
         error stop 1
      end if
   end subroutine require_equal

   subroutine compare_value(text_value, descriptor_value, field_name, channel)
      integer, intent(in) :: text_value, descriptor_value, channel
      character(len=*), intent(in) :: field_name
      if (text_value == descriptor_value) return
      write(error_unit, '(A)') 'Descriptor/text array mismatch: '//field_name
      write(error_unit, '(A,6(I0,1X))') 'block, CSF, shell/slot, channel, text, descriptor: ', &
                                       block_index, csf_index, shell_index, channel, text_value, descriptor_value
      error stop 1
   end subroutine compare_value
end program csf_descriptor_verify
