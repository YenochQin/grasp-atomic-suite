! Read the requested ASF indices for each symmetry block and record the largest index needed.
module block_level_selection
    implicit none (type, external)
    private
    public :: lodstate
contains
      subroutine lodstate(block_labels)
      use array_memory_management, only: dalloc
      use atomic_calculation_state, only: selected_level_indices => iccmin, selected_level_count => ncmin
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, &
         & highest_selected_level_by_block => ncmaxblk, symmetry_block_count
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use level_selection_parser, only: items
      implicit none (type, external)
      character (len = 8) :: block_labels(*)
      integer :: display_block_index, selection_index, block_index, block_csf_count, previous_selection_count, &
         & parse_status, highest_selected_level
      character :: selection_input*256

      write (stderr_unit, *) "There are ", symmetry_block_count, " blocks  ", &
         "(block   J/Parity   NCF):"
      write (stderr_unit, "( 4(I3, 1X, A5, I8, 5X) )") (display_block_index,block_labels(display_block_index)(1:5), &
         & csf_counts_by_block(display_block_index),display_block_index=1&
         ,symmetry_block_count)

      eigenvector_counts_by_block = 0
      highest_selected_level_by_block = 0

      write (stderr_unit, *)
      write (stderr_unit, *) "Enter ASF serial numbers for each block"

      selected_level_count = 0
  123 continue
      do block_index = 1, symmetry_block_count
         block_csf_count = csf_counts_by_block(block_index)
  234    continue
         write (stderr_unit, *) "Block ", block_index, "   ncf = ", block_csf_count, " id = ", block_labels(&
            block_index)(1:5)

         read (*, "(A)") selection_input
         write(734,"(a)") trim(selection_input)

         previous_selection_count = selected_level_count
         call items (selected_level_count, block_csf_count, selection_input, parse_status)
         if (selected_level_count == 0) call dalloc (selected_level_indices, "ICCMIN", "LODSTATE" )
         if (parse_status < 0) goto 234
         eigenvector_counts_by_block(block_index) = selected_level_count - previous_selection_count

         highest_selected_level = 0
         do selection_index = previous_selection_count + 1, selected_level_count
            highest_selected_level = max(highest_selected_level,selected_level_indices(selection_index))
         end do
         highest_selected_level_by_block(block_index) = highest_selected_level
      end do

      if (selected_level_count == 0) then
         write (stderr_unit, *) "At least one state should be selected"
         goto 123
      end if

      return
      end subroutine lodstate
end module block_level_selection
