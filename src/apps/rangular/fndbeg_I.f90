      module fndbeg_i
      implicit none (type, external)
      interface
      subroutine fndbeg (first_bra_csf_index, first_ket_csf_index, matrix_element_position, one_electron_coefficient_count, &
         & two_electron_counts_by_rank)
      implicit none (type, external)
      integer, intent(inout) :: first_bra_csf_index
      integer, intent(inout) :: first_ket_csf_index
      integer, intent(inout) :: matrix_element_position
      integer, intent(out) :: one_electron_coefficient_count
      integer, pointer, intent(in) :: two_electron_counts_by_rank(:)
      end subroutine fndbeg
      end interface
      end module fndbeg_i
