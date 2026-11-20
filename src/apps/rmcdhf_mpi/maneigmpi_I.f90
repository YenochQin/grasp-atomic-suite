      module maneigmpi_i
      implicit none (type, external)
      interface
      subroutine maneigmpi (dvdfirst,print_details, symmetry_block_index, ncfpat, ncminpat, nevecpat, total_csf_count)
      implicit none (type, external)
      logical, intent(in) ::  dvdfirst
      logical :: print_details
!VAST...Dummy argument LPRINT is not referenced in this routine.
      integer, intent(in) :: symmetry_block_index
      integer, intent(in) :: ncfpat
      integer, intent(in) :: ncminpat
      integer, intent(in) :: nevecpat
      integer :: total_csf_count
      end subroutine maneigmpi
      end interface
      end module maneigmpi_i
