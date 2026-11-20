      module outsda_i
      implicit none (type, external)
      interface
      subroutine outsda (print_details, nonzero_element_count, csf_count)
      implicit none (type, external)
      logical, intent(in) :: print_details
      integer, intent(in) :: nonzero_element_count
      integer, intent(in) :: csf_count
      end subroutine outsda
      end interface
      end module outsda_i
