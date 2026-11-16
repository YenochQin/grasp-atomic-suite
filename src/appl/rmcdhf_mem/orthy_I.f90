      module orthy_i
      implicit none (type, external)
      interface
      subroutine orthy (nw, jp, lsort)
      implicit none (type, external)
      integer, intent(in) :: nw
      integer, intent(in) :: jp
      logical, intent(in) :: lsort
      end subroutine orthy
      end interface
      end module orthy_i
