      module maxarr_i
      implicit none (type, external)
      interface
      subroutine maxarr (orbital_index)
      implicit none (type, external)
      integer, intent(out) :: orbital_index
      end subroutine maxarr
      end interface
      end module maxarr_i
