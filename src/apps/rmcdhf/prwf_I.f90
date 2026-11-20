      module prwf_i
      implicit none (type, external)
      interface
      subroutine prwf (orbital_index)
      implicit none (type, external)
      integer, intent(in) :: orbital_index
      end subroutine prwf
      end interface
      end module prwf_i
