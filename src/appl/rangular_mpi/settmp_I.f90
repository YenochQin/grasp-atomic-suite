      module settmp_i
      implicit none (type, external)
      interface
      subroutine settmp (nb, kmax, filehead)
      implicit none (type, external)
      integer, intent(in) :: nb
      integer, intent(in) :: kmax
      character (len = *), intent(in) :: filehead
      end subroutine
      end interface
      end module
