      module settmpgg_i
      implicit none (type, external)
      interface
      subroutine settmpgg (nb, k, filehead)
      implicit none (type, external)
      integer, intent(in) :: nb
      integer, intent(in) :: k
      character (len = *), intent(in) :: filehead
      end subroutine settmpgg
      end interface
      end module settmpgg_i
