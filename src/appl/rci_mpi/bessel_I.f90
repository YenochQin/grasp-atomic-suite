      module bessel_i
      implicit none (type, external)
      interface
      subroutine bessel (ia, ib, ik, iw, k)
      implicit none (type, external)
      integer, intent(in) :: ia
      integer, intent(in) :: ib
      integer, intent(in) :: ik
      integer, intent(in) :: iw
      integer, intent(in) :: k
      end subroutine
      end interface
      end module
