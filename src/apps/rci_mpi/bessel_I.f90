      module bessel_i
      implicit none (type, external)
      interface
      subroutine bessel (orbital_a, orbital_b, ik, iw, k)
      implicit none (type, external)
      integer, intent(in) :: orbital_a
      integer, intent(in) :: orbital_b
      integer, intent(in) :: ik
      integer, intent(in) :: iw
      integer, intent(in) :: k
      end subroutine bessel
      end interface
      end module bessel_i
