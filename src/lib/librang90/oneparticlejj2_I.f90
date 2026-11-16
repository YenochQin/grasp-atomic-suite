      module oneparticlejj2_i
      interface
!
      subroutine oneparticlejj2(ns,ka,ja,jb,coeff)
      use vast_kind_param, only:  double
      integer, intent(in) :: ns,ka,ja,jb
      real(double), intent(out) :: coeff
      end subroutine
      end interface
      end module
