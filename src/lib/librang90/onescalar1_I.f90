      module onescalar1_i
      interface
!
      subroutine onescalar1(ns,jja,jjb,ja,jb,coeff)
      use vast_kind_param, only:  double
      integer, intent(in) :: ns,jja,jjb,ja,jb
      real(double), intent(out) :: coeff
      end subroutine
      end interface
      end module
