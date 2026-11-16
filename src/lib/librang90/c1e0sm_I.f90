      module c1e0sm_i
      interface
!
      subroutine c1e0sm(q,qm,c,cm,a)
      use vast_kind_param, only:  double
      real(double), intent(in)  :: q, qm, c, cm
      real(double), intent(out) :: a
      end subroutine
      end interface
      end module
