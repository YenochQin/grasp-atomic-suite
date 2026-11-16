      module c0t5s_i
      interface
!
      subroutine c0t5s(q, qm, sm, c, cm, a)
      use vast_kind_param, only:  double
      real(double), intent(in)  :: q, qm, sm, c, cm
      real(double), intent(out) :: a
      end subroutine
      end interface
      end module
