      module cle0sm_i
      interface
!
      subroutine cle0sm(q,qm,s,c,cm,a)
      use vast_kind_param, only:  double
      real(double), intent(in)  :: q, qm, s, c, cm
      real(double), intent(out) :: a
      end subroutine
      end interface
      end module
