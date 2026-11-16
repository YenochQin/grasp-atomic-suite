      module a1jj_i
      interface
!
      subroutine a1jj (ik, bk, id, bd, qm1, a)
      use vast_kind_param, only:  double
      integer,      intent(in), dimension(7) :: ik, id
      real(double), intent(in), dimension(3) :: bk, bd
      real(double), intent(in)               :: qm1
      real(double), intent(out)              :: a
      end subroutine
      end interface
      end module
