      module perko1_i
      interface
!
      subroutine perko1(ja,bk,ik,bd,id)
      use vast_kind_param, only:  double
      integer,      intent(in)                :: ja
      integer,      intent(out), dimension(7) :: ik, id
      real(double), intent(out), dimension(3) :: bk, bd
      end subroutine
      end interface
      end module
