      module izas1_i
      interface
!
      integer function izas1(ib,qb,ik,qk)
      use vast_kind_param, only:  double
      integer,      intent(in) :: ib, ik
      real(double), intent(in) :: qb, qk
      end function
      end interface
      end module
