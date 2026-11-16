      module recop1_i
      interface
!
      subroutine recop1(ns,ja1,ka,ire,iat,recc)
      use vast_kind_param,  only:  double
      integer, intent(in)       :: ns, ja1, ka, ire
      integer, intent(out)      :: iat
      real(double), intent(out) :: recc
      end subroutine
      end interface
      end module
