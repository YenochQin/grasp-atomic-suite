      module recop2_i
      interface
!
      subroutine recop2(ns,ja1,ja2,k1,k2,ka,ire,iat,recc)
      use vast_kind_param,  only:  double
      integer, intent(in)       :: ns, ja1, ja2, k1, k2, ka, ire
      integer, intent(out)      :: iat
      real(double), intent(out) :: recc
      end subroutine
      end interface
      end module
