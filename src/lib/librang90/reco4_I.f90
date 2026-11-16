      module reco4_i
      interface
!
      subroutine reco4(ja1,ja2,ja3,ja4,k1,k2,k3,k4,ka,ire,iat,recc)
      use vast_kind_param,  only:  double
      integer, intent(in)       :: ja1,ja2,ja3,ja4,k1,k2,k3,k4,ka,ire
      integer, intent(out)      :: iat
      real(double), intent(out) :: recc
      end subroutine
      end interface
      end module
