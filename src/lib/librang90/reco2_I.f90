      module reco2_i
      interface
!
      subroutine reco2(ja1,ja2,ka,ire,iat,recc)
      use vast_kind_param,  only:  double
      integer, intent(in)       :: ja1,ja2,ka,ire
      integer, intent(out)      :: iat
      real(double), intent(out) :: recc
      end subroutine
      end interface
      end module
