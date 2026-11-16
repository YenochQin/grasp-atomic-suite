      module rec3_i
      interface
!
      subroutine rec3(ja1,ja2,ja3,k1,k2,k3,ire,iat,recc)
      use vast_kind_param,  only:  double
      integer, intent(in)       :: ja1,ja2,ja3,k1,k2,k3,ire
      integer, intent(out)      :: iat
      real(double), intent(out) :: recc
      end subroutine
      end interface
      end module
