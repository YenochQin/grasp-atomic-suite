      module diaga4_i
      interface
!
      subroutine diaga4(ja1,ja2,k1,k2,ka,ire,iat,recc)
      use vast_kind_param,  only:  double
      integer, intent(in)       :: ja1,ja2,k1,k2,ka,ire
!      INTEGER, INTENT(OUT)      :: IAT
      integer, intent(inout)      :: iat
      real(double), intent(out) :: recc
      end subroutine
      end interface
      end module
