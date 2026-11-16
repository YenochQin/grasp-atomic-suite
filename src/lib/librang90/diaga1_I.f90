      module diaga1_i
      interface
!
      subroutine diaga1(ja1,ka,ire,iat,recc)
      use vast_kind_param,  only:  double
      integer, intent(in)       :: ja1,ka,ire
!      INTEGER, INTENT(OUT)      :: IAT
      integer, intent(inout)      :: iat
      real(double), intent(out) :: recc
      end subroutine
      end interface
      end module
