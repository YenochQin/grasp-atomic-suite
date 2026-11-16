      module diaga5_i
      interface
!
      subroutine diaga5(npeelgg,ja1,ka,ire,iat,recc)
      use vast_kind_param,  only:  double
      integer, intent(in)       :: npeelgg,ja1,ka,ire
!      INTEGER, INTENT(OUT)      :: IAT
      integer, intent(inout)      :: iat
      real(double), intent(out) :: recc
      end subroutine
      end interface
      end module
