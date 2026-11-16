      module diaga3_i
      implicit none (type, external)
      interface
!
      subroutine diaga3(ja1,ja2,ka,ire,iat,recc)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in)       :: ja1,ja2,ka,ire
!      INTEGER, INTENT(OUT)      :: IAT
      integer, intent(inout)      :: iat
      real(real64), intent(out) :: recc
      end subroutine
      end interface
      end module
