      module reco4_i
      implicit none (type, external)
      interface
!
      subroutine reco4(ja1,ja2,ja3,ja4,k1,k2,k3,k4,ka,ire,iat,recc)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in)       :: ja1,ja2,ja3,ja4,k1,k2,k3,k4,ka,ire
      integer, intent(out)      :: iat
      real(real64), intent(out) :: recc
      end subroutine
      end interface
      end module
