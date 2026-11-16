      module recoonescalar_i
      implicit none (type, external)
      interface
!
      subroutine recoonescalar(ns,ja1,ja2,ja3,ja4,ka,iat)
      implicit none (type, external)
      integer, intent(in)  :: ns, ja1, ja2, ja3, ja4, ka
      integer, intent(out) :: iat
      end subroutine
      end interface
      end module
