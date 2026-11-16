      module recop00_i
      implicit none (type, external)
      interface
!
      subroutine recop00(ns,ja1,ja2,ka,iat)
      implicit none (type, external)
      integer, intent(in)  :: ns, ja1, ja2, ka
      integer, intent(out) :: iat
      end subroutine
      end interface
      end module
