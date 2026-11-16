      module ifnmnx_i
      implicit none (type, external)
      interface
      integer function ifnmnx (ivec, nel, imxmn)
      implicit none (type, external)
      integer, intent(in) :: nel
      integer, intent(in) :: ivec(nel)
      integer, intent(in) :: imxmn
!...This routine performs I/O.
      end function ifnmnx
      end interface
      end module ifnmnx_i
