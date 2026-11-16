      module setlag_i
      implicit none (type, external)
      interface
      subroutine setlag (eol)
      implicit none (type, external)
      logical :: eol
!...This routine performs I/O.
      end subroutine setlag
      end interface
      end module setlag_i
