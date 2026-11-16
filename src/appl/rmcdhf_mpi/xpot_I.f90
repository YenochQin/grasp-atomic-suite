      module xpot_i
      implicit none (type, external)
      interface
      subroutine xpot (j)
      implicit none (type, external)
      integer, intent(in) :: j
!...This routine performs I/O.
      end subroutine
      end interface
      end module
