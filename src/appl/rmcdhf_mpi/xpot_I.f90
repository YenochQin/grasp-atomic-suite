      module xpot_i
      implicit none (type, external)
      interface
      subroutine xpot (j)
      implicit none (type, external)
      integer, intent(in) :: j
!...This routine performs I/O.
      end subroutine xpot
      end interface
      end module xpot_i
