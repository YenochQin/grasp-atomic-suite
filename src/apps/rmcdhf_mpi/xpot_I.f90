      module xpot_i
      implicit none (type, external)
      interface
      subroutine xpot (orbital_index)
      implicit none (type, external)
      integer, intent(in) :: orbital_index
!...This routine performs I/O.
      end subroutine xpot
      end interface
      end module xpot_i
