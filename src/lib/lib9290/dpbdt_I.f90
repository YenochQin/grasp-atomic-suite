      module dpbdt_i
      implicit none (type, external)
      interface
      subroutine dpbdt (j)
      implicit none (type, external)
      integer, intent(in) :: j
!VAST.../GRID/ N(IN)
!VAST.../LIC13/ A(IN)
!VAST.../TATB/ TA(OUT), TB(OUT)
!VAST...Calls: PF, QF
      end subroutine
      end interface
      end module
