      module yzk_i
      implicit none (type, external)
      interface
      subroutine yzk (k, i, j)
      implicit none (type, external)
      integer, parameter :: klimit = 20
      integer, intent(in) :: k
      integer, intent(in) :: i
      integer, intent(in) :: j
!VAST.../CNC5/ CNC5C(IN)
!VAST.../DEF4/ ACCY(IN)
!VAST.../GRID/ R(IN), RP(IN), N(IN)
!VAST.../NCC/ C1(IN), C2(IN), C3(IN)
!VAST.../TATB/ MTP(OUT)
!VAST.../WAVE/ MF(IN)
!VAST...Calls: PF, QF
!...This routine performs I/O.
      end subroutine
      end interface
      end module
