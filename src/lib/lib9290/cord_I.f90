      module cord_i
      implicit none (type, external)
      interface
      subroutine cord (ja, jb, ja1, ipca, jb1)
      implicit none (type, external)
      integer :: ja
      integer :: jb
      integer, intent(in) :: ja1
      integer, intent(in) :: ipca
      integer, intent(in) :: jb1
!VAST.../ORB4/ NAK(IN)
!VAST.../M1/ NQ1(IN)
!VAST.../M3/ JLIST(IN), KLIST(IN)
!VAST...Calls: SPEAK, CLRX
      end subroutine
      end interface
      end module
