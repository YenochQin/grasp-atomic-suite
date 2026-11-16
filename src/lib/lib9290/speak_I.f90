      module speak_i
      implicit none (type, external)
      interface
      subroutine speak (ja, jb, ia1, ib1, ia2, ib2, k, x)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ja
      integer, intent(in) :: jb
      integer, intent(in) :: ia1
      integer, intent(in) :: ib1
      integer, intent(in) :: ia2
      integer, intent(in) :: ib2
      integer, intent(in) :: k
      real(real64), intent(in) :: x
!VAST.../BUFFER/ NBDIM(IN), NVCOEF(INOUT)
!VAST.../DEBUG/ IBUG1(IN)
!VAST.../ORB4/ NP(IN)
!VAST.../ORB10/ NH(IN)
!VAST...Calls: ALCBUF
!...This routine performs I/O.
      end subroutine
      end interface
      end module
