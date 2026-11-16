      module cslh_i
      implicit none (type, external)
      interface
      subroutine cslh (name, ncore, nblkin, idblk)
      implicit none (type, external)
      character (len = *), intent(inout) :: name
      integer, intent(out) :: ncore
      integer, intent(in) :: nblkin
      character (len = 8), intent(out) :: idblk(*)
!VAST.../HBLOCK/ NBLOCK(INOUT), PNCFBLK(INOUT)
!VAST.../ORB2/ NCFTOT(INOUT)
!VAST...Calls: ALLOC, SETCSLL, NCFBLK, RALLOC, LODCSH
!...This routine performs I/O.
      end subroutine
      end interface
      end module
