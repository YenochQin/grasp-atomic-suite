      module setcsll_i
      implicit none (type, external)
      interface
      subroutine setcsll (nunit, name, nblkin, nblock, ncfblk, ncftot, idblk)
      implicit none (type, external)
      integer, intent(in) :: nunit
      character (len = *), intent(inout) :: name
      integer, intent(in) :: nblkin
      integer, intent(out) :: nblock
      integer, intent(inout) :: ncfblk(*)
      integer, intent(out) :: ncftot
      character (len = 8), intent(out) :: idblk(*)
!VAST...Calls: OPENFL
!...This routine performs I/O.
      end subroutine
      end interface
      end module
