      module alclla_i
      implicit none (type, external)
      interface
      subroutine alclla(ibeg, ilab, ilast, ileft, iptcs, irigh, lblin, &
         lldim, imode)
      implicit none (type, external)
      integer, intent(inout) :: lldim
      integer, intent(in) :: imode
      integer, pointer :: ibeg(:), ilab(:), ilast(:), ileft(:), &
                                        iptcs(:), irigh(:), lblin(:)
!...This routine performs I/O.
      end subroutine
      end interface
      end module
