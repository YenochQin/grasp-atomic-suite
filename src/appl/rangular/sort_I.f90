      module sort_i
      implicit none (type, external)
      interface
      subroutine sort (nfile, ncoeff, ntgral, lprint, nb, fhead)
      implicit none (type, external)
      integer, intent(in)  :: nfile
      integer, intent(in)  :: ncoeff
      integer, intent(out) :: ntgral
      logical, intent(in)  :: lprint
      integer, intent(in)  :: nb
      character (len = *), intent(in) :: fhead
      end subroutine
      end interface
      end module
