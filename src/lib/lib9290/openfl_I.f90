      module openfl_i
      implicit none (type, external)
      interface
      subroutine openfl (nfile, filnam, rform, rstat, ierr)
      implicit none (type, external)
      integer, intent(in) :: nfile
      character (len = *), intent(in) :: filnam
      character (len = *), intent(in) :: rform
      character (len = *), intent(in) :: rstat
      integer, intent(out) :: ierr
!VAST.../IOUNIT/ ISTDE(IN)
!...This routine performs I/O.
      end subroutine
      end interface
      end module
