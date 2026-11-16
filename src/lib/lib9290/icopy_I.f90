      module icopy_i
      implicit none (type, external)
      interface
      subroutine icopy (n, ix, incx, iy, incy)
      implicit none (type, external)
      integer, intent(in) :: n
      integer, intent(in) :: ix(*)
      integer, intent(in) :: incx
!VAST...Dummy argument INCX is not referenced in this routine.
      integer, intent(out) :: iy(*)
      integer, intent(in) :: incy
!VAST...Dummy argument INCY is not referenced in this routine.
      end subroutine
      end interface
      end module
