      module pack_i
      implicit none (type, external)
      interface
      subroutine pack (iunpkd, isubsh, ipackd)
      use, intrinsic :: iso_fortran_env, only: int8
      implicit none (type, external)
      integer, intent(in) :: iunpkd
      integer, intent(in) :: isubsh
      integer(int8), intent(inout) :: ipackd(*)
!VAST.../IOUNIT/ ISTDE(IN)
!...This routine performs I/O.
      end subroutine
      end interface
      end module
