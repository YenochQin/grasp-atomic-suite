      module trsort_i
      implicit none (type, external)
      interface
      subroutine trsort (name, nfile, nfile2, lprint, jkp, iblki, iblkf)
      implicit none (type, external)
      character (len = 24), intent(in) :: name(2)
      integer, intent(in) :: nfile
      integer, intent(in) :: nfile2
      logical, intent(in) :: lprint
      integer :: jkp
      integer, intent(in) :: iblki
      integer, intent(in) :: iblkf
      end subroutine
      end interface
      end module
