      module angdata_i
      implicit none (type, external)
      interface
      subroutine angdata (name, avail, jkp, nfile2)
      implicit none (type, external)
      character (len = 24), intent(inout) :: name(2)
      logical, intent(out) :: avail
      integer :: jkp
      integer, intent(in) :: nfile2
      end subroutine
      end interface
      end module
