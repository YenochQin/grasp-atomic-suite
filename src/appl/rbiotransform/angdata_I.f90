      module angdata_i
      implicit none (type, external)
      interface
      subroutine angdata (name, avail, kamax)
      implicit none (type, external)
      character (len = 24), intent(inout) :: name
      logical, intent(out) :: avail
      integer, intent(in) :: kamax
      end subroutine
      end interface
      end module
