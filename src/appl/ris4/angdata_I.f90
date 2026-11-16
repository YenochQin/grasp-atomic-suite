      module angdata_i
      implicit none (type, external)
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      subroutine angdata (name,avail,whichone)
      implicit none (type, external)
      character (len = 24), intent(in) :: name
      logical, intent(out) :: avail
      integer, intent(in)  :: whichone
      end subroutine
      end interface
      end module
