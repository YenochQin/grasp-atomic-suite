      module angdata_i
      implicit none (type, external)
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      subroutine angdata (name,avail,whichone,status,cache_unit)
      implicit none (type, external)
      character (len = 24), intent(in) :: name
      logical, intent(out) :: avail
      integer, intent(in)  :: whichone
      integer, intent(out) :: status
      integer, intent(out) :: cache_unit
      end subroutine angdata
      end interface
      end module angdata_i
