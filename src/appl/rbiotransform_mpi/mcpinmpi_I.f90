      module mcpin_i
      interface
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      subroutine mcpin (name, startdir, ik, ntestg, inpci)
      character (len = 24), intent(in) :: name
      character(len=128), intent(in) :: startdir
      integer, intent(in) :: ik
      integer, intent(in) :: ntestg
      integer :: inpci
      end subroutine
      end interface
      end module
