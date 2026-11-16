      module mcpin_i
      implicit none (type, external)
      interface
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      subroutine mcpin (name, startdir, ik, ntestg, inpci)
      implicit none (type, external)
      character (len = 24), intent(in) :: name
      character(len=128), intent(in) :: startdir
      integer, intent(in) :: ik
      integer, intent(in) :: ntestg
      integer, intent(in) :: inpci
      end subroutine
      end interface
      end module
