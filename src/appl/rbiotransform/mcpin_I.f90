      module mcpin_i
      implicit none (type, external)
      interface
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      subroutine mcpin (name, ik, ntestg, inpci)
      implicit none (type, external)
      character (len = 24), intent(in) :: name
      integer, intent(in) :: ik
      integer, intent(in) :: ntestg
      integer, intent(in) :: inpci
      end subroutine
      end interface
      end module
