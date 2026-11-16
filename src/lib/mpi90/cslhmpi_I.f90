      module cslhmpi_i
      implicit none (type, external)
      interface
!
      subroutine cslhmpi (name, ncore, nblkin, idblk)
      implicit none (type, external)
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      character (len = *) :: name
      integer :: ncore
      integer, intent(in) :: nblkin
      character (len = 8) :: idblk(*)
      end subroutine
      end interface
      end module
