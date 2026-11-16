      module cslhmpi_i
      interface
!
      subroutine cslhmpi (name, ncore, nblkin, idblk)
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      character (len = *) :: name
      integer :: ncore
      integer, intent(in) :: nblkin
      character (len = 8), dimension(*) :: idblk
      end subroutine
      end interface
      end module
