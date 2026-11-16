      module lodrwfmpi_i
      implicit none (type, external)
      interface
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      subroutine lodrwfmpi (ierror)
      implicit none (type, external)
      integer , intent(out) :: ierror
      end subroutine
      end interface
      end module
