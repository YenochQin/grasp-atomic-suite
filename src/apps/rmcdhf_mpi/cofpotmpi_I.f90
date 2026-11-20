      module cofpotmpi_i
      implicit none (type, external)
      interface
      subroutine cofpotmpi (eol, j, npts)
      implicit none (type, external)
      logical :: eol
      integer :: j
      integer :: npts
      end subroutine cofpotmpi
      end interface
      end module cofpotmpi_i
