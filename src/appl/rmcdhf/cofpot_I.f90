      module cofpot_i
      implicit none (type, external)
      interface
      subroutine cofpot (eol, j, npts)
      implicit none (type, external)
      logical :: eol
      integer :: j
      integer :: npts
      end subroutine
      end interface
      end module
