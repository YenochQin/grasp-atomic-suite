      module improvmpi_i
      implicit none (type, external)
      interface
      subroutine improvmpi (eol, j, lsort, dampmx)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      logical, intent(in) :: eol
      integer, intent(in) :: j
      logical :: lsort
      real(real64), intent(inout) :: dampmx
      end subroutine
      end interface
      end module
