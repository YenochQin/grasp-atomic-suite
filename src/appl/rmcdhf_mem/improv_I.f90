      module improv_i
      implicit none (type, external)
      interface
      subroutine improv (eol, j, lsort, dampmx)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      logical, intent(in) :: eol
      integer, intent(in) :: j
      logical :: lsort
      real(real64), intent(inout) :: dampmx
      end subroutine improv
      end interface
      end module improv_i
