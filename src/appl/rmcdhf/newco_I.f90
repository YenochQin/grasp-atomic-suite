      module newco_i
      implicit none (type, external)
      interface
      subroutine newco (sum)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(out) :: sum
      end subroutine
      end interface
      end module
