      module bessj_i
      implicit none (type, external)
      interface
      subroutine bessj (w)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(in) :: w
      end subroutine
      end interface
      end module
