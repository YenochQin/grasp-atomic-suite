      module klamaq_i
      implicit none (type, external)
      interface
      subroutine klamaq (n, kappa, z, fzalfa)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: n
      integer, intent(in) :: kappa
      real(real64), intent(in) :: z
      real(real64), intent(out) :: fzalfa
      end subroutine klamaq
      end interface
      end module klamaq_i
