      module fzalf_i
      implicit none (type, external)
      interface
      real(real64) function fzalf (n, kappa, z)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: n
      integer, intent(in) :: kappa
      real(real64) :: z
      end function fzalf
      end interface
      end module fzalf_i
