      module fzalf_i
      implicit none (type, external)
      interface
      real(real64) function fzalf (principal_n, kappa, nuclear_charge)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: principal_n
      integer, intent(in) :: kappa
      real(real64) :: nuclear_charge
      end function fzalf
      end interface
      end module fzalf_i
