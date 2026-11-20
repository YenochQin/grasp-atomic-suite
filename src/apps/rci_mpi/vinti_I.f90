      module vinti_i
      implicit none (type, external)
      interface
      real(real64) function vinti (first_orbital, second_orbital)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: first_orbital
      integer, intent(in) :: second_orbital
      end function vinti
      end interface
      end module vinti_i
