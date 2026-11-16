      module vinti_i
      implicit none (type, external)
      interface
      real(real64) function vinti (j, k)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: j
      integer, intent(in) :: k
      end function
      end interface
      end module
