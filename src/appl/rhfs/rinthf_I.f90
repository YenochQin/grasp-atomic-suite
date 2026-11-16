      module rinthf_i
      implicit none (type, external)
      interface
      real(real64) function rinthf (i, j, k)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      end function rinthf
      end interface
      end module rinthf_i
