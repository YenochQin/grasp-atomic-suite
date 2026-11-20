      module rinthf_i
      implicit none (type, external)
      interface
      real(real64) function rinthf (first_orbital, second_orbital, k)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: first_orbital
      integer, intent(in) :: second_orbital
      integer, intent(in) :: k
      end function rinthf
      end interface
      end module rinthf_i
