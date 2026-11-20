      module rintii_i
      implicit none (type, external)
      interface
      real(real64) function rintii (first_orbital, second_orbital, k)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: first_orbital
      integer, intent(in) :: second_orbital
      integer, intent(in) :: k
      end function rintii
      end interface
      end module rintii_i
