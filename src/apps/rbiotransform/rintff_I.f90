      module rintff_i
      implicit none (type, external)
      interface
      real(real64) function rintff (first_orbital, second_orbital, k)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: first_orbital
      integer, intent(in) :: second_orbital
      integer, intent(in) :: k
      end function rintff
      end interface
      end module rintff_i
