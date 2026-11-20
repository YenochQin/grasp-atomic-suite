      module rintdens_i
      implicit none (type, external)
      interface
      real(real64) function rintdens (first_orbital, second_orbital)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer :: first_orbital
      integer :: second_orbital
      end function rintdens
      end interface
      end module rintdens_i
