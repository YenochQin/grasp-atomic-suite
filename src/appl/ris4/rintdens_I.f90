      module rintdens_i
      implicit none (type, external)
      interface
      real(real64) function rintdens (i, j)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer :: i
      integer :: j
      end function rintdens
      end interface
      end module rintdens_i
