      module systemfreemem_i
      implicit none (type, external)
      interface
!   Written by  G. Gaigalas               Vilnius, September 2021
      real(real64) function systemfreemem ()
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      end function
      end interface
      end module
