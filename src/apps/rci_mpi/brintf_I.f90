      module brintf_i
      implicit none (type, external)
      interface
      real(real64) function brintf (itype, orbital_a, orbital_b, orbital_c, orbital_d, k)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: itype
      integer, intent(in) :: orbital_a
      integer, intent(in) :: orbital_b
      integer, intent(in) :: orbital_c
      integer, intent(in) :: orbital_d
      integer, intent(in) :: k
      end function brintf
      end interface
      end module brintf_i
