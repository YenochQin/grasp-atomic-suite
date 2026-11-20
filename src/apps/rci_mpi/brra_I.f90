      module brra_i
      implicit none (type, external)
      interface
      real(real64) function brra (itype, orbital_a, orbital_c, orbital_b, orbital_d, k)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: itype
      integer, intent(in) :: orbital_a
      integer, intent(in) :: orbital_c
      integer, intent(in) :: orbital_b
      integer, intent(in) :: orbital_d
      integer :: k
      end function brra
      end interface
      end module brra_i
