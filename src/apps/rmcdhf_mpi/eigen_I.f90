      module eigen_i
      implicit none (type, external)
      interface
      real(real64) function eigen (orbital_index)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: orbital_index
      end function eigen
      end interface
      end module eigen_i
