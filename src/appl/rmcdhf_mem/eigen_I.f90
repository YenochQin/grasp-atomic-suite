      module eigen_i
      implicit none (type, external)
      interface
      real(real64) function eigen (j)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: j
      end function eigen
      end interface
      end module eigen_i
