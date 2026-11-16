      module systemmem_i
      implicit none (type, external)
      interface
      real(real64) function systemmem ()
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!...This routine performs I/O.
      end function
      end interface
      end module
