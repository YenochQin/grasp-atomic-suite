      module arctan_i
      implicit none (type, external)
      interface
      real(real64) function arctan (arg1, arg2)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(in) :: arg1
      real(real64), intent(in) :: arg2
!VAST.../DEF9/ PI(IN)
      end function
      end interface
      end module
