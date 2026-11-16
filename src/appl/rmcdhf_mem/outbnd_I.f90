      module outbnd_i
      implicit none (type, external)
      interface
      logical function outbnd (etry)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(in) :: etry
      end function
      end interface
      end module
