      module inprod_i
      implicit none (type, external)
      interface
      real(real64) function inprod (a, b, ndim)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(in) :: a(*)
      real(real64), intent(in) :: b(*)
      integer, intent(in) :: ndim
      end function inprod
      end interface
      end module inprod_i
