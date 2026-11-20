      module copvec_i
      implicit none (type, external)
      interface
      subroutine copvec (from, to, ndim)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ndim
      real(real64), intent(in) :: from(ndim)
      real(real64), intent(out) :: to(ndim)
      end subroutine copvec
      end interface
      end module copvec_i
