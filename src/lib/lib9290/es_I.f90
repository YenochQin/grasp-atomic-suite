      module es_i
      implicit none (type, external)
      interface
      subroutine es (f, s2f, s3f)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(in) :: f
      real(real64), intent(out) :: s2f
      real(real64), intent(out) :: s3f
      end subroutine
      end interface
      end module
