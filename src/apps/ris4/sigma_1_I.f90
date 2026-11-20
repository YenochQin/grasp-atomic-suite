      module sigma_1_i
      implicit none (type, external)
      interface
      subroutine sigma_1 (ipar,i1,i2,apart)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ipar, i1, i2
      real(real64), intent(out) :: apart
      end subroutine sigma_1
      end interface
      end module sigma_1_i
