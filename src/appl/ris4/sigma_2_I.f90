      module sigma_2_i
      implicit none (type, external)
      interface
      subroutine sigma_2 (ipar,i1,i2,apart)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ipar, i1, i2
      real(real64), intent(out) :: apart
      end subroutine sigma_2
      end interface
      end module
