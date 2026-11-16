      module sigma_1_i
      interface
      subroutine sigma_1 (ipar,i1,i2,apart)
      use vast_kind_param, only: double
      integer, intent(in) :: ipar, i1, i2
      real(double), intent(out) :: apart
      end subroutine sigma_1
      end interface
      end module
