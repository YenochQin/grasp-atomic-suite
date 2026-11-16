      module sigma_2_i
      interface
      subroutine sigma_2 (ipar,i1,i2,apart)
      use vast_kind_param, only: double
      integer, intent(in) :: ipar, i1, i2
      real(double), intent(out) :: apart
      end subroutine sigma_2
      end interface
      end module
