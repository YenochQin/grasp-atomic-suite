      module tiinig_i
      implicit none (type, external)
      interface
      subroutine tiinig (input_mixing_coefficients, csf_count, level_count, i, l, const, transformed_mixing_coefficients, ntestg)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: csf_count
      integer, intent(in) :: level_count
      real(real64), intent(in) :: input_mixing_coefficients(csf_count, level_count)
      integer, intent(in) :: i
      integer, intent(in) :: l
      real(real64), intent(in) :: const
      real(real64), intent(inout) :: transformed_mixing_coefficients(csf_count, level_count)
      integer, intent(in) :: ntestg
      end subroutine tiinig
      end interface
      end module tiinig_i
