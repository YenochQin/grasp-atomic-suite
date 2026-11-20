      module citrag_i
      implicit none (type, external)
      interface
      subroutine citrag (input_mixing_coefficients, csf_count, level_count, l, subshell_count, t, nin, ntestg, &
         & transformed_mixing_coefficients, scr)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: csf_count
      integer, intent(in) :: level_count
      real(real64) :: input_mixing_coefficients(csf_count, level_count)
      integer, intent(in) :: l
      integer, intent(in) :: subshell_count
      real(real64), intent(in) :: t(subshell_count, subshell_count)
      integer, intent(in) :: nin
      integer, intent(in) :: ntestg
      real(real64) :: transformed_mixing_coefficients(csf_count, level_count)
      real(real64) :: scr(csf_count, level_count)
      end subroutine citrag
      end interface
      end module citrag_i
