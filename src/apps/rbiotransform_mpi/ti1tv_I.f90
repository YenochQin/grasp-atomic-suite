      module ti1tv_i
      implicit none (type, external)
      interface
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      subroutine ti1tv(input_mixing_coefficients,csf_count,level_count,i,l,t,subshell_count,transformed_mixing_coefficients,ntestg)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: csf_count, level_count, subshell_count, i, l, ntestg
      real(real64), intent(in) :: input_mixing_coefficients(csf_count, level_count)
      real(real64), intent(in) :: t(subshell_count)
      real(real64), intent(out) :: transformed_mixing_coefficients(csf_count, level_count)
      end subroutine ti1tv
      end interface
      end module ti1tv_i
