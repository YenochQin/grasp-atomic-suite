      module intrpqf_i
      implicit none (type, external)
      interface
      subroutine intrpqf (pa, qa, ma, ra, j, dnorm)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(in) :: pa(*)
      real(real64), intent(in) :: qa(*)
      integer, intent(in) :: ma
      real(real64), intent(in) :: ra(*)
      integer, intent(in) :: j
      real(real64), intent(out) :: dnorm
      end subroutine intrpqf
      end interface
      end module intrpqf_i
