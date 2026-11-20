      module spme_i
      implicit none (type, external)
      interface
      subroutine spme (first_orbital, second_orbital, hcoul, hbab, hmag)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: first_orbital
      integer, intent(in) :: second_orbital
      real(real64), intent(out) :: hcoul
      real(real64), intent(out) :: hbab
      real(real64), intent(out) :: hmag
      end subroutine spme
      end interface
      end module spme_i
