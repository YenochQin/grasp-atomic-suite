      module auxblk_i
      implicit none (type, external)
      interface
      subroutine auxblk (max_orbital_two_j, atwinv)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: max_orbital_two_j
      real(real64), intent(out) :: atwinv
      end subroutine auxblk
      end interface
      end module auxblk_i
