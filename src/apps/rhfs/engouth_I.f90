      module engouth_i
      implicit none (type, external)
      interface
      subroutine engouth (hamiltonian_energy_shift, e, jtot, ipar, level_indices, level_count, mode)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(in) :: hamiltonian_energy_shift
      integer, intent(in) :: level_count
      real(real64), intent(in) :: e(level_count)
      integer, intent(in) :: jtot(level_count)
      integer, intent(in) :: ipar(level_count)
      integer, intent(in) :: level_indices(level_count)
      integer, intent(in) :: mode
      end subroutine engouth
      end interface
      end module engouth_i
