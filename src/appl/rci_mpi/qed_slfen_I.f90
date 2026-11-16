      module qed_slfen_i
      implicit none (type, external)
      interface
      subroutine qed_slfen (slfint)
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: max_orbitals
      implicit none (type, external)
      real(real64), intent(out) :: slfint(max_orbitals)
      end subroutine
      end interface
      end module
