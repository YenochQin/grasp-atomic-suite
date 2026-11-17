      module qed_i
      implicit none (type, external)
      interface
      subroutine qed (jstate, slfint, ucf)
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: max_orbitals
      implicit none (type, external)
      integer, intent(in) :: jstate
      real(real64), intent(out) :: slfint(max_orbitals)
      real(real64), intent(out) :: ucf(:)
      end subroutine qed
      end interface
      end module qed_i
