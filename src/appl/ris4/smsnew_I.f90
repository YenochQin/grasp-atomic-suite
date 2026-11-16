      module smsnew_i
      implicit none (type, external)
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      subroutine smsnew (doit,vint,vint2)
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: max_orbitals
      implicit none (type, external)
      real(real64), intent(in) :: vint(max_orbitals, max_orbitals), vint2(max_orbitals, max_orbitals)
      integer, intent(in) :: doit
      end subroutine smsnew
      end interface
      end module smsnew_i
