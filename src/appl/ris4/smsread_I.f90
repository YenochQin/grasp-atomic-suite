      module smsread_i
      implicit none (type, external)
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      subroutine smsread (vint,vint2)
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: max_orbitals
      implicit none (type, external)
      real(real64), intent(in) :: vint(max_orbitals, max_orbitals), vint2(max_orbitals, max_orbitals)
      end subroutine smsread
      end interface
      end module smsread_i
