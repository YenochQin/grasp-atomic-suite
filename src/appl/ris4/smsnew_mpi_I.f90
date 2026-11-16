      module smsnew_mpi_i
      implicit none (type, external)
      interface
      subroutine smsnew_mpi (doit,vint,vint2)
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def, only: max_orbitals
      implicit none (type, external)
      integer, intent(in) :: doit
      real(real64), intent(in) :: vint(max_orbitals, max_orbitals), vint2(max_orbitals, max_orbitals)
      end subroutine smsnew_mpi
      end interface
      end module smsnew_mpi_i
