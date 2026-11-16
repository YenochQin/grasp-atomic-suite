      module smsnew_mpi_i
      interface
      subroutine smsnew_mpi (doit,vint,vint2)
      use vast_kind_param, only: double
      use parameter_def, only: max_orbitals
      integer, intent(in) :: doit
      real(double), dimension(max_orbitals,max_orbitals), intent(in) :: vint, vint2
      end subroutine
      end interface
      end module
