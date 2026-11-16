      module smsnew_mpi_i
      interface
      subroutine smsnew_mpi (doit,vint,vint2)
      use vast_kind_param, only: double
      use parameter_def, only: nnnw
      integer, intent(in) :: doit
      real(double), dimension(nnnw,nnnw), intent(in) :: vint, vint2
      end subroutine
      end interface
      end module
