      module smsread_mpi_i
      interface
      subroutine smsread_mpi (vint,vint2)
      use vast_kind_param, only: double
      use parameter_def, only: max_orbitals
      real(double), dimension(max_orbitals,max_orbitals), intent(in) :: vint, vint2
      end subroutine
      end interface
      end module
