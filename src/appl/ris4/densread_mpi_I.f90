      module densread_mpi_i
      interface
      subroutine densread_mpi (dint1,dint2,dint3,dint4,dint5,dint6,dint7)
      use vast_kind_param, only: double
      use parameter_def, only: nnnw
      real(double), dimension(nnnw,nnnw), intent(in) :: dint1, dint2,  &
                                          dint3, dint4, dint5, dint6,  &
                                          dint7
      end subroutine
      end interface
      end module
