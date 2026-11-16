      module densnew_mpi_i
      interface
      subroutine densnew_mpi (doit,dint1,dint2,dint3,dint4,dint5,dint6,dint7)
      use vast_kind_param, only: double
      use parameter_def, only: nnnw
      integer, intent(in) :: doit
      real(double), dimension(nnnw,nnnw), intent(in) :: dint1, dint2, &
                                                        dint3, dint4, &
                                                        dint5, dint6, &
                                                        dint7
      end subroutine
      end interface
      end module
