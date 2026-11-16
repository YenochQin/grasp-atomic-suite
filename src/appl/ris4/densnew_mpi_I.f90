      module densnew_mpi_i
      interface
      subroutine densnew_mpi (doit,dint1,dint2,dint3,dint4,dint5,dint6,dint7)
      use vast_kind_param, only: double
      use parameter_def, only: max_orbitals
      integer, intent(in) :: doit
      real(double), dimension(max_orbitals,max_orbitals), intent(in) :: dint1, dint2, &
                                                        dint3, dint4, &
                                                        dint5, dint6, &
                                                        dint7
      end subroutine
      end interface
      end module
