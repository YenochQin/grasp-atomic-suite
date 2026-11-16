      module densread_i
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      subroutine densread (dint1,dint2,dint3,                          &
                           dint4,dint5,dint6,                          &
                           dint7)
      use vast_kind_param,  only: double
      use parameter_def,    only: max_orbitals, max_radial_points
      use prnt_c,           only: nvec
      real(double), dimension(max_orbitals,max_orbitals), intent(in) :: dint1, dint2,  &
                                          dint3, dint4, dint5, dint6,  &
                                          dint7
      end subroutine
      end interface
      end module
